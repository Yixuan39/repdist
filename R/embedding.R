# Everything here is plain R. The torch code is inst/python/repdist_embed.py,
# run inside the pinned basilisk environment defined in R/basilisk.R.

#' Read and validate protein sequences from a FASTA file
#'
#' @param path Path to an amino-acid FASTA file.
#' @return Named character vector, names taken from the first
#'   whitespace-delimited token of each header. Sequences are uppercased and
#'   terminal stops are removed. Invalid, empty, duplicate-named, and
#'   model-incompatible records are rejected.
#' @examples
#' fa <- tempfile(fileext = ".fasta")
#' writeLines(c(">p1 first protein", "ACDEFGHIKL", ">p2", "MNPQRSTVWY"), fa)
#' read_fasta(fa)
#' @export
read_fasta <- function(path) {
    aa <- withCallingHandlers(
        Biostrings::readAAStringSet(path),
        warning = function(w) stop(conditionMessage(w), call. = FALSE)
    )
    names(aa) <- sub("\\s.*", "", names(aa))
    .check_sequences(aa)
}

# Sample and protein names are how counts, embeddings and similarity matrices
# are matched, so every entry point that takes them checks them here.
.check_ids <- function(ids, what) {
    if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) ||
        any(!nzchar(trimws(ids)))) {
        stop(what, " names must be unique, non-missing and non-empty.",
            call. = FALSE
        )
    }
}

# Validate once for every public entry point that accepts protein sequences,
# or a FASTA path. Returns the cleaned named character vector.
.check_sequences <- function(seqs) {
    if (is.character(seqs) && length(seqs) == 1L && file.exists(seqs)) {
        return(read_fasta(seqs))
    }
    ids <- names(seqs) # as.character() drops them, and AAStringSet in
    seqs <- as.character(seqs) # is as valid an input here as a character vector
    if (!length(seqs)) {
        stop("`seqs` must contain at least one sequence.", call. = FALSE)
    }
    .check_ids(ids, "`seqs`")
    if (anyNA(seqs)) {
        stop("`seqs` must not contain missing values.", call. = FALSE)
    }

    seqs <- toupper(sub("\\*$", "", seqs))
    if (any(empty <- !nzchar(seqs))) {
        stop("`seqs` contains empty sequence(s): ",
            paste(ids[empty], collapse = ", "),
            call. = FALSE
        )
    }
    # Gaps, internal stops and stray symbols reach the tokenizer as <unk> and
    # embed to noise without complaint; an aligned FASTA is the usual way in.
    if (any(off <- grepl("[^ACDEFGHIKLMNPQRSTVWYBXZUO]", seqs))) {
        stop("model-incompatible symbol(s): ", paste(ids[off], collapse = ", "),
            call. = FALSE
        )
    }
    stats::setNames(seqs, ids)
}

#' Registered protein embedding models
#'
#' Head/backbone pairings live in `inst/extdata/models.json`, keyed by
#' model; register another there without touching R or Python. These keys are
#' exactly the models [embed_proteins()]
#' supports: unregistered HuggingFace repos are rejected, since the pinned
#' Python environment fixes which architectures can be loaded at all.
#'
#' @return Named list keyed by model, each entry carrying `description`
#'   (what the embedding measures and when to pick this model), `head`
#'   (architecture, or `"generic"`), `backbone` (HuggingFace repo),
#'   and `training_max_length` -- recorded for reference, never enforced.
#'   `backbone_revision` pins the tokenizer and model to an immutable
#'   commit. A `tmvec1` head carries `head_repo` and its immutable
#'   `head_revision`; a head with no HuggingFace release instead carries
#'   the `url`, `sha256`, and `config` needed to fetch and
#'   rebuild it.
#' @examples
#' str(known_models())
#' @export
known_models <- function() {
    jsonlite::fromJSON(
        system.file("extdata", "models.json",
            package = "repdist", mustWork = TRUE),
        simplifyVector = FALSE
    )
}

# Weights for heads with no HuggingFace release, from the `url`/`sha256` the
# registry gives them. Hash-checked on every load because they are deserialised
# as torch weights; a mismatch also covers a truncated or partial download.
.cached_weights <- function(spec, model) {
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    path <- BiocFileCache::bfcrpath(bfc, spec$url)
    if (digest::digest(path, algo = "sha256", file = TRUE) != spec$sha256) {
        BiocFileCache::bfcremove(
            bfc, BiocFileCache::bfcquery(
                bfc, spec$url, "rname", exact = TRUE)$rid
        )
        stop("Cached weights for model `", model,
            "` failed their SHA-256 check. ",
            "They have been dropped from the cache; ",
            "retry to download them again.",
            call. = FALSE
        )
    }
    path
}

# The registry entry for `model`, with the weights its head needs: nothing for
# a plain encoder, a HuggingFace repo for a released head, a hash-checked
# download for one with no release.
.model_spec <- function(model) {
    if (!is.character(model) || length(model) != 1L || is.na(model) ||
        !nzchar(model)) {
        stop("`model` must be one non-missing string from known_models().",
            call. = FALSE
        )
    }
    registry <- known_models()
    spec <- registry[[model]]
    if (is.null(spec)) {
        stop("Unknown model `", model, "`. Registered models are: ",
            paste(names(registry), collapse = ", "),
            ". Unregistered HuggingFace repos are not supported -- the pinned ",
            "Python environment fixes which architectures load; ",
            "add an entry to ",
            "inst/extdata/models.json to register one.",
            call. = FALSE
        )
    }
    weights <- switch(spec$head,
        generic = "",
        tmvec1 = spec$head_repo,
        `tmvec1-large` = .cached_weights(spec, model),
        stop("Unknown head type `", spec$head, "` for model `", model, "`.",
            call. = FALSE
        )
    )
    if (is.null(weights)) {
        stop("Registry entry `", model, "` has head `", spec$head,
            "` but no `head_repo` naming its weights.",
            call. = FALSE
        )
    }
    list(spec = spec, weights = weights)
}

# Runs in the basilisk process, so Python is touched only here.
.embed_worker <- function(seqs, spec, weights, device, batch_size) {
    py <- reticulate::import_from_path(
        "repdist_embed",
        path = system.file("python", package = "repdist", mustWork = TRUE)
    )
    py$embed(seqs, spec, weights, device, batch_size)
}

#' Embed protein sequences with a HuggingFace protein language model
#'
#' Dispatches on `model`. A TM-Vec head regresses the TM-score of a
#' protein pair, so cosine similarity between two rows is a predicted
#' TM-score; the generic presets are mean-pooled `AutoModel` hidden
#' states, with no such calibration claim. `basilisk` provisions the
#' pinned Python environment on first use, and model weights are cached
#' (HuggingFace models under `~/.cache/huggingface`, the large TM-Vec
#' head under `BiocFileCache`). HuggingFace weights and tokenizers are
#' pinned to commit revisions, and Figshare weights are checked against a
#' SHA-256 digest. Floating-point differences can still occur across devices
#' and numerical libraries.
#'
#' Every sequence is embedded whole and in input order: nothing is truncated,
#' subset, or silently dropped, and no length limit is imposed -- a sequence
#' past a model's `training_max_length` is embedded, just outside the
#' range it was fitted on. A sequence too large for device memory even in a
#' batch of its own is an error naming the offender. Batches are capped at
#' `batch_size` sequences (and an internal padded-token limit); a batch
#' that still exceeds device memory is split and retried automatically.
#'
#' @param seqs Named character vector of amino-acid sequences, **or** a
#'   path to a FASTA file (see [read_fasta()]). Names
#'   must be unique, non-missing and non-empty. A terminal stop `*` is
#'   removed; gaps, whitespace, internal stops, non-standard letters like
#'   `J`, empty sequences, and missing values are rejected.
#' @param model A key from [known_models()].
#'   Unregistered HuggingFace repos are rejected: the pinned Python
#'   environment fixes which architectures can be loaded at all, so a repo
#'   outside the registry has no supported-model guarantee behind it. Register
#'   one in `inst/extdata/models.json` instead. `"tmvec"` (the
#'   default) and `"tmvec-300"` measure structural similarity -- their
#'   cosine is a predicted TM-score; the `"esm2-*"` models (by parameter
#'   count: `-8m`, `-150m`, `-650m`) measure general sequence
#'   features instead, for questions where structure is not the axis of
#'   interest. Each registry entry's `description` says when to pick it.
#' @param device `"auto"` picks CUDA, then MPS, then CPU. An explicit
#'   `"mps"`, `"cuda"` or `"cpu"` errors if unavailable rather
#'   than falling back. A device too small for the catalog is not an error:
#'   batches that exhaust device memory are split and retried, and a single
#'   sequence that still does not fit is embedded on the CPU with a warning.
#' @param batch_size Maximum number of sequences per model batch, a positive
#'   integer. A throughput knob, not a correctness one: lower it to reduce
#'   peak memory when sharing a device with other processes.
#' @return An `AAStringSet` of the sequences, named by accession, with one
#'   unit-norm embedding row per protein in the matrix `mcols(x)$embedding`.
#'   To use embeddings from elsewhere, build the same object:
#'   `x <- Biostrings::AAStringSet(seqs)` then
#'   `S4Vectors::mcols(x)$embedding <- Z`.
#' @usage embed_proteins(
#'     seqs,
#'     model = "tmvec",
#'     device = c("auto", "mps", "cuda", "cpu"),
#'     batch_size = 16L
#' )
#' @references Hamamsy et al. (2024). Protein remote homology detection and
#'   structural alignment using deep learning. Nature Biotechnology.
#'   \doi{10.1038/s41587-023-01917-2}.
#' @examples
#' if (interactive()) {
#'     # esm2-8m is 8M parameters; the default downloads several GB.
#'     seqs <- c(
#'         p1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
#'         p2 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVE"
#'     )
#'     x <- embed_proteins(seqs, model = "esm2-8m")
#'     dim(S4Vectors::mcols(x)$embedding)
#' }
#' @export
embed_proteins <- function(
    seqs, model = "tmvec",
    device = c("auto", "mps", "cuda", "cpu"), batch_size = 16L
) {
    device <- match.arg(device)
    if (!is.numeric(batch_size) || length(batch_size) != 1L ||
        !is.finite(batch_size) || batch_size < 1 ||
        batch_size != floor(batch_size) || batch_size > .Machine$integer.max) {
        stop("`batch_size` must be a single positive integer.", call. = FALSE)
    }
    resolved <- .model_spec(model)
    seqs <- .check_sequences(seqs)
    ids <- names(seqs)

    proc <- basilisk::basiliskStart(repdist_env)
    on.exit(basilisk::basiliskStop(proc), add = TRUE)
    Z <- basilisk::basiliskRun(
        proc, .embed_worker,
        as.list(unname(seqs)),
        as.character(jsonlite::toJSON(resolved$spec, auto_unbox = TRUE)),
        resolved$weights, device, as.integer(batch_size)
    )

    if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != length(ids) ||
        ncol(Z) < 1L || any(!is.finite(Z)) ||
        any(abs(rowSums(Z^2) - 1) > 1e-6)) {
        stop("Embedding backend returned invalid or non-unit rows.",
            call. = FALSE)
    }
    dimnames(Z) <- list(ids, paste0("d", seq_len(ncol(Z))))
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- Z
    x
}

#' Read embeddings computed outside R
#'
#' Reads a `build-db` database (`.npz`) from the tmvec command-line tool
#' (<https://github.com/valentynbez/tmvec>), or an HDF5 file in the layout
#' of `tmvec embed` -- one group per protein holding datasets `seq` and
#' `emb` -- which is the way in for embeddings from any other model. A
#' per-residue `emb` matrix is averaged over residues; rows are L2-normalised.
#'
#' @param path Path to a `.npz` database or an HDF5 file.
#' @param seqs For a `.npz`, which stores no sequences: a FASTA path or named
#'   character vector covering its proteins. Defaults to the FASTA the
#'   database was built from. Ignored for HDF5, which stores them.
#' @return As [embed_proteins()].
#' @examples
#' if (interactive()) {
#'     x <- read_embeddings("db/proteins.npz", seqs = "proteins.fasta")
#' }
#' @export
read_embeddings <- function(path, seqs = NULL) {
    if (grepl("\\.npz$", path, ignore.case = TRUE)) {
        Z <- .read_npy(path, "embeddings")
        rownames(Z) <- .read_npy(path, "headers")
        if (is.null(seqs)) {
            seqs <- .read_npy(path, "input_fasta")
            if (!file.exists(seqs)) {
                stop("The database's FASTA, ", seqs, ", is missing; pass ",
                    "`seqs`.", call. = FALSE)
            }
        }
    } else {
        h5 <- rhdf5::h5dump(path)
        Z <- do.call(rbind, lapply(h5, function(g) rowMeans(as.matrix(g$emb))))
        seqs <- vapply(h5, function(g) as.character(g$seq), "")
    }
    .check_ids(rownames(Z), "Embedded protein")
    seqs <- .check_sequences(seqs)
    missing <- setdiff(rownames(Z), names(seqs))
    if (length(missing)) {
        stop(length(missing), " embedded protein(s) have no sequence, e.g. ",
            missing[[1]], ".", call. = FALSE)
    }
    x <- Biostrings::AAStringSet(seqs[rownames(Z)])
    S4Vectors::mcols(x)$embedding <- Z / sqrt(rowSums(Z^2))
    x
}

# One array from a NumPy .npz archive, as numpy.savez writes it: C order,
# little-endian float or unicode.
.read_npy <- function(path, name) {
    con <- unz(path, paste0(name, ".npy"), "rb")
    on.exit(close(con))
    magic <- readBin(con, "raw", 8L) # \x93NUMPY, then the format version
    size <- if (magic[[7]] == as.raw(1L)) 2L else 4L
    header <- rawToChar(readBin(con, "raw",
        readBin(con, "integer", size = size, endian = "little")))
    descr <- sub(".*'descr': *'([^']*)'.*", "\\1", header)
    dims <- sub(".*'shape': *\\(([^)]*)\\).*", "\\1", header)
    dims <- as.integer(strsplit(dims, ", *")[[1]])
    n <- prod(dims)
    width <- as.integer(substring(descr, 3L))
    x <- switch(substr(descr, 2L, 2L),
        f = readBin(con, "double", n, size = width, endian = "little"),
        U = vapply(split(
            readBin(con, "integer", n * width, size = 4L, endian = "little"),
            rep(seq_len(n), each = width)
        ), function(u) intToUtf8(u[u > 0L]), ""),
        stop("Unsupported array type ", descr, " in ", path, ".",
            call. = FALSE)
    )
    if (length(dims) == 2L) matrix(x, dims[[1]], byrow = TRUE) else unname(x)
}

# The unit-norm embedding matrix of an embed_proteins() result, rows named by
# accession.
.embeddings <- function(x) {
    Z <- if (inherits(x, "AAStringSet")) S4Vectors::mcols(x)$embedding
    if (!is.matrix(Z) || !is.numeric(Z) || nrow(Z) != length(x)) {
        stop("`x` must be an AAStringSet with an `embedding` matrix in ",
            "mcols(), as returned by embed_proteins().", call. = FALSE)
    }
    .check_ids(names(x), "Protein")
    rownames(Z) <- names(x)
    Z <- Z / sqrt(rowSums(Z^2))
    stopifnot(all(is.finite(Z)))
    Z
}
