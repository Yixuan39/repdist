# Everything here is plain R. The torch code is inst/python/repdist_embed.py,
# run inside the pinned basilisk environment defined in R/basilisk.R.

# Validate once for every public entry point that accepts protein sequences.
.repdist_qc_sequences <- function(seqs, required_names = NULL) {
  if (!length(seqs))
    stop("`seqs` must contain at least one sequence.", call. = FALSE)
  ids <- names(seqs)
  if (is.null(ids) || anyDuplicated(ids) || any(!nzchar(ids)))
    stop("`seqs` names must be unique and non-empty.", call. = FALSE)
  missing <- setdiff(required_names, ids)
  if (length(missing))
    stop("`seqs` is missing ", length(missing), " protein(s), e.g. ",
         paste(utils::head(missing, 5), collapse = ", "), call. = FALSE)
  if (anyNA(seqs))
    stop("`seqs` must not contain missing values.", call. = FALSE)

  aa <- Biostrings::AAStringSet(sub("\\*$", "", as.character(seqs)))
  names(aa) <- ids
  empty <- Biostrings::width(aa) == 0L
  if (any(empty))
    stop("`seqs` contains empty sequence(s): ",
         paste(ids[empty], collapse = ", "), call. = FALSE)
  bad <- rowSums(Biostrings::alphabetFrequency(aa)[,
    c("J", "*", "-", "+", "."), drop = FALSE]) > 0L
  if (any(bad))
    stop("model-incompatible symbol(s): ", paste(ids[bad], collapse = ", "),
         call. = FALSE)
  aa
}

#' Read and validate protein sequences from a FASTA file
#'
#' @param path Path to an amino-acid FASTA file.
#' @return Named character vector, names taken from the first whitespace-delimited
#'   token of each header. Sequences are uppercased and terminal stops are
#'   removed. Invalid, empty, duplicate-named, and model-incompatible records
#'   are rejected.
#' @examples
#' fa <- tempfile(fileext = ".fasta")
#' writeLines(c(">p1 first protein", "ACDEFGHIKL", ">p2", "MNPQRSTVWY"), fa)
#' read_fasta(fa)
#' @export
read_fasta <- function(path) {
  aa <- withCallingHandlers(
    Biostrings::readAAStringSet(path),
    warning = function(w) stop(conditionMessage(w), call. = FALSE))
  names(aa) <- sub("\\s.*", "", names(aa))
  as.character(.repdist_qc_sequences(aa))
}

# Either the caller named an array that must exist, or we take the first
# conventional name the archive happens to carry.
.repdist_npz_key <- function(key, keys, defaults) {
  if (is.null(key)) {
    key <- intersect(defaults, keys)[1L]
    if (is.na(key))
      stop("NPZ archive must contain an array named one of: ",
           paste0("`", defaults, "`", collapse = ", "), ".", call. = FALSE)
  } else if (!key %in% keys) {
    stop("NPZ archive has no `", key, "` array; available arrays: ",
         paste(keys, collapse = ", "), call. = FALSE)
  }
  key
}

.repdist_read_npz_worker <- function(path, embedding_key, id_key) {
  np <- reticulate::import("numpy", convert = FALSE)
  archive <- np$load(path, allow_pickle = FALSE)
  on.exit(archive$close(), add = TRUE)

  keys <- unlist(reticulate::py_to_r(archive$files), use.names = FALSE)
  embedding_key <- .repdist_npz_key(embedding_key, keys,
                                    c("embeddings", "embedding"))
  id_key <- .repdist_npz_key(id_key, keys, c("ids", "headers", "accession"))

  list(
    embeddings = reticulate::py_to_r(archive$get(embedding_key)),
    ids = unlist(reticulate::py_to_r(
      archive$get(id_key)$astype("str")$tolist()), use.names = FALSE)
  )
}

#' Read pre-computed protein embeddings
#'
#' Imports a NumPy `.npz` archive produced outside R and returns the same matrix
#' shape as [embed_proteins()]: proteins in rows and embedding dimensions in
#' columns. By default, `embeddings` or `embedding` identifies the
#' two-dimensional numeric array, and `ids`, `headers`, or `accession`
#' identifies the matching one-dimensional protein identifier array.
#'
#' NumPy object arrays are deliberately rejected because loading them requires
#' Python pickle. Save identifiers as a Unicode array, for example
#' `numpy.asarray(ids, dtype=str)`.
#'
#' @param path Path to a NumPy `.npz` archive.
#' @param embedding_key Name of the embedding array. `NULL` detects
#'   `embeddings`, then `embedding`.
#' @param id_key Name of the protein-identifier array. `NULL` detects `ids`,
#'   `headers`, then `accession`.
#' @return Numeric matrix with protein identifiers as row names.
#' @examples
#' if (interactive()) {
#'     Z <- read_embeddings("protein_embeddings.npz")
#'     sample_repdist(counts, Z)
#' }
#' @export
read_embeddings <- function(path, embedding_key = NULL, id_key = NULL) {
  if (length(path) != 1L || is.na(path) || !file.exists(path))
    stop("`path` must name one existing file.", call. = FALSE)
  if (tolower(tools::file_ext(path)) != "npz")
    stop("`path` must be a NumPy `.npz` archive.", call. = FALSE)
  if (!is.null(embedding_key) &&
      (length(embedding_key) != 1L || is.na(embedding_key) ||
       !nzchar(embedding_key)))
    stop("`embedding_key` must be NULL or one non-empty string.",
         call. = FALSE)
  if (!is.null(id_key) && (length(id_key) != 1L || is.na(id_key) ||
                           !nzchar(id_key)))
    stop("`id_key` must be NULL or one non-empty string.", call. = FALSE)

  proc <- basilisk::basiliskStart(.repdist_io_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  out <- basilisk::basiliskRun(
    proc, .repdist_read_npz_worker, normalizePath(path), embedding_key, id_key)

  Z <- out$embeddings
  ids <- as.character(out$ids)
  if (!is.matrix(Z) || !is.numeric(Z) || any(dim(Z) == 0L))
    stop("The embedding array must be non-empty, two-dimensional, and numeric.",
         call. = FALSE)
  if (length(ids) != nrow(Z))
    stop("The identifier array must have one entry per embedding row.",
         call. = FALSE)
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Embedding identifiers must be unique, non-missing, and non-empty.",
         call. = FALSE)
  if (!all(is.finite(Z)))
    stop("Embeddings must contain only finite values.", call. = FALSE)

  dimnames(Z) <- list(ids, paste0("d", seq_len(ncol(Z))))
  Z
}

#' Registered protein embedding models
#'
#' A trained head and the backbone it was trained on are separate repos, and the
#' head's own `config.json` does not record the pairing. Those pairings live in
#' `inst/extdata/models.json`, keyed by model; register another there without
#' touching R or Python.
#'
#' @return Named list keyed by model, each entry carrying `head` (architecture,
#'   or `"generic"`), `backbone` (HuggingFace repo), `training_max_length`
#'   (the range used to train the head) and, when the architecture
#'   has one, `context_max_length` (an architectural limit). A
#'   head with no HuggingFace release also carries the `url` its weights come
#'   from, the `sha256` they are verified against, and the `config` to rebuild
#'   the architecture with, since a bare checkpoint ships none.
#' @examples
#' str(known_models())
#' @export
known_models <- function() jsonlite::fromJSON(
  system.file("extdata", "models.json", package = "repdist", mustWork = TRUE),
  simplifyVector = FALSE)

# Weights for heads with no HuggingFace release, from the `url`/`sha256` the
# registry gives them. Hash-checked on every load because they are deserialised
# as torch weights; a mismatch also covers a truncated or partial download.
.repdist_cached_weights <- function(spec, model) {
  bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
  path <- BiocFileCache::bfcrpath(bfc, spec$url)
  if (digest::digest(path, algo = "sha256", file = TRUE) != spec$sha256) {
    BiocFileCache::bfcremove(
      bfc, BiocFileCache::bfcquery(bfc, spec$url, "rname", exact = TRUE)$rid)
    stop("Cached weights for model `", model, "` failed their SHA-256 check. ",
         "They have been dropped from the cache; retry to download them again.",
         call. = FALSE)
  }
  path
}

# Runs in the basilisk process, so Python is touched only here.
.repdist_embed_worker <- function(seqs, spec, weights, device, memory_fraction) {
  py <- reticulate::import_from_path(
    "repdist_embed",
    path = system.file("python", package = "repdist", mustWork = TRUE))
  py$embed(seqs, spec, weights, device, memory_fraction)
}

#' Embed protein sequences with a HuggingFace protein language model
#'
#' Dispatches on `model`. A TM-Vec entry in [known_models()] pairs its registered
#' head with its registered backbone; that head regresses the TM-score of a
#' protein pair, so cosine similarity between two rows is a predicted TM-score
#' (see `docs/ground_metric_models.md`). Generic presets and unregistered
#' HuggingFace repos are mean-pooled `AutoModel` hidden states instead, with no
#' such calibration claim.
#'
#' `basilisk` provisions the pinned Python environment on first use; no manual
#' setup is needed. Weights are cached outside it -- HuggingFace models under
#' `~/.cache/huggingface`, the large TM-Vec head under `BiocFileCache`.
#'
#' For TM-Vec, available memory is measured after its weights reach the selected
#' device. A conservative estimate of ProtT5 and head attention memory then
#' reduces each batch. A sequence that cannot fit alone, or exceeds a backbone's
#' architectural context, is skipped with a warning; sequences are never
#' truncated.
#' The model weights themselves must fit before this estimate can be made.
#' Linux container limits (cgroups), Windows physical memory, and macOS unified
#' memory are handled separately; CUDA reports its own device memory.
#' Generic HuggingFace encoders retain their own architecture and memory rules.
#'
#' @param seqs Named character vector of amino-acid sequences, **or** a path to a
#'   FASTA file (see [read_fasta()]). Names must be unique and non-empty. A
#'   terminal stop `*` is removed; gaps, whitespace, internal stops, non-standard
#'   letters like `J`, empty sequences, and missing values are rejected.
#' @param model A key from [known_models()] -- `"tmvec-swissmodel-large"` (the
#'   default, trained through 1000 residues), `"scikit-bio/tmvec-swissmodel"` (300
#'   residues, for in-range sensitivity analyses), `"esm2"`, `"esm2-small"` -- or
#'   any HuggingFace encoder repo, embedded generically.
#' @param device `"auto"` picks CUDA, then MPS, then CPU. An explicit `"mps"`,
#'   `"cuda"` or `"cpu"` errors if unavailable rather than falling back.
#' @param memory_fraction Fraction of currently available device memory that
#'   TM-Vec may use for estimated temporary tensors. Must be greater than zero
#'   and no greater than one. Lower it when sharing a device with other
#'   processes.
#' @return Numeric matrix with one unit-norm row per embedded sequence. Skipped
#'   sequences are omitted and named in a warning.
#' @examples
#' if (interactive()) {
#'     # esm2-small is 8M parameters; the default downloads several GB.
#'     seqs <- c(p1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
#'               p2 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVE")
#'     Z <- embed_proteins(seqs, model = "esm2-small")
#'     dim(Z)
#' }
#' @export
embed_proteins <- function(
    seqs, model = "tmvec-swissmodel-large",
    device = c("auto", "mps", "cuda", "cpu"), memory_fraction = 0.8) {
  device <- match.arg(device)
  valid_memory <- length(memory_fraction) == 1L && is.finite(memory_fraction) &&
    memory_fraction > 0 && memory_fraction <= 1
  if (!valid_memory) stop("`memory_fraction` must be one number in (0, 1].", call. = FALSE)
  if (is.character(seqs) && length(seqs) == 1L && file.exists(seqs)) {
    seqs <- read_fasta(seqs)
  } else {
    seqs <- as.character(.repdist_qc_sequences(seqs))
  }
  ids <- names(seqs)

  spec <- known_models()[[model]]
  if (is.null(spec)) spec <- list(head = "generic", backbone = model)
  weights <- switch(spec$head,
    generic        = "",
    tmvec1         = model,
    `tmvec1-large` = .repdist_cached_weights(spec, model),
    stop("Unknown head type `", spec$head, "` for model `", model, "`.",
         call. = FALSE))

  proc <- basilisk::basiliskStart(.repdist_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  result <- basilisk::basiliskRun(proc, .repdist_embed_worker,
    as.list(unname(seqs)),
    as.character(jsonlite::toJSON(spec, auto_unbox = TRUE)),
    weights, device, memory_fraction)

  keep <- as.integer(result$kept) + 1L
  skipped <- setdiff(seq_along(seqs), keep)
  if (length(skipped))
    warning(sprintf("Skipped %d sequence(s) over the model/memory limit: %s",
      length(skipped), paste(ids[skipped], collapse = ", ")), call. = FALSE)
  dimnames(result$embeddings) <- list(
    ids[keep], paste0("d", seq_len(ncol(result$embeddings))))
  result$embeddings
}
