# Everything here is plain R. The torch code is inst/python/repdist_embed.py,
# run inside the pinned basilisk environment defined in R/basilisk.R.

#' Read and validate protein sequences from a FASTA file
#'
#' @param path Path to an amino-acid FASTA file. At most 46341 sequences;
#'   dereplicate a larger catalog before reading it.
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
  .repdist_qc_sequences(aa)
}

#' Registered protein embedding models
#'
#' Head/backbone pairings live in `inst/extdata/models.json`, keyed by model;
#' register another there without touching R or Python. These keys are exactly
#' the models [embed_proteins()] supports: unregistered HuggingFace repos are
#' rejected, since the pinned Python environment fixes which architectures can
#' be loaded at all.
#'
#' @return Named list keyed by model, each entry carrying `description` (what
#'   the embedding measures and when to pick this model), `head` (architecture,
#'   or `"generic"`), `backbone` (HuggingFace repo), and `training_max_length`
#'   -- recorded for reference, never enforced. A `tmvec1` head
#'   carries `head_repo`, the HuggingFace repo its weights live in; a head with
#'   no HuggingFace release instead carries the `url`, `sha256`, and `config`
#'   needed to fetch and rebuild it.
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
.repdist_embed_worker <- function(seqs, spec, weights, device, batch_size) {
  py <- reticulate::import_from_path(
    "repdist_embed",
    path = system.file("python", package = "repdist", mustWork = TRUE))
  py$embed(seqs, spec, weights, device, batch_size)
}

#' Embed protein sequences with a HuggingFace protein language model
#'
#' Dispatches on `model`. A TM-Vec head regresses the TM-score of a protein
#' pair, so cosine similarity between two rows is a predicted TM-score;
#' the generic presets are mean-pooled `AutoModel` hidden states, with no such
#' calibration claim. `basilisk` provisions the pinned Python environment on
#' first use, and model weights are cached (HuggingFace models under
#' `~/.cache/huggingface`, the large TM-Vec head under `BiocFileCache`).
#'
#' Every sequence is embedded whole and in input order: nothing is truncated,
#' subset, or silently dropped, and no length limit is imposed -- a sequence
#' past a model's `training_max_length` is embedded, just outside the range it
#' was fitted on. A sequence too large for device memory even in a batch of its
#' own is an error naming the offender. Batches are capped at `batch_size`
#' sequences (and an internal padded-token limit); a batch that still exceeds
#' device memory is split and retried automatically.
#'
#' @param seqs Named character vector of amino-acid sequences, **or** a path to a
#'   FASTA file (see [read_fasta()]). Names must be unique and non-empty. A
#'   terminal stop `*` is removed; gaps, whitespace, internal stops, non-standard
#'   letters like `J`, empty sequences, and missing values are rejected. At most
#'   46341 sequences -- dereplicate a larger catalog before embedding it.
#' @param model A key from [known_models()]. Unregistered HuggingFace repos are
#'   rejected: the pinned Python environment fixes which architectures can be
#'   loaded at all, so a repo outside the registry has no supported-model
#'   guarantee behind it. Register one in `inst/extdata/models.json` instead.
#'   `"tmvec"` (the default) and `"tmvec-300"` measure structural similarity --
#'   their cosine is a predicted TM-score; the
#'   `"esm2-*"` models (by parameter count: `-8m`, `-150m`, `-650m`) measure
#'   general sequence features instead, for questions where structure is not
#'   the axis of interest. Each registry entry's `description` says when to
#'   pick it.
#' @param device `"auto"` picks CUDA, then MPS, then CPU. An explicit `"mps"`,
#'   `"cuda"` or `"cpu"` errors if unavailable rather than falling back. A
#'   device too small for the catalog is not an error: batches that exhaust
#'   device memory are split and retried, and a single sequence that still does
#'   not fit is embedded on the CPU with a warning.
#' @param batch_size Maximum number of sequences per model batch. A throughput
#'   knob, not a correctness one: lower it to reduce peak memory when sharing
#'   a device with other processes.
#' @return Numeric matrix with one unit-norm row per input sequence, in input
#'   order. Rows are L2-normalised, so the inner product of two rows is their
#'   cosine similarity.
#' @examples
#' if (interactive()) {
#'     # esm2-8m is 8M parameters; the default downloads several GB.
#'     seqs <- c(p1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
#'               p2 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVE")
#'     Z <- embed_proteins(seqs, model = "esm2-8m")
#'     dim(Z)
#' }
#' @export
embed_proteins <- function(
    seqs, model = "tmvec",
    device = c("auto", "mps", "cuda", "cpu"), batch_size = 16L) {
  device <- match.arg(device)
  if (length(batch_size) != 1L || !is.finite(batch_size) || batch_size < 1)
    stop("`batch_size` must be one number >= 1.", call. = FALSE)
  if (is.character(seqs) && length(seqs) == 1L && file.exists(seqs)) {
    seqs <- read_fasta(seqs)
  } else {
    seqs <- .repdist_qc_sequences(seqs)
  }
  ids <- names(seqs)

  registry <- known_models()
  spec <- registry[[model]]
  if (is.null(spec))
    stop("Unknown model `", model, "`. Registered models are: ",
         paste(names(registry), collapse = ", "),
         ". Unregistered HuggingFace repos are not supported -- the pinned ",
         "Python environment fixes which architectures load; add an entry to ",
         "inst/extdata/models.json to register one.", call. = FALSE)
  weights <- switch(spec$head,
    generic        = "",
    tmvec1         = spec$head_repo,
    `tmvec1-large` = .repdist_cached_weights(spec, model),
    stop("Unknown head type `", spec$head, "` for model `", model, "`.",
         call. = FALSE))
  if (is.null(weights))
    stop("Registry entry `", model, "` has head `", spec$head,
         "` but no `head_repo` naming its weights.", call. = FALSE)

  proc <- basilisk::basiliskStart(repdist_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  Z <- basilisk::basiliskRun(proc, .repdist_embed_worker,
    as.list(unname(seqs)),
    as.character(jsonlite::toJSON(spec, auto_unbox = TRUE)),
    weights, device, as.integer(batch_size))

  dimnames(Z) <- list(ids, paste0("d", seq_len(ncol(Z))))
  Z
}
