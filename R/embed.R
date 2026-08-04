# Two embedding paths:
# - TM-Vec1 (below): a custom trained head bolted onto a frozen ProtT5 backbone,
#   verified against a published sanity check -- see the docstring on embed_proteins().
# - generic (bottom): plain AutoModel + mean pooling, for any other HuggingFace
#   encoder repo. No calibration claim, just a mean-pooled last hidden state.
.repdist_embed_py <- r"(
import re, json, torch
from functools import lru_cache
from torch import nn
from safetensors.torch import load_file
from huggingface_hub import hf_hub_download
from transformers import T5Tokenizer, T5EncoderModel, AutoTokenizer, AutoModel

class _RepdistTMVec1Head(nn.Module):
    def __init__(self, cfg):
        super().__init__()
        layer = nn.TransformerEncoderLayer(
            d_model=cfg["d_model"], nhead=cfg["nhead"],
            dim_feedforward=cfg["dim_feedforward"], dropout=cfg["dropout"],
            activation=cfg["activation"], batch_first=True)
        # enable_nested_tensor's fast path calls an op MPS doesn't implement
        # (aten::_nested_tensor_from_mask_left_aligned) -- off regardless of device.
        self.encoder = nn.TransformerEncoder(layer, num_layers=cfg["num_layers"],
                                              enable_nested_tensor=False)
        self.mlp = nn.Linear(cfg["d_model"], cfg["out_dim"])

    def forward(self, x, pad):
        x = self.encoder(x, src_key_padding_mask=pad)
        x = x.masked_fill(pad.unsqueeze(-1), 0.0)
        pooled = x.sum(1) / torch.logical_not(pad).sum(1, keepdim=True)
        return self.mlp(pooled)

@lru_cache(maxsize=None)
def _repdist_tmvec1_models(head_repo, backbone_repo, device):
    cfg = json.load(open(hf_hub_download(head_repo, "config.json")))
    head = _RepdistTMVec1Head(cfg)
    state = load_file(hf_hub_download(head_repo, "model.safetensors"))
    head.load_state_dict(state, strict=True)
    head.eval().to(device)
    tok = T5Tokenizer.from_pretrained(backbone_repo, do_lower_case=False)
    backbone = T5EncoderModel.from_pretrained(backbone_repo).eval().to(device)
    return tok, backbone, head

# ProtT5 expects space-separated residues, rare amino acids mapped to X.
def _repdist_tmvec1(seqs, head_repo, backbone_repo, device):
    tok, backbone, head = _repdist_tmvec1_models(head_repo, backbone_repo, device)
    prepped = [" ".join(re.sub(r"[UZOB]", "X", s)) for s in seqs]
    enc = tok(prepped, padding=True, return_tensors="pt")
    ids, attn = enc["input_ids"].to(device), enc["attention_mask"].to(device)
    with torch.no_grad():
        h = backbone(input_ids=ids, attention_mask=attn).last_hidden_state
        lens = attn.sum(1) - 1   # tokenizer appends </s>; residues are [0, n)
        L = int(lens.max())
        x = torch.zeros(h.shape[0], L, h.shape[2], device=h.device, dtype=torch.float32)
        pad = torch.ones(h.shape[0], L, dtype=torch.bool, device=h.device)
        for i, n in enumerate(lens.tolist()):
            x[i, :n] = h[i, :n].float()
            pad[i, :n] = False
        return head(x, pad).cpu().numpy()

@lru_cache(maxsize=None)
def _repdist_generic_models(model_repo, device):
    tok = AutoTokenizer.from_pretrained(model_repo)
    model = AutoModel.from_pretrained(model_repo).eval().to(device)
    return tok, model

# Mean-pools the residue tokens -- a naive default that works for any encoder,
# not a calibrated ground metric like the TM-Vec1 head above. Pads and the
# tokenizer's own markers (<cls>/<eos>/<sep>) are excluded: they carry no residue
# and averaging them in makes short sequences look alike.
def _repdist_generic(seqs, model_repo, device):
    tok, model = _repdist_generic_models(model_repo, device)
    enc = tok(list(seqs), padding=True, return_tensors="pt",
              return_special_tokens_mask=True)
    special = enc.pop("special_tokens_mask").to(device)
    enc = {k: v.to(device) for k, v in enc.items()}
    with torch.no_grad():
        h = model(**enc).last_hidden_state
    mask = (enc["attention_mask"] * (1 - special)).unsqueeze(-1).float()
    return ((h * mask).sum(1) / mask.sum(1).clamp(min=1)).cpu().numpy()
)"

#' Read protein sequences from a FASTA file
#'
#' @param path Path to an amino-acid FASTA file.
#' @return Named character vector, names taken from the first whitespace-delimited
#'   token of each header.
#' @export
read_fasta <- function(path) {
  aa <- Biostrings::readAAStringSet(path)
  stats::setNames(as.character(aa), sub("\\s.*", "", names(aa)))
}

#' Models `embed_proteins()` assembles from more than one HuggingFace repo
#'
#' A trained head and the frozen backbone it was trained on are two separate
#' repos, and the head's own `config.json` does not record which backbone it
#' expects. That pairing lives in `inst/extdata/models.json` -- add an entry
#' there to teach [embed_proteins()] a new head, no code change needed.
#'
#' @return Named list, one entry per known head repo, each with `head` (the
#'   architecture to build) and `backbone` (the repo to run underneath it).
#' @export
known_models <- function() {
  path <- system.file("extdata", "models.json", package = "repdist")
  # system.file() only resolves once repdist is installed; this repo's own
  # analyses source() R/*.R directly from the package root instead (see the
  # analysis/*.Rmd setup chunks), so fall back to that same relative layout.
  if (!nzchar(path)) path <- "inst/extdata/models.json"
  jsonlite::fromJSON(path, simplifyVector = FALSE)
}

#' Embed protein sequences with a HuggingFace protein language model
#'
#' Dispatches on `model`. A repo listed in [known_models()] is assembled with
#' its registered head over its registered backbone -- by default TM-Vec1,
#' whose head regresses the TM-score of protein pairs, so cosine similarity
#' between two vectors is a predicted TM-score (see
#' `docs/ground_metric_models.md`).
#'
#' Any other repo is embedded directly: loaded with `AutoModel`/`AutoTokenizer`
#' and represented by the mean-pooled last hidden state. No calibration claim
#' -- unlike TM-Vec1, cosine similarity there isn't a named, verified quantity.
#'
#' Downloads are cached under `~/.cache/huggingface`.
#'
#' Needs `torch`, `transformers`, `sentencepiece`, `safetensors` and
#' `huggingface_hub` importable from the Python reticulate is using (see
#' `reticulate::py_config()`). Errors with the exact `pip install` command if
#' any are missing -- install them yourself rather than relying on this
#' function to provision an environment.
#'
#' @param seqs Named character vector of amino-acid sequences, **or** a path to a
#'   FASTA file (see [read_fasta()]).
#' @param model HuggingFace repo. Either a head listed in [known_models()] (the
#'   default TM-Vec1 head was trained on chains up to 300 residues; longer
#'   sequences still embed but are out of the training range), or any encoder
#'   repo, embedded generically.
#' @param device "mps", "cuda", or "cpu".
#' @param batch_size Maximum sequences per forward pass. The real limit is
#'   `token_budget`; this caps the short-sequence end.
#' @param token_budget Maximum batch x padded-length product per forward pass.
#'   A hardware knob, not a model property -- lower it if the device runs out
#'   of memory, raise it to use a bigger card more fully.
#' @return Numeric matrix, one row per sequence, rownames = names(seqs).
#' @export
embed_proteins <- function(seqs,
                           model = "scikit-bio/tmvec-swissmodel",
                           device = c("mps", "cuda", "cpu"),
                           batch_size = 16,
                           token_budget = 16384L) {
  device <- match.arg(device)
  if (is.character(seqs) && length(seqs) == 1L && file.exists(seqs))
    seqs <- read_fasta(seqs)
  stopifnot(!is.null(names(seqs)), !anyDuplicated(names(seqs)))

  Sys.setenv(KMP_DUPLICATE_LIB_OK = "TRUE")  # macOS: duplicate OpenMP runtimes abort at import

  needed <- c("torch", "transformers", "sentencepiece", "safetensors", "huggingface_hub")
  missing <- needed[!vapply(needed, reticulate::py_module_available, logical(1))]
  if (length(missing) > 0L)
    stop("embed_proteins() needs these Python packages, not found in ",
         reticulate::py_config()$python, ":\n  pip install ",
         paste(missing, collapse = " "), call. = FALSE)

  main <- reticulate::import_main(convert = FALSE)
  if (!reticulate::py_has_attr(main, "_repdist_generic"))
    reticulate::py_run_string(.repdist_embed_py)
  torch <- reticulate::import("torch")
  torch$set_grad_enabled(FALSE)  # simpler than bridging no_grad()'s context manager
  on.exit(torch$set_grad_enabled(TRUE), add = TRUE)

  prepped <- unname(toupper(seqs))
  lens <- nchar(seqs)
  out <- NULL

  spec <- known_models()[[model]]
  if (is.null(spec)) {
    embed_fn <- reticulate::py$`_repdist_generic`
    extra_args <- list(model, device)
  } else {
    stopifnot(spec$head == "tmvec1")   # only head built by .repdist_embed_py so far
    embed_fn <- reticulate::py$`_repdist_tmvec1`
    extra_args <- list(model, spec$backbone, device)
  }

  # Batch by token budget, not by count. Padding makes every sequence in a batch as
  # long as the longest, and attention is quadratic in that length, so a fixed count
  # varies ~100x in peak memory between the shortest and longest batches -- which is
  # what exhausts MPS on a length-sorted universe. `batch_size` caps the short end.
  chunks <- list()
  i <- 1L
  while (i <= length(seqs)) {
    j <- i
    while (j < length(seqs) && (j - i + 2L) * max(lens[i:(j + 1L)]) <= token_budget &&
           (j - i + 2L) <= batch_size) j <- j + 1L
    chunks[[length(chunks) + 1L]] <- i:j
    i <- j + 1L
  }

  for (idx in chunks) {
    arr <- reticulate::py_to_r(do.call(embed_fn, c(list(as.list(prepped[idx])), extra_args)))
    if (is.null(out)) out <- matrix(NA_real_, nrow = length(seqs), ncol = ncol(arr))
    out[idx, ] <- arr
    # torch's device allocators cache freed blocks; over a whole universe that cache
    # grows past the limit and allocation fails even though nothing is live.
    if (device == "mps") torch$mps$empty_cache() else if (device == "cuda") torch$cuda$empty_cache()
  }

  rownames(out) <- names(seqs)
  colnames(out) <- paste0("d", seq_len(ncol(out)))
  out
}
