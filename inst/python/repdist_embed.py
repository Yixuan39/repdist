"""Protein-embedding backends for the repdist R package.

Called from R through reticulate inside a basilisk environment; `embed()` is the
only entry point. Two paths, both returning one L2-normalised row per sequence:
`tmvec1`, a trained TM-Vec head over a frozen ProtT5 backbone whose cosine
similarity is a predicted TM-score, and `generic`, mean-pooled AutoModel hidden
states for the registry's plain encoder models, with no such calibration claim.
"""

import json
import os
import re
import warnings
from functools import lru_cache

# Must precede the torch import: on macOS a second OpenMP runtime in the process
# aborts at load time. setdefault, so an explicit user setting still wins.
os.environ.setdefault("KMP_DUPLICATE_LIB_OK", "TRUE")

import numpy as np  # noqa: E402
import torch  # noqa: E402
from huggingface_hub import hf_hub_download  # noqa: E402
from safetensors.torch import load_file  # noqa: E402
from torch import nn  # noqa: E402
from transformers import (  # noqa: E402
    AutoModel, AutoTokenizer, T5EncoderModel, T5Tokenizer)

def resolve_device(requested):
    """Map "auto" onto the best available accelerator, or check an explicit one."""
    mps = getattr(torch.backends, "mps", None)
    available = {"cuda": torch.cuda.is_available(),
                 "mps": mps is not None and mps.is_available(),
                 "cpu": True}
    if requested == "auto":
        return next(name for name, ok in available.items() if ok)
    if not available[requested]:
        raise ValueError(
            "requested device '{}' is not available (torch reports {})".format(
                requested,
                ", ".join("{} = {}".format(*kv) for kv in available.items())))
    return requested


class _TMVecHead(nn.Module):
    """TM-Vec's transformer head over frozen per-residue backbone states."""

    def __init__(self, cfg):
        super().__init__()
        layer = nn.TransformerEncoderLayer(
            d_model=cfg["d_model"], nhead=cfg["nhead"],
            dim_feedforward=cfg["dim_feedforward"], dropout=cfg["dropout"],
            activation=cfg["activation"], batch_first=True)
        # enable_nested_tensor's fast path calls an op MPS does not implement
        # (aten::_nested_tensor_from_mask_left_aligned) -- off on every device.
        self.encoder = nn.TransformerEncoder(
            layer, num_layers=cfg["num_layers"], enable_nested_tensor=False)
        self.mlp = nn.Linear(cfg["d_model"], cfg["out_dim"])

    def forward(self, x, pad):
        x = self.encoder(x, src_key_padding_mask=pad)
        x = x.masked_fill(pad.unsqueeze(-1), 0.0)
        return self.mlp(x.sum(1) / torch.logical_not(pad).sum(1, keepdim=True))


# Cached so one embed() call loads the multi-GB backbone once, not per batch.

@lru_cache(maxsize=None)
def _load_prott5(backbone, device):
    tok = T5Tokenizer.from_pretrained(backbone, do_lower_case=False)
    return tok, T5EncoderModel.from_pretrained(backbone).eval().to(device)


@lru_cache(maxsize=None)
def _load_generic(repo, device):
    return (AutoTokenizer.from_pretrained(repo),
            AutoModel.from_pretrained(repo).eval().to(device))


@lru_cache(maxsize=None)
def _load_head(weights, device, cfg_json):
    """An empty cfg_json means `weights` is a repo carrying its own config.json."""
    if cfg_json:
        cfg = json.loads(cfg_json)
        state = torch.load(weights, map_location="cpu", weights_only=True)["state_dict"]
    else:
        cfg = json.load(open(hf_hub_download(weights, "config.json")))
        state = load_file(hf_hub_download(weights, "model.safetensors"))
    head = _TMVecHead(cfg)
    # strict=True checks every shape cfg implies; only nhead and activation,
    # which leave no shape behind, would pass silently if the registry were wrong.
    head.load_state_dict(state, strict=True)
    return head.eval().to(device)


def _run_tmvec(seqs, device, backbone, head):
    tok, body = _load_prott5(backbone, device)
    # ProtT5 wants space-separated residues, rare amino acids mapped to X.
    enc = tok([" ".join(re.sub(r"[UZOB]", "X", s)) for s in seqs],
              padding=True, return_tensors="pt")
    ids, attn = enc["input_ids"].to(device), enc["attention_mask"].to(device)
    h = body(input_ids=ids, attention_mask=attn).last_hidden_state
    lens = (attn.sum(1) - 1).tolist()   # tokenizer appends </s>; residues are [0, n)
    width = int(max(lens))
    x = torch.zeros(len(lens), width, h.shape[2], device=device)
    pad = torch.ones(len(lens), width, dtype=torch.bool, device=device)
    for i, n in enumerate(lens):
        x[i, :n] = h[i, :n].float()
        pad[i, :n] = False
    return head(x, pad)


def _run_generic(seqs, device, tok, model):
    enc = tok(list(seqs), padding=True, return_tensors="pt",
              return_special_tokens_mask=True)
    special = enc.pop("special_tokens_mask").to(device)
    enc = {k: v.to(device) for k, v in enc.items()}
    h = model(**enc).last_hidden_state
    # Pads and the tokenizer's own markers (<cls>/<eos>/<sep>) hold no residue;
    # averaging them in would make short sequences look alike.
    mask = (enc["attention_mask"] * (1 - special)).unsqueeze(-1).float()
    return (h * mask).sum(1) / mask.sum(1).clamp(min=1)


def _batches(lengths, batch_size):
    """Yield (start, stop) index pairs partitioning the sequences in order.

    Padding stretches every sequence in a batch to the longest one, so batches
    are capped by padded token count as well as by `batch_size`. A batch that
    still exceeds device memory is split and retried by the caller, so these
    caps are throughput tuning, not correctness.
    """
    if not lengths:
        return
    start, longest = 0, 0
    for i, n in enumerate(lengths):
        widest = max(longest, n)
        if i > start and (i - start + 1 > batch_size
                          or (i - start + 1) * widest > 16384):
            yield start, i
            start, longest = i, n
        else:
            longest = widest
    yield start, len(lengths)


def _is_oom(error):
    return (isinstance(error, MemoryError)
            or "out of memory" in str(error).lower()
            or "allocate memory" in str(error).lower())


@torch.no_grad()
def embed(seqs, spec_json, weights, device, batch_size):
    """Embed a list of amino-acid strings into one L2-normalised row each.

    `spec_json` is one entry of the model registry (inst/extdata/models.json);
    the R side rejects anything not in it. Its `head` picks the path and says
    what `weights` is: unused for "generic",
    a HuggingFace head repo for "tmvec1", or a local Lightning checkpoint
    (already hash-verified on the R side) for "tmvec1-large" -- which, carrying
    no config.json of its own, takes its architecture from the registry's
    `config`. Rows come back in input order, one per sequence: a sequence that
    cannot be embedded raises rather than being dropped.
    """
    spec = json.loads(spec_json)
    head, backbone = spec["head"], spec["backbone"]
    seqs = list(seqs)
    batch_size = int(batch_size)
    device = resolve_device(device)

    def make_run(dev):
        """The embedding closure for one device. Loads are lru_cached, so
        asking for "cpu" after "cuda" costs a second copy of the backbone in
        host RAM -- only paid if the OOM fallback below actually fires."""
        if head == "generic":
            tok, model = _load_generic(backbone, dev)
            return lambda batch: _run_generic(batch, dev, tok, model)
        if head in ("tmvec1", "tmvec1-large"):
            cfg_json = json.dumps(spec["config"]) if "config" in spec else ""
            top = _load_head(weights, dev, cfg_json)
            _load_prott5(backbone, dev)
            return lambda batch: _run_tmvec(batch, dev, backbone, top)
        raise ValueError("unknown head type '{}'".format(head))

    run = make_run(device)

    lengths = [len(s) for s in seqs]

    def run_batch(batch):
        try:
            return [run(batch).cpu().numpy()]
        except (RuntimeError, MemoryError) as error:
            if not _is_oom(error):
                raise
            if device != "cpu":
                getattr(torch, device).empty_cache()
            if len(batch) == 1:
                if device == "cpu":
                    raise MemoryError(
                        "a {}-residue sequence does not fit in host memory "
                        "even in a batch of 1; use a smaller model".format(
                            len(batch[0]))) from error
                # Splitting has nothing left to give, but host RAM is usually
                # far larger than device memory, so run this one sequence on
                # the CPU rather than failing the whole catalog.
                warnings.warn(
                    "a {}-residue sequence does not fit in {} memory; "
                    "embedding it on the CPU instead, which is much slower"
                    .format(len(batch[0]), device))
                return [make_run("cpu")(batch).cpu().numpy()]
            middle = len(batch) // 2
            return run_batch(batch[:middle]) + run_batch(batch[middle:])

    pieces = []
    for start, stop in _batches(lengths, batch_size):
        pieces.extend(run_batch(seqs[start:stop]))
        # torch's device allocators cache freed blocks; over a whole protein
        # universe that cache outgrows the device even though nothing is live.
        if device != "cpu":
            getattr(torch, device).empty_cache()

    # Raw TM-Vec norms vary with protein length, so only the normalised inner
    # product is the calibrated cosine the head predicts a TM-score from
    # (see docs/ground_metric_models.md).
    out = np.concatenate(pieces)
    return out / np.linalg.norm(out, axis=1, keepdims=True)


if __name__ == "__main__":
    # `python -m repdist_embed`: the only logic checkable without a model.
    assert list(_batches([10] * 17, 16)) == [(0, 16), (16, 17)]
    assert list(_batches([10] * 5, 2)) == [(0, 2), (2, 4), (4, 5)]
    assert list(_batches([9000, 9000], 16)) == [(0, 1), (1, 2)]
    assert list(_batches([50000], 16)) == [(0, 1)]
    assert list(_batches([], 16)) == []
    assert _is_oom(MemoryError())
    assert _is_oom(RuntimeError("CUDA out of memory."))
    assert not _is_oom(RuntimeError("shape mismatch"))
    print("ok")
