"""Protein-embedding backends for the repdist R package.

Called from R through reticulate inside a basilisk environment; `embed()` is the
only entry point. Two paths, both returning one L2-normalised row per sequence:
`tmvec1`, a trained TM-Vec head over a frozen ProtT5 backbone whose cosine
similarity is a predicted TM-score, and `generic`, mean-pooled AutoModel hidden
states for any other HuggingFace encoder, with no such calibration claim.
"""

import json
import os
import re
import subprocess
import sys
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
        self.d_model = cfg["d_model"]
        self.nhead = cfg["nhead"]
        self.dim_feedforward = cfg["dim_feedforward"]
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


def _cgroup_remaining(limit, current):
    """Parse a Linux cgroup memory limit and current use."""
    if limit == "max":
        return None
    try:
        limit, current = int(limit), int(current)
    except ValueError:
        return None
    # Cgroup v1 represents "unlimited" with a value close to signed LONG_MAX.
    if limit >= 1 << 60:
        return None
    return max(0, limit - current)


def _linux_cgroup_available():
    """Return the current process's cgroup allowance, when constrained."""
    paths = {
        ("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
        ("/sys/fs/cgroup/memory/memory.limit_in_bytes",
         "/sys/fs/cgroup/memory/memory.usage_in_bytes"),
    }
    try:
        with open("/proc/self/cgroup", encoding="ascii") as handle:
            for line in handle:
                _, controllers, group = line.rstrip().split(":", 2)
                group = group.lstrip("/")
                if not controllers:
                    base = os.path.join("/sys/fs/cgroup", group)
                    paths.add((os.path.join(base, "memory.max"),
                               os.path.join(base, "memory.current")))
                elif "memory" in controllers.split(","):
                    base = os.path.join("/sys/fs/cgroup/memory", group)
                    paths.add((os.path.join(base, "memory.limit_in_bytes"),
                               os.path.join(base, "memory.usage_in_bytes")))
    except (OSError, ValueError):
        pass

    remaining = []
    for limit_path, current_path in paths:
        try:
            with open(limit_path, encoding="ascii") as handle:
                limit = handle.read().strip()
            with open(current_path, encoding="ascii") as handle:
                current = handle.read().strip()
        except OSError:
            continue
        value = _cgroup_remaining(limit, current)
        if value is not None:
            remaining.append(value)
    return min(remaining) if remaining else None


def _cpu_available_memory():
    """Return available physical memory using only platform facilities."""
    if sys.platform.startswith("linux"):
        available = None
        try:
            with open("/proc/meminfo", encoding="ascii") as handle:
                fields = dict(line.split(":", 1) for line in handle
                              if ":" in line)
            if "MemAvailable" in fields:
                available = int(fields["MemAvailable"].split()[0]) * 1024
        except (OSError, ValueError):
            pass
        cgroup = _linux_cgroup_available()
        if cgroup is not None:
            available = cgroup if available is None else min(available, cgroup)
        if available is not None:
            return available

    if os.name == "nt":
        import ctypes

        class MemoryStatus(ctypes.Structure):
            _fields_ = [
                ("length", ctypes.c_ulong),
                ("load", ctypes.c_ulong),
                ("total_physical", ctypes.c_ulonglong),
                ("available_physical", ctypes.c_ulonglong),
                ("total_page_file", ctypes.c_ulonglong),
                ("available_page_file", ctypes.c_ulonglong),
                ("total_virtual", ctypes.c_ulonglong),
                ("available_virtual", ctypes.c_ulonglong),
                ("available_extended_virtual", ctypes.c_ulonglong),
            ]

        status = MemoryStatus()
        status.length = ctypes.sizeof(status)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(
                ctypes.byref(status)):
            raise OSError("GlobalMemoryStatusEx failed")
        return int(status.available_physical)

    total = os.sysconf("SC_PAGE_SIZE") * os.sysconf("SC_PHYS_PAGES")
    if sys.platform == "darwin":
        try:
            report = subprocess.check_output(
                ["memory_pressure", "-Q"], text=True,
                stderr=subprocess.DEVNULL)
            percent = re.search(r"free percentage:\s*(\d+)%", report)
            if percent:
                return total * int(percent.group(1)) // 100
        except (OSError, subprocess.SubprocessError):
            pass
    # ponytail: half of physical RAM is the portable fallback; add a platform
    # query if a supported build host proves this too conservative.
    return total // 2


def _available_device_memory(device):
    """Return bytes available after model weights have reached the device."""
    if device == "cuda":
        torch.cuda.synchronize()
        torch.cuda.empty_cache()
        return int(torch.cuda.mem_get_info()[0])
    if device == "mps":
        torch.mps.synchronize()
        torch.mps.empty_cache()
        device_free = max(0, int(torch.mps.recommended_max_memory()
                                 - torch.mps.driver_allocated_memory()))
        # MPS uses unified memory, so applications outside this process matter.
        return min(device_free, _cpu_available_memory())
    return _cpu_available_memory()


def _tmvec_working_bytes(batch, length, body, head):
    """Conservative peak-memory estimate for one TM-Vec forward pass."""
    batch, length = int(batch), int(length)
    cfg = body.config
    inner = cfg.num_heads * cfg.d_kv

    # T5 keeps one relative-position matrix and materializes attention scores
    # plus softmax weights. Linear terms cover q/k/v, feed-forward states, and
    # hidden outputs. Layers run serially under no_grad(), so layer count does
    # not multiply peak activation memory.
    body_quadratic = (2 * batch + 1) * cfg.num_heads * length * length
    body_linear = batch * length * (
        4 * inner + 2 * cfg.d_ff + 4 * cfg.d_model)

    head_quadratic = 2 * batch * head.nhead * length * length
    head_linear = batch * length * (
        4 * head.d_model + 2 * head.dim_feedforward)
    bridge = 2 * batch * length * head.d_model

    body_bytes = next(body.parameters()).element_size()
    head_bytes = next(head.parameters()).element_size()
    estimated = (body_bytes * (body_quadratic + body_linear)
                 + head_bytes * (head_quadratic + head_linear + bridge))
    return int(1.25 * estimated + 256 * 1024**2)


def _batches(lengths, memory_budget=None, memory_estimate=None):
    """Yield (start, stop) index pairs partitioning the sequences in order.

    Padding stretches every sequence in a batch to the longest one. Small
    internal count/token caps prevent oversized generic-model batches; when
    supplied, `memory_estimate` adds the device-specific quadratic limit. These
    caps never reject a sequence that can run alone.
    """
    if not lengths:
        return
    start, longest = 0, 0
    for i, n in enumerate(lengths):
        widest = max(longest, n)
        count = i - start + 1
        over_memory = (memory_estimate is not None
                       and memory_estimate(count, widest) > memory_budget)
        if i > start and (i - start + 1 > 16
                          or (i - start + 1) * widest > 16384
                          or over_memory):
            yield start, i
            start, longest = i, n
        else:
            longest = widest
    yield start, len(lengths)


@torch.no_grad()
def embed(seqs, spec_json, weights, device, memory_fraction):
    """Embed a list of amino-acid strings into a float64 matrix, one row each.

    `spec_json` is one entry of the model registry (inst/extdata/models/). Its
    `head` picks the path and says what `weights` is: unused for "generic", a
    HuggingFace head repo for "tmvec1", or a local Lightning checkpoint (already
    hash-verified on the R side) for "tmvec1-large" -- which, carrying no
    config.json of its own, takes its architecture from the registry's `config`.
    """
    spec = json.loads(spec_json)
    head, backbone = spec["head"], spec["backbone"]
    seqs = list(seqs)
    device = resolve_device(device)
    if head == "generic":
        tok, model = _load_generic(backbone, device)
        run = lambda batch: _run_generic(batch, device, tok, model)
        out_dim = model.config.hidden_size
        memory_budget = memory_estimate = None
    elif head in ("tmvec1", "tmvec1-large"):
        cfg_json = json.dumps(spec["config"]) if "config" in spec else ""
        top = _load_head(weights, device, cfg_json)
        _, body = _load_prott5(backbone, device)
        run = lambda batch: _run_tmvec(batch, device, backbone, top)
        out_dim = top.mlp.out_features
        available = _available_device_memory(device)
        memory_budget = int(float(memory_fraction) * available)
        memory_estimate = lambda count, width: _tmvec_working_bytes(
            count, width, body, top)
    else:
        raise ValueError("unknown head type '{}'".format(head))

    lengths = [len(s) for s in seqs]
    context = spec.get("context_max_length")
    if head == "generic" and context is None:
        limits = (getattr(model.config, "max_position_embeddings", None),
                  getattr(tok, "model_max_length", None))
        limits = [n for n in limits if isinstance(n, int) and n < 10**9]
        context = (min(limits) - tok.num_special_tokens_to_add()
                   if limits else None)
    eligible = [i for i, length in enumerate(lengths)
                if (context is None or length <= context)
                and (memory_estimate is None
                     or memory_estimate(1, length) <= memory_budget)]
    pieces, kept = [], []

    def run_indices(indices):
        try:
            pieces.append(run([seqs[i] for i in indices]).cpu().numpy())
            kept.extend(indices)
        except (RuntimeError, MemoryError) as error:
            oom = (isinstance(error, MemoryError)
                   or "out of memory" in str(error).lower()
                   or "allocate memory" in str(error).lower())
            if not oom:
                raise
            if device != "cpu":
                getattr(torch, device).empty_cache()
            if len(indices) == 1:
                return
            middle = len(indices) // 2
            run_indices(indices[:middle])
            run_indices(indices[middle:])

    eligible_lengths = [lengths[i] for i in eligible]
    for start, stop in _batches(eligible_lengths, memory_budget, memory_estimate):
        run_indices(eligible[start:stop])
        # torch's device allocators cache freed blocks; over a whole protein
        # universe that cache outgrows the device even though nothing is live.
        if device != "cpu":
            getattr(torch, device).empty_cache()

    # Raw TM-Vec norms vary with protein length, so only the normalised inner
    # product is the calibrated cosine the head predicts a TM-score from
    # (see docs/ground_metric_models.md).
    out = np.concatenate(pieces) if pieces else np.empty((0, out_dim))
    if len(out):
        out = out / np.linalg.norm(out, axis=1, keepdims=True)
    return {"embeddings": out, "kept": np.asarray(kept, dtype=np.int64)}


if __name__ == "__main__":
    # `python -m repdist_embed`: the only logic checkable without a model.
    assert _cgroup_remaining("max", "1") is None
    assert _cgroup_remaining(str(1 << 62), "1") is None
    assert _cgroup_remaining("200", "75") == 125
    assert list(_batches([10] * 17)) == [(0, 16), (16, 17)]
    assert list(_batches([9000, 9000])) == [(0, 1), (1, 2)]
    assert list(_batches([50000])) == [(0, 1)]
    assert list(_batches([])) == []
    square = lambda batch, width: batch * width * width
    assert list(_batches([10, 10, 10], 250, square)) == [
        (0, 2), (2, 3)]
    print("ok")
