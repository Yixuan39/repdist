# repdist 0.99.0

- `sample_repdist()` now only uses a Euclidean ground metric, while
  `seq_repdist()` and `repdist_bin()` keep `"cosine"`. MMD is a distance only
  under a positive-semi-definite kernel, and a Gaussian RBF guarantees that
  only on a true metric; `1 - cosine` is not one, so the old default produced
  an indefinite kernel and `sample_repdist()` aborted on most inputs where
  samples shared few accessions. Nothing is lost: `embed_proteins()` returns
  unit-norm rows, for which `||a - b|| = sqrt(2 * (1 - cos))`, so the two
  metrics rank protein pairs identically. The cosine option was removed from
  `sample_repdist()` and remains available only where the TM-score scale is
  valid, such as `repdist_bin()`'s `cluster_threshold`.
- `sample_repdist()` rejects a precomputed non-Euclidean `dist` rather than
  letting an indefinite kernel pass for a distance.
- `mmd_matrix()` scales its negative-MMD tolerance to the kernel rather than to
  the MMD matrix itself. The old tolerance collapsed onto rounding noise when
  every sample had the same composition, rejecting a valid kernel.
- `sample_repdist()` now computes RBF-MMD only. OT, transport matrices, core
  controls, and parallel imports were removed.
- `embed_proteins()` embeds complete proteins without truncation. Its single
  `memory_fraction` control sizes work to available Linux, Windows, macOS,
  CUDA, or MPS memory and warns when a sequence cannot fit or exceeds a true
  model context limit.
- `read_fasta()` now performs the shared Biostrings-based protein QC used by
  sequence-consuming functions. It normalizes terminal stops and rejects
  missing, empty, duplicate-named, invalid, or model-incompatible sequences.
- `read_embeddings()` imports pre-computed NumPy `.npz` embeddings using either
  `ids` or the official TM-Vec `headers` identifier convention.
- `length_abundance_profile()` moved to `R/QC.R`.
- `repdist_bin()` now uses complete linkage so every within-bin predicted
  TM-score meets the requested threshold. It returns the original similarity
  tree, structural bins, the tree pruned to medoid representatives, and a
  bin-named `AAStringSet` of representative sequences when sequences are given.
- `repdist_tree()` was removed; tree construction is part of `repdist_bin()`.
