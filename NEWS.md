# repdist 0.99.0

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
