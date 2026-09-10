# repdist 0.99.0

Initial Bioconductor submission.

- `bin_proteins()` cuts a complete-linkage tree and nothing else. The
  optional MCL, DBSCAN*/DBSCAN and HDBSCAN methods have been removed: on the
  EC benchmark none of them, nor an emergent SOM, improved annotation purity
  at a matched bin count, and each gave up the all-pairs similarity floor and
  the recuttable tree. The `method`, `min_pts`, `border_points`, `inflation`
  and `mcl_bin` arguments and the `$method`, `$parameters` and `$noise` return
  fields are gone; `$tree` is always an `hclust`.
- `embed_proteins()` embeds every sequence whole -- never truncated, never
  silently dropped -- with a registered TM-Vec head or any HuggingFace
  encoder, through a pinned basilisk Python environment. Batches are capped
  at `batch_size` sequences and split automatically on out-of-memory; a
  sequence that cannot be embedded at all (over the model's architectural
  context, or too large for the device even alone) is an error naming the
  offender.
- `read_fasta()` reads amino-acid FASTA with the shared protein QC:
  terminal stops normalized; missing, empty, duplicate-named, or
  model-incompatible sequences rejected.
- `repdist_matrix()` builds the protein-by-protein `dist` everything else
  consumes, condensed so a square matrix is never materialized. `"cosine"`
  is `1 - predicted TM-score` for TM-Vec embeddings, computed as
  `||a - b||^2 / 2` on unit-normed rows through `stats::dist()`, so no extra
  dependency is loaded; `"euclidean"` is the true metric a kernel needs.
- `sample_repdist()` computes RBF-MMD between samples under a
  median-heuristic bandwidth, on a euclidean ground metric only: MMD is a
  distance only under a positive-semi-definite kernel, and a Gaussian RBF on
  `1 - cosine` can be indefinite. On unit-norm embeddings nothing is lost --
  `||a - b|| = sqrt(2 * (1 - cos))`, so the two scales rank pairs
  identically. The kernel is accumulated directly into the sample Gram matrix
  from two-dimensional protein blocks, and sparse count matrices stay sparse,
  so no allocation grows quadratically with the protein count.
  The bandwidth used is returned as `attr(, "sigma")`; pass it back to hold it
  fixed across studies. With the median default the kernel is nearly linear in TM-score,
  so the distance is close to the euclidean distance between mean
  embeddings -- on the shipped data, also where separation peaks.
- `bin_proteins()` cuts a complete-linkage tree over a `repdist_matrix()`
  distance at `1 - min_sim`, so every within-bin pair meets the requested
  similarity. It returns the cluster list, the similarity matrix between bin
  representatives (medoids), the representative accessions, and the full
  `hclust` tree for recutting at other thresholds.
- Three QC plots: `plot_length_profile()` shows the fraction of the catalog
  and of total abundance beyond each length, for checking a model's training
  range against the data; `plot_similarity_profile()` histograms each
  protein's nearest-neighbour similarity, exposing near-duplicates,
  saturation, and would-be singleton bins; `plot_bin_profile()` sweeps the
  binning threshold over the already-computed tree, showing bin count and
  largest-bin share at every `min_sim`.
- Vignette `repdist` walks the whole workflow on real data (Blakeley-Ruiz
  et al., ISME J 2025 mouse diet metaproteome): MMD beta diversity against a
  Bray-Curtis baseline with a vegan PERMANOVA, structural bins at a predicted
  TM-score of 0.7, bin-level differential abundance with DESeq2, and the
  functional annotations of the most responsive bin. The example data ships
  as `inst/extdata/diet_metaproteome.rds`.
- Vignette `simulation` builds the case the workflow exists for: two groups
  whose difference is a functional shift spread across non-overlapping
  homologs (10 CATH FunFams, shipped as `inst/extdata/funfam_universe.rds`).
  RBF-MMD, weighted UniFrac on the structure tree, and Bray-Curtis on TM 0.7
  bins all detect the groups; protein-level Bray-Curtis does not.
