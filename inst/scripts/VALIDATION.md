# Maintainer validation

Before submission, on Bioconductor devel:

* `R CMD build .`, then `R CMD check` on the tarball.
* `BiocCheck::BiocCheck(new-package = TRUE)` and
  `BiocCheck::BiocCheckGitClone()` on a clean checkout.

`StagedInstall: no` is required by basilisk, which cannot have staged
installation paths baked into its Python environments; R CMD check's NOTE about
it is expected.

The test suite covers the numerical kernel, dense/sparse agreement, identifier
alignment, argument rejection and the mocked R-to-Python dispatch. It does not
download models. For a real inference check, install the package and run
`Rscript inst/scripts/check_embedding.R` (ESM2-8m on CPU: unit norms, row
order, repeated sequences, batch-size invariance). Set `REPDIST_CHECK_MODEL` to
check another registered preset.

The 0.99.0 numerical fixes, model pins, added regression checks and data
provenance were AI-assisted; disclose this in the submission issue under
Bioconductor's AI and third-party code policy.
