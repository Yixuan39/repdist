# Build the simulation vignette's protein universe from the FunFam benchmark.
#
# data/funfam/beta_universe.rds carries 10 CATH FunFams x 150 domains with
# TM-Vec embeddings (see the simulation_study branch for how it was built).
# The vignette only needs the membership table and the embeddings, rounded to
# 4 decimal places -- that moves a cosine similarity by at most 2e-4, four
# orders of magnitude below any threshold used, and halves the shipped size.
#
# Run from the repository root:  Rscript inst/scripts/prepare_simulation_data.R
# Output: inst/extdata/funfam_universe.rds

u <- readRDS(file.path("data", "funfam", "beta_universe.rds"))

out <- list(
  members    = u$members[, c("protein", "family", "label")],
  embeddings = round(u$embeddings, 4),
  provenance = u$provenance[c("source", "accessed", "embedding_model")])

path <- file.path("inst", "extdata", "funfam_universe.rds")
saveRDS(out, path, compress = "xz")
message(sprintf("wrote %s (%.2f MB): %d proteins x %d dims", path,
                file.size(path) / 1048576, nrow(out$embeddings),
                ncol(out$embeddings)))
