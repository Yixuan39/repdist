# Writes inst/extdata/simulate_ec_similarity.rds: protein_similarity() of
# every data(simulate_ec) protein, for the clustering vignette. All-pairs
# alignment of the 1,160 proteins takes about 20 minutes on one core (a few
# minutes on several); shipping the result keeps the vignette fast. The gate
# is deterministic regardless of n_cores (see test-binning.R), so this can be
# regenerated with any core count and get the same matrix. See
# DATA_PROVENANCE.md.
#
# Usage: Rscript inst/scripts/make_simulate_ec_similarity.R [n_cores]

library(repdist)
data(simulate_ec, package = "repdist")
n_cores <- as.integer(commandArgs(TRUE)[1])
if (is.na(n_cores)) n_cores <- 1L

S <- protein_similarity(simulate_ec$proteins, n_cores = n_cores)
saveRDS(Matrix::Matrix(S, sparse = TRUE),
    "inst/extdata/simulate_ec_similarity.rds", compress = "xz")
