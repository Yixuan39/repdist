# How much of a representation ground metric is homology, and how much is taxonomy?
#
# The claim a referee will attack is that the PLM metric is a taxonomic metric wearing
# a functional hat. That is answerable by one regression on the upper triangle of the
# U x U cost matrix, before any sample or abundance enters -- no group labels, no
# PERMANOVA, no data beyond the protein universe itself.
#
# Unique R^2 is what matters. Family difference and taxonomic distance are correlated
# by biology (orthologs ARE homologs), so marginal R^2 overstates both.

upper <- function(M) M[upper.tri(M)]

r2 <- function(y, X) {
  f <- stats::lm.fit(cbind(1, X), y)
  1 - sum(f$residuals^2) / sum((y - mean(y))^2)
}

decompose <- function(C, predictors) {
  y <- upper(C)
  P <- lapply(predictors, upper)
  full <- r2(y, do.call(cbind, P))
  unique_share <- vapply(names(P), function(k) {
    rest <- P[setdiff(names(P), k)]
    full - if (length(rest)) r2(y, do.call(cbind, rest)) else 0
  }, numeric(1))
  list(full = full, unique = unique_share)
}

test_that("a known family > taxonomy > length mixture is recovered", {
  # Embeddings with a KNOWN mixture. A decomposition that cannot recover the injected
  # ordering is not fit to answer the question on real embeddings. The length term is a
  # fixed direction scaled by log length -- the mechanism method_notes.md §2 warns about.
  set.seed(0)
  U <- 300L; d <- 64L; n_fam <- 12L; n_tax <- 8L
  fam <- sample.int(n_fam, U, replace = TRUE)
  tax <- sample.int(n_tax, U, replace = TRUE)
  len <- runif(U, 80, 900)

  Z <- 1.00 * matrix(rnorm(n_fam * d), n_fam, d)[fam, ] +
       0.45 * matrix(rnorm(n_tax * d), n_tax, d)[tax, ] +
       0.20 * outer(log(len), rnorm(d)) +
       0.30 * matrix(rnorm(U * d), U, d)
  Z <- Z / sqrt(rowSums(Z^2))
  C <- as.matrix(dist(Z))

  preds <- list(
    # 1[family differs] is the ground metric an annotation-based method uses, so its
    # share is "how much of this metric a Pfam/eggNOG table already has".
    family_differs = outer(fam, fam, "!=") * 1,
    taxon_differs  = outer(tax, tax, "!=") * 1,
    log_length_gap = abs(outer(log(len), log(len), "-"))
  )

  d3 <- decompose(C, preds)
  expect_gt(d3$full, 0.30)
  expect_gt(d3$unique[["family_differs"]], d3$unique[["taxon_differs"]])
  expect_gt(d3$unique[["taxon_differs"]], d3$unique[["log_length_gap"]])
  # The point of the exercise: taxonomy is separately detectable, so it is a
  # measurable competitor rather than an unmeasurable confound.
  expect_gt(d3$unique[["taxon_differs"]], 0.005)
})
