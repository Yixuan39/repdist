# MMD is a distance only when the kernel is positive semi-definite. A Gaussian
# RBF guarantees that on a true metric (Schoenberg), and `1 - cosine` is not
# one -- it violates the triangle inequality. So sample_repdist() defaults to
# "euclidean" while seq_repdist()/repdist_bin() default to "cosine", which is
# the calibrated TM-score scale cluster_threshold is expressed in.
#
# These embeddings are the shape that broke the old cosine default: four
# structural families on disjoint dimension blocks, and samples that share no
# accession at all -- the scenario the method exists for.

sim_universe <- function(seed, n_fam = 4L, per_fam = 12L, d = 32L) {
  set.seed(seed)
  centres <- matrix(0.3, n_fam, d)
  for (f in seq_len(n_fam)) centres[f, ((f - 1L) * 8L + 1L):(f * 8L)] <- 1.5
  family <- rep(seq_len(n_fam), each = per_fam)
  Z <- centres[family, ] +
    matrix(abs(rnorm(n_fam * per_fam * d, 0, 0.2)), n_fam * per_fam, d)
  Z <- Z / sqrt(rowSums(Z^2))            # unit-norm, as embed_proteins() returns
  rownames(Z) <- sprintf("prot%02d", seq_len(nrow(Z)))

  counts <- matrix(0L, 8L, nrow(Z),
                   dimnames = list(sprintf("s%d", 1:8), rownames(Z)))
  for (i in seq_len(8L)) {
    fams <- if (i <= 4L) c(1L, 2L) else c(3L, 4L)
    k <- if (i <= 4L) i else i - 4L
    idx <- unlist(lapply(fams,
      function(f) (f - 1L) * per_fam + ((k - 1L) * 3L + 1:3)))
    counts[i, idx] <- rpois(length(idx), 30) + 10L
  }
  list(Z = Z, counts = counts)
}

test_that("the default ground metric yields a PSD kernel on disjoint samples", {
  # The old "cosine" default failed roughly 9 runs in 10 here.
  for (seed in 1:12) {
    x <- sim_universe(seed)
    expect_silent(D <- sample_repdist(x$counts, x$Z))
    expect_s3_class(D, "dist")
    expect_true(all(is.finite(D)))
    expect_true(all(D >= 0))
  }
})

test_that("sample_repdist defaults to euclidean, seq_repdist to cosine", {
  x <- sim_universe(1)
  expect_equal(sample_repdist(x$counts, x$Z),
               sample_repdist(x$counts, x$Z, distance = "euclidean"))
  expect_equal(seq_repdist(x$Z), seq_repdist(x$Z, distance = "cosine"))
})

test_that("an RBF kernel on a euclidean ground metric is PSD", {
  Z <- sim_universe(1)$Z
  K <- rbf_kernel(as.matrix(seq_repdist(Z, "euclidean")))
  ev <- eigen(K, symmetric = TRUE, only.values = TRUE)$values
  expect_gt(min(ev) / max(ev), -1e-8)
})

test_that("cosine and euclidean rank protein pairs identically when unit-norm", {
  # Why the split costs nothing: ||a - b|| = sqrt(2 * (1 - cos)) for unit-norm
  # rows, so the two ground metrics are monotone transforms of one another.
  Z <- sim_universe(1)$Z
  dc <- seq_repdist(Z, "cosine")
  de <- seq_repdist(Z, "euclidean")
  expect_equal(as.vector(de), sqrt(2 * as.vector(dc)), tolerance = 1e-8)
  expect_equal(cor(as.vector(dc), as.vector(de), method = "spearman"), 1,
               tolerance = 1e-8)
})

test_that("mmd_matrix rejects an indefinite kernel rather than returning junk", {
  K <- matrix(c(1, 0.9, 0.1, 0.9, 1, 0.9, 0.1, 0.9, 1), 3, 3)
  expect_lt(min(eigen(K, symmetric = TRUE, only.values = TRUE)$values), 0)
  # An indefinite kernel only misbehaves along its negative directions, which is
  # why this failure is data-dependent and can hide in testing. p - q here is
  # proportional to (1, -2, 1), one such direction, so MMD^2 goes negative.
  P <- rbind(s1 = c(0.5, 0, 0.5), s2 = c(0, 1, 0))
  expect_lt(drop(crossprod(P[1, ] - P[2, ], K %*% (P[1, ] - P[2, ]))), 0)
  expect_error(mmd_matrix(P, K), "positive-semidefinite")
})
