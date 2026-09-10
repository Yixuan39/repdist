set.seed(0)

Z <- rbind(
  matrix(rnorm(10 * 4, 0, 0.05), 10, 4),
  matrix(rnorm(10 * 4, 5, 0.05), 10, 4),
  matrix(rnorm(10 * 4, 10, 0.05), 10, 4)
)
rownames(Z) <- paste0("p", 1:30)
counts <- matrix(
  rpois(6 * 30, 20) + 1L, 6, 30,
  dimnames = list(paste0("s", 1:6), rownames(Z)))

test_that("sample_repdist takes the embedding matrix, not a precomputed dist", {
  expect_error(sample_repdist(counts, stats::dist(Z)),
               "embedding matrix")
  expect_error(sample_repdist(counts, repdist_matrix(Z / sqrt(rowSums(Z^2)))),
               "embedding matrix")
})

test_that("the blocked kernel matches the dense one", {
  G <- as.matrix(stats::dist(Z))
  sigma <- stats::median(G[upper.tri(G)])
  P <- counts / rowSums(counts)
  dense <- repdist:::mmd_matrix(tcrossprod(
    P %*% repdist:::rbf_kernel(G, sigma), P))
  D <- sample_repdist(counts, Z)
  expect_equal(attr(D, "sigma"), sigma)
  expect_equal(as.numeric(D), as.numeric(stats::as.dist(dense)))
  # more proteins than one block, so the loop and the triangle offsets run
  set.seed(4)
  big <- matrix(rnorm(2500 * 3), 2500, 3, dimnames = list(paste0("q", 1:2500), NULL))
  cb <- matrix(rpois(4 * 2500, 1), 4, 2500, dimnames = list(letters[1:4], rownames(big)))
  Gb <- as.matrix(stats::dist(big[colnames(cb)[colSums(cb) > 0], ]))
  Db <- sample_repdist(cb, big)
  # Past 2000 proteins the bandwidth is the median over a subsample, so it
  # tracks the all-pairs median rather than matching it exactly.
  expect_equal(attr(Db, "sigma"), stats::median(Gb[upper.tri(Gb)]),
               tolerance = 0.01)
  Pb <- cb[, colSums(cb) > 0] / rowSums(cb)
  expect_equal(as.numeric(Db), as.numeric(stats::as.dist(
    repdist:::mmd_matrix(tcrossprod(
      Pb %*% repdist:::rbf_kernel(Gb, attr(Db, "sigma")), Pb)))))
})

test_that("sample_repdist preserves sparse counts", {
  sparse_counts <- counts
  sparse_counts[(row(sparse_counts) + col(sparse_counts)) %% 4 != 0] <- 0L
  sparse <- Matrix::Matrix(sparse_counts, sparse = TRUE)
  expect_equal(sample_repdist(sparse, Z), sample_repdist(sparse_counts, Z),
               tolerance = 1e-7)
  expect_equal(sample_repdist(sparse, Z, weighted = FALSE),
               sample_repdist(sparse_counts, Z, weighted = FALSE),
               tolerance = 1e-7)
})

test_that("kernel block and bandwidth sample sizes are configurable", {
  cc <- counts[1:3, 1:7]
  cc[, 7] <- 0L
  Zs <- Z[7:1, ]   # embedding order must be aligned after dropping column 7
  idx <- round(seq(1, 6, length.out = 3))
  expected_sigma <- stats::median(stats::dist(Z[idx, ]))
  for (weighted in c(TRUE, FALSE)) {
    expected <- sample_repdist(cc, Zs, weighted = weighted, sigma = expected_sigma)
    for (size in c(1, 4, 20)) {
      actual <- sample_repdist(cc, Zs, weighted = weighted,
                               block_size = size, sigma_max_proteins = 3)
      expect_equal(actual, expected, tolerance = 1e-7)
      expect_equal(sample_repdist(Matrix::Matrix(cc, sparse = TRUE), Zs,
                                   weighted = weighted, block_size = size,
                                   sigma_max_proteins = 3), actual,
                   tolerance = 1e-7)
    }
  }
  for (bad in list(0, -1, 1.5, NA, Inf, c(1, 2), "2", NULL)) {
    expect_error(sample_repdist(cc, Zs, block_size = bad), "block_size")
    expect_error(sample_repdist(cc, Zs, sigma_max_proteins = bad), "sigma_max_proteins")
  }
  expect_error(sample_repdist(cc, Zs, sigma_max_proteins = 1), "sigma_max_proteins")
})

test_that("sample_repdist returns a dist that adonis2 accepts", {
  skip_if_not_installed("vegan")
  D <- sample_repdist(counts, Z)
  expect_s3_class(D, "dist")
  expect_equal(attr(D, "Size"), nrow(counts))
  expect_equal(labels(D), rownames(counts))

  md <- data.frame(group = rep(c("A", "B"), each = 3))
  a <- vegan::adonis2(D ~ group, data = md, permutations = 99)
  expect_true(is.finite(a$F[1]))
})

test_that("RBF-MMD is symmetric, non-negative, and zero on the diagonal", {
  for (weighted in c(TRUE, FALSE)) {
    D <- as.matrix(sample_repdist(counts, Z, weighted = weighted))
    expect_true(isSymmetric(D))
    expect_true(all(D >= 0))
    expect_equal(diag(D), setNames(rep(0, nrow(D)), rownames(D)))
  }
})

test_that("identical composition gives zero distance", {
  cc <- rbind(a = c(10, 0, 0), b = c(10, 0, 0), d = c(0, 10, 0))
  colnames(cc) <- rownames(Z)[1:3]
  D <- as.matrix(sample_repdist(cc, Z[1:3, ]))
  expect_equal(D["a", "b"], 0)
  expect_gt(D["a", "d"], 0)
})

test_that("weighted = FALSE uses support rather than abundance", {
  Zb <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(-1, 0))
  cc <- rbind(a = c(98, 1, 1), b = c(1, 1, 98))
  colnames(cc) <- rownames(Zb)

  expect_gt(as.numeric(sample_repdist(cc, Zb)), 0.1)
  expect_equal(as.numeric(sample_repdist(cc, Zb, weighted = FALSE)), 0)
})

test_that("weighted = FALSE is exactly sign(abundance)", {
  cc <- rbind(a = c(90, 5, 5, 0), b = c(1, 1, 1, 97))
  colnames(cc) <- rownames(Z)[1:4]
  Zs <- Z[1:4, ]

  unweighted <- sample_repdist(cc, Zs, weighted = FALSE)
  expect_equal(unweighted, sample_repdist(sign(cc), Zs))
  expect_false(isTRUE(all.equal(
    as.numeric(unweighted), as.numeric(sample_repdist(cc, Zs)))))
})

test_that("RBF-MMD sees mass rearranged at a fixed weighted mean", {
  Zb <- rbind(p1 = c(-1, 0), p2 = c(1, 0), p3 = c(0, 0))
  cc <- rbind(a = c(50, 50, 0), b = c(0, 0, 100))
  colnames(cc) <- rownames(Zb)
  P <- cc / rowSums(cc)
  expect_equal(sqrt(sum((P[1, ] %*% Zb - P[2, ] %*% Zb)^2)), 0)
  expect_gt(as.numeric(sample_repdist(cc, Zb)), 0)
})

test_that("rarefaction shrinks depth-driven distance at fixed composition", {
  skip_if_not_installed("vegan")
  set.seed(1)
  p <- c(rep(0.06, 10), rep(0.03, 10), rep(0.01, 10))
  cc <- rbind(deep = rmultinom(1, 60000, p)[, 1],
              shallow = rmultinom(1, 1500, p)[, 1])
  colnames(cc) <- rownames(Z)

  raw <- as.numeric(sample_repdist(cc, Z))
  rarefied <- suppressWarnings(as.numeric(
    sample_repdist(vegan::rrarefy(cc, min(rowSums(cc))), Z)))
  expect_lt(rarefied, raw)
})

test_that("missing or unnamed representations fail loudly", {
  expect_error(sample_repdist(counts, Z[1:5, ]))
  G <- repdist_matrix(Z / sqrt(rowSums(Z^2)))
  attr(G, "Labels") <- NULL
  expect_error(sample_repdist(counts, G))
})

test_that("sample_repdist rejects malformed counts", {
  bad_names <- counts
  colnames(bad_names) <- NULL
  expect_error(sample_repdist(bad_names, Z))

  negative <- counts
  negative[1, 1] <- -1
  expect_error(sample_repdist(negative, Z))

  nonfinite <- counts
  nonfinite[1, 1] <- NA
  expect_error(sample_repdist(nonfinite, Z))

  empty_row <- counts
  empty_row[1, ] <- 0
  expect_error(sample_repdist(empty_row, Z))
  expect_error(sample_repdist(counts, Z, weighted = NA))
})

test_that("globally zero-count proteins are dropped before the ground distance", {
  extra <- rbind(Z, junk = c(NA, NA, NA, NA))
  cc <- cbind(counts, junk = 0L)
  expect_equal(sample_repdist(cc, extra), sample_repdist(counts, Z))
})

test_that("one retained protein gives zero sample distance", {
  cc <- rbind(a = 10, b = 20)
  colnames(cc) <- "p1"
  expect_equal(as.numeric(sample_repdist(cc, Z[1, , drop = FALSE])), 0)
  expect_equal(as.numeric(sample_repdist(cc, 1e10 * Z[1, , drop = FALSE])), 0)
})

test_that("repdist_matrix refuses more proteins than usedist can index", {
  # 46341 is the largest n with n * (n - 1) under .Machine$integer.max
  expect_identical(46341L * 46340L, 2147441940L)
  expect_true(is.na(suppressWarnings(46342L * 46341L)))
  expect_error(repdist_matrix(matrix(1, 46342L, 1L)), "46341")
})

test_that("a degenerate bandwidth is an error, not a NaN kernel", {
  same <- Z[rep(1, nrow(Z)), ]          # every protein identical
  rownames(same) <- rownames(Z)
  expect_error(sample_repdist(counts, same), "identical embeddings")
  near <- same
  near[, 1] <- near[, 1] + seq_len(nrow(near)) * 1e-10
  expect_error(sample_repdist(counts, near), "too small relative to")
  expect_error(sample_repdist(counts, Z, sigma = 0), "positive")
  expect_error(sample_repdist(counts, Z, sigma = -1), "positive")
})

test_that("mmd_matrix matches the closed form on a linear kernel", {
  K <- Z %*% t(Z)
  P <- counts / rowSums(counts)
  slow <- repdist:::mmd_matrix(tcrossprod(P %*% K, P))
  fast <- as.matrix(stats::dist(P %*% Z))
  dimnames(slow) <- dimnames(fast)
  expect_equal(fast, slow)
})

test_that("mmd_matrix rejects materially negative squared distances", {
  K <- matrix(c(1, 2, 2, 1), 2, 2)
  P <- rbind(a = c(1, 0), b = c(0, 1))
  expect_error(repdist:::mmd_matrix(tcrossprod(P %*% K, P)),
               "Negative squared MMD")
  # Roundoff is clipped, relative to the Gram matrix's scale.
  for (scale in c(1, 1e-6, 1e6)) {
    rounded <- matrix(c(1, 1 + 1e-10, 1 + 1e-10, 1), 2) * scale
    expect_equal(repdist:::mmd_matrix(rounded), matrix(0, 2, 2))
    invalid <- matrix(c(1, 1 + 1e-6, 1 + 1e-6, 1), 2) * scale
    expect_error(repdist:::mmd_matrix(invalid), "Negative squared MMD")
  }
})

test_that("repdist_matrix computes inner-product distances on unit-norm rows", {
  # More than one block, with a partial final block and named rows.
  set.seed(14)
  Z <- matrix(rnorm(1003 * 4), 1003, 4,
              dimnames = list(paste0("p", 1:1003), NULL))
  Z <- Z / sqrt(rowSums(Z^2))
  D <- repdist_matrix(Z)
  expect_equal(as.vector(D), as.vector(stats::dist(Z)^2 / 2))
  expect_equal(labels(D), rownames(Z))
  expect_identical(attr(D, "method"), "cosine")
  expect_s3_class(D, "dist")
  expect_length(repdist_matrix(Z[1, , drop = FALSE]), 0)
  expect_length(repdist_matrix(Z[FALSE, , drop = FALSE]), 0)
  expect_equal(as.vector(repdist_matrix(rbind(c(1, 0), c(1, 0),
                                             c(0, 1), c(-1, 0)))),
               c(0, 1, 2, 1, 2, 1))
  # Float32-sized norm error is accepted, with distances clipped to [0, 2].
  rounded <- rbind(c(1 + 1e-7, 0), c(1 + 1e-7, 0), c(-1 - 1e-7, 0))
  expect_equal(as.vector(repdist_matrix(rounded)), c(0, 2, 2))
  expect_error(repdist_matrix(2 * Z), "L2-normalised")
  expect_error(repdist_matrix(rbind(Z[1:3, ], dead = 0)), "L2-normalised")
  expect_error(repdist_matrix(matrix(NA_real_, 2, 2)), "is.finite")
  small <- Z[1:7, ]
  for (size in c(1, 3, 7, 20))
    expect_equal(repdist_matrix(small, block_size = size), repdist_matrix(small))
  for (size in list(0, -1, 1.5, NA, Inf, c(1, 2), "2", NULL))
    expect_error(repdist_matrix(small, block_size = size), "block_size")
})

test_that("on unit-norm embeddings the euclidean RBF is a kernel in cosine distance", {
  set.seed(13)
  Z <- matrix(rnorm(30 * 8), 30, 8, dimnames = list(paste0("p", 1:30), NULL))
  Z <- Z / sqrt(rowSums(Z^2))
  Ceuc <- as.matrix(stats::dist(Z))
  Ccos <- as.matrix(repdist_matrix(Z))
  expect_equal(Ceuc^2, 2 * Ccos, tolerance = 1e-6)
  sigma <- stats::median(Ceuc[upper.tri(Ceuc)])
  expect_equal(repdist:::rbf_kernel(Ceuc, sigma), exp(-Ccos / sigma^2), tolerance = 1e-6)
})

test_that("sample_repdist reports its bandwidth and holds it fixed when passed back", {
  D <- sample_repdist(counts, Z)
  G <- as.matrix(stats::dist(Z))
  expect_equal(attr(D, "sigma"), stats::median(G[upper.tri(G)]))
  expect_equal(sample_repdist(counts, Z, sigma = attr(D, "sigma")), D)

  # The default moves with the retained proteins; an explicit sigma does not.
  cc <- counts; cc[, 11:30] <- 0                       # six samples on one cluster
  cc7 <- rbind(cc, s7 = c(rep(0, 10), rep(5, 20)))     # a seventh brings in the rest
  d6 <- as.matrix(sample_repdist(cc, Z))
  expect_false(isTRUE(all.equal(d6, as.matrix(sample_repdist(cc7, Z))[1:6, 1:6])))
  fixed <- sample_repdist(cc7, Z, sigma = attr(sample_repdist(cc, Z), "sigma"))
  expect_equal(as.matrix(fixed)[1:6, 1:6], d6)
})

test_that("malformed sigma and duplicated protein names are rejected", {
  expect_error(sample_repdist(counts, Z, sigma = c(1, 2)), "single finite")
  expect_error(sample_repdist(counts, Z, sigma = NA), "single finite")
  expect_error(sample_repdist(counts, rbind(Z, Z[1, , drop = FALSE])))
})

test_that("the returned dist is a plain stats dist that base all.equal accepts", {
  # proxy's S3 methods on "dist" make all.equal() error on any dist object in
  # the session; repdist_matrix returns a plain dist without loading proxy.
  expect_true(isTRUE(all.equal(sample_repdist(counts, Z), sample_repdist(counts, Z))))
})
