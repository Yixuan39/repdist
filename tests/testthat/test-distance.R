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

test_that("sample_repdist accepts embeddings or a precomputed ground distance", {
  expect_equal(
    sample_repdist(counts, Z),
    sample_repdist(counts, seq_repdist(Z)))
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

  for (distance in c("cosine", "euclidean")) {
    expect_gt(as.numeric(sample_repdist(cc, Zb, distance)), 0.1)
    expect_equal(
      as.numeric(sample_repdist(cc, Zb, distance, weighted = FALSE)), 0)
  }
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
  expect_gt(as.numeric(sample_repdist(cc, Zb, distance = "euclidean")), 0)
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
  G <- seq_repdist(Z)
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

  C <- as.matrix(seq_repdist(Z))
  C <- rbind(cbind(C, junk = 1), junk = c(rep(1, ncol(C)), 0))
  expect_equal(
    sample_repdist(cc, stats::as.dist(C)),
    sample_repdist(counts, seq_repdist(Z)))
})

test_that("one retained protein gives zero sample distance", {
  cc <- rbind(a = 10, b = 20)
  colnames(cc) <- "p1"
  expect_equal(as.numeric(sample_repdist(cc, Z[1, , drop = FALSE])), 0)
})

test_that("rbf_kernel handles a zero bandwidth without NaN", {
  C <- matrix(0, 3, 3)
  expect_equal(rbf_kernel(C), matrix(1, 3, 3))

  C2 <- rbind(c(0, 0, 2), c(0, 0, 2), c(2, 2, 0))
  expect_equal(
    rbf_kernel(C2, sigma = 0),
    matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 1), 3, 3))
  expect_error(rbf_kernel(C2, sigma = -1), "non-negative")
})

test_that("mmd_matrix matches the closed form on a linear kernel", {
  K <- Z %*% t(Z)
  P <- counts / rowSums(counts)
  slow <- mmd_matrix(P, K)
  fast <- as.matrix(stats::dist(P %*% Z))
  dimnames(slow) <- dimnames(fast)
  expect_equal(fast, slow)
})

test_that("mmd_matrix rejects a kernel that is not positive semidefinite", {
  K <- matrix(c(1, 2, 2, 1), 2, 2)
  P <- rbind(a = c(1, 0), b = c(0, 1))
  expect_error(mmd_matrix(P, K), "positive-semidefinite")
})

test_that("seq_repdist computes cosine and Euclidean distances", {
  Zn <- Z / sqrt(rowSums(Z^2))
  expect_equal(as.vector(seq_repdist(Z, "cosine")),
               as.vector(stats::as.dist(1 - tcrossprod(Zn))))
  expect_equal(as.vector(seq_repdist(Z, "euclidean")),
               as.vector(stats::dist(Z)))

  zero <- rbind(Z[1:3, ], dead = 0)
  expect_error(seq_repdist(zero, "cosine"), "zero-length")
})

test_that("cosine blocks produce the standard condensed ordering", {
  set.seed(1)
  x <- matrix(rnorm(7 * 4), 7, 4,
              dimnames = list(paste0("x", 1:7), NULL))
  x <- x / sqrt(rowSums(x^2))
  observed <- .repdist_cosine_dist(x, block_bytes = 2 * 8 * nrow(x))
  expected <- stats::as.dist(1 - tcrossprod(x))

  expect_equal(as.matrix(observed), as.matrix(expected))
  expect_identical(attr(observed, "Labels"), rownames(x))
})

test_that("condensed distances can be subset and reordered directly", {
  full <- seq_repdist(Z)
  keep <- c("p4", "p1", "p3")
  observed <- .repdist_subset_dist(full, keep)

  expect_equal(as.matrix(observed),
               as.matrix(full)[keep, keep], ignore_attr = TRUE)
  expect_identical(attr(observed, "Labels"), keep)
  expect_error(.repdist_subset_dist(full, "missing"), "missing requested")
})
