set.seed(0)

base <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Z <- do.call(rbind, lapply(1:3, function(k) {
  matrix(base[k, ], 5, 3, byrow = TRUE) + matrix(rnorm(15, 0, 0.03), 5, 3)
}))
Z <- Z / sqrt(rowSums(Z^2))
rownames(Z) <- paste0("p", 1:15)
family <- rep(1:3, each = 5)
D <- repdist_matrix(Z)

test_that("bin_proteins returns clusters, similarity, representatives, tree", {
  out <- bin_proteins(D, min_sim = 0.5)
  expect_named(out, c("clusters", "similarity", "representatives", "tree",
                      "min_sim"))
  expect_setequal(unlist(out$clusters), rownames(Z))
  expect_named(out$representatives, names(out$clusters))
  expect_s3_class(out$tree, "hclust")
  expect_equal(dim(out$similarity), c(3, 3))
  expect_equal(dimnames(out$similarity),
               list(names(out$clusters), names(out$clusters)))

  # bins are named largest first
  expect_false(is.unsorted(rev(lengths(out$clusters))))
})

test_that("the partition is exactly a complete-linkage cut of the distance", {
  out <- bin_proteins(D, min_sim = 0.5)
  # compared by co-membership: bin labels are ranked by size while cutree
  # numbers by order of appearance -- same partition, different names
  membership <- setNames(rep(names(out$clusters), lengths(out$clusters)),
                         unlist(out$clusters))[rownames(Z)]
  expected <- stats::cutree(
    fastcluster::hclust(D, method = "complete"), h = 1 - 0.5)[rownames(Z)]
  expect_equal(outer(membership, membership, "=="),
               outer(expected, expected, "=="), ignore_attr = TRUE)
})

test_that("complete linkage enforces the all-pairs similarity threshold", {
  out <- bin_proteins(D, min_sim = 0.5)
  expect_equal(lengths(out$clusters), c(bin1 = 5L, bin2 = 5L, bin3 = 5L))
  membership <- setNames(rep(names(out$clusters), lengths(out$clusters)),
                         unlist(out$clusters))
  for (fam in 1:3)
    expect_length(unique(membership[rownames(Z)[family == fam]]), 1L)

  similarity <- 1 - as.matrix(D)
  minimum <- vapply(out$clusters, function(members) {
    if (length(members) < 2L) return(1)
    min(similarity[members, members][upper.tri(diag(length(members)))])
  }, numeric(1))
  expect_true(all(minimum >= 0.5))
})

test_that("a lower threshold merges families and a higher one splits them", {
  # two families 45 degrees apart: similarity 0.71, so 0.5 merges, 0.9 splits
  W <- rbind(a1 = c(1, 0), a2 = c(1, 0), b1 = c(1, 1), b2 = c(1, 1))
  W <- W / sqrt(rowSums(W^2))
  expect_length(bin_proteins(repdist_matrix(W), min_sim = 0.5)$clusters, 1)
  expect_length(bin_proteins(repdist_matrix(W), min_sim = 0.9)$clusters, 2)

  expect_gt(length(bin_proteins(D, min_sim = 0.999)$clusters), 3)
})

test_that("the returned tree reproduces any other threshold via cutree", {
  out <- bin_proteins(D, min_sim = 0.5)
  again <- bin_proteins(D, min_sim = 0.999)
  expect_equal(length(again$clusters),
               length(unique(stats::cutree(out$tree, h = 1 - 0.999))))
})

test_that("min_sim outside (0, 1] is rejected rather than silently cut", {
  expect_error(bin_proteins(D, min_sim = -0.5))
  expect_error(bin_proteins(D, min_sim = 1.5))
  expect_error(bin_proteins(D, min_sim = 0))
  expect_error(bin_proteins(D, min_sim = c(0.5, 0.7)))
})

test_that("a euclidean or non-finite distance is rejected", {
  expect_error(bin_proteins(stats::dist(Z)), "euclidean")
  bad <- D
  bad[1] <- NA
  expect_error(bin_proteins(bad), "non-finite")
})

test_that("a square distance matrix is accepted", {
  expect_equal(bin_proteins(as.matrix(D), min_sim = 0.5)$clusters,
               bin_proteins(D, min_sim = 0.5)$clusters)
})

test_that("representatives are within-bin medoids and similarity matches", {
  out <- bin_proteins(D, min_sim = 0.5)
  C <- as.matrix(D)
  for (bin in names(out$clusters)) {
    members <- out$clusters[[bin]]
    expect_equal(out$representatives[[bin]],
                 members[which.min(colSums(C[members, members, drop = FALSE]))])
  }
  reps <- out$representatives
  expect_equal(out$similarity, 1 - C[reps, reps], ignore_attr = TRUE)
})

test_that("the condensed distance is never expanded to a square matrix", {
  gm <- repdist_matrix(Z)
  class(gm) <- c("no_square_dist", class(gm))
  as.matrix.no_square_dist <- function(...) stop("expanded to a square matrix")

  expect_no_error(bin_proteins(gm))
})

test_that("bin_proteins reports the threshold it cut at", {
  Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
  Z <- Z / sqrt(rowSums(Z^2))
  expect_equal(bin_proteins(repdist_matrix(Z), min_sim = 0.6)$min_sim, 0.6)
})

test_that("invalid distances fail before clustering", {
  bad <- as.matrix(D)
  bad[1, 2] <- 1.5
  expect_error(bin_proteins(bad), "symmetric")
  expect_error(bin_proteins(bad[, -1]), "symmetric")
  bad <- D
  bad[1] <- -.1
  expect_error(bin_proteins(bad), "\\[0, 2\\]")
  attr(bad, "Labels")[1] <- NA
  expect_error(bin_proteins(bad))
})
