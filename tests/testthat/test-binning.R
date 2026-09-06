set.seed(0)

base <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Z <- do.call(rbind, lapply(1:3, function(k) {
  matrix(base[k, ], 5, 3, byrow = TRUE) + matrix(rnorm(15, 0, 0.03), 5, 3)
}))
rownames(Z) <- paste0("p", 1:15)
family <- rep(1:3, each = 5)
D <- repdist_matrix(Z)

test_that("bin_proteins returns clusters, similarity, representatives, tree", {
  out <- bin_proteins(D, min_sim = 0.5)
  expect_named(out, c("clusters", "similarity", "representatives", "tree",
                      "min_sim", "method", "parameters", "noise"))
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
  expect_error(bin_proteins(repdist_matrix(Z, "euclidean")), "euclidean")
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
  expect_equal(bin_proteins(repdist_matrix(Z), min_sim = 0.6)$min_sim, 0.6)
})

test_that("density methods preserve noise as separate bins and support QC", {
  skip_if_not_installed("dbscan")
  W <- rbind(Z, noise1 = c(-1, 0, 0), noise2 = c(0, -1, 0))
  gm <- repdist_matrix(W)
  for (method in c("dbscan", "hdbscan")) {
    out <- bin_proteins(gm, method = method, min_pts = 3)
    expect_null(out$tree)
    expect_identical(out$method, method)
    expect_setequal(unlist(out$clusters), rownames(W))
    expect_length(unlist(out$clusters), nrow(W))
    expect_setequal(out$noise, c("noise1", "noise2"))
    expect_equal(unname(lengths(out$clusters)), c(5L, 5L, 5L, 1L, 1L))
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      expect_s3_class(plot_bin_similarity(out, "bin1", gm), "ggplot")
      expect_error(plot_bin_profile(out), "no tree")
    }
    if (method == "hdbscan") expect_null(out$min_sim)
  }
  all_noise <- bin_proteins(gm, method = "dbscan", min_pts = 100)
  expect_setequal(all_noise$noise, rownames(W))
  expect_true(all(lengths(all_noise$clusters) == 1L))
  all_noise <- bin_proteins(gm, method = "hdbscan", min_pts = 100)
  expect_setequal(all_noise$noise, rownames(W))
  expect_true(all(lengths(all_noise$clusters) == 1L))
  expect_error(bin_proteins(gm, method = "dbscan", min_pts = 2.5))
  expect_error(bin_proteins(gm, method = "dbscan", border_points = NA))
  expect_warning(bin_proteins(gm, method = "hdbscan", min_sim = 0.5),
                  "ignored")
})

test_that("DBSCAN radius is not interpreted as a complete-linkage floor", {
  skip_if_not_installed("dbscan")
  # A chain with a weak endpoint pair: DBSCAN connects it; hclust cannot.
  gm <- as.dist(matrix(c(0, .2, .4, .2, 0, .2, .4, .2, 0), 3,
                        dimnames = list(letters[1:3], letters[1:3])))
  out <- bin_proteins(gm, method = "dbscan", min_pts = 2)
  expect_length(out$clusters, 1L)
  expect_length(bin_proteins(gm)$clusters, 2L)
  core <- bin_proteins(gm, method = "dbscan", min_pts = 3)
  border <- bin_proteins(gm, method = "dbscan", min_pts = 3,
                         border_points = TRUE)
  expect_setequal(core$noise, c("a", "c"))
  expect_length(border$clusters, 1L)
})

test_that("MCL separates weakly connected groups and preserves isolates and IDs", {
  skip_if(!nzchar(Sys.which("mcl")), "MCL executable unavailable")
  ids <- c(paste0("protein ", 1:6), "isolate\twith whitespace")
  S <- diag(7)
  S[1:3, 1:3] <- S[4:6, 4:6] <- .95
  S[3, 4] <- S[4, 3] <- .51
  diag(S) <- 1
  dimnames(S) <- list(ids, ids)
  gm <- as.dist(1 - S)
  out <- bin_proteins(gm, min_sim = 0, method = "mcl")
  expect_null(out$tree)
  expect_identical(out$method, "mcl")
  expect_length(out$noise, 0L)
  expect_equal(out$clusters, list(bin1 = ids[1:3], bin2 = ids[4:6], bin3 = ids[7]))
  expect_true(all(out$representatives %in% ids))
  expect_equal(unname(out$similarity), unname(S[out$representatives,
                                              out$representatives]))
  if (requireNamespace("ggplot2", quietly = TRUE))
    expect_s3_class(plot_bin_similarity(out, "bin1", gm), "ggplot")
  empty <- bin_proteins(gm, min_sim = 1, method = "mcl")
  expect_true(all(lengths(empty$clusters) == 1L))
  expect_error(bin_proteins(gm, method = "mcl", inflation = 1))
  expect_error(bin_proteins(gm, method = "mcl", mcl_bin = "repdist-no-mcl"),
                "executable not found")
})

test_that("invalid distances fail before any clustering backend is called", {
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

test_that("failed or incomplete MCL output cannot silently lose proteins", {
  skip_on_os("windows")
  executable <- tempfile("mcl stub ")
  writeLines(c("#!/bin/sh", "exit 7"), executable)
  Sys.chmod(executable, "0755")
  expect_error(bin_proteins(D, method = "mcl", mcl_bin = executable), "MCL failed")
  writeLines(c("#!/bin/sh", 'while [ "$1" != "-o" ]; do shift; done',
               'echo 1 > "$2"'), executable)
  expect_error(bin_proteins(D, method = "mcl", mcl_bin = executable),
                "exactly once")
  unlink(executable)
})
