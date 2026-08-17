set.seed(0)

base <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Z <- do.call(rbind, lapply(1:3, function(k) {
  matrix(base[k, ], 5, 3, byrow = TRUE) + matrix(rnorm(15, 0, 0.03), 5, 3)
}))
rownames(Z) <- paste0("p", 1:15)
family <- rep(1:3, each = 5)
counts <- matrix(
  rpois(4 * 15, 20) + 1L, 4, 15,
  dimnames = list(paste0("s", 1:4), rownames(Z)))

test_that("repdist_bin returns trees, bins, and representative sequences", {
  out <- repdist_bin(counts, Z)
  expect_named(out, c("tree", "bins", "binned_tree", "representative_sequences"),
               ignore.order = FALSE)
  expect_s3_class(out$tree, "phylo")
  expect_s3_class(out$binned_tree, "phylo")
  expect_setequal(out$tree$tip.label, colnames(counts))
  expect_setequal(out$binned_tree$tip.label, colnames(out$bins$counts))
  expect_null(out$representative_sequences)

  expected <- ape::as.phylo(fastcluster::hclust(seq_repdist(Z), method = "complete"))
  observed_C <- ape::cophenetic.phylo(out$tree)[rownames(Z), rownames(Z)]
  expected_C <- ape::cophenetic.phylo(expected)[rownames(Z), rownames(Z)]
  expect_equal(observed_C, expected_C)

  expected_binned <- ape::keep.tip(out$tree, out$bins$reps$protein)
  expected_binned$tip.label <- out$bins$reps$bin[
    match(expected_binned$tip.label, out$bins$reps$protein)]
  expect_equal(out$binned_tree, expected_binned)
})

test_that("complete-linkage bins enforce the requested all-pairs TM threshold", {
  bins <- repdist_bin(counts, Z, cluster_threshold = 0.5)$bins
  expect_equal(length(unique(bins$bin[family == 1])), 1)
  expect_equal(length(unique(bins$bin[family == 2])), 1)
  expect_equal(length(unique(bins$bin[family == 3])), 1)
  expect_equal(length(unique(bins$bin)), 3)

  similarity <- 1 - as.matrix(seq_repdist(Z))
  minimum <- vapply(split(names(bins$bin), bins$bin), function(members) {
    if (length(members) < 2L) return(1)
    min(similarity[members, members][upper.tri(similarity[members, members])])
  }, numeric(1))
  expect_true(all(minimum >= 0.5))
})

test_that("binning conserves total abundance per sample", {
  bins <- repdist_bin(counts, Z)$bins
  expect_equal(unname(rowSums(bins$counts)), unname(rowSums(counts)))
  expect_equal(sum(bins$counts), sum(counts))
})

test_that("a lower threshold merges families and a higher one splits them", {
  low <- repdist_bin(counts, Z, cluster_threshold = -0.5)$bins
  high <- repdist_bin(counts, Z, cluster_threshold = 0.999)$bins
  expect_equal(length(unique(low$bin)), 1)
  expect_gt(length(unique(high$bin)), 3)
})

test_that("missing embeddings fail loudly", {
  expect_error(repdist_bin(counts, Z[1:10, ]))
})

test_that("a single protein produces valid one-tip trees and one bin", {
  cc <- counts[, 1, drop = FALSE]
  out <- repdist_bin(cc, Z)
  expect_equal(out$bins$bin, c(p1 = "bin1"))
  expect_equal(unname(out$bins$counts[, 1]), unname(cc[, 1]))
  expect_equal(out$tree$tip.label, "p1")
  expect_equal(out$binned_tree$tip.label, "bin1")
})

test_that("malformed counts or embeddings are rejected", {
  negative <- counts
  negative[1, 1] <- -1
  expect_error(repdist_bin(negative, Z))

  nonfinite <- counts
  nonfinite[1, 1] <- NA
  expect_error(repdist_bin(nonfinite, Z))

  bad_Z <- Z
  bad_Z[1, 1] <- NA
  expect_error(repdist_bin(counts, bad_Z))

  unnamed <- counts
  colnames(unnamed) <- NULL
  expect_error(repdist_bin(unnamed, Z))
})

test_that("a threshold vector shares one full tree and returns named cuts", {
  out <- repdist_bin(counts, Z, cluster_threshold = c(0.5, 0.999))
  expect_named(out$bins, c("0.5", "0.999"))
  expect_named(out$binned_tree, c("0.5", "0.999"))
  expect_named(out$representative_sequences, c("0.5", "0.999"))
  expect_equal(out$bins[["0.5"]], repdist_bin(counts, Z, 0.5)$bins)
  expect_equal(out$bins[["0.999"]], repdist_bin(counts, Z, 0.999)$bins)
  expect_equal(out$tree, repdist_bin(counts, Z, 0.5)$tree)

  for (threshold in names(out$bins))
    expect_setequal(out$binned_tree[[threshold]]$tip.label,
                    colnames(out$bins[[threshold]]$counts))
})

test_that("bin_dist feeds sample_repdist at bin resolution", {
  bins <- repdist_bin(counts, Z, cluster_threshold = 0.5)$bins
  C <- as.matrix(seq_repdist(Z))

  expect_s3_class(bins$bin_dist, "dist")
  expect_setequal(labels(bins$bin_dist), colnames(bins$counts))
  expect_equal(as.matrix(bins$bin_dist),
               C[bins$reps$protein, bins$reps$protein], ignore_attr = TRUE)
  expect_equal(attr(sample_repdist(bins$counts, bins$bin_dist), "Size"),
               nrow(counts))
})

test_that("precomputed distances are not expanded for binning", {
  gm <- seq_repdist(Z)
  class(gm) <- c("no_square_dist", class(gm))
  as.matrix.no_square_dist <- function(...) stop("expanded to a square matrix")

  expect_no_error(repdist_bin(counts, gm))
})

test_that("representatives are within-bin medoids with named sequences", {
  seqs <- stats::setNames(rep(c("ACDE", "FGHI", "KLMN"), each = 5), rownames(Z))
  out <- repdist_bin(counts, Z, seqs = seqs)
  bins <- out$bins
  C <- as.matrix(seq_repdist(Z))

  expect_setequal(bins$reps$bin, unique(bins$bin))
  expect_equal(sum(bins$reps$size), ncol(counts))
  expect_equal(unname(bins$bin[bins$reps$protein]), bins$reps$bin)
  for (bin in bins$reps$bin) {
    members <- names(bins$bin)[bins$bin == bin]
    expect_equal(bins$reps$protein[bins$reps$bin == bin],
                 members[which.min(colSums(C[members, members, drop = FALSE]))])
  }
  expect_s4_class(out$representative_sequences, "AAStringSet")
  expect_equal(names(out$representative_sequences), bins$reps$bin)
  expect_equal(unname(as.character(out$representative_sequences)),
               unname(seqs[bins$reps$protein]))
  expect_error(repdist_bin(counts, Z, seqs = c(nope = "X")), "missing")
})
