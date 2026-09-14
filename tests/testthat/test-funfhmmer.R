set.seed(0)

# three tight families, far apart
base <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Z <- do.call(rbind, lapply(1:3, function(k) {
  matrix(base[k, ], 5, 3, byrow = TRUE) + matrix(rnorm(15, 0, 0.03), 5, 3)
}))
Z <- Z / sqrt(rowSums(Z^2))
rownames(Z) <- paste0("p", 1:15)
family <- rep(1:3, each = 5)
D <- repdist_matrix(Z)

test_that("funfhmmer returns a full partition shaped like bin_proteins", {
  out <- funfhmmer(D)
  expect_named(out, c("clusters", "similarity", "representatives", "tree",
                      "separation", "auto_merge", "test"))
  expect_equal(out$test, "distance")
  expect_setequal(unlist(out$clusters), rownames(Z))
  expect_equal(sum(lengths(out$clusters)), 15L)   # no protein in two families
  expect_named(out$representatives, names(out$clusters))
  expect_s3_class(out$tree, "hclust")
  expect_equal(dim(out$similarity),
               rep(length(out$clusters), 2))
  expect_false(is.unsorted(rev(lengths(out$clusters))))
})

test_that("the three planted families are never merged with each other", {
  membership <- setNames(rep(names(funfhmmer(D)$clusters),
                             lengths(funfhmmer(D)$clusters)),
                         unlist(funfhmmer(D)$clusters))
  for (fam in 1:3)
    for (other in setdiff(1:3, fam))
      expect_length(intersect(membership[rownames(Z)[family == fam]],
                              membership[rownames(Z)[family == other]]), 0L)
})

test_that("separation trades family count against size", {
  expect_gte(length(funfhmmer(D, separation = 0.5)$clusters),
             length(funfhmmer(D, separation = 5)$clusters))
  expect_length(funfhmmer(D, separation = 1e6, min_sim = NULL)$clusters, 1L)
})

test_that("auto_merge merges the low end of the tree without testing", {
  expect_gte(length(funfhmmer(D, auto_merge = 0)$clusters),
             length(funfhmmer(D, auto_merge = 1)$clusters))
  expect_length(funfhmmer(D, auto_merge = 1, min_sim = NULL)$clusters, 1L)
})

test_that("the cut is adaptive, not one height for the whole tree", {
  # An evenly spread triple 10 degrees apart, and two tight pairs only 5
  # degrees apart from each other. The triple merges across a wider gap than
  # the one that splits the pairs, so no single cutree height reproduces the
  # partition -- it can only come from testing each branch on its own scale.
  angle <- c(w1 = 0, w2 = 10, w3 = 20, a1 = 90, a2 = 91, b1 = 96, b2 = 97)
  W <- cbind(cos(angle * pi / 180), sin(angle * pi / 180))
  rownames(W) <- names(angle)
  out <- funfhmmer(repdist_matrix(W), separation = 3, auto_merge = 0,
                   min_sim = NULL)

  membership <- setNames(rep(names(out$clusters), lengths(out$clusters)),
                         unlist(out$clusters))
  expect_length(unique(membership[c("w1", "w2", "w3")]), 1L)
  expect_length(unique(membership[c("a1", "a2", "b1", "b2")]), 2L)

  # no height of the same tree gives this partition
  co <- outer(membership, membership, "==")[rownames(W), rownames(W)]
  for (h in c(0, out$tree$height)) {
    cut <- stats::cutree(out$tree, h = h)[rownames(W)]
    expect_false(identical(unname(outer(cut, cut, "==")), unname(co)))
  }
})

test_that("a supplied tree is traversed rather than a fresh one", {
  own <- bin_proteins(D, min_sim = 0.5)$tree
  expect_identical(funfhmmer(D, tree = own)$tree, own)
  expect_error(funfhmmer(D, tree = fastcluster::hclust(stats::dist(Z[1:5, ]))),
               "same 15 proteins")
})

test_that("representatives are within-family medoids", {
  out <- funfhmmer(D)
  C <- as.matrix(D)
  for (fam in names(out$clusters)) {
    members <- out$clusters[[fam]]
    expect_equal(out$representatives[[fam]],
                 members[which.min(colSums(C[members, members, drop = FALSE]))])
  }
})

test_that("invalid input is rejected the way bin_proteins rejects it", {
  expect_error(funfhmmer(stats::dist(Z)), "euclidean")
  expect_error(funfhmmer(D, separation = 0))
  expect_error(funfhmmer(D, separation = c(1, 2)))
  expect_error(funfhmmer(D, auto_merge = 1.5))
  bad <- D
  bad[1] <- NA
  expect_error(funfhmmer(bad), "non-finite")
})

test_that("min_sim rails the result to a refinement of bin_proteins", {
  for (s in c(0.3, 0.5, 0.8)) {
    bins <- bin_proteins(D, min_sim = s)$clusters
    fams <- funfhmmer(D, min_sim = s, separation = 0.1, auto_merge = 0)$clusters
    # every family sits inside one bin, however hard the test pushes to split
    for (members in fams)
      expect_length(Filter(function(b) all(members %in% b), bins), 1L)
  }
})

test_that("no family breaks the all-pairs similarity floor", {
  similarity <- 1 - as.matrix(D)
  for (s in c(0.3, 0.5, 0.8)) {
    fams <- funfhmmer(D, min_sim = s, separation = 1e6, auto_merge = 1)$clusters
    worst <- vapply(fams, function(members) {
      if (length(members) < 2L) return(1)
      min(similarity[members, members][upper.tri(diag(length(members)))])
    }, numeric(1))
    expect_true(all(worst >= s))
  }
})

test_that("the rail overrides both the test and the auto-merge head", {
  # auto_merge = 1 says take every merge untested; min_sim still refuses
  expect_gt(length(funfhmmer(D, min_sim = 0.999, auto_merge = 1)$clusters),
            length(funfhmmer(D, min_sim = 0.3, auto_merge = 1)$clusters))
  # and with the rail off, auto_merge = 1 really does take everything
  expect_length(funfhmmer(D, min_sim = NULL, auto_merge = 1)$clusters, 1L)
})
