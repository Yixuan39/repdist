test_that("duplicate, missing, or empty names are rejected", {
  expect_error(embed_proteins(c("ACDE", "KLMN")), "names")
  expect_error(embed_proteins(c(p1 = "ACDE", p1 = "KLMN")), "unique")
  expect_error(embed_proteins(stats::setNames(c("ACDE", "KLMN"), c("p1", ""))),
               "names")
})

test_that("an oversized catalog is rejected before any embedding work", {
  big <- stats::setNames(rep("ACDE", 46342L), paste0("p", seq_len(46342L)))
  expect_error(embed_proteins(big), "46341")
  expect_error(embed_proteins(big), "redundancy")
})

test_that("NA and empty sequences are rejected", {
  expect_error(embed_proteins(stats::setNames(character(), character())),
               "at least one")
  expect_error(embed_proteins(c(p1 = "ACDE", p2 = NA)), "missing")
  expect_error(embed_proteins(c(p1 = "ACDE", p2 = "")), "empty")
  expect_error(embed_proteins(c(p1 = "*")), "empty")
})

test_that("gaps, internal stops and stray symbols are rejected", {
  expect_error(embed_proteins(c(p1 = "AC-DE")), "model-incompatible")
  expect_error(embed_proteins(c(p1 = "AC DE")))
  expect_error(embed_proteins(c(p1 = "AC*DE")), "model-incompatible")
  expect_error(embed_proteins(c(p1 = "ACJDE")), "model-incompatible")
})

test_that("training ranges are recorded but never enforced", {
  models <- known_models()
  expect_equal(models$tmvec$training_max_length, 1000)
  expect_equal(models$`esm2-650m`$training_max_length, 1022)
  # no length limit at inference: a sequence past the training range is embedded
  expect_true(all(vapply(models, function(m) is.null(m$context_max_length),
                         logical(1))))
})

test_that("batch_size is validated before model loading", {
  seqs <- c(p1 = "ACDE")
  expect_error(embed_proteins(seqs, batch_size = 0), ">= 1")
  expect_error(embed_proteins(seqs, batch_size = c(4, 8)), ">= 1")
  expect_error(embed_proteins(seqs, batch_size = NA), ">= 1")
})

test_that("an unregistered model is rejected before Python starts", {
  seqs <- c(p1 = "ACDE")
  # the pinned basilisk env fixes which architectures load, so a repo outside
  # the registry has nothing behind it; erroring beats a wrong embedding
  expect_error(embed_proteins(seqs, model = "facebook/esm2_t6_8M_UR50D"),
               "Unknown model")
  expect_error(embed_proteins(seqs, model = "tmvec2"), "tmvec, tmvec-300")
})

test_that("plot_similarity_profile reports each protein's nearest neighbour", {
  skip_if_not_installed("ggplot2")
  Z <- rbind(p1 = c(1, 0), p2 = c(0.999, 0.04), p3 = c(0, 1))
  Z <- Z / sqrt(rowSums(Z^2))
  p <- plot_similarity_profile(repdist_matrix(Z), min_sim = 0.7)
  expect_s3_class(p, "ggplot")
  # p1 and p2 are mutual nearest neighbours; p3's nearest is whichever of the
  # two it is more similar to
  expect_equal(p$data$similarity[1:2], rep(sum(Z[1, ] * Z[2, ]), 2))
  expect_equal(p$data$similarity[3], max(sum(Z[3, ] * Z[1, ]), sum(Z[3, ] * Z[2, ])))
  expect_error(plot_similarity_profile(repdist_matrix(Z, "euclidean")), "euclidean")
})

test_that("plot_bin_profile sweeps the tree monotonically", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  Z <- matrix(rnorm(60), 30, 2, dimnames = list(paste0("p", 1:30), NULL))
  Z <- Z / sqrt(rowSums(Z^2))
  b <- bin_proteins(repdist_matrix(Z), min_sim = 0.7)
  p <- plot_bin_profile(b, min_sim = 0.7)
  expect_s3_class(p, "ggplot")
  d <- split(p$data, p$data$curve)
  # raising the threshold can only split bins, never merge them
  expect_true(all(diff(d$bins$value) >= 0))
  expect_true(all(diff(d$`largest bin`$value) <= 0))
  # and every split can only strand more proteins on their own
  expect_true(all(diff(d$singletons$value) >= 0))
  expect_equal(d$singletons$value[[100]], 100)   # min_sim = 1: all singletons
  # the bare hclust tree is accepted too; anything else is not
  expect_s3_class(plot_bin_profile(b$tree), "ggplot")
  expect_error(plot_bin_profile(list(tree = 1)), "hclust")
})

test_that("plot_bin_similarity shows one bin's pairwise similarities", {
  skip_if_not_installed("ggplot2")
  Z <- rbind(p1 = c(1, 0), p2 = c(0.99, 0.14), p3 = c(0.97, 0.24), p4 = c(0, 1))
  Z <- Z / sqrt(rowSums(Z^2))
  D <- repdist_matrix(Z)
  b <- bin_proteins(D, min_sim = 0.9)

  p <- plot_bin_similarity(b, "bin1", D)
  expect_s3_class(p, "ggplot")
  members <- b$clusters$bin1
  expect_equal(nrow(p$data), length(members)^2)
  # every off-diagonal cell is the ground metric's similarity, and the block
  # honours the threshold bin_proteins() promised
  expect_setequal(as.character(p$data$protein), members)
  expect_equal(max(p$data$similarity), 1)
  expect_gte(min(p$data$similarity), 0.9)
  expect_equal(sort(p$data$similarity),
               sort(as.numeric(1 - as.matrix(usedist::dist_subset(D, members)))))

  expect_error(plot_bin_similarity(b, 1L, D), "must name one")   # names only
  expect_error(plot_bin_similarity(b, "bin99", D), "must name one")
  expect_error(plot_bin_similarity(b, "bin1", repdist_matrix(Z, "euclidean")),
               "euclidean")
  expect_error(plot_bin_similarity(list(tree = 1), "bin1", D), "bin_proteins")
  singleton <- names(which(lengths(b$clusters) == 1L))[[1]]
  expect_error(plot_bin_similarity(b, singleton, D), "no pair")
})
