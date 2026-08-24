test_that("duplicate, missing, or empty names are rejected", {
  expect_error(embed_proteins(c("ACDE", "KLMN")), "names")
  expect_error(embed_proteins(c(p1 = "ACDE", p1 = "KLMN")), "unique")
  expect_error(embed_proteins(stats::setNames(c("ACDE", "KLMN"), c("p1", ""))),
               "names")
})

test_that("NA and empty sequences are rejected", {
  expect_error(embed_proteins(stats::setNames(character(), character())),
               "at least one")
  expect_error(embed_proteins(c(p1 = "ACDE", p2 = NA)), "missing")
  expect_error(embed_proteins(c(p1 = "ACDE", p2 = "")), "empty")
  expect_error(embed_proteins(c(p1 = "*")), "empty")
})

test_that("Biostrings QC rejects symbols unsupported by the embedding models", {
  expect_error(embed_proteins(c(p1 = "AC-DE")), "model-incompatible")
  expect_error(embed_proteins(c(p1 = "AC DE")))
  expect_error(embed_proteins(c(p1 = "AC*DE")), "model-incompatible")
  expect_error(embed_proteins(c(p1 = "ACJDE")), "model-incompatible")
})

test_that("registry keeps architectural limits separate from training ranges", {
  models <- known_models()
  expect_null(models[["tmvec-swissmodel-large"]]$context_max_length)
  expect_equal(models[["tmvec-swissmodel-large"]]$training_max_length, 1000)
  expect_equal(models$esm2$context_max_length, 1022)
})

test_that("the embedding memory control rejects invalid values before model loading", {
  seqs <- c(p1 = "ACDE")
  expect_error(embed_proteins(seqs, memory_fraction = 0), "in \\(0, 1\\]")
  expect_error(embed_proteins(seqs, memory_fraction = 1.1), "in \\(0, 1\\]")
})

test_that("length_abundance_profile returns two step curves that start at/near 100%", {
  skip_if_not_installed("ggplot2")
  set.seed(0)
  n <- 200
  lens <- stats::setNames(round(exp(rnorm(n, log(300), 0.6))) + 20L, paste0("p", 1:n))
  cc <- matrix(rpois(3 * n, 5) + 1L, 3, n, dimnames = list(paste0("s", 1:3), names(lens)))

  p <- length_abundance_profile(cc, lens, max_length = 1022)
  expect_s3_class(p, "ggplot")
  expect_setequal(unique(p$data$curve), c("proteins", "abundance"))
  expect_true(all(p$data$pct >= 0 & p$data$pct <= 100))

  by_curve <- split(p$data, p$data$curve)
  # at the shortest observed length, everything else is still "longer than it"
  expect_equal(by_curve$proteins$pct[which.min(by_curve$proteins$length)], 100 * (1 - 1 / n))
  # both curves must be non-increasing in length (survival functions)
  expect_true(all(diff(by_curve$proteins$pct) <= 1e-8))
  expect_true(all(diff(by_curve$abundance$pct) <= 1e-8))
})

test_that("length_abundance_profile rejects a sample with zero total abundance", {
  skip_if_not_installed("ggplot2")
  n <- 20
  lens <- stats::setNames(round(exp(rnorm(n, log(300), 0.4))) + 20L, paste0("p", 1:n))
  cc <- matrix(rpois(3 * n, 5) + 1L, 3, n, dimnames = list(paste0("s", 1:3), names(lens)))
  cc[1, ] <- 0L   # relative abundance for this sample is 0/0
  expect_error(length_abundance_profile(cc, lens))
})

test_that("length_abundance_profile is invariant to per-sample sequencing depth", {
  # the whole point of using relative, not raw, abundance: a sample sequenced
  # deeper but with the same underlying profile must not shift the curve
  skip_if_not_installed("ggplot2")
  set.seed(1)
  n <- 150
  lens <- stats::setNames(round(exp(rnorm(n, log(300), 0.5))) + 20L, paste0("p", 1:n))
  shallow <- matrix(rpois(4 * n, 4) + 1L, 4, n, dimnames = list(paste0("s", 1:4), names(lens)))
  deep <- shallow
  deep[1, ] <- deep[1, ] * 20L   # one sample at 20x depth, identical relative profile

  pct_shallow <- length_abundance_profile(shallow, lens)$data$pct
  pct_deep    <- length_abundance_profile(deep, lens)$data$pct
  expect_equal(pct_shallow, pct_deep, tolerance = 1e-8)
})

test_that("the abundance curve tracks where abundance actually concentrates, not protein count", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  n <- 150
  # two well-separated length clusters, short and long, equal protein count
  lens <- stats::setNames(c(rep(150L, n / 2), rep(5000L, n / 2)), paste0("p", 1:n))
  cc <- matrix(rpois(4 * n, 4) + 1L, 4, n, dimnames = list(paste0("s", 1:4), names(lens)))
  cc[, lens == 5000L] <- cc[, lens == 5000L] + 500L  # nearly all abundance in the long cluster

  d <- length_abundance_profile(cc, lens)$data
  proteins <- d[d$curve == "proteins", ]
  abund    <- d[d$curve == "abundance", ]
  # equal-sized clusters: right at the boundary, half the catalog remains by count ...
  boundary <- proteins$pct[which.min(abs(proteins$length - 150))]
  expect_equal(boundary, 50, tolerance = 5)
  # ... but almost none of the abundance is short, so almost all of it remains
  boundary_abund <- abund$pct[which.min(abs(abund$length - 150))]
  expect_gt(boundary_abund, 90)
})

test_that("length_abundance_profile accepts sequences in place of lengths", {
  skip_if_not_installed("ggplot2")
  seqs <- stats::setNames(vapply(1:30, function(i) strrep("A", 50 + i * 20), ""), paste0("p", 1:30))
  cc <- matrix(rpois(2 * 30, 5) + 1L, 2, 30, dimnames = list(paste0("s", 1:2), names(seqs)))
  expect_s3_class(length_abundance_profile(cc, seqs), "ggplot")
})
