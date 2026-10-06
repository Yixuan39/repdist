# The EC3 simulation design: group A upregulates one EC3 through private
# accessions, group B scatters across the others, and no accession is
# upregulated in two samples. test_file() runs it from this directory:
# testthat::test_file("validation/test-simulation.R").
library(testthat)
sim_env <- new.env()
sys.source("simulate_ec3.R", envir = sim_env)
u <- readRDS("microbial_ec_universe.rds")

test_that("upregulated accessions are disjoint and group-specific", {
    set.seed(1)
    target <- u$members$ec3[[1]]
    d <- sim_env$design_upregulation(u$members, target, 10L, 5L)
    ec3 <- stats::setNames(u$members$ec3, u$members$protein)
    all_up <- unlist(d$up)
    expect_length(d$up, 20L)
    expect_true(all(lengths(d$up) == 5L))
    expect_false(anyDuplicated(all_up) > 0)
    a <- d$group == "A"
    expect_true(all(ec3[unlist(d$up[a])] == target))
    expect_false(any(ec3[unlist(d$up[!a])] == target))
    # group B: distinct EC3s within a sample, spread evenly across the group
    for (s in d$up[!a]) expect_false(anyDuplicated(ec3[s]) > 0)
    n_other <- length(unique(ec3)) - 1L
    expect_lte(max(table(ec3[unlist(d$up[!a])])), ceiling(50 / n_other))
})

test_that("counts keep every protein in every sample and fold only where up", {
    set.seed(2)
    d <- sim_env$design_upregulation(u$members, u$members$ec3[[1]], 10L, 5L)
    sim <- sim_env$simulate_counts(u$members$protein, d)
    expect_identical(dim(sim$counts), c(20L, nrow(u$members)))
    expect_true(all(sim$counts == round(sim$counts)))
    expect_true(all(sim$truth$baseline > 0))
    expect_identical(sim$truth$up, d$up)
    expect_named(sim$truth$size_factor, names(d$up))
})

test_that("both groups move the same baseline abundance", {
    set.seed(3)
    d <- sim_env$design_upregulation(u$members, u$members$ec3[[1]], 10L, 5L)
    base <- sim_env$simulate_counts(u$members$protein, d)$truth$baseline
    for (k in 1:10) {
        expect_equal(sort(unname(base[d$up[[k]]])),
            sort(unname(base[d$up[[10 + k]]])))
    }
})

test_that("private target effects average to fourfold EC3 ground truth", {
    set.seed(3)
    target <- u$members$ec3[[1]]
    d <- sim_env$design_upregulation(u$members, target)
    sim <- sim_env$simulate_counts(u$members$protein, d)
    base <- sim$truth$baseline
    target_total <- sum(base[u$members$protein[u$members$ec3 == target]])
    expected <- vapply(d$up[d$group == "A"], function(ids) {
        1 + (sim$truth$fold - 1) * sum(base[ids]) / target_total
    }, numeric(1))
    expect_equal(mean(expected), 4)
})

test_that("validation records failed fits as NA with their error", {
    template <- c(recall = 0.5, discoveries = 3)
    failed <- sim_env$record_fit(stop("dispersion fit failed"), template)
    expect_identical(failed$scores, template * NA_real_)
    expect_identical(failed$error, "dispersion fit failed")
    ok <- sim_env$record_fit(template, template)
    expect_identical(ok$scores, template)
    expect_true(is.na(ok$error))
})
