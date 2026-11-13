# Runs embed_proteins() against a real registered model through the
# basilisk environment, rather than the fabricated AAStringSet/embedding
# fixtures the other test files use. It is the only test that actually
# downloads model weights and provisions the pinned Python environment
# (basilisk + a standalone Python via uv), so it is off by default: everyday
# `devtools::test()` / `R CMD check` runs never touch the network or a
# multi-GB cache. It is skipped whenever the REPDIST_CHECK_MODEL environment
# variable is unset or empty.
#
# To run it, set REPDIST_CHECK_MODEL to any key from known_models() (e.g.
# "esm2-8m", the smallest registered model, to keep the download small) and
# run the tests from a shell, not from within an R session that has already
# loaded reticulate/basilisk for something else:
#   REPDIST_CHECK_MODEL=esm2-8m Rscript -e 'devtools::test()'
# or, scoped to just this file:
#   REPDIST_CHECK_MODEL=esm2-8m Rscript -e 'testthat::test_file("tests/testthat/test-embed-model.R")'
# The first run downloads the Python environment and model weights (up to a
# few GB for the default "tmvec" model); later runs reuse the cache.
test_that("a registered model embeds in order, unit-norm, batch-invariant", {
    model <- Sys.getenv("REPDIST_CHECK_MODEL")
    skip_if(!nzchar(model), "REPDIST_CHECK_MODEL is not set")
    s <- c(
        p1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
        p2 = "ACDEFGHIKL", p3 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    )
    z1 <- S4Vectors::mcols(embed_proteins(s, model, "cpu", 1))$embedding
    z3 <- S4Vectors::mcols(embed_proteins(s, model, "cpu", 3))$embedding
    expect_identical(rownames(z1), names(s))
    expect_equal(unname(rowSums(z1^2)), rep(1, 3), tolerance = 1e-6)
    expect_equal(z1[1, ], z1[3, ], tolerance = 1e-6)
    expect_equal(z1, z3, tolerance = 1e-5)
})
