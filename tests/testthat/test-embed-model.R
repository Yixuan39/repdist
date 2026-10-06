# Real inference through the basilisk environment. Off by default: it
# downloads model weights. Run with, e.g.,
#   REPDIST_CHECK_MODEL=esm2-8m Rscript -e 'testthat::test_local()'
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
