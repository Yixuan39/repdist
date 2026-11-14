# The only test that runs a real model: embed_proteins() through the basilisk
# Python environment, which downloads the environment and the model weights
# (up to a few GB) on first use. Every other test uses the fabricated
# embeddings of helper-repdist.R. It is skipped unless the REPDIST_CHECK_MODEL
# environment variable names a model from known_models(), set in the shell
# that runs the tests, e.g. with the smallest model:
#   REPDIST_CHECK_MODEL=esm2-8m Rscript -e 'devtools::test()'
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
