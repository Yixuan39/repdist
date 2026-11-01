test_that("simulate_ec is aligned and matches its design", {
    data("simulate_ec", package = "repdist", envir = environment())
    x <- simulate_ec
    ids <- names(x$proteins)
    expect_identical(colnames(x$counts), ids)
    expect_identical(x$annotation$protein, ids)
    expect_identical(rownames(x$counts), x$samples$sample)
    expect_identical(nrow(S4Vectors::mcols(x$proteins)$embedding), length(ids))
    expect_type(x$counts, "integer")
    enr <- x$annotation[!is.na(x$annotation$enriched_in), ]
    # each enriched protein sits in its sample's class, one 50% cluster each
    expect_identical(enr$ec3,
        x$samples$enriched_ec3[match(enr$enriched_in, x$samples$sample)])
    expect_false(anyDuplicated(enr$mmseqs50) > 0)
})
