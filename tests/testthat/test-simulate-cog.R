test_that("simulate_cog is aligned and matches its design", {
    data("simulate_cog", package = "repdist", envir = environment())
    x <- simulate_cog
    ids <- names(x$proteins)
    expect_identical(x$annotation$protein, ids)
    expect_identical(nrow(S4Vectors::mcols(x$proteins)$embedding), length(ids))
    expect_identical(colnames(x$counts), ids)
    expect_identical(rownames(x$counts), x$samples$sample)
    expect_type(x$counts, "integer")
    grpe <- x$annotation[x$annotation$role == "GrpE", ]
    # one GrpE protein per group-A sample, each in its own 30% cluster
    expect_setequal(grpe$enriched_in, x$samples$sample[x$samples$group == "A"])
    expect_false(anyDuplicated(grpe$mmseqs30) > 0)
    # a GrpE protein is counted only in its own sample
    own <- outer(x$samples$grpe_protein, grpe$protein, "==")
    own[is.na(own)] <- FALSE
    expect_true(all(x$counts[, grpe$protein][!own] == 0L))
    expect_true(all(x$counts[, grpe$protein][own] > 0L))
})
