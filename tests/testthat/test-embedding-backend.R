test_that("R dispatch preserves sequences and validates backend output", {
    s <- c(p1 = "ACDE", p2 = "KLMN")
    result <- diag(2)
    testthat::local_mocked_bindings(
        basiliskStart = function(...) NULL,
        basiliskStop = function(...) NULL,
        basiliskRun = function(proc, fun, seqs, spec, weights, device, batch) {
            expect_identical(seqs, as.list(unname(s)))
            expect_identical(device, "cpu")
            expect_identical(batch, 2L)
            expect_match(
                jsonlite::fromJSON(spec)$backbone_revision,
                "^[0-9a-f]{40}$"
            )
            result
        }, .package = "basilisk"
    )
    z <- embed_proteins(s, "esm2-8m", "cpu", 2)
    expect_identical(rownames(z), names(s))
    for (bad in list(matrix(0, 2, 2), matrix(NA_real_, 2, 2), diag(3))) {
        result <- bad
        expect_error(embed_proteins(s, "esm2-8m", "cpu", 2), "backend")
    }
})

test_that("all remote model artifacts have immutable revisions", {
    for (model in known_models()) {
        expect_match(model$backbone_revision, "^[0-9a-f]{40}$")
        if (!is.null(model$head_repo)) {
            expect_match(model$head_revision, "^[0-9a-f]{40}$")
        }
    }
})
