test_that("duplicate, missing, or empty names are rejected", {
    expect_error(embed_proteins(c("ACDE", "KLMN")), "names")
    expect_error(embed_proteins(c(p1 = "ACDE", p1 = "KLMN")), "unique")
    expect_error(
        embed_proteins(stats::setNames(c("ACDE", "KLMN"), c("p1", ""))),
        "names"
    )
    expect_error(embed_proteins(stats::setNames("ACDE", NA_character_)), "names")
    expect_error(embed_proteins(stats::setNames("ACDE", " ")), "names")
})

test_that("NA and empty sequences are rejected", {
    expect_error(
        embed_proteins(stats::setNames(character(), character())),
        "at least one"
    )
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
    expect_true(all(vapply(
        models, function(m) is.null(m$context_max_length),
        logical(1)
    )))
})

test_that("batch_size is validated before model loading", {
    seqs <- c(p1 = "ACDE")
    expect_error(embed_proteins(seqs, batch_size = 0), "batch_size")
    expect_error(embed_proteins(seqs, batch_size = c(4, 8)), "batch_size")
    expect_error(embed_proteins(seqs, batch_size = NA), "batch_size")
    for (bad in list(1.5, 2^31, TRUE, "2", NULL)) {
        expect_error(embed_proteins(seqs, batch_size = bad), "batch_size")
    }
    for (bad in list(3, NA_character_, character(), c("tmvec", "esm2-8m"))) {
        expect_error(embed_proteins(seqs, model = bad), "model")
    }
})

test_that("an unregistered model is rejected before Python starts", {
    seqs <- c(p1 = "ACDE")
    # the pinned basilisk env fixes which architectures load, so a repo outside
    # the registry has nothing behind it; erroring beats a wrong embedding
    expect_error(
        embed_proteins(seqs, model = "facebook/esm2_t6_8M_UR50D"),
        "Unknown model"
    )
    expect_error(embed_proteins(seqs, model = "tmvec2"), "tmvec, tmvec-300")
})
