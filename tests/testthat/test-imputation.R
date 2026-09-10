imputation_example <- function() {
    n <- 96L
    ids <- sprintf("p%03d", seq_len(n))
    Z <- outer(seq_len(n), seq_len(7), function(i, j)
        sin(i * j * sqrt(2)) + cos(i * j * sqrt(3)))
    Z <- Z / sqrt(rowSums(Z^2))
    rownames(Z) <- ids
    labels <- setNames(c("A", "B", "C")[1L +
        floor(abs(sin(seq_len(72) * sqrt(5))) * 2.99)], ids[seq_len(72)])
    list(D = repdist_matrix(Z), labels = labels,
         groups = setNames(rep(seq_len(n / 2), each = 2), ids))
}

test_that("KNN retains full coverage and traceable original-label evidence", {
    x <- imputation_example()
    original <- x$labels
    out <- impute_labels(x$D, x$labels, groups = x$groups)
    expect_identical(x$labels, original)
    expect_identical(out$protein, setdiff(attr(x$D, "Labels"), names(x$labels)))
    expect_false(anyNA(out$predicted_label))
    expect_identical(unname(x$labels[out$supporting_protein]),
                     out$predicted_label)
    expect_equal(out$supporting_similarity, 1 - usedist::dist_get(
        x$D, out$protein, out$supporting_protein))
    expect_true(all(out$n_support == 1L & out$n_neighbors == 1L))
    expect_named(attr(out, "models"), c("knn", "score"))
    expect_identical(attr(out, "models")$knn$labels, original)
    expect_equal(impute_labels(as.matrix(x$D), factor(x$labels),
        groups = x$groups)$confidence, out$confidence)

    opposite <- repdist_matrix(rbind(a = c(1, 0), q = c(-1, 0)))
    expect_warning(far <- impute_labels(opposite, c(a = "A")),
                   "confidence is unavailable")
    expect_identical(far$predicted_label, "A")
    expect_identical(far$supporting_protein, "a")
    expect_equal(far$supporting_similarity, -1)
    expect_true(is.na(far$confidence))
    expect_false(far$has_alternative)
    expect_equal(far$margin, 0)
})

test_that("KNN includes exact kth-distance ties and resolves votes by evidence", {
    D <- repdist_matrix(matrix(1, 5, 1,
        dimnames = list(c("a", "b", "c", "d", "q"), NULL)))
    labels <- c(a = "B", b = "B", c = "A", d = "A")
    expect_warning(out <- impute_labels(D, labels), "confidence is unavailable")
    expect_identical(out$predicted_label, "B")
    expect_identical(out$supporting_protein, "a")
    expect_identical(out$n_support, 2L)
    expect_identical(out$n_neighbors, 4L)
    expect_equal(out$agreement, .5)
    expect_equal(out$margin, 0)

    Z <- rbind(a = c(1, 0), b = c(.8, .6), c = c(.6, .8),
               d = c(0, 1), q = c(1, 0))
    D <- repdist_matrix(Z)
    expect_warning(one <- impute_labels(D, labels, k = 1),
                   "confidence is unavailable")
    expect_warning(three <- impute_labels(D, c(a = "A", b = "B", c = "B",
        d = "C"), k = 3), "confidence is unavailable")
    expect_identical(one$predicted_label, "B")
    expect_identical(three$predicted_label, "B")
    expect_identical(three$supporting_protein, "b")
    expect_identical(three$n_neighbors, 3L)
    expect_identical(three$n_support, 2L)
    expect_equal(three$margin, -.2)
})

test_that("calibration masks groups and matches references for source and target", {
    x <- imputation_example()
    out <- impute_labels(x$D, x$labels, groups = x$groups)
    cal <- attr(out, "calibration")
    fold <- setNames(cal$fold, cal$protein)
    expect_true(all(vapply(split(cal$fold, x$groups[cal$protein]),
        function(f) length(unique(f)) == 1L, logical(1))))
    expect_true(all(fold[cal$supporting_protein] != cal$fold))
    expect_false(any(cal$protein == cal$supporting_protein))
    expect_true(all(is.na(cal$correct[!cal$annotated])))
    expect_identical(cal$correct[cal$annotated], as.numeric(
        cal$predicted_label[cal$annotated] ==
        x$labels[cal$protein[cal$annotated]]))
    for (i in seq_len(nrow(cal))) {
        reference <- names(x$labels)[fold[names(x$labels)] != cal$fold[[i]]]
        d <- usedist::dist_get(x$D, cal$protein[[i]], reference)
        closest <- reference[order(d, reference, method = "radix")[[1L]]]
        expect_identical(cal$supporting_protein[[i]], closest)
        expect_equal(cal$label_count[[i]], sum(x$labels[reference] ==
            cal$predicted_label[[i]]))
    }
    expect_equal(cal$local_similarity, unname(vapply(cal$protein, function(id)
        1 - min(usedist::dist_get(x$D, id, setdiff(cal$protein, id))), numeric(1))))
})

test_that("automatic groups include unannotated single-linkage bridges", {
    D <- matrix(c(0, .04, .08, 1, .04, 0, .04, 1,
                  .08, .04, 0, 1, 1, 1, 1, 0), 4, 4,
        dimnames = list(c("a", "bridge", "b", "q"),
                        c("a", "bridge", "b", "q")))
    expect_warning(out <- impute_labels(D, c(a = "A", b = "B")),
                   "confidence is unavailable")
    cal <- attr(out, "calibration")
    expect_equal(length(unique(cal$group[1:3])), 1L)
    expect_equal(length(unique(cal$fold[1:3])), 1L)
    expect_true(all(is.na(cal$predicted_label[1:3])))
})

test_that("insufficient outcomes retain labels without invented confidence", {
    x <- imputation_example()
    unique_labels <- setNames(names(x$labels), names(x$labels))
    expect_warning(out <- impute_labels(x$D, unique_labels, groups = x$groups),
                   "confidence is unavailable")
    expect_true(all(is.na(out$confidence)))
    expect_true(all(out$predicted_label %in% unique_labels))
    cal <- attr(out, "calibration")
    expect_true(all(cal$correct[cal$annotated] == 0))
    expect_true(all(is.na(cal$correct[!cal$annotated])))

    one_group <- setNames(rep("all", length(x$groups)), names(x$groups))
    expect_warning(out <- impute_labels(x$D, x$labels, groups = one_group),
                   "Fewer than three annotated groups")
    expect_false(anyNA(out$predicted_label))
    expect_true(all(is.na(out$confidence)))
    expect_true(all(is.na(attr(out, "calibration")$supporting_protein)))
})

test_that("scoring preserves labels, random state, and protein-order invariance", {
    skip_if_not_installed("ranger")
    x <- imputation_example()
    baseline <- impute_labels(x$D, x$labels, groups = x$groups)
    label_columns <- c("protein", "predicted_label", "supporting_protein",
        "supporting_similarity", "n_support", "n_neighbors", "margin",
        "label_count", "local_similarity", "agreement", "has_alternative")
    set.seed(31)
    rng <- .Random.seed
    out <- impute_labels(x$D, x$labels, groups = x$groups)
    expect_identical(.Random.seed, rng)
    expect_identical(attr(out, "models")$score$status, "ok")
    expect_true(all(is.finite(out$confidence)))
    expect_true(all(out$confidence >= 0 & out$confidence <= 1))
    expect_equal(as.matrix(out[, label_columns]),
                 as.matrix(baseline[, label_columns]))
    reverse <- rev(seq_len(attr(x$D, "Size")))
    again <- impute_labels(as.matrix(x$D)[reverse, reverse],
        rev(x$labels), groups = rev(x$groups))
    at <- match(out$protein, again$protein)
    expect_equal(out$confidence, again$confidence[at], tolerance = 1e-10)
    expect_identical(out$predicted_label, again$predicted_label[at])
    expect_identical(out$score_extrapolation, again$score_extrapolation[at])
})

test_that("scoring features exclude has_alternative", {
    skip_if_not_installed("ranger")
    features <- c("similarity", "margin", "logSupport", "localSimilarity",
                  "agreement")
    x <- imputation_example()
    out <- impute_labels(x$D, x$labels, groups = x$groups)
    expect_identical(names(.imputeFeatures(out)), features)
    expect_identical(attr(out, "models")$score$columns, features)
    expect_true(all(out$has_alternative))
    expect_named(attr(out, "diagnostics"), c("nAnnotated", "nGroups", "nFolds",
        "nCorrect", "nIncorrect", "nMissingScores", "featureRanges"))
    expect_identical(rownames(attr(out, "diagnostics")$featureRanges), NULL)
    expect_identical(colnames(attr(out, "diagnostics")$featureRanges), features)
})

test_that("the forest fits, predicts in [0, 1] and reports importance", {
    skip_if_not_installed("ranger")
    x <- data.frame(similarity = seq_len(40) / 40,
                    margin = rev(seq_len(40)) / 40,
                    row.names = sprintf("p%02d", seq_len(40)))
    model <- .imputeFit(x, rep(0:1, 20), 1)
    expect_identical(model$status, "ok")
    p <- .imputePredict(model, x)
    expect_true(all(is.finite(p) & p >= 0 & p <= 1))
    expect_named(ranger::importance(model$fit), c("similarity", "margin"))
    expect_identical(.imputePredict(list(status = "no"), x), rep(NA_real_, 40))
    expect_identical(.imputeFit(x[1:10, ], rep(0:1, 5), 1)$status,
                     "Insufficient usable outcomes.")
    expect_identical(.imputeFit(x, rep(0:1, c(38, 2)), 1)$status,
                     "Insufficient usable outcomes.")
})

test_that("RF fit errors leave KNN labels available", {
    skip_if_not_installed("ranger")
    x <- imputation_example()
    testthat::local_mocked_bindings(
        ranger = function(...) stop("Scoring engine failed."), .package = "ranger")
    expect_warning(out <- impute_labels(x$D, x$labels, groups = x$groups),
        "Scoring engine failed")
    expect_false(anyNA(out$predicted_label))
    expect_true(all(is.na(out$confidence)))
})

test_that("input validation rejects malformed distances, labels and settings", {
    Z <- rbind(a = c(1, 0), b = c(0, 1), q = c(1, 0))
    D <- repdist_matrix(Z)
    expect_equal(nrow(impute_labels(D, c(a = "A", b = "B", q = "A"))), 0)
    expect_equal(nrow(impute_labels(repdist_matrix(Z[1, , drop = FALSE]),
        c(a = "A"))), 0)
    expect_error(impute_labels(D, c(a = NA_character_, b = " ")), "at least one")
    for (bad in list(c("A", "B"), c(a = "A", a = "B"), c(other = "A"),
        setNames("A", ""), setNames("A", NA_character_), c(a = 1),
        list(a = "A"), matrix("A", dimnames = list("a", "label"))))
        expect_error(impute_labels(D, bad), "labels")
    for (name in c("k", "nFolds", "seed")) {
        for (bad in list(0, -1, 1.5, NA, Inf, "2", NULL, c(1, 2))) {
            args <- list(D = D, labels = c(a = "A"))
            args[name] <- list(bad)
            expect_error(do.call(impute_labels, args), name)
        }
    }
    expect_error(impute_labels(D, c(a = "A"), nFolds = 2), "at least 3")
    for (bad in list(-1.1, 1.1, NA, Inf, c(.9, .95)))
        expect_error(impute_labels(D, c(a = "A"), groupSimilarity = bad),
                     "groupSimilarity")
    for (bad in list(c(a = 1), c(a = 1, b = 1, q = NA),
        c(a = 1, a = 2, q = 3), c(a = "", b = "b", q = "q"),
        seq_len(3), c(a = 1, b = 2, other = 3)))
        expect_error(impute_labels(D, c(a = "A"), groups = bad), "groups")
    expect_error(impute_labels(stats::dist(Z), c(a = "A")), "euclidean")
    for (value in c(NA, Inf, -.1, 2.1)) {
        bad <- D
        bad[1] <- value
        expect_error(impute_labels(bad, c(a = "A")), "finite|\\[0, 2\\]")
    }
    bad <- as.matrix(D)
    bad[1, 2] <- .5
    expect_error(impute_labels(bad, c(a = "A")), "symmetric")
    bad <- D
    attr(bad, "Labels")[1] <- "b"
    expect_error(impute_labels(bad, c(a = "A")), "unique")
    bad <- D[-1L]
    class(bad) <- "dist"
    attributes(bad) <- attributes(D)
    expect_error(impute_labels(bad, c(a = "A")), "valid dist")
    expect_error(impute_labels(unname(as.matrix(D)), c(a = "A")), "named")
})
