# Run from the package root: Rscript inst/scripts/annotation_imputation.R
# Case study 2: nested TM95-grouped evaluation of 1NN correctness scores,
# plus an annotation-evidence check inside the unannotated population.
code <- c("R/imputation.R", "inst/scripts/annotation_imputation.R")
code_md5 <- tools::md5sum(code)
source("R/qc.R")
source("R/distances.R")
source("R/imputation.R")
stopifnot(requireNamespace("ranger", quietly = TRUE))

input <- c("results/case_study2_annotation_clustering/comparison.rds",
           "results/case_study2_annotation_clustering/annotations.csv")
input_md5 <- tools::md5sum(input)
out_dir <- "results/annotation_imputation_rf"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
comparison <- readRDS(input[[1L]])
D <- comparison$G
ids <- attr(D, "Labels")
labels <- setNames(comparison$labels$consensus_KO, comparison$labels$protein)[ids]
annotated <- !is.na(labels) & nzchar(trimws(labels))
labels[!annotated] <- NA_character_
rm(comparison)
stopifnot(identical(names(labels), ids), !anyNA(ids), !anyDuplicated(ids))
groups <- stats::cutree(fastcluster::hclust(D, method = "single"), h = 1 - .95)

# Components include unannotated bridges. Balance annotation counts without
# using label identities, and keep every component within one outer fold.
seed <- 20260909L
members <- split(seq_along(ids), groups)
counts <- vapply(members, function(at) sum(annotated[at]), integer(1))
key <- vapply(members, function(at) sort(ids[at], method = "radix")[[1L]], "")
order_groups <- order(-counts, digest::digest2int(key, seed = seed), key,
                      method = "radix")
fold <- integer(length(ids))
load <- integer(5L)
for (g in order_groups) {
    f <- which.min(load)
    fold[members[[g]]] <- f
    load[[f]] <- load[[f]] + counts[[g]]
}
stopifnot(all(load > 0L), all(vapply(split(fold, groups), function(x)
    length(unique(x)) == 1L, logical(1))))
write.csv(data.frame(protein = ids, annotated, label = unname(labels),
    group = unname(groups[ids]), outer_fold = fold, seed),
    file.path(out_dir, "splits.csv"), row.names = FALSE)
write.csv(data.frame(proteins = length(ids), annotated = sum(annotated),
    unannotated = sum(!annotated), labels = length(unique(labels[annotated])),
    singleton_labels = sum(table(labels[annotated]) == 1L),
    TM95_groups = length(unique(groups)),
    annotated_TM95_groups = length(unique(groups[annotated])),
    largest_TM95_group = max(lengths(members))),
    file.path(out_dir, "cohort.csv"), row.names = FALSE)

auc <- function(score, y) {
    if (length(unique(y)) < 2L) return(NA_real_)
    r <- rank(score)
    (sum(r[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) /
        (sum(y == 1) * sum(y == 0))
}
score_metrics <- function(x, outcome = "correct") {
    y <- x[[outcome]]
    keep <- is.finite(x$confidence) & !is.na(y)
    x <- x[keep, , drop = FALSE]
    y <- y[keep]
    if (!nrow(x)) return(data.frame(n_scored = 0L, accuracy = NA_real_,
        mean_score = NA_real_, brier = NA_real_, log_loss = NA_real_,
        auc = NA_real_))
    p <- pmin(pmax(x$confidence, 1e-6), 1 - 1e-6)
    data.frame(n_scored = nrow(x), accuracy = mean(y),
        mean_score = mean(x$confidence), brier = mean((x$confidence - y)^2),
        log_loss = mean(-y * log(p) - (1 - y) * log1p(-p)),
        auc = auc(x$confidence, y))
}

predictions <- runtimes <- diagnostics <- list()
for (f in seq_len(5L)) {
    test <- which(annotated & fold == f)
    train <- which(annotated & fold != f)
    truth <- unname(labels[test])
    reference <- labels[train]
    message("Outer fold ", f, "/5: ", length(train), " reference proteins; ",
            length(test), " hidden labels")
    elapsed <- system.time(result <- impute_labels(D, reference, k = 1L,
        groups = groups, nFolds = 5L, seed = seed))[["elapsed"]]
    models <- attr(result, "models")
    cal <- attr(result, "calibration")
    diagnostics[[paste0("fold", f)]] <- attr(result, "diagnostics")
    stopifnot(setequal(models$knn$training_proteins, names(reference)),
        !any(models$knn$training_proteins %in% ids[test]),
        !any(cal$annotated[match(ids[test], cal$protein)]),
        all(is.na(cal$correct[match(ids[test], cal$protein)])))
    at <- match(ids[test], result$protein)
    stopifnot(!anyNA(at))
    result <- result[at, , drop = FALSE]
    stopifnot(all(result$supporting_protein %in% names(reference)),
        all(result$predicted_label == reference[result$supporting_protein]),
        all(is.finite(result$confidence)),
        all(result$confidence >= 0 & result$confidence <= 1))
    base <- cbind(data.frame(outer_fold = f), result, truth,
        correct = as.integer(result$predicted_label == truth),
        label_seen = truth %in% reference)
    predictions[[length(predictions) + 1L]] <- cbind(data.frame(model = "rf"), base)
    # Baselines: the inner out-of-fold mean correctness, and raw similarity
    # clipped to [0, 1]. Neither is fitted on outer-test outcomes.
    constant <- base
    constant$confidence <- mean(cal$correct[cal$annotated], na.rm = TRUE)
    predictions[[length(predictions) + 1L]] <-
        cbind(data.frame(model = "constant"), constant)
    similarity <- base
    similarity$confidence <- pmin(1, pmax(0, similarity$supporting_similarity))
    predictions[[length(predictions) + 1L]] <-
        cbind(data.frame(model = "similarity"), similarity)
    runtimes[[length(runtimes) + 1L]] <- data.frame(outer_fold = f,
        n_reference = length(train), n_test = length(test),
        elapsed_seconds = elapsed)
    print(cbind(data.frame(outer_fold = f), score_metrics(base)),
          row.names = FALSE)
}
predictions <- do.call(rbind, predictions)
write.csv(predictions, file.path(out_dir, "test_predictions.csv"), row.names = FALSE)
fold_metrics <- do.call(rbind, lapply(split(predictions,
    interaction(predictions$model, predictions$outer_fold, drop = TRUE)), function(x)
        cbind(x[1L, c("model", "outer_fold")], score_metrics(x))))
write.csv(fold_metrics, file.path(out_dir, "fold_metrics.csv"), row.names = FALSE)
summary <- do.call(rbind, lapply(split(predictions, predictions$model), function(x)
    cbind(data.frame(model = x$model[[1L]], n_test = nrow(x),
        n_labels_unseen = sum(!x$label_seen),
        seen_label_accuracy = mean(x$correct[x$label_seen]),
        n_extrapolation = sum(x$score_extrapolation)), score_metrics(x))))
write.csv(summary, file.path(out_dir, "test_summary.csv"), row.names = FALSE)

# Fixed bins and thresholds describe score reliability; none selects the fit.
calibration <- do.call(rbind, lapply(split(predictions, predictions$model), function(x) {
    x$bin <- cut(x$confidence, seq(0, 1, .1), include.lowest = TRUE)
    do.call(rbind, lapply(split(x, x$bin, drop = TRUE), function(z)
        data.frame(model = z$model[[1L]], bin = as.character(z$bin[[1L]]),
            n = nrow(z), mean_score = mean(z$confidence),
            observed_accuracy = mean(z$correct))))
}))
write.csv(calibration, file.path(out_dir, "calibration.csv"), row.names = FALSE)
predictions$similarity_bin <- cut(predictions$supporting_similarity,
    c(-Inf, .5, .6, .7, .8, .9, .95, Inf), right = FALSE)
strata <- do.call(rbind, lapply(split(predictions, interaction(predictions$model,
    predictions$similarity_bin, drop = TRUE)), function(x)
        cbind(x[1L, c("model", "similarity_bin")], score_metrics(x))))
write.csv(strata, file.path(out_dir, "similarity_strata.csv"), row.names = FALSE)
thresholds <- do.call(rbind, lapply(split(predictions, predictions$model), function(x)
    do.call(rbind, lapply(c(.5, .8, .9), function(threshold) {
        chosen <- x$confidence >= threshold
        data.frame(model = x$model[[1L]], threshold, selected = sum(chosen),
            coverage = mean(chosen), accuracy = if (any(chosen))
                mean(x$correct[chosen]) else NA_real_)
    }))))
write.csv(thresholds, file.path(out_dir, "thresholds.csv"), row.names = FALSE)

message("Full catalog fit")
elapsed <- system.time(full <- impute_labels(D, labels[annotated], k = 1L,
    groups = groups, nFolds = 5L, seed = seed))[["elapsed"]]
stopifnot(identical(full$protein, ids[!annotated]), all(is.finite(full$confidence)))
diagnostics$full <- attr(full, "diagnostics")
runtimes[[length(runtimes) + 1L]] <- data.frame(outer_fold = 0L,
    n_reference = sum(annotated), n_test = sum(!annotated),
    elapsed_seconds = elapsed)
write.csv(attr(full, "calibration"), file.path(out_dir, "full_calibration.csv"),
          row.names = FALSE)
write.csv(full, file.path(out_dir, "unannotated_predictions.csv"), row.names = FALSE)
write.csv(do.call(rbind, runtimes), file.path(out_dir, "runtimes.csv"), row.names = FALSE)
importance <- ranger::importance(attr(full, "models")$score$fit)
write.csv(data.frame(feature = names(importance),
    permutation_importance = unname(importance)),
    file.path(out_dir, "feature_importance.csv"), row.names = FALSE)
saveRDS(attr(full, "models"), file.path(out_dir, "models.rds"))

# ---------------------------------------------------------------------------
# Annotation-evidence stratification.
#
# Consensus KO requires agreement between pipelines, so proteins whose tools
# disagreed, or where only one tool spoke, are formally unannotated. They sit
# in the target population and still carry a checkable call. Matching the
# copied KO against those calls is a calibration check on unannotated
# proteins, which grouped calibration on annotated proteins cannot provide.
#
# Circularity: EggNOG is homology-based and contributed to the donor's own
# consensus KO, so agreement with a query's EggNOG call is partly homology
# confirming homology. Every quantity is therefore reported twice: against
# all calls, and against non-EggNOG calls only.
# ---------------------------------------------------------------------------
tools <- c("MantisKO", "MicrobeKO", "EggNOG_KO")
independent <- setdiff(tools, "EggNOG_KO")
ann <- read.csv(input[[2L]], stringsAsFactors = FALSE)
ann <- ann[!duplicated(ann$protein), ]
evidence <- merge(full, ann[, c("protein", tools)], by = "protein")
stopifnot(nrow(evidence) == nrow(full))
matches <- function(x, cols) vapply(seq_len(nrow(x)), function(i)
    as.integer(x$predicted_label[[i]] %in%
        stats::na.omit(unlist(x[i, cols, drop = FALSE]))), integer(1))
evidence$n_calls <- rowSums(!is.na(evidence[, tools]))
evidence$n_independent <- rowSums(!is.na(evidence[, independent]))
evidence <- evidence[evidence$n_calls >= 1L, , drop = FALSE]
evidence$match_any <- matches(evidence, tools)
evidence$match_independent <- matches(evidence, independent)
evidence$depth <- ifelse(evidence$n_calls == 1L, "1 call",
    paste(evidence$n_calls, "conflicting calls"))
write.csv(evidence[, c("protein", "predicted_label", "supporting_protein",
    "supporting_similarity", "confidence", "score_extrapolation", tools,
    "n_calls", "n_independent", "depth", "match_any", "match_independent")],
    file.path(out_dir, "evidence_predictions.csv"), row.names = FALSE)

evidence_metrics <- function(x, outcome, subset = rep(TRUE, nrow(x))) {
    x <- x[subset, , drop = FALSE]
    do.call(rbind, lapply(split(x, x$depth), function(z)
        cbind(data.frame(criterion = outcome, depth = z$depth[[1L]]),
              score_metrics(z, outcome))))
}
evidence_summary <- rbind(
    cbind(data.frame(criterion = "any call", depth = "all"),
          score_metrics(evidence, "match_any")),
    evidence_metrics(evidence, "match_any"),
    cbind(data.frame(criterion = "non-EggNOG call", depth = "all"),
          score_metrics(evidence[evidence$n_independent > 0L, ],
                        "match_independent")),
    evidence_metrics(evidence, "match_independent", evidence$n_independent > 0L))
evidence_summary$criterion[evidence_summary$criterion == "match_any"] <- "any call"
evidence_summary$criterion[evidence_summary$criterion == "match_independent"] <-
    "non-EggNOG call"
write.csv(evidence_summary, file.path(out_dir, "evidence_summary.csv"),
          row.names = FALSE)

evidence$bin <- cut(evidence$confidence, seq(0, 1, .2), include.lowest = TRUE)
evidence_calibration <- do.call(rbind, lapply(c("match_any", "match_independent"),
    function(outcome) {
        x <- if (outcome == "match_any") evidence else
            evidence[evidence$n_independent > 0L, ]
        do.call(rbind, lapply(split(x, x$bin, drop = TRUE), function(z)
            data.frame(criterion = outcome, bin = as.character(z$bin[[1L]]),
                n = nrow(z), mean_score = mean(z$confidence),
                observed_match = mean(z[[outcome]]))))
    }))
write.csv(evidence_calibration, file.path(out_dir, "evidence_calibration.csv"),
          row.names = FALSE)

png(file.path(out_dir, "calibration.png"), width = 1400, height = 700, res = 140)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
x <- calibration[calibration$model == "rf", ]
plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted probability",
     ylab = "Observed fraction correct", main = "Held-out annotated proteins")
abline(0, 1, col = "grey70", lty = 2)
lines(x$mean_score, x$observed_accuracy, type = "b", pch = 16, col = "#0072B2")
plot(c(0, 1), c(0, 1), type = "n", xlab = "Mean predicted probability",
     ylab = "Observed fraction matching a pipeline call",
     main = "Unannotated proteins with pipeline evidence")
abline(0, 1, col = "grey70", lty = 2)
colors <- c(match_any = "#D55E00", match_independent = "#009E73")
for (outcome in names(colors)) {
    z <- evidence_calibration[evidence_calibration$criterion == outcome, ]
    lines(z$mean_score, z$observed_match, type = "b", pch = 16,
          col = colors[[outcome]])
}
legend("topleft", legend = c("any call", "non-EggNOG call"), col = colors,
       pch = 16, lty = 1, bty = "n", cex = .8)
invisible(dev.off())

saveRDS(diagnostics, file.path(out_dir, "diagnostics.rds"))
stopifnot(identical(code_md5, tools::md5sum(code)),
          identical(input_md5, tools::md5sum(input)))
write.csv(data.frame(file = input, md5 = unname(input_md5)),
    file.path(out_dir, "inputs.csv"), row.names = FALSE)
write.csv(data.frame(file = code, md5 = unname(code_md5)),
    file.path(out_dir, "code_inputs.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))
print(summary, row.names = FALSE, digits = 4)
cat("\nFeature importance:\n"); print(round(sort(importance, decreasing = TRUE), 5))
cat("\nAnnotation-evidence check:\n")
print(evidence_summary, row.names = FALSE, digits = 4)
