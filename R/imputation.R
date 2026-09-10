#' Propagate protein labels and score their reliability
#'
#' KNN assigns every unannotated protein a label. A random forest estimates
#' agreement with original annotations, using grouped held-out predictions.
#' Scoring never changes the labels or rejects a prediction.
#'
#' @details
#' Within each calibration fold, annotated and unannotated queries use the
#' same references: original annotations outside that fold. Labels absent from
#' references count as errors. Final KNN predictions use all original labels.
#'
#' Scoring features are supporting similarity, its difference from the nearest
#' alternative-label similarity, log-transformed candidate-label support,
#' nearest-any-protein similarity (excluding self), and top-five label
#' agreement. The scorer is a `ranger` probability forest with 500 trees,
#' minimum node size 20, permutation importance and one thread. No
#' hyperparameters are selected on queries.
#'
#' A forest cannot predict beyond the feature range it was trained on, so
#' scores plateau rather than extrapolate above the calibration ceiling.
#' Grouped calibration caps supporting similarity below the grouping
#' threshold, so queries with closer references than any calibration record
#' receive a conservative plateau score and a `score_extrapolation` flag.
#' That flag checks marginal feature ranges only.
#'
#' Scores require three annotated groups, 20 held-out records, five correct
#' and five incorrect outcomes. Insufficient data or failed fits produce NA
#' confidence with a warning, while labels are still returned. Confidence is a
#' model estimate, not a guaranteed probability of biological correctness.
#' Grouped calibration does not validate full-reference predictions on actual
#' unknowns.
#'
#' @param D Named protein dist from [repdist_matrix()], or a named symmetric
#'   numeric matrix of cosine distances in `[0, 2]`, with zero diagonal.
#' @param labels Named character vector or factor of original labels. Names
#'   must uniquely match proteins in D. Missing IDs, NA and blank values are
#'   unannotated. Each value is one categorical label, not a multi-label set.
#' @param k Single positive integer neighbour count, default 1. All exact
#'   kth-distance ties are included. Vote ties favour the nearest supporting
#'   reference, then protein ID in bytewise order.
#' @param groups Optional named vector of non-missing group IDs for every
#'   protein in D. Groups cannot cross calibration folds. Otherwise groups are
#'   single-linkage components at groupSimilarity, including unannotated bridge
#'   proteins. For TM-Vec, the default corresponds to predicted TM95, not a
#'   sequence-identity or experimentally measured homology threshold.
#' @param groupSimilarity Similarity threshold in `[-1, 1]`, default 0.95,
#'   for automatic groups. Ignored when groups is supplied.
#' @param nFolds Integer number of calibration folds, at least 3, default 5.
#'   Reduced to the number of annotated groups when necessary.
#' @param seed Positive integer for deterministic folds and the forest,
#'   default 1. Does not modify R's random-number state.
#' @return A data frame, in D order, with one row per unannotated protein:
#'   protein, predicted_label, supporting_protein, confidence,
#'   supporting_similarity, n_support, n_neighbors, margin, label_count,
#'   local_similarity, agreement, has_alternative, and score_extrapolation.
#'   Missing alternatives have zero margin and has_alternative = FALSE;
#'   has_alternative is reported for interpreting margin and is not a
#'   scoring feature.
#'
#'   Attribute "models" contains two blocks: knn (settings and original
#'   reference labels) and score (fit and status).
#'   Attribute "calibration" contains fold assignments, held-out predictions
#'   and correctness.
#'   Attribute "diagnostics" records calibration sizes and source feature
#'   ranges. These are training diagnostics, not independent performance
#'   estimates. All-annotated input returns zero rows and no fitted models.
#' @examples
#' Z <- matrix(sin(seq_len(360)), nrow = 60)
#' Z <- Z / sqrt(rowSums(Z^2))
#' rownames(Z) <- paste0("p", seq_len(nrow(Z)))
#' labels <- setNames(rep(c("A", "B"), 25), rownames(Z)[seq_len(50)])
#' groups <- setNames(rownames(Z), rownames(Z))
#' result <- impute_labels(repdist_matrix(Z), labels, groups = groups)
#' head(result)
#' attr(result, "diagnostics")
#' @importFrom utils globalVariables
#' @export
impute_labels <- function(D, labels, k = 1L, groups = NULL,
                          groupSimilarity = 0.95, nFolds = 5L, seed = 1L) {
    D <- .imputeValidate(D, labels, k, groupSimilarity, nFolds, seed)
    ids <- attr(D, "Labels")
    annotation <- rep(NA_character_, length(ids))
    annotation[match(names(labels), ids)] <- as.character(labels)
    known <- which(!is.na(annotation) & nzchar(trimws(annotation)))
    if (!length(known))
        stop("Supply at least one annotated protein.", call. = FALSE)
    known <- known[order(ids[known], method = "radix")]
    query <- setdiff(seq_along(ids), known)
    if (!length(query))
        return(.imputeEmpty(character()))
    if (!requireNamespace("ranger", quietly = TRUE))
        stop("Label scoring requires the ranger package.", call. = FALSE)
    groups <- .imputeGroups(D, groups, groupSimilarity)
    fold <- .imputeFolds(groups, ids, known, nFolds, seed)
    local <- vapply(seq_along(ids), function(i) {
        1 - min(usedist::dist_get(D, i, setdiff(seq_along(ids), i)))
    }, numeric(1))
    out <- .imputeKnn(D, annotation, known, query, k, local)
    calibration <- .imputeCalibration(D, annotation, known, fold, k, local)
    calibration$group <- unname(groups)
    source <- calibration[calibration$annotated, , drop = FALSE]
    x <- .imputeFeatures(source)
    score <- .imputeFit(x, source$correct, seed)
    if (length(unique(fold[known])) < 3L)
        score <- list(status = "Fewer than three annotated groups.")
    if (score$status == "ok") {
        out$confidence <- .imputePredict(score, .imputeFeatures(out))
        ranges <- vapply(x, range, numeric(2))
        out$score_extrapolation <- .imputeExtrapolation(out, ranges)
        invalid <- !is.finite(out$confidence) |
            out$confidence < 0 | out$confidence > 1
        if (any(invalid)) {
            out$confidence[invalid] <- NA_real_
            warning("Confidence is unavailable for ", sum(invalid),
                    " non-estimable predictions; labels are retained.",
                    call. = FALSE)
        }
    } else {
        ranges <- NULL
        warning("Labels were propagated but confidence is unavailable: ",
                score$status, call. = FALSE)
    }
    attr(out, "models") <- list(
        knn = list(k = k, training_proteins = ids[known],
                   labels = stats::setNames(annotation[known], ids[known])),
        score = score)
    attr(out, "calibration") <- calibration
    attr(out, "diagnostics") <- list(
        nAnnotated = length(known), nGroups = length(unique(groups[known])),
        nFolds = length(unique(fold[known])),
        nCorrect = sum(source$correct == 1, na.rm = TRUE),
        nIncorrect = sum(source$correct == 0, na.rm = TRUE),
        nMissingScores = sum(is.na(out$confidence)), featureRanges = ranges)
    out
}

.imputeValidate <- function(D, labels, k, groupSimilarity, nFolds, seed) {
    if (!inherits(D, "dist")) {
        if (!is.matrix(D) || !is.numeric(D) || !isSymmetric(D) ||
            anyNA(D) || any(diag(D) != 0) || is.null(rownames(D)) ||
            !identical(rownames(D), colnames(D)))
            stop("D must be a named symmetric numeric distance matrix ",
                 "with matching dimnames and zero diagonal.", call. = FALSE)
        D <- stats::as.dist(D)
    }
    n <- attr(D, "Size")
    ids <- attr(D, "Labels")
    if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
        n != floor(n) || !is.numeric(D) || length(D) != n * (n - 1) / 2 ||
        !is.character(ids) || length(ids) != n || anyNA(ids) ||
        any(!nzchar(trimws(ids))) || anyDuplicated(ids))
        stop("D must be a valid dist with unique non-empty protein names.",
             call. = FALSE)
    .repdist_check_n(n, "proteins")
    if (identical(attr(D, "method"), "euclidean") ||
        any(!is.finite(D)) || any(D < 0 | D > 2))
        stop("D must contain finite cosine distances in [0, 2], ",
             "not euclidean distances.", call. = FALSE)
    if (!(is.character(labels) || is.factor(labels)) ||
        !is.null(dim(labels)) || is.null(names(labels)) ||
        anyNA(names(labels)) || any(!nzchar(trimws(names(labels)))) ||
        anyDuplicated(names(labels)) || any(!names(labels) %in% ids))
        stop("labels must have unique non-empty protein names in D.",
             call. = FALSE)
    for (name in c("k", "nFolds", "seed")) {
        value <- get(name)
        if (!is.numeric(value) || length(value) != 1L ||
            !is.finite(value) || value < 1 || value != floor(value) ||
            value > .Machine$integer.max)
            stop(name, " must be a single positive integer.", call. = FALSE)
    }
    if (nFolds < 3L)
        stop("nFolds must be at least 3.", call. = FALSE)
    if (!is.numeric(groupSimilarity) || length(groupSimilarity) != 1L ||
        !is.finite(groupSimilarity) || abs(groupSimilarity) > 1)
        stop("groupSimilarity must be in [-1, 1].", call. = FALSE)
    D
}

.imputeGroups <- function(D, groups, threshold) {
    ids <- attr(D, "Labels")
    if (is.null(groups)) {
        ## Sort before clustering so input order cannot change fold identity.
        sorted <- usedist::dist_subset(D, sort(ids, method = "radix"))
        tree <- fastcluster::hclust(sorted, method = "single")
        groups <- stats::cutree(tree, h = 1 - threshold)
    }
    if (!is.atomic(groups) || !is.null(dim(groups)) ||
        is.null(names(groups)) || anyNA(names(groups)) ||
        anyDuplicated(names(groups)) || !setequal(names(groups), ids) ||
        anyNA(groups) || any(!nzchar(trimws(as.character(groups)))))
        stop("groups must name every protein exactly once with a ",
             "non-missing group ID.", call. = FALSE)
    groups <- as.character(groups[ids])
    members <- split(ids, groups)
    representative <- vapply(members, function(x)
        sort(x, method = "radix")[[1L]], character(1))
    stats::setNames(unname(representative[groups]), ids)
}

.imputeFolds <- function(groups, ids, known, nFolds, seed) {
    members <- split(seq_along(ids), groups)
    count <- vapply(members, function(at) sum(at %in% known), integer(1))
    nFolds <- min(nFolds, sum(count > 0))
    priority <- order(-count,
        digest::digest2int(names(members), seed = as.integer(seed)),
        names(members), method = "radix")
    fold <- integer(length(ids))
    annotatedLoad <- totalLoad <- numeric(nFolds)
    for (i in priority) {
        f <- order(annotatedLoad, totalLoad, seq_len(nFolds))[[1L]]
        if (count[[i]] == 0L)
            f <- which.min(totalLoad)
        fold[members[[i]]] <- f
        annotatedLoad[[f]] <- annotatedLoad[[f]] + count[[i]]
        totalLoad[[f]] <- totalLoad[[f]] + length(members[[i]])
    }
    fold
}

.imputeEmpty <- function(ids) {
    n <- length(ids)
    data.frame(protein = ids, predicted_label = rep(NA_character_, n),
        supporting_protein = rep(NA_character_, n),
        confidence = rep(NA_real_, n), supporting_similarity = rep(NA_real_, n),
        n_support = integer(n), n_neighbors = integer(n),
        margin = rep(NA_real_, n), label_count = integer(n),
        local_similarity = rep(NA_real_, n), agreement = rep(NA_real_, n),
        has_alternative = rep(FALSE, n),
        score_extrapolation = rep(NA, n))
}

.imputeKnn <- function(D, annotation, train, query, k, local) {
    ids <- attr(D, "Labels")
    train <- train[order(ids[train], method = "radix")]
    out <- .imputeEmpty(ids[query])
    out$local_similarity <- local[query]
    if (!length(train))
        return(out)
    reference <- annotation[train]
    counts <- table(reference)
    ## ponytail: exact O(n_query * n_reference) scans; use indexed neighbours
    ## when reference catalogs outgrow the existing condensed distance input.
    for (i in seq_along(query)) {
        d <- usedist::dist_get(D, query[[i]], train)
        at <- order(d, ids[train], method = "radix")
        voters <- at[d[at] <= d[at[[min(k, length(at))]]]]
        candidates <- unique(reference[voters])
        votes <- tabulate(match(reference[voters], candidates))
        winner <- which.max(votes)
        label <- candidates[[winner]]
        support <- voters[match(label, reference[voters])]
        other <- reference != label
        out$predicted_label[[i]] <- label
        out$supporting_protein[[i]] <- ids[train[support]]
        out$supporting_similarity[[i]] <- 1 - d[[support]]
        out$n_support[[i]] <- votes[[winner]]
        out$n_neighbors[[i]] <- length(voters)
        out$has_alternative[[i]] <- any(other)
        out$margin[[i]] <- if (any(other)) min(d[other]) - d[[support]] else 0
        out$label_count[[i]] <- unname(counts[[label]])
        nearby <- at[d[at] <= d[at[[min(5L, length(at))]]]]
        out$agreement[[i]] <- mean(reference[nearby] == label)
    }
    out
}

.imputeCalibration <- function(D, annotation, known, fold, k, local) {
    calibration <- .imputeEmpty(attr(D, "Labels"))
    for (f in sort(unique(fold))) {
        query <- which(fold == f)
        train <- known[fold[known] != f]
        calibration[query, ] <- .imputeKnn(D, annotation, train, query,
                                          k, local)
    }
    calibration$annotated <- seq_along(annotation) %in% known
    calibration$correct <- NA_real_
    calibration$correct[known] <- as.numeric(
        calibration$predicted_label[known] == annotation[known])
    calibration$fold <- fold
    calibration
}

## has_alternative is constant on any multi-label catalog, so it is reported
## for interpreting margin but never fitted.
.imputeFeatures <- function(x) {
    data.frame(similarity = x$supporting_similarity, margin = x$margin,
        logSupport = log1p(x$label_count), localSimilarity = x$local_similarity,
        agreement = x$agreement, row.names = x$protein)
}

.imputeFit <- function(x, y, seed) {
    if (nrow(x) < 20L || anyNA(y) || any(!is.finite(as.matrix(x))) ||
        min(sum(y == 0), sum(y == 1)) < 5L)
        return(list(status = "Insufficient usable outcomes."))
    at <- order(rownames(x), method = "radix")
    x <- x[at, , drop = FALSE]
    y <- y[at]
    varying <- vapply(x, function(v) length(unique(v)) > 1L, logical(1))
    columns <- names(x)[varying]
    if (!length(columns))
        return(list(status = "No varying scoring features."))
    data <- x[, columns, drop = FALSE]
    data$response <- factor(y, levels = c(0, 1))
    fit <- tryCatch(withr::with_preserve_seed(ranger::ranger(
        dependent.variable.name = "response", data = data, probability = TRUE,
        num.trees = 500L, min.node.size = 20L, num.threads = 1L,
        importance = "permutation", seed = as.integer(seed))),
        error = function(e) e)
    if (inherits(fit, "error"))
        return(list(status = conditionMessage(fit)))
    list(status = "ok", fit = fit, columns = columns)
}

.imputePredict <- function(model, x) {
    if (model$status != "ok")
        return(rep(NA_real_, nrow(x)))
    tryCatch(withr::with_preserve_seed(as.numeric(stats::predict(
        model$fit, data = x, num.threads = 1L)$predictions[, "1"])),
        error = function(e) rep(NA_real_, nrow(x)))
}

.imputeExtrapolation <- function(out, ranges) {
    x <- as.matrix(.imputeFeatures(out))
    rowSums(sweep(x, 2L, ranges[1L, ], "<") |
            sweep(x, 2L, ranges[2L, ], ">")) > 0L
}
