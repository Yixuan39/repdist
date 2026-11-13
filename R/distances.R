# Relative abundances `P` over the proteins seen in any sample, and their
# unit-norm embeddings `Z` in the same order.
.sample_inputs <- function(counts, x, weighted) {
    counts <- as.matrix(counts)
    Z <- .embeddings(x)
    .check_ids(rownames(counts), "Sample")
    .check_ids(colnames(counts), "Protein")
    stopifnot(
        is.numeric(counts), all(is.finite(counts)), all(counts >= 0),
        all(rowSums(counts) > 0), all(colnames(counts) %in% rownames(Z)),
        is.logical(weighted), length(weighted) == 1L, !is.na(weighted)
    )
    if (!weighted) counts <- sign(counts)
    counts <- counts[, colSums(counts) > 0, drop = FALSE]
    list(P = counts / rowSums(counts), Z = Z[colnames(counts), , drop = FALSE])
}

#' RBF-MMD distance between samples
#'
#' Treats each sample as an abundance-weighted distribution over protein
#' embeddings and returns the maximum mean discrepancy between those
#' distributions under a Gaussian kernel on Euclidean distance, which is
#' positive semidefinite. On the unit-norm embeddings the kernel is
#' `exp(-(1 - cosine similarity) / sigma^2)`; for TM-Vec, cosine similarity is
#' the predicted TM-score.
#'
#' Proteins with zero abundance in every sample are dropped first. Row
#' normalisation removes total-count scaling, but not the effect of depth on
#' which proteins were seen.
#'
#' @section Estimator:
#' This is the exact MMD between the observed relative-abundance
#' distributions, including within-sample diagonal terms. As an estimate
#' of a population MMD it is the biased plug-in form. The unbiased
#' U-statistic estimates squared population MMD and can be negative, so it
#' is not used as a sample distance for ordination or PERMANOVA.
#'
#' @section Memory:
#' Zero-abundance proteins are removed across all samples in `counts`,
#' not separately for each sample pair. The retained protein kernel is a
#' dense matrix, so memory grows quadratically with the number of proteins
#' detected in at least one supplied sample.
#'
#' @section Bandwidth:
#' `sigma = NULL` uses the median pairwise Euclidean distance between the
#' retained proteins. It therefore moves when a new sample brings in new
#' proteins; the value used is returned as `attr(, "sigma")`, so pass it back
#' to hold it fixed. For TM-Vec embeddings, `sigma = sqrt(1 - t)` puts the
#' kernel's `1/e` point at predicted TM-score `t`.
#'
#' @param counts Numeric matrix of counts or relative abundances, or a sparse
#'   `Matrix`, with samples in rows and proteins in columns. Every row must
#'   have a positive total. Sample and protein names must be unique,
#'   non-missing and non-empty.
#' @param x Proteins with embeddings, from
#'   [embed_proteins()] or
#'   [read_embeddings()],
#'   covering `colnames(counts)`.
#' @param weighted Use abundance as the weight (default). `FALSE` gives every
#'   detected protein equal weight, the presence/absence counterpart.
#' @param sigma RBF bandwidth; see the Bandwidth section.
#' @return A `dist` object for `vegan::adonis2()`, with the bandwidth as
#'   `attr(, "sigma")`.
#' @references Gretton et al. (2012). A Kernel Two-Sample Test. JMLR 13,
#'   723-773.
#' @examples
#' x <- Biostrings::AAStringSet(c(p1 = "MKV", p2 = "MAL", p3 = "MKL"))
#' S4Vectors::mcols(x)$embedding <- rbind(c(1, 0), c(0, 1), c(1, 1))
#' counts <- rbind(s1 = c(p1 = 8, p2 = 1, p3 = 1), s2 = c(1, 8, 1))
#' D <- sample_mmd(counts, x)
#' attr(D, "sigma")
#' @export
sample_mmd <- function(counts, x, weighted = TRUE, sigma = NULL) {
    input <- .sample_inputs(counts, x, weighted)
    d <- stats::dist(input$Z)
    if (is.null(sigma)) {
        # One protein: the kernel is 1 whatever the bandwidth.
        sigma <- if (length(d)) stats::median(d) else 1
        if (sigma == 0) {
            stop("The median protein distance is zero. Check for identical ",
                "embeddings, or pass a positive `sigma`.", call. = FALSE)
        }
    }
    if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) ||
        sigma <= 0) {
        stop("`sigma` must be a single finite, positive number.", call. = FALSE)
    }
    K <- exp(-0.5 * (as.matrix(d) / sigma)^2)
    G <- tcrossprod(input$P %*% K, input$P)
    g <- diag(G)
    # The kernel is PSD, so a negative squared distance is only rounding.
    D <- sqrt(pmax(outer(g, g, "+") - 2 * G, 0))
    structure(stats::as.dist(D), sigma = sigma)
}

utils::globalVariables(c("PCo1", "PCo2", "group"))

#' PCoA of samples on their most variable protein regions
#'
#' Principal coordinates of [sample_mmd()] on the proteins whose structural
#' region varies most across samples, as DESeq2's `plotPCA()` plots the most
#' variable genes. A protein's region is every protein within predicted
#' TM-score `min_sim` of it; its abundance in a sample is the summed relative
#' abundance of those proteins. Proteins are ranked by the coefficient of
#' variation of their region's abundance across samples, and the `top` are
#' kept. Regions present in every sample at a steady share, such as
#' translation and glycolysis, rank last. Ranking by variance instead would
#' keep them: being the most abundant, they also vary the most.
#'
#' The selection ignores sample labels and is for display only; test on the
#' full distance from [sample_mmd()].
#'
#' @inheritParams sample_mmd
#' @param group Optional sample labels to colour by, in the row order of
#'   `counts` or named by sample. Each group is outlined by its convex hull.
#' @param top Proteins to keep: a proportion in `(0, 1]` (default 0.9) or a
#'   whole number above 1.
#' @param min_sim Predicted TM-score that defines a region (default 0.9).
#' @param sigma RBF bandwidth. `NULL` uses the median distance over all
#'   proteins, the value `sample_mmd(counts, x)` uses, so the kernel is the
#'   one of the full distance.
#' @return A ggplot object whose `data` holds the coordinates. Needs the
#'   `ggplot2` package.
#' @examples
#' data(simulate_ec, package = "repdist")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_pcoa(simulate_ec$counts, simulate_ec$proteins,
#'         group = simulate_ec$samples$group)
#' }
#' @export
plot_pcoa <- function(counts, x, group = NULL, top = 0.9, min_sim = 0.9,
    weighted = TRUE, sigma = NULL) {
    input <- .sample_inputs(counts, x, weighted)
    P <- input$P
    m <- ncol(P)
    if (!is.numeric(top) || length(top) != 1L || !is.finite(top) ||
        top <= 0 || (top > 1 && (top != round(top) || top > m))) {
        stop("`top` must be a proportion in (0, 1] or a whole number of ",
            "proteins up to ", m, ".", call. = FALSE)
    }
    if (!is.numeric(min_sim) || length(min_sim) != 1L || is.na(min_sim) ||
        min_sim < 0 || min_sim > 1) {
        stop("`min_sim` must be a single number in [0, 1].", call. = FALSE)
    }
    if (nrow(P) < 3L) stop("A PCoA needs at least 3 samples.", call. = FALSE)
    if (!is.null(group)) {
        if (!is.null(names(group))) group <- group[rownames(P)]
        if (length(group) != nrow(P) || anyNA(group)) {
            stop("`group` must give one label per sample.", call. = FALSE)
        }
    }
    cs <- tcrossprod(input$Z) # cosine; for TM-Vec, predicted TM-score
    if (is.null(sigma)) {
        # sample_mmd()'s default: on unit-norm rows, distance^2 = 2 (1 - cosine)
        sigma <- if (m > 1L) {
            stats::median(sqrt(pmax(2 - 2 * cs[upper.tri(cs)], 0)))
        } else 1
    }
    region <- P %*% (cs >= min_sim)
    rm(cs)
    cv <- apply(region, 2, stats::sd) / colMeans(region)
    k <- if (top <= 1) ceiling(top * m) else top
    keep <- rank(-cv, ties.method = "first") <= k
    left <- rowSums(P[, keep, drop = FALSE])
    if (any(left == 0)) {
        stop("No protein kept for sample(s) ",
            paste(rownames(P)[left == 0], collapse = ", "), "; raise `top`.",
            call. = FALSE)
    }
    e <- stats::cmdscale(sample_mmd(P[, keep, drop = FALSE], x, weighted,
        sigma), k = 2, eig = TRUE)
    pct <- 100 * e$eig[seq_len(2L)] / sum(e$eig[e$eig > 0])
    df <- data.frame(sample = rownames(P), PCo1 = e$points[, 1],
        PCo2 = e$points[, 2])
    if (!is.null(group)) df$group <- factor(group)
    p <- ggplot2::ggplot(df, ggplot2::aes(PCo1, PCo2)) +
        ggplot2::labs(x = sprintf("PCo1 (%.0f%%)", pct[1]),
            y = sprintf("PCo2 (%.0f%%)", pct[2]),
            subtitle = sprintf(paste("MMD on the %d most variable of %d",
                "proteins (%.0f%% of abundance)"), k, m, 100 * mean(left))) +
        ggplot2::theme_bw()
    if (is.null(group)) return(p + ggplot2::geom_point())
    hull <- do.call(rbind, lapply(split(df, df$group),
        function(g) g[grDevices::chull(g$PCo1, g$PCo2), ]))
    p + ggplot2::geom_polygon(data = hull,
            ggplot2::aes(colour = group, group = group), fill = NA) +
        ggplot2::geom_point(ggplot2::aes(colour = group)) +
        ggplot2::labs(colour = NULL)
}
