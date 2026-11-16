# Relative abundances `P` over the proteins seen in any sample, and their
# unit-norm embeddings `Z` in the same order.
.sample_inputs <- function(counts, x, weighted) {
    counts <- as.matrix(counts)
    Z <- .embeddings(x)
    .check_ids(rownames(counts), "Sample")
    .check_ids(colnames(counts), "Protein")
    stopifnot(
        is.numeric(counts), all(is.finite(counts)), all(counts >= 0),
        is.logical(weighted), length(weighted) == 1L, !is.na(weighted)
    )
    missing <- setdiff(colnames(counts), rownames(Z))
    if (length(missing)) {
        stop(length(missing), " of ", ncol(counts), " protein(s) in `counts` ",
            "have no embedding in `x`, e.g. ", missing[[1]], ". `counts` ",
            "needs samples in rows and proteins in columns.", call. = FALSE)
    }
    if (any(empty <- rowSums(counts) == 0)) {
        stop("Sample(s) with no counts: ",
            paste(rownames(counts)[empty], collapse = ", "), ".", call. = FALSE)
    }
    if (!weighted) counts <- sign(counts)
    counts <- counts[, colSums(counts) > 0, drop = FALSE]
    list(P = counts / rowSums(counts), Z = Z[colnames(counts), , drop = FALSE])
}

# sample_mmd()'s default bandwidth: the median of the pairwise Euclidean
# distances `d` between proteins. One protein: the kernel is 1 whatever it is.
.median_sigma <- function(d) {
    sigma <- if (length(d)) stats::median(d) else 1
    if (sigma == 0) {
        stop("The median protein distance is zero. Check for identical ",
            "embeddings, or pass a positive `sigma`.", call. = FALSE)
    }
    sigma
}

#' MMD distance between samples
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
#'   [embed_proteins()], covering
#'   `colnames(counts)`.
#' @param weighted Use abundance as the weight (default). `FALSE` gives every
#'   detected protein equal weight, the presence/absence counterpart.
#' @param sigma Gaussian kernel bandwidth; see the Bandwidth section.
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
    if (is.null(sigma)) sigma <- .median_sigma(d)
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

#' PCoA of samples without their most abundant protein regions
#'
#' Principal coordinates of [sample_mmd()] after
#' dropping the `drop` most abundant protein regions. A region is a
#' connected group of proteins within predicted TM-score `min_sim` of each
#' other, so its kernels nearly coincide. The most abundant regions are
#' usually housekeeping proteins such as elongation factor Tu and GAPDH,
#' present in every sample; much of their share varies with things other
#' than the design, and that variation can crowd the plot.
#'
#' The selection ignores sample labels and is for display only; test on the
#' full distance from [sample_mmd()].
#'
#' @inheritParams sample_mmd
#' @param group Optional sample labels to colour by, in the row order of
#'   `counts` or named by sample. Each group is outlined by its convex hull.
#' @param drop Number of most abundant regions to drop (default 2). `0`
#'   plots the full distance.
#' @param min_sim Predicted TM-score that links proteins into a region
#'   (default 0.9).
#' @param sigma Gaussian kernel bandwidth. `NULL` uses the median distance
#'   over all proteins, the value `sample_mmd(counts, x)` uses, so the kernel
#'   is the one of the full distance.
#' @return A ggplot object whose `data` holds the coordinates. A message
#'   names the most abundant protein of each dropped region and the region's
#'   size. Needs the `ggplot2` package.
#' @examples
#' data(simulate_cog, package = "repdist")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_pcoa(simulate_cog$counts, simulate_cog$proteins,
#'         group = simulate_cog$samples$group)
#' }
#' @export
plot_pcoa <- function(counts, x, group = NULL, drop = 2, min_sim = 0.9,
    weighted = TRUE, sigma = NULL) {
    input <- .sample_inputs(counts, x, weighted)
    P <- input$P
    m <- ncol(P)
    if (!is.numeric(drop) || length(drop) != 1L || !is.finite(drop) ||
        drop < 0 || drop != round(drop)) {
        stop("`drop` must be a whole number of regions.", call. = FALSE)
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
        # on unit-norm rows, distance^2 = 2 (1 - cosine)
        sigma <- .median_sigma(sqrt(pmax(2 - 2 * cs[upper.tri(cs)], 0)))
    }
    ij <- which(cs >= min_sim, arr.ind = TRUE)
    rm(cs)
    ij <- ij[ij[, 1] < ij[, 2], , drop = FALSE]
    region <- .components(ij[, 1], ij[, 2], m)
    mass <- colSums(P)
    n_region <- max(region)
    if (drop >= n_region) {
        stop("`drop` must be below the ", n_region, " regions.", call. = FALSE)
    }
    gone <- order(-tapply(mass, region, sum))[seq_len(drop)]
    keep <- !region %in% gone
    left <- rowSums(P[, keep, drop = FALSE])
    if (any(left == 0)) {
        stop("Sample(s) ", paste(rownames(P)[left == 0], collapse = ", "),
            " hold only dropped regions; lower `drop`.", call. = FALSE)
    }
    e <- stats::cmdscale(sample_mmd(P[, keep, drop = FALSE], x, weighted,
        sigma), k = 2, eig = TRUE)
    if (ncol(e$points) < 2L) {
        stop("The samples differ along fewer than two principal coordinates.",
            call. = FALSE)
    }
    pct <- 100 * e$eig[seq_len(2L)] / sum(e$eig[e$eig > 0])
    df <- data.frame(sample = rownames(P), PCo1 = e$points[, 1],
        PCo2 = e$points[, 2])
    if (!is.null(group)) df$group <- factor(group)
    if (drop) {
        # name each dropped region by its most abundant protein
        top <- vapply(gone, function(r) {
            j <- which(region == r)
            sprintf("%s (%d protein%s)", colnames(P)[j][which.max(mass[j])],
                length(j), if (length(j) > 1L) "s" else "")
        }, character(1))
        message(sprintf("Dropped %d region%s, %.0f%% of abundance: %s.", drop,
            if (drop > 1) "s" else "", 100 * (1 - mean(left)),
            paste(top, collapse = ", ")))
    }
    p <- ggplot2::ggplot(df, ggplot2::aes(PCo1, PCo2)) +
        ggplot2::labs(x = sprintf("PCo1 (%.0f%%)", pct[1]),
            y = sprintf("PCo2 (%.0f%%)", pct[2])) +
        ggplot2::theme_bw()
    if (is.null(group)) return(p + ggplot2::geom_point())
    hull <- do.call(rbind, lapply(split(df, df$group),
        function(g) g[grDevices::chull(g$PCo1, g$PCo2), ]))
    p + ggplot2::geom_polygon(data = hull,
            ggplot2::aes(colour = group, group = group), fill = NA) +
        ggplot2::geom_point(ggplot2::aes(colour = group)) +
        ggplot2::labs(colour = NULL)
}
