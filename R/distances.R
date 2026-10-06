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
