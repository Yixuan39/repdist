# `block_size` and `sigma_max_proteins` are single positive integers; the
# second needs at least a pair of proteins to take a median distance over.
.repdist_check_count_arg <- function(x, what, min = 1) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < min ||
        x != floor(x)) {
        stop(what, " must be a single integer >= ", min, ".", call. = FALSE)
    }
}

#' Protein distance matrix from embeddings
#'
#' Builds the protein-by-protein `dist` that
#' [bin_proteins()] and the QC plots consume.
#' Compute it once and reuse it: at real proteome sizes it is the expensive
#' object (n(n-1)/2 doubles, roughly 400 MB at 10,000 proteins). Capped at
#' 46341 proteins, the largest set this package can index; dereplicate the
#' sequences before embedding if you have more.
#'
#' Cosine distance is `1 - cosine similarity` -- for TM-Vec embeddings,
#' `1 - predicted TM-score`, the scale
#' [bin_proteins()]'s `min_sim` is read
#' against. [sample_repdist()] uses a Gaussian
#' kernel on Euclidean distance and builds its own ground distance from the
#' embeddings, in blocks, and never reads this one. On the unit-norm rows
#' [embed_proteins()] returns, the two are
#' monotone transforms of each other (`||a - b|| = sqrt(2 * (1 - cos))`)
#' and rank protein pairs identically.
#'
#' @param embeddings Numeric matrix, one L2-normalised row per protein, as
#'   returned by [embed_proteins()]. Computes
#'   `1 - inner product` in blocks without renormalising rows or
#'   allocating a full square matrix.
#' @param block_size Positive integer, the number of protein columns per
#'   inner-product block. Smaller blocks reduce temporary memory use.
#' @return A `dist` object over proteins.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' Z <- Z / sqrt(rowSums(Z^2))
#' repdist_matrix(Z)
#' @export
repdist_matrix <- function(embeddings, block_size = 1000L) {
    .repdist_check_count_arg(block_size, "`block_size`")
    Z <- as.matrix(embeddings)
    stopifnot(is.numeric(Z), all(is.finite(Z)))
    if (!is.null(rownames(Z))) .repdist_check_ids(rownames(Z), "`embeddings`")
    .repdist_check_n(nrow(Z), "proteins")
    # Allow float32 rounding from the embedding backend.
    if (any(abs(rowSums(Z^2) - 1) > 1e-6)) {
        stop("`embeddings` must have L2-normalised rows (unit norm).",
            call. = FALSE
        )
    }
    n <- nrow(Z)
    d <- numeric(n * (n - 1) / 2)
    columns <- seq_len(max(0L, n - 1L))
    for (j in split(columns, ceiling(columns / block_size))) {
        S <- tcrossprod(Z, Z[j, , drop = FALSE])
        for (k in seq_along(j)) {
            col <- j[[k]]
            i <- seq.int(col + 1L, n)
            # Lower triangle in stats::dist order; clip fp overshoot.
            d[n * (col - 1) - col * (col - 1) / 2 + i - col] <-
                pmax(0, pmin(2, 1 - S[i, k]))
        }
    }
    # Wrap the packed lower triangle as a standard dist; no square allocation.
    structure(d,
        Size = n, Labels = rownames(Z), Diag = FALSE, Upper = FALSE,
        method = "cosine", class = "dist"
    )
}

# Coerce and validate sample_repdist()'s two matrices together: a sparse
# `counts` stays sparse, and `embeddings` must be the matrix, never a `dist`.
.repdist_sample_inputs <- function(counts, embeddings, weighted) {
    sparse <- inherits(counts, "sparseMatrix")
    if (sparse && !inherits(counts, "dsparseMatrix")) {
        stop("`counts` must be numeric.", call. = FALSE)
    }
    if (!sparse) counts <- as.matrix(counts)
    if (inherits(embeddings, "dist")) {
        stop("`embeddings` must be the embedding matrix, not a `dist`: the ",
            "ground distance is computed in blocks and never held whole.",
            call. = FALSE
        )
    }
    Z <- as.matrix(embeddings)
    .repdist_check_sample_inputs(
        counts, Z, if (sparse) counts@x else counts, weighted
    )
    list(counts = counts, Z = Z, sparse = sparse)
}

.repdist_check_sample_inputs <- function(counts, Z, count_values, weighted) {
    .repdist_check_ids(rownames(counts), "Sample")
    .repdist_check_ids(colnames(counts), "Protein")
    .repdist_check_ids(rownames(Z), "`embeddings`")
    stopifnot(
        nrow(counts) > 0L, ncol(Z) > 0L,
        all(colnames(counts) %in% rownames(Z)),
        is.numeric(count_values), all(is.finite(count_values)),
        all(count_values >= 0),
        all(Matrix::rowSums(counts) > 0),
        is.logical(weighted), length(weighted) == 1L, !is.na(weighted)
    )
}

# Median pairwise distance over a deterministic subsample.
# ponytail: capped at sigma_max_proteins; pass `sigma` when the bandwidth must
# be fixed across catalogs or protein orderings.
.repdist_median_sigma <- function(Z, sigma_max_proteins) {
    n <- nrow(Z)
    # One protein: the kernel is the 1 x 1 matrix 1 whatever the bandwidth.
    if (n < 2L) return(1)
    idx <- round(seq(1, n, length.out = min(n, sigma_max_proteins)))
    sigma <- stats::median(stats::dist(Z[idx, , drop = FALSE]))
    if (sigma == 0) {
        stop("The median protein distance is zero. Check for duplicate or ",
            "identical embeddings, or pass a positive `sigma`.",
            call. = FALSE
        )
    }
    sigma
}

# Sample inner products G = P K P', one kernel tile at a time; K is never held
# in full. mmd_matrix() then converts G to distances, so these stay separate.
.repdist_gram <- function(P, Z, sq, blocks, sigma) {
    G <- matrix(0, nrow(P), nrow(P))
    for (u in seq_along(blocks)) {
        i <- blocks[[u]]
        Pi <- P[, i, drop = FALSE]
        for (v in seq_len(u)) {
            j <- blocks[[v]]
            Gij <- as.matrix(Matrix::tcrossprod(
                Pi %*% rbf_kernel(.repdist_dist_block(Z, sq, j, i), sigma),
                P[, j, drop = FALSE]
            ))
            # Off-diagonal tiles have a transposed partner; add both once.
            G <- G + if (u == v) Gij else Gij + t(Gij)
        }
    }
    G
}

# sample_repdist()'s Gaussian kernel still uses Euclidean distances, through
# ||a - b||^2 = |a|^2 + |b|^2 - 2 a.b. `sq` is rowSums(Z^2), computed once.
.repdist_dist_block <- function(Z, sq, j, i) {
    scale <- outer(sq[i], sq[j], "+")
    D2 <- scale -
        2 * tcrossprod(Z[i, , drop = FALSE], Z[j, , drop = FALSE])
    # The Gram identity subtracts comparable squared norms, so a separation
    # below sqrt(eps) * scale is lost in the cancellation -- and clipping it to
    # zero erases a genuine distance when `sigma` is smaller still. Recompute
    # exactly there; that also leaves every D2 non-negative, so sqrt() needs
    # no clamp. sample_repdist() centres Z first, which keeps `scale` (and so
    # this tolerance) as small as the embeddings allow.
    for (k in seq_along(j)) {
        close <- which(D2[, k] <= sqrt(.Machine$double.eps) * scale[, k])
        if (length(close)) {
            delta <- sweep(Z[i[close], , drop = FALSE], 2, Z[j[k], ], "-")
            D2[close, k] <- rowSums(delta^2)
        }
    }
    sqrt(D2)
}

#' RBF-MMD distance between samples
#'
#' Treats each sample as an abundance-weighted distribution over protein
#' embeddings and returns the maximum mean discrepancy between those
#' distributions under a median-heuristic RBF kernel. The ground metric is
#' always Euclidean, whose Gaussian RBF is positive semidefinite. On unit-norm
#' embeddings it equals `exp(-(1 - cosine similarity) / sigma^2)`; for
#' TM-Vec, cosine similarity is the predicted TM-score. Squaring cosine
#' distance in the Gaussian formula instead would give a different kernel
#' without this guarantee.
#'
#' Proteins with zero abundance in every sample are dropped before the kernel
#' or bandwidth is computed. Row normalization removes total-count scaling,
#' but sampling effects on composition and richness still need assessment.
#' Rarefaction is one possible preprocessing choice, with information loss.
#'
#' @section Bandwidth:
#' With `sigma = NULL`, the bandwidth is chosen automatically as the median
#' pairwise Euclidean distance between retained proteins, without abundance
#' weighting or sample-group labels. At most `sigma_max_proteins` rows are
#' selected at evenly spaced positions (all rows when below that limit).
#' This median heuristic is a common default, not a guarantee of optimal test
#' power. Smaller bandwidths emphasise local differences; larger ones smooth
#' over them. For unit-norm TM-Vec embeddings, `sigma = sqrt(1 - t)` puts the
#' kernel's `1/e` point at predicted TM-score `t`.
#'
#' The default depends on which proteins are retained, so the same two samples
#' sit at a slightly different distance when a third sample brings in new
#' proteins. The bandwidth actually used is returned as
#' `attr(, "sigma")`; pass it back as `sigma` to hold it fixed
#' across runs or studies.
#'
#' @section Memory:
#' The full protein kernel is never stored: tiles of
#' `block_size`-by-`block_size` are accumulated
#' into the samples-by-samples Gram matrix. Sparse `Matrix` counts stay
#' sparse. Automatic bandwidth selection holds pairwise distances for at most
#' `sigma_max_proteins` proteins, using quadratic memory in that limit;
#' lowering it reduces memory use at the cost of a smaller subsample.
#'
#' @param counts Numeric matrix of counts or relative abundances, including a
#'   sparse `Matrix`, with samples in rows and proteins in columns. Every
#'   row must have a finite positive total. Sample and protein names must be
#'   unique, non-missing and non-empty.
#' @param embeddings Numeric matrix, one row per protein with rownames
#'   matching `colnames(counts)`.
#' @param weighted Use abundance as the weight (default). `FALSE`
#'   replaces each sample's abundances with `sign(abundance)`, so every
#'   expressed protein carries equal weight -- the presence/absence
#'   counterpart, as unweighted UniFrac stands to weighted UniFrac.
#' @param sigma RBF bandwidth. Defaults to the median heuristic; see the
#'   Bandwidth section.
#' @param block_size Positive integer, the maximum number of proteins per side
#'   of a kernel block. Smaller blocks reduce temporary memory use.
#' @param sigma_max_proteins Integer at least 2, the maximum number of
#'   proteins used to estimate the median bandwidth. Only used when
#'   `sigma = NULL`.
#' @return A `dist` object, for direct use with `vegan::adonis2()`,
#'   with the bandwidth used as `attr(, "sigma")`.
#' @usage sample_repdist(
#'     counts,
#'     embeddings,
#'     weighted = TRUE,
#'     sigma = NULL,
#'     block_size = 1000L,
#'     sigma_max_proteins = 2000L
#' )
#' @importClassesFrom Matrix sparseMatrix dsparseMatrix
#' @references Garreau, Jitkrittum and Kanagawa, *Large sample analysis of the
#'   median heuristic*. <https://arxiv.org/abs/1707.07269>
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' counts <- rbind(s1 = c(8, 1, 1), s2 = c(1, 8, 1))
#' colnames(counts) <- rownames(Z)
#' D <- sample_repdist(counts, Z)
#' attr(D, "sigma")
#' @export
sample_repdist <- function(
    counts, embeddings, weighted = TRUE, sigma = NULL,
    block_size = 1000L, sigma_max_proteins = 2000L) {
    .repdist_check_count_arg(block_size, "`block_size`")
    .repdist_check_count_arg(sigma_max_proteins, "`sigma_max_proteins`", 2)
    input <- .repdist_sample_inputs(counts, embeddings, weighted)
    counts <- input$counts
    Z <- input$Z

    if (!weighted) counts <- sign(counts)

    # Remove proteins absent from every sample before estimating the bandwidth.
    counts <- counts[, Matrix::colSums(counts) > 0, drop = FALSE]
    totals <- Matrix::rowSums(counts)
    # Each row of P is a probability distribution; diagonal scaling keeps
    # sparse counts sparse. With weighted = FALSE this is uniform over support.
    P <- if (input$sparse) {
        Matrix::Diagonal(x = 1 / totals) %*% counts
    } else {
        counts / totals
    }
    # Align embedding rows with P's protein columns before matrix
    # multiplication.
    Z <- Z[colnames(counts), , drop = FALSE]
    stopifnot(is.numeric(Z), all(is.finite(Z)))

    # A common translation leaves Euclidean distances unchanged and removes
    # any large shared offset before the Gram identity, which loses relative
    # precision against big squared norms.
    Z <- sweep(Z, 2, colMeans(Z), "-")

    n <- nrow(Z)
    sq <- rowSums(Z^2)
    blocks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
    if (is.null(sigma)) sigma <- .repdist_median_sigma(Z, sigma_max_proteins)
    if (!is.numeric(sigma) || length(sigma) != 1L ||
        !is.finite(sigma) || sigma <= 0) {
        stop("`sigma` must be a single finite, positive number.", call. = FALSE)
    }

    D <- mmd_matrix(.repdist_gram(P, Z, sq, blocks, sigma))
    dimnames(D) <- list(rownames(counts), rownames(counts))
    structure(stats::as.dist(D), sigma = sigma)
}

#' RBF kernel from ground-metric distances
#'
#' @param C Distance matrix; [sample_repdist()]
#'   passes protein-block tiles.
#' @param sigma Bandwidth, a single positive number.
#' @return `exp(-C^2 / (2 sigma^2))`, the same shape as `C`.
#' @noRd
rbf_kernel <- function(C, sigma) exp(-0.5 * (C / sigma)^2)

#' Abundance-weighted MMD distance matrix
#'
#' @param G Samples-by-samples Gram matrix `P %*% K %*% t(P)`.
#' @return A samples x samples matrix of MMD distances.
#' @noRd
mmd_matrix <- function(G) {
    self_term <- diag(G)
    M2 <- outer(self_term, self_term, "+") - 2 * G
    # Tolerance scaled to the kernel, not to M2: when samples share a
    # composition M2 is exactly 0, and a tolerance relative to max(abs(M2))
    # collapses onto the rounding noise it should absorb.
    tol <- 1e-8 * max(self_term)
    # Detect invalid squared distances; this is not a full PSD test of G or K.
    if (any(M2 < -tol)) {
        stop("Negative squared MMD beyond floating-point tolerance; ",
            "check the kernel and numerical stability.",
            call. = FALSE
        )
    }
    diag(M2) <- 0
    sqrt(pmax(M2, 0))
}
