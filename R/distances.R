#' Protein distance matrix from embeddings
#'
#' Builds the protein-by-protein `dist` that [bin_proteins()] and the QC
#' plots consume. Compute it once and reuse it: at real proteome
#' sizes it is the expensive object (n(n-1)/2 doubles, roughly 400 MB at
#' 10,000 proteins). Capped at 46341 proteins, the largest set this package
#' can index; dereplicate the sequences before embedding if you have more.
#'
#' Cosine distance is `1 - cosine similarity` -- for TM-Vec embeddings, `1 -
#' predicted TM-score`, the scale [bin_proteins()]'s `min_sim` is read
#' against. [sample_repdist()] uses a Gaussian kernel on Euclidean distance
#' and builds its own ground distance from
#' the embeddings, in blocks, and never reads this one. On the unit-norm rows
#' [embed_proteins()] returns, the two are monotone transforms of each other
#' (`||a - b|| = sqrt(2 * (1 - cos))`) and rank protein pairs identically.
#'
#' @param embeddings Numeric matrix, one L2-normalised row per protein, as
#'   returned by [embed_proteins()]. Computes `1 - inner product` in blocks
#'   without renormalising rows or allocating a full square matrix.
#' @param block_size Positive integer, the number of protein columns per
#'   inner-product block. Smaller blocks reduce temporary memory use.
#' @return A `dist` object over proteins.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' Z <- Z / sqrt(rowSums(Z^2))
#' repdist_matrix(Z)
#' @export
repdist_matrix <- function(embeddings, block_size = 1000L) {
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1 || block_size != floor(block_size))
    stop("`block_size` must be a single positive integer.", call. = FALSE)
  Z <- as.matrix(embeddings)
  stopifnot(is.numeric(Z), all(is.finite(Z)))
  .repdist_check_n(nrow(Z), "proteins")
  # Allow float32 rounding from the embedding backend.
  if (any(abs(rowSums(Z^2) - 1) > 1e-6))
    stop("`embeddings` must have L2-normalised rows (unit norm).",
         call. = FALSE)
  n <- nrow(Z)
  d <- numeric(n * (n - 1) / 2)
  columns <- seq_len(max(0L, n - 1L))
  for (j in split(columns, ceiling(columns / block_size))) {
    S <- tcrossprod(Z, Z[j, , drop = FALSE])
    for (k in seq_along(j)) {
      col <- j[[k]]
      i <- seq.int(col + 1L, n)
      # Lower triangle in stats::dist order; clip floating-point overshoot.
      d[n * (col - 1) - col * (col - 1) / 2 + i - col] <-
        pmax(0, pmin(2, 1 - S[i, k]))
    }
  }
  # Wrap the packed lower triangle as a standard dist; no square allocation.
  structure(d, Size = n, Labels = rownames(Z), Diag = FALSE, Upper = FALSE,
            method = "cosine", class = "dist")
}

# sample_repdist()'s Gaussian kernel still uses Euclidean distances, through
# ||a - b||^2 = |a|^2 + |b|^2 - 2 a.b. `sq` is rowSums(Z^2), computed once.
.repdist_dist_block <- function(Z, sq, j, i) {
  D2 <- outer(sq[i], sq[j], "+") -
    2 * tcrossprod(Z[i, , drop = FALSE], Z[j, , drop = FALSE])
  sqrt(pmax(D2, 0))   # rounding can push a coincident pair a hair below zero
}

#' RBF-MMD distance between samples
#'
#' Treats each sample as an abundance-weighted distribution over protein
#' embeddings and returns the maximum mean discrepancy between those
#' distributions under a median-heuristic RBF kernel. The ground metric is
#' always Euclidean, whose Gaussian RBF is positive semidefinite. On unit-norm
#' embeddings it equals `exp(-(1 - cosine similarity) / sigma^2)`; for TM-Vec,
#' cosine similarity is the predicted TM-score. Squaring cosine distance in
#' the Gaussian formula instead would give a different kernel without this
#' guarantee.
#'
#' Proteins with zero abundance in every sample are dropped before the kernel
#' or bandwidth is computed. Rarefy `counts` first (e.g. `vegan::rrarefy()`)
#' if samples weren't collected at a common depth.
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
#' proteins. The bandwidth actually used is returned as `attr(, "sigma")`;
#' pass it back as `sigma` to hold it fixed across runs or studies.
#'
#' @section Memory:
#' The full protein kernel is never stored: `block_size`-by-`block_size` tiles
#' are accumulated into the samples-by-samples Gram matrix. Sparse `Matrix`
#' counts stay sparse. Automatic bandwidth selection holds pairwise distances
#' for at most `sigma_max_proteins` proteins, using quadratic memory in that
#' limit; lowering it reduces memory use at the cost of a smaller subsample.
#'
#' @param counts Numeric matrix of counts or relative abundances, including a
#'   sparse `Matrix`, with samples in rows and proteins in columns. Every row
#'   must have a positive total.
#' @param embeddings Numeric matrix, one row per protein with rownames matching
#'   `colnames(counts)`.
#' @param weighted Use abundance as the weight (default). `FALSE` replaces each
#'   sample's abundances with `sign(abundance)`, so every expressed protein
#'   carries equal weight -- the presence/absence counterpart, as unweighted
#'   UniFrac stands to weighted UniFrac.
#' @param sigma RBF bandwidth. Defaults to the median heuristic; see the
#'   Bandwidth section.
#' @param block_size Positive integer, the maximum number of proteins per side
#'   of a kernel block. Smaller blocks reduce temporary memory use.
#' @param sigma_max_proteins Integer at least 2, the maximum number of proteins
#'   used to estimate the median bandwidth. Only used when `sigma = NULL`.
#' @return A `dist` object, for direct use with `vegan::adonis2()`, with the
#'   bandwidth used as `attr(, "sigma")`.
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
sample_repdist <- function(counts,
                           embeddings,
                           weighted = TRUE,
                           sigma = NULL,
                           block_size = 1000L,
                           sigma_max_proteins = 2000L) {
  if (!is.numeric(block_size) || length(block_size) != 1L ||
      !is.finite(block_size) || block_size < 1 || block_size != floor(block_size))
    stop("`block_size` must be a single positive integer.", call. = FALSE)
  if (!is.numeric(sigma_max_proteins) || length(sigma_max_proteins) != 1L ||
      !is.finite(sigma_max_proteins) || sigma_max_proteins < 2 ||
      sigma_max_proteins != floor(sigma_max_proteins))
    stop("`sigma_max_proteins` must be a single integer >= 2.", call. = FALSE)
  sparse_counts <- inherits(counts, "sparseMatrix")
  if (sparse_counts && !inherits(counts, "dsparseMatrix"))
    stop("`counts` must be numeric.", call. = FALSE)
  if (!sparse_counts) counts <- as.matrix(counts)
  count_values <- if (sparse_counts) counts@x else counts
  if (inherits(embeddings, "dist"))
    stop("`embeddings` must be the embedding matrix, not a `dist`: the ground ",
         "distance is computed in blocks and never held whole.", call. = FALSE)
  Z <- as.matrix(embeddings)

  stopifnot(
    !is.null(rownames(counts)), !is.null(colnames(counts)),
    !is.null(rownames(Z)), !anyDuplicated(rownames(Z)),
    all(colnames(counts) %in% rownames(Z)), !anyDuplicated(colnames(counts)),
    is.numeric(count_values), all(is.finite(count_values)),
    all(count_values >= 0),
    all(Matrix::rowSums(counts) > 0),
    is.logical(weighted), length(weighted) == 1L, !is.na(weighted)
  )

  if (!weighted) counts <- sign(counts)

  # Remove proteins absent from every sample before estimating the bandwidth.
  counts <- counts[, Matrix::colSums(counts) > 0, drop = FALSE]
  totals <- Matrix::rowSums(counts)
  # Each row of P is a probability distribution; diagonal scaling keeps
  # sparse counts sparse. With weighted = FALSE this is uniform over support.
  P <- if (sparse_counts) Matrix::Diagonal(x = 1 / totals) %*% counts else
    counts / totals
  # Align embedding rows with P's protein columns before matrix multiplication.
  Z <- Z[colnames(counts), , drop = FALSE]
  stopifnot(is.numeric(Z), all(is.finite(Z)))

  n <- nrow(Z)
  sq <- rowSums(Z^2)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / block_size))
  if (is.null(sigma)) {
    # ponytail: deterministic subsample capped at sigma_max_proteins; pass sigma
    # when the bandwidth must be fixed across catalogs or protein orderings.
    idx <- round(seq(1, n, length.out = min(n, sigma_max_proteins)))
    # One protein: the kernel is the 1 x 1 matrix 1 whatever the bandwidth.
    sigma <- if (n > 1L) stats::median(stats::dist(Z[idx, , drop = FALSE])) else 1
    # Conservative squared-distance cutoff: the Gram identity subtracts large
    # squared norms, losing relative precision when separations are tiny.
    distance_tol <- sqrt(.Machine$double.eps) * max(sq)
    if (n > 1L && sigma^2 <= distance_tol)
      stop("The median protein distance is zero or too small relative to ",
           "embedding norms for stable kernel computation. Check for duplicate ",
           "or near-identical embeddings, or pass a larger `sigma`.", call. = FALSE)
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("`sigma` must be a single finite, positive number.", call. = FALSE)

  # Accumulate sample inner products G = P K P'; K is never held in full.
  # mmd_matrix() then converts G to distances, so these are separate steps.
  G <- matrix(0, nrow(P), nrow(P))
  for (u in seq_along(blocks)) {
    i <- blocks[[u]]
    Pi <- P[, i, drop = FALSE]
    for (v in seq_len(u)) {
      j <- blocks[[v]]
      Gij <- as.matrix(Matrix::tcrossprod(
        Pi %*% rbf_kernel(.repdist_dist_block(Z, sq, j, i), sigma),
        P[, j, drop = FALSE]))
      # Off-diagonal kernel tiles have a transposed partner; add both once.
      G <- G + if (u == v) Gij else Gij + t(Gij)
    }
  }
  D <- mmd_matrix(G)
  dimnames(D) <- list(rownames(counts), rownames(counts))
  structure(stats::as.dist(D), sigma = sigma)
}

#' RBF kernel from ground-metric distances
#'
#' @param C Distance matrix; [sample_repdist()] passes protein-block tiles.
#' @param sigma Bandwidth, a single positive number.
#' @return `exp(-C^2 / (2 sigma^2))`, the same shape as `C`.
#' @noRd
rbf_kernel <- function(C, sigma) exp(-(C^2) / (2 * sigma^2))

#' Abundance-weighted MMD distance matrix
#'
#' @param G Samples-by-samples Gram matrix `P %*% K %*% t(P)`.
#' @return A samples x samples matrix of MMD distances.
#' @noRd
mmd_matrix <- function(G) {
  self_term <- diag(G)
  M2 <- outer(self_term, self_term, "+") - 2 * G
  # Tolerance scaled to the kernel, not to M2: when samples share a composition
  # M2 is exactly 0, and a tolerance relative to max(abs(M2)) collapses onto
  # the rounding noise it should absorb.
  tol <- 1e-8 * max(self_term)
  # Detect invalid squared distances; this is not a full PSD test of G or K.
  if (any(M2 < -tol))
    stop("Negative squared MMD beyond floating-point tolerance; ",
         "check the kernel and numerical stability.", call. = FALSE)
  diag(M2) <- 0
  sqrt(pmax(M2, 0))
}
