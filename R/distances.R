#' Protein distance matrix from embeddings
#'
#' Builds the protein-by-protein `dist` that [bin_proteins()] and the QC
#' plots consume. Compute it once and reuse it: at real proteome
#' sizes it is the expensive object (n(n-1)/2 doubles, roughly 400 MB at
#' 10,000 proteins). Capped at 46341 proteins, the largest set this package
#' can index; dereplicate the sequences before embedding if you have more.
#'
#' `"cosine"` is `1 - cosine similarity` -- for TM-Vec embeddings, `1 -
#' predicted TM-score`, the scale [bin_proteins()]'s `min_sim` is read
#' against. It is not a true metric, so it cannot feed an RBF kernel;
#' [sample_repdist()] therefore builds its own euclidean ground distance from
#' the embeddings, in blocks, and never reads this one. On the unit-norm rows
#' [embed_proteins()] returns, the two are monotone transforms of each other
#' (`||a - b|| = sqrt(2 * (1 - cos))`) and rank protein pairs identically.
#'
#' @param embeddings Numeric matrix, one row per protein.
#' @param distance `"cosine"` (default) or `"euclidean"`.
#' @return A `dist` object over proteins.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' repdist_matrix(Z)
#' @export
repdist_matrix <- function(embeddings, distance = c("cosine", "euclidean")) {
  distance <- match.arg(distance)
  Z <- as.matrix(embeddings)
  stopifnot(is.numeric(Z), all(is.finite(Z)))
  .repdist_check_n(nrow(Z), "proteins")
  if (distance == "euclidean") return(stats::dist(Z))

  norm <- sqrt(rowSums(Z^2))
  if (any(dead <- norm == 0))
    stop("Cosine distance is undefined for zero-length embedding rows: ",
         paste(utils::head(rownames(Z)[dead], 5), collapse = ", "),
         call. = FALSE)
  # 1 - cos = ||a - b||^2 / 2 once rows are unit-norm, so stats::dist yields
  # the condensed cosine distance without a square matrix or another package.
  d <- stats::dist(Z / norm)^2 / 2
  attr(d, "method") <- "cosine"
  d
}

# Euclidean distances between protein blocks `i` and `j`, through
# ||a - b||^2 = |a|^2 + |b|^2 - 2 a.b. `sq` is rowSums(Z^2), computed once.
.repdist_dist_block <- function(Z, sq, j, i = seq_len(nrow(Z))) {
  D2 <- outer(sq[i], sq[j], "+") -
    2 * tcrossprod(Z[i, , drop = FALSE], Z[j, , drop = FALSE])
  sqrt(pmax(D2, 0))   # rounding can push a coincident pair a hair below zero
}

# The median heuristic on an evenly spaced subsample. Visiting all n(n-1)/2
# pairs costs n^2 * 8 bytes (7 GB at 30,000 proteins) and was the only
# quadratic allocation left in sample_repdist(); .repdist_max_median_n rows
# bound it at ~2 million pairs, far more than a bandwidth is sensitive to.
# ponytail: a stride, not a random draw, so the default sigma is reproducible
# without a seed. Pass `sigma` if you need it fixed across protein sets.
.repdist_max_median_n <- 2000L

.repdist_median_dist <- function(Z, sq) {
  n <- nrow(Z)
  if (n > .repdist_max_median_n) {
    i <- round(seq(1, n, length.out = .repdist_max_median_n))
    Z <- Z[i, , drop = FALSE]
    sq <- sq[i]
  }
  d <- .repdist_dist_block(Z, sq, seq_len(nrow(Z)))
  stats::median(d[lower.tri(d)])
}

#' RBF-MMD distance between samples
#'
#' Treats each sample as an abundance-weighted distribution over protein
#' embeddings and returns the maximum mean discrepancy between those
#' distributions under a median-heuristic RBF kernel. The ground metric is
#' always euclidean: MMD is a distance only under a positive-semi-definite
#' kernel, a Gaussian RBF guarantees that only on a true metric, and `1 -
#' cosine` is not one. On unit-norm embeddings nothing is lost -- the
#' euclidean RBF equals `exp(-(1 - TM) / sigma^2)`.
#'
#' Proteins with zero abundance in every sample are dropped before the kernel
#' or bandwidth is computed. Rarefy `counts` first (e.g. `vegan::rrarefy()`)
#' if samples weren't collected at a common depth.
#'
#' @section Bandwidth:
#' The default `sigma` is the median distance between the retained proteins,
#' taken over an evenly spaced subsample of 2000 of them once there are more.
#' TM-Vec distances are bounded (every pair lies within about 1.3), so with
#' that bandwidth the exponent stays small and the kernel is nearly linear in
#' TM-score: the result is then close to the euclidean distance between the
#' samples' abundance-weighted mean embeddings. On the data shipped with the
#' package that is also where PERMANOVA separation peaks; a smaller `sigma`
#' only loses power. If you want the kernel to discount pairs below a chosen
#' structural similarity, `sigma = sqrt(1 - t)` puts its `1/e` point at
#' TM-score `t`.
#'
#' The default depends on which proteins are retained, so the same two samples
#' sit at a slightly different distance when a third sample brings in new
#' proteins. The bandwidth actually used is returned as `attr(, "sigma")`;
#' pass it back as `sigma` to hold it fixed across runs or studies.
#'
#' @section Memory:
#' Nothing quadratic in the protein count is held. The protein kernel is read
#' in 1000-by-1000 blocks and accumulated directly into the samples-by-samples
#' Gram matrix; sparse `Matrix` counts stay sparse. The default bandwidth reads
#' its median off at most 2000 of the retained proteins.
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
#' @return A `dist` object, for direct use with `vegan::adonis2()`, with the
#'   bandwidth used as `attr(, "sigma")`.
#' @importClassesFrom Matrix sparseMatrix dsparseMatrix
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
                           sigma = NULL) {
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

  counts <- counts[, Matrix::colSums(counts) > 0, drop = FALSE]
  totals <- Matrix::rowSums(counts)
  P <- if (sparse_counts) Matrix::Diagonal(x = 1 / totals) %*% counts else
    counts / totals
  Z <- Z[colnames(counts), , drop = FALSE]
  stopifnot(is.numeric(Z), all(is.finite(Z)))

  n <- nrow(Z)
  sq <- rowSums(Z^2)
  blocks <- split(seq_len(n), ceiling(seq_len(n) / 1000L))
  if (is.null(sigma)) {
    # One protein: the kernel is the 1 x 1 matrix 1 whatever the bandwidth.
    sigma <- if (n > 1L) .repdist_median_dist(Z, sq) else 1
    # A median at rounding level (the Gram identity resolves no finer) means
    # at least half of the sampled pairs coincide, and the kernel would be noise.
    if (sigma^2 <= sqrt(.Machine$double.eps) * max(sq))
      stop("The median protein distance is 0: at least half of the sampled ",
           "protein pairs have identical embeddings. Dereplicate the ",
           "sequences, or pass `sigma`.", call. = FALSE)
  }
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("`sigma` must be a single finite, positive number.", call. = FALSE)

  G <- matrix(0, nrow(P), nrow(P))   # the samples x samples Gram matrix P K P'
  for (u in seq_along(blocks)) {
    i <- blocks[[u]]
    Pi <- P[, i, drop = FALSE]
    for (v in seq_len(u)) {
      j <- blocks[[v]]
      Gij <- as.matrix(Matrix::tcrossprod(
        Pi %*% rbf_kernel(.repdist_dist_block(Z, sq, j, i), sigma),
        P[, j, drop = FALSE]))
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
  if (any(M2 < -tol))
    stop("Negative squared MMD beyond floating-point tolerance -- `K` is not a ",
         "valid positive-semidefinite kernel.", call. = FALSE)
  diag(M2) <- 0
  sqrt(pmax(M2, 0))
}
