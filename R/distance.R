#' Weighted-representation distance between samples
#'
#' Treats each sample as its abundance-weighted mean over protein
#' representations, and returns an MMD-family distance between samples,
#' ready for `vegan::adonis2()`.
#'
#' `"mmd"` is the Euclidean distance between weighted mean representations
#' (linear kernel). `"rbf_mmd"` swaps in a median-heuristic RBF kernel, which
#' also picks up higher-order differences in how mass is spread out, not just
#' the mean.
#'
#' Rarefy `counts` first (e.g. `vegan::rrarefy()`) if samples weren't
#' collected at a common depth.
#'
#' @param counts Integer matrix of counts, or relative abundances, samples in
#'   rows, proteins in columns.
#' @param embeddings Numeric matrix, one row per protein, rownames matching
#'   `colnames(counts)`.
#' @param method "mmd" (linear MMD) or "rbf_mmd".
#' @param sigma RBF bandwidth. Defaults to the median heuristic.
#' @return A `dist` object.
#' @export
weightedRep <- function(counts, embeddings,
                        method = c("mmd", "rbf_mmd"),
                        sigma = NULL) {
  method <- match.arg(method)
  counts <- as.matrix(counts)

  stopifnot(all(colnames(counts) %in% rownames(embeddings)))
  Z <- embeddings[colnames(counts), , drop = FALSE]

  P <- counts / rowSums(counts)   # the abundance weights, one row per sample

  K <- switch(method,
    mmd     = Z %*% t(Z),
    rbf_mmd = rbf_kernel(as.matrix(stats::dist(Z)), sigma)
  )
  D <- mmd_matrix(P, K)
  dimnames(D) <- list(rownames(counts), rownames(counts))
  stats::as.dist(D)
}

#' RBF Gram matrix from a ground-metric distance matrix
#'
#' @param C Protein x protein distance matrix.
#' @param sigma Bandwidth; median off-diagonal distance if NULL.
#' @export
rbf_kernel <- function(C, sigma = NULL) {
  if (is.null(sigma)) sigma <- stats::median(C[upper.tri(C)])
  exp(-(C^2) / (2 * sigma^2))
}

#' Abundance-weighted MMD distance matrix
#'
#' @param P Relative abundance matrix, samples in rows.
#' @param K Protein x protein kernel matrix.
#' @export
mmd_matrix <- function(P, K) {
  self_term <- rowSums((P %*% K) * P)
  M2 <- outer(self_term, self_term, "+") - 2 * (P %*% K %*% t(P))
  diag(M2) <- 0
  sqrt(pmax(M2, 0))
}
