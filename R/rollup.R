#' Roll up protein abundance into structural bins
#'
#' Collapses `counts` from protein-level to structural-bin level by clustering
#' `embeddings` and summing abundance within each bin -- the structural
#' analogue of collapsing an ASV table to a coarser rank. Point this at a
#' per-bin differential-abundance test (corncob, ALDEx2, ANCOM-BM, ...)
#' instead of comparing raw proteins one at a time.
#'
#' Cosine similarity between TM-Vec vectors is the predicted TM-score, so
#' `1 - cluster_threshold` is a distance cutoff. Clustering is complete
#' linkage, so every pair within a bin meets the threshold, not just nearest
#' neighbors. Default 0.5 is TM-score's conventional same-fold boundary (Xu &
#' Zhang, Bioinformatics 2010).
#'
#' @param counts Integer matrix of counts, or relative abundances, samples in
#'   rows, proteins in columns.
#' @param embeddings Numeric matrix, one row per protein, rownames matching
#'   `colnames(counts)`. Must be cosine-calibrated to a similarity score for
#'   `cluster_threshold` to mean anything.
#' @param method Distance method passed to `proxy::dist()`. Only `"cosine"`
#'   (the default) is calibrated against `cluster_threshold` as a TM-score
#'   cutoff -- changing it changes what the threshold means.
#' @param cluster_threshold Similarity score defining the bin boundary
#'   (default 0.5, TM-score's same-fold threshold).
#' @return List with `counts` (samples x bins matrix, same row order as
#'   input) and `bin` (named character vector, protein -> bin id; `bin1` is
#'   the highest-total-abundance bin, `bin2` the next, ...).
#' @export
rollupStructural <- function(counts, embeddings, method = "cosine", cluster_threshold = 0.5) {
  counts <- as.matrix(counts)
  stopifnot(all(colnames(counts) %in% rownames(embeddings)))
  Z <- embeddings[colnames(counts), , drop = FALSE]
  D <- proxy::dist(Z, method)

  clust <- stats::cutree(fastcluster::hclust(D, method = "complete"), h = 1 - cluster_threshold)
  totals <- tapply(colSums(counts), clust, sum)
  rank_id <- rank(-totals, ties.method = "first")
  bin <- stats::setNames(paste0("bin", rank_id[as.character(clust)]), rownames(Z))

  rolled <- t(base::rowsum(t(counts), group = bin[colnames(counts)]))
  list(counts = rolled, bin = bin)
}
