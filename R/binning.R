# Medoid of `members`: the one closest to all the others, read straight out of
# the condensed dist so the full matrix is never expanded.
.repdist_medoid <- function(members, D) {
  k <- length(members)
  if (k == 1L) return(members)
  idx <- match(members, attr(D, "Labels"))
  totals <- numeric(k)
  for (j in seq_len(k - 1L)) {
    at <- seq.int(j + 1L, k)
    values <- usedist::dist_get(D, idx[[j]], idx[at])
    totals[[j]] <- totals[[j]] + sum(values)
    totals[at] <- totals[at] + values
  }
  members[[which.min(totals)]]
}

#' Cluster proteins into structural bins
#'
#' Cuts a complete-linkage tree over a [repdist_matrix()] cosine distance at
#' `1 - min_sim`, so every within-bin pair has similarity at least `min_sim`.
#' Bins are named `bin1`, `bin2`, ... from largest to smallest.
#'
#' Structural similarity is not shared function: compare annotation agreement
#' together with non-singleton coverage and fragmentation before transferring
#' annotations to uncharacterised members.
#'
#' Collapse a sample-by-protein counts table to bin level with the returned
#' membership, e.g. `t(rowsum(t(counts), bin[colnames(counts)]))` where `bin`
#' maps protein to bin name.
#'
#' @param D Protein `dist` from [repdist_matrix()], or a square distance matrix
#'   on the same `1 - similarity` scale. A euclidean `dist` is rejected,
#'   because `1 - min_sim` would be a meaningless cut height on it.
#' @param min_sim All-pairs similarity threshold in `(0, 1]` (default 0.7).
#' @return A list:
#'   \describe{
#'     \item{`clusters`}{Named list, one character vector of member accessions
#'       per bin.}
#'     \item{`similarity`}{Bin-by-bin similarity matrix between the
#'       representatives.}
#'     \item{`representatives`}{Bin-named character vector of medoid
#'       accessions -- each bin's member closest to all the others.}
#'     \item{`tree`}{The complete-linkage [stats::hclust] tree over all
#'       proteins. Cut it at `1 - s` with [stats::cutree()] to explore other
#'       thresholds without reclustering.}
#'     \item{`min_sim`}{The similarity parameter the tree was cut at.}
#'   }
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
#' Z <- Z / sqrt(rowSums(Z^2))
#' bin_proteins(repdist_matrix(Z), min_sim = 0.7)
#' @export
bin_proteins <- function(D, min_sim = 0.7) {
  if (!inherits(D, "dist")) {
    if (!is.matrix(D) || !is.numeric(D) || nrow(D) != ncol(D) ||
        !isSymmetric(D) || any(diag(D) != 0))
      stop("`D` must be a symmetric numeric distance matrix with zero diagonal.",
           call. = FALSE)
    D <- stats::as.dist(D)
  }
  labels <- attr(D, "Labels")
  n <- attr(D, "Size")
  stopifnot(length(n) == 1L, is.finite(n), n >= 2L, n == floor(n),
            is.numeric(D), length(D) == n * (n - 1) / 2,
            length(labels) == n, !anyNA(labels), all(nzchar(labels)),
            !anyDuplicated(labels))
  .repdist_check_n(n, "proteins")
  stopifnot(is.numeric(min_sim), length(min_sim) == 1L,
            is.finite(min_sim), min_sim > 0, min_sim <= 1)
  if (identical(attr(D, "method"), "euclidean"))
    stop("`D` must be on the 1 - similarity scale, not euclidean; ",
         "use repdist_matrix(Z).", call. = FALSE)
  # sum(), not all(is.finite()): the latter allocates a logical the size of the
  # ground metric. Distances are non-negative, so any NA/NaN/Inf reaches the sum.
  if (!is.finite(sum(D)))
    stop("`D` contains non-finite values.", call. = FALSE)
  if (min(D) < 0 || max(D) > 2 + 1e-8)
    stop("`D` must contain cosine distances in [0, 2].", call. = FALSE)

  tree <- fastcluster::hclust(D, method = "complete")
  clusters <- unname(split(labels, stats::cutree(tree, h = 1 - min_sim)))
  # cluster names are assigned by the size of clusters.
  clusters <- clusters[order(-lengths(clusters))]
  names(clusters) <- paste0("bin", seq_along(clusters))
  reps <- vapply(clusters, .repdist_medoid, character(1), D = D)
  similarity <- 1 - as.matrix(usedist::dist_subset(D, unname(reps)))
  dimnames(similarity) <- list(names(reps), names(reps))
  list(clusters = clusters,
       similarity = similarity,
       representatives = reps,
       tree = tree,
       min_sim = min_sim)
}
