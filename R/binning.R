.repdist_medoid <- function(members, gm) {
  k <- length(members)
  if (k == 1L) return(members)
  idx <- match(members, attr(gm, "Labels"))
  totals <- numeric(k)
  for (j in seq_len(k - 1L)) {
    at <- seq.int(j + 1L, k)
    values <- .repdist_dist_values(gm, idx[[j]], idx[at])
    totals[[j]] <- totals[[j]] + sum(values)
    totals[at] <- totals[at] + values
  }
  members[[which.min(totals)]]
}

# One dendrogram cut. `gm` is the condensed ground metric shared across cuts.
.bin_result <- function(clust, counts, gm) {
  keep <- names(clust)
  totals <- tapply(colSums(counts), clust, sum)
  rank_id <- rank(-totals, ties.method = "first")
  bin <- stats::setNames(paste0("bin", rank_id[as.character(clust)]), keep)

  # split() and rowsum() both order groups lexically (bin1, bin10, bin2), so
  # everything below is put back in bin order -- counts, reps and bin_dist line
  # up positionally, not just by name.
  members <- split(keep, bin)
  members <- members[order(as.integer(sub("^bin", "", names(members))))]
  reps <- vapply(members, .repdist_medoid, character(1), gm = gm)

  rolled <- t(base::rowsum(t(counts), group = bin[colnames(counts)]))
  rolled <- rolled[, names(members), drop = FALSE]

  rep_df <- data.frame(
    bin = names(reps), protein = unname(reps), size = lengths(members))
  list(counts = rolled, bin = bin,
       bin_dist = .repdist_subset_dist(gm, unname(reps), names(reps)),
       reps = rep_df)
}

# Cluster proteins from a condensed ground distance. Complete linkage cuts a
# dendrogram at a fixed height, so one global threshold applies everywhere;
# the density methods let each cluster keep its own, which matters when
# within-family similarity varies (it ranges 0.47-0.99 across CATH FunFams).
.repdist_cluster_labels <- function(gm, keep, method, cluster_threshold,
                                    min_size, hc = NULL) {
  if (length(keep) < 2L) return(stats::setNames(1L, keep))
  cl <- switch(
    method,
    complete = stats::cutree(
      if (is.null(hc)) fastcluster::hclust(gm, method = "complete") else hc,
      h = 1 - cluster_threshold),
    dbscan = dbscan::dbscan(gm, eps = 1 - cluster_threshold,
                            minPts = min_size)$cluster,
    hdbscan = dbscan::hdbscan(gm, minPts = min_size)$cluster
  )
  # Density methods mark unclustered points 0. A protein that belongs with
  # nothing is its own bin, not a dropped observation.
  cl[cl == 0L] <- max(cl) + seq_len(sum(cl == 0L))
  stats::setNames(as.integer(cl), keep)
}

#' Cluster proteins by structural similarity
#'
#' Returns bin membership without collapsing an abundance table. This is the
#' clustering step of [repdist_bin()], exposed for benchmarking a partition
#' against a reference grouping.
#'
#' @param embeddings Numeric matrix, one row per protein, or a precomputed
#'   protein `dist` from [seq_repdist()].
#' @param cluster_threshold Similarity score defining the bin boundary. Used by
#'   `"complete"` and `"dbscan"`; ignored by `"hdbscan"`.
#' @param method `"complete"` (complete-linkage hierarchical clustering at one
#'   global cut), `"dbscan"` (density-based at one global radius), or
#'   `"hdbscan"` (density-based, each cluster keeping its own radius).
#' @param min_size Minimum points forming a dense region, for the two density
#'   methods.
#' @return A named integer vector of cluster ids, one per protein.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
#' repdist_cluster(Z)
#' @export
repdist_cluster <- function(embeddings, cluster_threshold = 0.5,
                            method = c("complete", "dbscan", "hdbscan"),
                            min_size = 5L) {
  method <- match.arg(method)
  proteins <- .repdist_labels(embeddings)
  stopifnot(!is.null(proteins), length(cluster_threshold) == 1L,
            is.finite(cluster_threshold), min_size >= 1L)
  .repdist_cluster_labels(
    .repdist_ground(embeddings, proteins), proteins, method,
    cluster_threshold, min_size)
}

#' Cluster proteins into structural bins
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
#' Pass a vector of thresholds to compare several cuts. They share one pairwise
#' distance and one dendrogram. Complete linkage defines bin membership so every
#' protein pair in a bin meets the requested TM-score threshold.
#'
#' `bin_dist` is the cosine distance between bin representatives. For
#' bin-resolution MMD, index the original embedding matrix by
#' `out$bins$reps$protein`, rename those rows with `out$bins$reps$bin`, and pass
#' them with `out$bins$counts` to [sample_repdist()].
#'
#' Cost: the ground distance contains n(n-1)/2 doubles (roughly 400 MB at
#' 10,000 proteins). It remains condensed through clustering and medoid
#' selection. Pass a precomputed [seq_repdist()] as `embeddings` to build it
#' once across calls.
#'
#' @param counts Integer matrix of counts, or relative abundances, samples in
#'   rows, proteins in columns.
#' @param embeddings Numeric matrix, one row per protein, rownames matching
#'   `colnames(counts)` -- or a precomputed protein `dist` from [seq_repdist()].
#'   Must be cosine-calibrated to a similarity score for `cluster_threshold` to
#'   mean anything.
#' @param cluster_threshold Similarity score defining the bin boundary
#'   (default 0.5, TM-score's same-fold threshold). A vector cuts the same
#'   dendrogram at each value. Ignored when `method = "hdbscan"`.
#' @param method Clustering method, see [repdist_cluster()]. Default
#'   `"complete"`.
#' @param min_size Minimum points forming a dense region, for the two density
#'   methods.
#' @param seqs Optional named character vector or [Biostrings::AAStringSet] of
#'   protein sequences. It must contain every protein in `counts`.
#' @return A list with two elements:
#'   \describe{
#'     \item{`bins`}{For a scalar threshold, a list containing `counts` (samples
#'       x bins), `bin` (protein to bin membership), `bin_dist` (distance between
#'       representative proteins), and `reps` (bin medoids and sizes). For a
#'       threshold vector, a named list of these objects.}
#'     \item{`representative_sequences`}{A bin-named
#'       [Biostrings::AAStringSet] containing each medoid sequence, or `NULL` if
#'       `seqs` was not supplied. For a threshold vector, a named list of these
#'       objects.}
#'   }
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
#' counts <- rbind(s1 = c(5, 3, 0), s2 = c(0, 2, 7))
#' colnames(counts) <- rownames(Z)
#' repdist_bin(counts, Z)
#' @export
repdist_bin <- function(counts, embeddings, cluster_threshold = 0.5,
                        method = c("complete", "dbscan", "hdbscan"),
                        min_size = 5L, seqs = NULL) {
  method <- match.arg(method)
  counts <- as.matrix(counts)
  proteins <- .repdist_labels(embeddings)
  stopifnot(
    !is.null(rownames(counts)), !is.null(colnames(counts)),
    !is.null(proteins), all(colnames(counts) %in% proteins),
    ncol(counts) >= 1L, is.numeric(counts), all(is.finite(counts)), all(counts >= 0),
    length(cluster_threshold) >= 1L, all(is.finite(cluster_threshold)),
    length(min_size) == 1L, min_size >= 1L
  )
  keep <- colnames(counts)
  if (!is.null(seqs)) seqs <- .repdist_qc_sequences(seqs, keep)
  gm <- .repdist_ground(embeddings, keep)
  # sum(), not all(is.finite()): the latter allocates a logical vector the size
  # of the ground metric (8 GB at 64k proteins). Distances are non-negative, so
  # nothing cancels and any NA/NaN/Inf propagates to the sum.
  if (!is.finite(sum(gm)))
    stop("The ground distance contains non-finite values.", call. = FALSE)

  lone <- length(keep) < 2L   # hclust needs two objects; one protein is one bin
  # fastcluster allocates its own working copy of `gm`, so at 60k+ proteins the
  # two together approach R's vector limit. Drop the block-loop garbage from
  # seq_repdist() first -- without this the peak includes transients that are
  # dead but not yet collected. ponytail: explicit gc() earns its place only
  # here, where the object handed to C is multiple GB.
  if (length(keep) > 20000L) gc(full = TRUE)
  # One dendrogram is shared across cuts; the density methods are re-run per
  # threshold, which is cheap next to the distance itself.
  hc <- if (!lone && method == "complete")
    fastcluster::hclust(gm, method = "complete")

  out <- lapply(cluster_threshold, function(th) {
    x <- .bin_result(
      .repdist_cluster_labels(gm, keep, method, th, min_size, hc),
      counts, gm)
    list(bins = x,
         representative_sequences = if (!is.null(seqs))
           stats::setNames(seqs[x$reps$protein], x$reps$bin))
  })
  if (length(out) == 1L) return(out[[1L]])
  names(out) <- as.character(cluster_threshold)
  list(bins = lapply(out, `[[`, "bins"),
       representative_sequences = lapply(out, `[[`, "representative_sequences"))
}
