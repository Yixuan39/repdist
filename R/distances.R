#' Protein ground metric from embeddings
#'
#' The one distance matrix everything else is built on: [sample_repdist()] and
#' [repdist_bin()] consume it. At real proteome
#' sizes it is the expensive object (n(n-1)/2 doubles -- roughly 400 MB at
#' 10,000 proteins), so compute it once and pass it around rather than letting
#' each function rebuild it. Cosine distances are written directly into this
#' condensed representation in blocks; a square protein-by-protein matrix is
#' never materialized.
#'
#' `"cosine"` is `1 - cosine similarity`. For TM-Vec embeddings cosine
#' similarity is a *predicted TM-score*, so this is `1 - TM-score` -- the same
#' calibrated scale [repdist_bin()]'s `cluster_threshold` is expressed in,
#' and the reason it is this function's default and the one [repdist_bin()]
#' uses. Caveat worth knowing: `1 - cosine` is a dissimilarity, not a metric --
#' it violates the triangle inequality (0, 45, 90 degrees apart gives
#' 0.293 + 0.293 < 1). That is harmless for the complete-linkage clustering
#' [repdist_bin()] does.
#'
#' It is **not** harmless for [sample_repdist()]. A Gaussian RBF is
#' positive-semi-definite only on a conditionally-negative-definite metric
#' (Schoenberg), so `rbf_kernel(1 - cosine)` can be indefinite -- on real
#' TM-Vec embeddings its smallest eigenvalue runs to `-5e-1` -- and MMD under
#' an indefinite kernel is not a distance. That is why [sample_repdist()] always
#' uses `"euclidean"`, while this function and [repdist_bin()] default to
#' `"cosine"` where that scale is valid.
#'
#' Nothing is given up by that split. `embed_proteins()` returns unit-norm
#' rows, and for unit-norm vectors `||a - b|| = sqrt(2 * (1 - cos))` -- the two
#' ground metrics are monotone transforms of each other (rank correlation 1.0
#' to eight decimals on real TM-Vec embeddings), so they induce the same
#' ordering of protein pairs. Euclidean simply expresses it on a scale that is
#' a true metric. The equivalence holds *only* for unit-norm rows: on an
#' embedding matrix with varying row norms, `"euclidean"` also picks up
#' magnitude, which `"cosine"` normalizes away.
#'
#' @param embeddings Numeric matrix, one row per protein.
#' @param distance "cosine" (1 - predicted TM-score; not a metric) or
#'   "euclidean" (a true metric, and the one to use for a kernel).
#' @return A `dist` object over proteins.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' seq_repdist(Z)
#' @export
seq_repdist <- function(embeddings, distance = c("cosine", "euclidean")) {
  distance <- match.arg(distance)
  Z <- as.matrix(embeddings)
  stopifnot(is.numeric(Z), all(is.finite(Z)))
  if (distance == "euclidean") return(stats::dist(Z))

  norms <- sqrt(rowSums(Z^2))
  if (any(norms == 0))
    stop("Cosine distance is undefined for zero-length embedding rows: ",
         paste(utils::head(rownames(Z)[norms == 0], 5), collapse = ", "),
         call. = FALSE)
  .repdist_cosine_dist(Z / norms)
}

.repdist_cosine_dist <- function(Z, block_bytes = 256 * 1024^2) {
  n <- nrow(Z)
  out <- numeric(as.double(n) * (n - 1) / 2)
  if (n >= 2L) {
    block_cols <- max(1L, min(n - 1L,
      as.integer(floor(block_bytes / (8 * n)))))
    cursor <- 1
    for (start in seq.int(1L, n - 1L, by = block_cols)) {
      cols <- seq.int(start, min(n - 1L, start + block_cols - 1L))
      similarities <- tcrossprod(
        Z[seq.int(start + 1L, n), , drop = FALSE],
        Z[cols, , drop = FALSE])
      for (k in seq_along(cols)) {
        len <- n - cols[[k]]
        at <- seq.int(cursor, length.out = len)
        out[at] <- 1 - similarities[seq.int(k, nrow(similarities)), k]
        cursor <- cursor + len
      }
    }
  }
  structure(out, Size = n, Labels = rownames(Z), Diag = FALSE, Upper = FALSE,
            method = "cosine", class = "dist")
}

.repdist_dist_values <- function(x, i, j) {
  n <- as.double(attr(x, "Size"))
  lo <- pmin(as.double(i), as.double(j))
  hi <- pmax(as.double(i), as.double(j))
  out <- numeric(length(lo))
  keep <- lo != hi
  at <- n * (lo[keep] - 1) - lo[keep] * (lo[keep] - 1) / 2 +
    hi[keep] - lo[keep]
  out[keep] <- unclass(x)[at]
  out
}

# Protein labels, however the ground metric arrived.
.repdist_labels <- function(x)
  if (inherits(x, "dist")) attr(x, "Labels") else rownames(x)

# The ground metric over `keep`, built from embeddings or reused from a `dist`.
.repdist_ground <- function(x, keep, distance = "cosine") {
  if (!inherits(x, "dist")) return(seq_repdist(x[keep, , drop = FALSE], distance))
  if (identical(keep, attr(x, "Labels"))) x else .repdist_subset_dist(x, keep)
}

.repdist_subset_dist <- function(x, members, labels = members) {
  idx <- match(members, attr(x, "Labels"))
  if (anyNA(idx))
    stop("Distance object is missing requested protein(s).", call. = FALSE)
  k <- length(idx)
  out <- numeric(as.double(k) * (k - 1) / 2)
  cursor <- 1
  if (k >= 2L) {
    for (j in seq_len(k - 1L)) {
      values <- .repdist_dist_values(x, idx[[j]], idx[seq.int(j + 1L, k)])
      at <- seq.int(cursor, length.out = length(values))
      out[at] <- values
      cursor <- cursor + length(values)
    }
  }
  structure(out, Size = k, Labels = labels, Diag = FALSE, Upper = FALSE,
            method = attr(x, "method"), class = "dist")
}

#' RBF-MMD distance between samples
#'
#' Treats each sample as a distribution over protein representations, and
#' returns their maximum mean discrepancy under a median-heuristic RBF
#' kernel -- the CMMD form. It picks up higher-order differences in how mass is
#' spread out, not just where its mean sits. The protein x protein Gram matrix
#' is materialized once; pass a precomputed [seq_repdist()] object to reuse a
#' ground distance already in memory.
#'
#' Proteins with zero abundance in every sample are dropped before the weights,
#' kernel, or bandwidth are computed -- they carry no signal and would
#' otherwise be free to distort the median-heuristic bandwidth.
#'
#' The ground metric is always Euclidean here. MMD is only a distance when the
#' kernel is positive semi-definite, and a Gaussian RBF on `1 - cosine` can be
#' indefinite. Nothing is lost on unit-norm embeddings, because their
#' Euclidean RBF *is* a kernel in
#' cosine distance: `d_euclid^2 = 2 * d_cosine` exactly, so
#' `exp(-d_euclid^2 / 2 sigma^2) = exp(-d_cosine / sigma^2)`. Choosing
#' `"euclidean"` is therefore not a retreat from predicted TM-score -- it is
#' `exp(-(1 - TM) / sigma^2)`, with the exponent linear in TM-distance rather
#' than squared.
#'
#' Rarefy `counts` first (e.g. `vegan::rrarefy()`) if samples weren't
#' collected at a common depth.
#'
#' @param counts Integer matrix of counts, or relative abundances, samples in
#'   rows, proteins in columns. Every row must have a positive total.
#' @param embeddings Numeric matrix, one row per protein with rownames matching
#'   `colnames(counts)`, or a named protein [stats::dist()] object such as the
#'   result of `seq_repdist(Z, "euclidean")`. Non-Euclidean precomputed
#'   distances are rejected.
#' @param weighted Use abundance as the weight (default). `FALSE` replaces each
#'   sample's abundances with `sign(abundance)`, so every expressed protein
#'   carries equal weight -- the presence/absence counterpart, standing to
#'   `weighted = TRUE` as unweighted UniFrac stands to weighted UniFrac.
#' @param sigma RBF bandwidth. Defaults to the median heuristic.
#' @return A `dist` object.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(1, 1))
#' counts <- rbind(s1 = c(8, 1, 1), s2 = c(1, 8, 1))
#' colnames(counts) <- rownames(Z)
#' sample_repdist(counts, Z)
#' @export
sample_repdist <- function(counts,
                           embeddings,
                           weighted = TRUE,
                           sigma = NULL) {
  counts <- as.matrix(counts)
  proteins <- .repdist_labels(embeddings)

  stopifnot(
    !is.null(rownames(counts)), !is.null(colnames(counts)),
    !is.null(proteins), all(colnames(counts) %in% proteins),
    is.numeric(counts), all(is.finite(counts)), all(counts >= 0),
    all(rowSums(counts) > 0),
    is.logical(weighted), length(weighted) == 1L, !is.na(weighted)
  )

  is_dist <- inherits(embeddings, "dist")
  if (is_dist && !identical(attr(embeddings, "method"), "euclidean"))
    stop("`embeddings` must be a matrix or a Euclidean `dist`; cosine and ",
         "unknown precomputed distances do not define a valid RBF-MMD.",
         call. = FALSE)

  if (!weighted) counts <- sign(counts)

  counts <- counts[, colSums(counts) > 0, drop = FALSE]
  keep <- colnames(counts)
  P <- counts / rowSums(counts)   # the weights, one row per sample

  G <- as.matrix(.repdist_ground(embeddings, keep, "euclidean"))
  D <- mmd_matrix(P, rbf_kernel(G, sigma))
  dimnames(D) <- list(rownames(counts), rownames(counts))
  stats::as.dist(D)
}

#' RBF Gram matrix from a ground-metric distance matrix
#'
#' @param C Protein x protein distance matrix.
#' @param sigma Bandwidth; median off-diagonal distance if NULL.
#' @return A protein x protein RBF Gram matrix.
#' @examples
#' C <- as.matrix(stats::dist(matrix(c(0, 0, 1, 0, 0, 1), ncol = 2,
#'                                   byrow = TRUE)))
#' rbf_kernel(C)
#' @export
rbf_kernel <- function(C, sigma = NULL) {
  if (nrow(C) < 2L)
    return(matrix(1, nrow(C), ncol(C), dimnames = dimnames(C)))
  if (is.null(sigma)) sigma <- stats::median(C[upper.tri(C)])
  if (!is.finite(sigma) || sigma < 0)
    stop("`sigma` must be a finite, non-negative number.", call. = FALSE)
  if (sigma == 0)
    # Median off-diagonal distance of exactly 0 means at least half the proteins
    # are exact duplicates in embedding space -- e.g. the diagonal always. The
    # division this would otherwise do is 0/0 there, not a real degeneracy: the
    # sigma -> 0 limit of a Gaussian kernel is the indicator of exact coincidence.
    return(matrix(as.numeric(C == 0), nrow(C), ncol(C), dimnames = dimnames(C)))
  exp(-(C^2) / (2 * sigma^2))
}

#' Abundance-weighted MMD distance matrix
#'
#' @param P Relative abundance matrix, samples in rows.
#' @param K Protein x protein kernel matrix.
#' @return A samples x samples matrix of MMD distances.
#' @examples
#' P <- rbind(s1 = c(0.8, 0.2), s2 = c(0.1, 0.9))
#' K <- matrix(c(1, 0.2, 0.2, 1), 2)
#' mmd_matrix(P, K)
#' @export
mmd_matrix <- function(P, K) {
  PK <- P %*% K   # the samples x proteins product, reused for both terms
  self_term <- rowSums(PK * P)
  M2 <- outer(self_term, self_term, "+") - 2 * tcrossprod(PK, P)
  # Scale the tolerance to the kernel (the self terms), not to M2 itself. When
  # every sample has the same composition M2 is exactly 0 in exact arithmetic,
  # so a tolerance relative to max(abs(M2)) collapses onto the rounding noise it
  # is meant to absorb and rejects a perfectly valid kernel.
  tol <- 1e-8 * max(self_term)
  if (any(M2 < -tol))
    stop("Negative squared MMD beyond floating-point tolerance -- `K` is not a ",
         "valid positive-semidefinite kernel.", call. = FALSE)
  diag(M2) <- 0
  sqrt(pmax(M2, 0))
}

#' Optimal-transport distance between samples
#'
#' Earth mover's distance between samples, where each sample is the distribution
#' of its relative abundance over proteins and the cost of moving abundance from
#' one protein to another is supplied as a ground metric. Where Bray-Curtis can
#' only ask whether two samples hold the *same* protein, this asks how far apart
#' the proteins they hold actually are. Under a ground metric that puts every
#' pair of distinct proteins exactly 1 apart, it *is* Bray-Curtis, to machine
#' precision; under one that measures function, samples built from different
#' accessions doing the same job come out close together.
#'
#' @section Scale:
#' The result inherits the scale of `D`. Supply a ground metric already on 0-1 --
#' `1 - GO similarity`, `1 - TM-score` and `1 - cosine` all are -- and the output
#' is directly comparable with Bray-Curtis: two samples sharing no proteins, all
#' of them maximally distant, score exactly 1.
#'
#' It reaches 1 only when *every* cross-sample protein pair is maximally distant.
#' Samples that share no accession but whose proteins are merely different score
#' lower, because transport routes abundance through whichever pairs are
#' cheapest. That is the purpose of a ground metric rather than a defect: two
#' samples with no accession in common but the same functions are genuinely not
#' maximally far apart, and this is the sense in which Bray-Curtis, which calls
#' them 1, is wrong.
#'
#' `D` is used as given, deliberately: dividing it by its largest entry would
#' make the answer depend on which proteins happen to be present, stretching a
#' set of functionally identical proteins 0.02 apart across the full range to
#' score a maximal 1.
#'
#' @section Exact transport:
#' Solved with \pkg{T4transport}'s `wassersteinD()`. Each sample pair is its own
#' transport problem, sized to the two samples' supports rather than the whole
#' ground metric, and at the sizes a beta-diversity comparison usually sees --
#' tens to a few hundred distinct proteins per sample -- that is already fast:
#' under a second for thousands of sample pairs.
#'
#' @param counts Sample-by-protein abundance matrix, samples in rows. Rows are
#'   converted to relative abundance, so samples need not share a depth.
#' @param D Protein ground metric, a `dist` or a square matrix, on whatever scale
#'   the answer should be in. Names must match `colnames(counts)`. Proteins
#'   absent from every sample are dropped.
#' @return A `dist` object over the samples, for direct use with
#'   `vegan::adonis2()` or [stats::cmdscale()].
#' @seealso [sample_repdist()], the MMD counterpart.
#'
#' @examples
#' counts <- rbind(s1 = c(10, 0, 0), s2 = c(0, 10, 0), s3 = c(0, 0, 10))
#' colnames(counts) <- c("p1", "p2", "p3")
#' # p1 and p2 are near-identical; p3 is unrelated to both
#' D <- matrix(c(0, 0.1, 1, 0.1, 0, 1, 1, 1, 0), 3,
#'             dimnames = list(colnames(counts), colnames(counts)))
#' sample_repdist_OT(counts, D)
#'
#' @export
sample_repdist_OT <- function(counts, D) {
  counts <- as.matrix(counts)
  if (inherits(D, "dist")) D <- as.matrix(D)
  stopifnot(!is.null(colnames(counts)), nrow(counts) >= 2L,
            is.numeric(counts), all(counts >= 0), all(rowSums(counts) > 0))

  keep <- colnames(counts)[colSums(counts) > 0]
  gone <- setdiff(keep, rownames(D))
  if (length(gone))
    stop(sprintf("%d protein(s) missing from the ground metric, e.g. %s",
                 length(gone), paste(utils::head(gone, 3L), collapse = ", ")),
         call. = FALSE)

  D <- D[keep, keep, drop = FALSE]
  if (anyNA(D) || min(D) < -1e-8)
    stop("`D` must be complete and non-negative", call. = FALSE)
  D[D < 0] <- 0                 # cosine round-off, not a real negative

  P <- counts[, keep, drop = FALSE]
  P <- P / rowSums(P)

  # Supports are found once per sample rather than once per pair: every pair
  # reuses them, and each transport problem is only as big as the two supports.
  supp <- lapply(seq_len(nrow(P)), function(i) which(P[i, ] > 0))
  pairs <- utils::combn(nrow(P), 2L)
  vals <- apply(pairs, 2L, function(ij) {
    a <- supp[[ij[1L]]]; b <- supp[[ij[2L]]]
    T4transport::wassersteinD(D[a, b, drop = FALSE], p = 1,
                              wx = P[ij[1L], a], wy = P[ij[2L], b])$distance
  })

  out <- matrix(0, nrow(P), nrow(P), dimnames = list(rownames(P), rownames(P)))
  out[t(pairs)] <- vals
  stats::as.dist(out + t(out))
}
