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
#' Cluster a [repdist_matrix()] cosine distance. The default complete-linkage
#' tree is cut at `1 - min_sim`, so every within-bin pair has similarity at
#' least `min_sim`. Optional graph and density methods explore groups for
#' annotation; none of these methods guarantees shared function. Bins are named
#' `bin1`, `bin2`, ... from largest to smallest.
#'
#' MCL clusters an undirected graph weighted by `1 - D`, keeping positive edges
#' with similarity at least `min_sim`. Self-loops retain isolated proteins;
#' MCL applies its default loop weighting.
#' Inflation controls granularity; increasing it usually yields smaller groups.
#' `min_sim = 0` keeps all positive edges but can make the graph very large.
#' MCL needs the external `mcl` executable, available from
#' <https://micans.org/mcl/>. Edges are streamed from the condensed distance.
#'
#' DBSCAN uses neighbourhood radius `1 - min_sim`; this is still a global
#' parameter, not an all-pairs bound. By default, border points are noise
#' (DBSCAN*) to avoid arbitrary assignment between dense groups. HDBSCAN
#' selects stable density clusters across radii and ignores `min_sim`.
#' Density methods need the optional `dbscan` R package. Every noise protein
#' becomes its own singleton bin and is also listed in `noise`, so abundance
#' is preserved without treating unrelated noise as one functional group.
#' Compare annotation agreement together with non-singleton coverage and
#' fragmentation before transferring annotations to uncharacterised members.
#'
#' Collapse a sample-by-protein counts table to bin level with the returned
#' membership, e.g. `t(rowsum(t(counts), bin[colnames(counts)]))` where `bin`
#' maps protein to bin name.
#'
#' @param D Protein `dist` from [repdist_matrix()], or a square distance matrix
#'   on the same `1 - similarity` scale. A euclidean `dist` is rejected,
#'   because `1 - min_sim` would be a meaningless cut height on it.
#' @param min_sim Similarity threshold (default 0.7): all-pairs minimum for
#'   `"hclust"`, edge minimum for `"mcl"`, neighbourhood minimum for
#'   `"dbscan"`. In `(0, 1]`, or `[0, 1]` for MCL. Ignored by HDBSCAN.
#' @param method `"hclust"` (default), `"mcl"`, `"dbscan"`, or `"hdbscan"`.
#' @param min_pts Integer at least 2, including the point itself, defining a
#'   dense neighbourhood for DBSCAN/HDBSCAN (default 5).
#' @param border_points Include DBSCAN border points (default `FALSE`, DBSCAN*).
#' @param inflation MCL inflation, a finite number greater than 1 (default 2).
#' @param mcl_bin Name or path of the MCL executable (default `"mcl"`).
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
#'       thresholds without reclustering. `NULL` for other methods.}
#'     \item{`min_sim`}{The similarity parameter; `NULL` for HDBSCAN. Only
#'       hclust guarantees this minimum within bins.}
#'     \item{`method`}{The selected clustering method.}
#'     \item{`parameters`}{Additional parameters used by the selected method.}
#'     \item{`noise`}{Accessions rejected by a density method, each retained
#'       in its own singleton bin. Empty for hclust and MCL.}
#'   }
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
#' bin_proteins(repdist_matrix(Z), min_sim = 0.7)
#' @export
bin_proteins <- function(D, min_sim = 0.7,
                         method = c("hclust", "mcl", "dbscan", "hdbscan"),
                         min_pts = 5L, border_points = FALSE,
                         inflation = 2, mcl_bin = "mcl") {
  method <- match.arg(method)
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
  if (method != "hdbscan")
    stopifnot(is.numeric(min_sim), length(min_sim) == 1L,
              is.finite(min_sim), min_sim <= 1,
              if (method == "mcl") min_sim >= 0 else min_sim > 0)
  else if (!missing(min_sim))
    warning("`min_sim` is ignored by HDBSCAN.", call. = FALSE)
  if (identical(attr(D, "method"), "euclidean"))
    stop("`D` must be on the 1 - similarity scale, not euclidean; ",
         "use repdist_matrix(Z, \"cosine\").", call. = FALSE)
  # sum(), not all(is.finite()): the latter allocates a logical the size of the
  # ground metric. Distances are non-negative, so any NA/NaN/Inf reaches the sum.
  if (!is.finite(sum(D)))
    stop("`D` contains non-finite values.", call. = FALSE)
  if (min(D) < 0 || max(D) > 2 + 1e-8)
    stop("`D` must contain cosine distances in [0, 2].", call. = FALSE)

  tree <- NULL
  noise <- character()
  parameters <- list()
  if (method == "hclust") {
    tree <- fastcluster::hclust(D, method = "complete")
    cl <- stats::cutree(tree, h = 1 - min_sim)
  } else if (method == "mcl") {
    stopifnot(is.numeric(inflation), length(inflation) == 1L,
              is.finite(inflation), inflation > 1)
    cl <- .repdist_mcl(D, min_sim, inflation, mcl_bin)
    parameters <- list(inflation = inflation, mcl_bin = mcl_bin)
  } else {
    stopifnot(is.numeric(min_pts), length(min_pts) == 1L,
              is.finite(min_pts), min_pts >= 2, min_pts <= .Machine$integer.max,
              min_pts == floor(min_pts), is.logical(border_points),
              length(border_points) == 1L, !is.na(border_points))
    if (!requireNamespace("dbscan", quietly = TRUE))
      stop("This method needs the dbscan R package.", call. = FALSE)
    parameters <- list(min_pts = min_pts)
    if (method == "dbscan") {
      cl <- dbscan::dbscan(D, eps = 1 - min_sim, minPts = min_pts,
                           borderPoints = border_points)$cluster
      parameters <- c(parameters, list(border_points = border_points))
    } else {
      cl <- if (min_pts > n) integer(n) else
        dbscan::hdbscan(D, minPts = min_pts)$cluster
      min_sim <- NULL
    }
    noise <- labels[cl == 0L]
    cl[cl == 0L] <- max(cl) + seq_along(noise)
  }
  clusters <- unname(split(labels, cl))
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
       min_sim = min_sim,
       method = method,
       parameters = parameters,
       noise = noise)
}

# Use the installed MCL implementation, with integer node IDs so arbitrary
# accessions cannot change the ABC file format or the command being executed.
.repdist_mcl <- function(D, min_sim, inflation, mcl_bin) {
  if (!is.character(mcl_bin) || length(mcl_bin) != 1L || is.na(mcl_bin) ||
      !nzchar(mcl_bin) || !nzchar(Sys.which(mcl_bin)))
    stop("MCL executable not found; install mcl and/or set `mcl_bin`.",
         call. = FALSE)
  work <- tempfile("repdist-mcl-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE), add = TRUE)
  input <- file.path(work, "graph.abc")
  output <- file.path(work, "clusters.txt")
  log <- file.path(work, "mcl.log")
  con <- file(input, "wt")
  tryCatch({
    n <- attr(D, "Size")
    writeLines(paste(seq_len(n), seq_len(n), 1, sep = "\t"), con)
    offset <- 0
    # ponytail: scans all O(n^2) distances; accept a sparse neighbour graph
    # directly if catalogs outgrow the existing condensed-distance API.
    for (j in seq_len(n - 1L)) {
      at <- seq_len(n - j)
      similarity <- 1 - D[offset + at]
      keep <- which(similarity > 0 & similarity >= min_sim)
      if (length(keep))
        writeLines(paste(j, j + keep, sprintf("%.17g", similarity[keep]),
                         sep = "\t"), con)
      offset <- offset + length(at)
    }
  }, finally = close(con))
  status <- system2(mcl_bin, c(shQuote(input), "--abc", "-I",
                              format(inflation, scientific = FALSE),
                              "-o", shQuote(output)), stdout = log, stderr = log)
  if (status != 0L || !file.exists(output))
    stop("MCL failed: ", paste(utils::tail(readLines(log, warn = FALSE), 5L),
                               collapse = "\n"), call. = FALSE)
  groups <- strsplit(readLines(output, warn = FALSE), "[[:space:]]+")
  ids <- suppressWarnings(as.integer(unlist(groups)))
  if (length(ids) != n || anyNA(ids) || anyDuplicated(ids) ||
      !setequal(ids, seq_len(n)))
    stop("MCL output does not partition every input protein exactly once.",
         call. = FALSE)
  cl <- integer(n)
  cl[ids] <- rep(seq_along(groups), lengths(groups))
  cl
}
