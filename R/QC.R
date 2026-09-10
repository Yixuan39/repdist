utils::globalVariables(c("similarity", "min_sim", "value", "curve", "metric",
                         "protein", "partner"))

# usedist indexes condensed distances with integer arithmetic, so n * (n - 1)
# must stay under .Machine$integer.max; 46341 is the exact maximum. Past it
# usedist::dist_get() returns NA instead of a distance, silently. Called from
# .repdist_qc_sequences() below, so an oversized catalog fails before hours of
# GPU work, and from repdist_matrix() for embeddings that arrived elsewhere.
.repdist_max_proteins <- 46341L

.repdist_check_n <- function(n, what) {
  if (n > .repdist_max_proteins)
    stop(n, " ", what, " exceeds the ", .repdist_max_proteins,
         " this package can index. Reduce sequence redundancy first -- ",
         "cluster with MMseqs2 or diamond and keep one representative per ",
         "cluster.", call. = FALSE)
}

# Validate once for every public entry point that accepts protein sequences.
# Returns the cleaned named character vector the embedders consume.
.repdist_qc_sequences <- function(seqs) {
  ids <- names(seqs)          # as.character() drops them, and AAStringSet in
  seqs <- as.character(seqs)  # is as valid an input here as a character vector
  if (!length(seqs))
    stop("`seqs` must contain at least one sequence.", call. = FALSE)
  .repdist_check_n(length(seqs), "sequences")
  if (is.null(ids) || anyDuplicated(ids) || any(!nzchar(ids)))
    stop("`seqs` names must be unique and non-empty.", call. = FALSE)
  if (anyNA(seqs))
    stop("`seqs` must not contain missing values.", call. = FALSE)

  seqs <- toupper(sub("\\*$", "", seqs))
  if (any(empty <- !nzchar(seqs)))
    stop("`seqs` contains empty sequence(s): ",
         paste(ids[empty], collapse = ", "), call. = FALSE)
  # Gaps, internal stops and stray symbols reach the tokenizer as <unk> and
  # embed to noise without complaint; an aligned FASTA is the usual way in.
  if (any(off <- grepl("[^ACDEFGHIKLMNPQRSTVWYBXZUO]", seqs)))
    stop("model-incompatible symbol(s): ", paste(ids[off], collapse = ", "),
         call. = FALSE)
  stats::setNames(seqs, ids)
}


#' Plot the distribution of nearest-neighbour protein similarity
#'
#' Histogram of each protein's similarity to its nearest neighbour in the
#' ground metric -- for TM-Vec embeddings, a predicted TM-score. A pile-up
#' at 1 means near-duplicate proteins, or an embedding model too saturated to
#' resolve similarity thresholds; every protein left of a candidate `min_sim`
#' becomes a singleton bin in [bin_proteins()].
#'
#' @param D Protein `dist` from [repdist_matrix()], or a square matrix on the
#'   same `1 - similarity` scale. A euclidean `dist` is rejected.
#' @param min_sim Optional binning threshold to mark with a vertical line.
#' @return A ggplot object. Needs the `ggplot2` package.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   Z <- matrix(rnorm(200), 50, 4, dimnames = list(paste0("p", 1:50), NULL))
#'   Z <- Z / sqrt(rowSums(Z^2))
#'   plot_similarity_profile(repdist_matrix(Z), min_sim = 0.7)
#' }
#' @export
plot_similarity_profile <- function(D, min_sim = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_similarity_profile() needs the ggplot2 package.", call. = FALSE)
  if (!inherits(D, "dist")) D <- stats::as.dist(D)
  if (identical(attr(D, "method"), "euclidean"))
    stop("`D` must be on the 1 - similarity scale, not euclidean; ",
         "use repdist_matrix(Z).", call. = FALSE)
  stopifnot(attr(D, "Size") >= 2L)

  # ponytail: expands the condensed dist to a square matrix; switch to a
  # blocked scan (see .repdist_medoid) if QC at bin_proteins()'s ~75k ceiling
  # ever matters.
  M <- as.matrix(D)
  diag(M) <- NA
  df <- data.frame(similarity = 1 - apply(M, 1, min, na.rm = TRUE))

  p <- ggplot2::ggplot(df, ggplot2::aes(similarity)) +
    ggplot2::geom_histogram(binwidth = 0.02, boundary = 0) +
    ggplot2::labs(x = "nearest-neighbour similarity", y = "proteins") +
    ggplot2::theme_minimal()
  if (!is.null(min_sim))
    p <- p + ggplot2::geom_vline(
      xintercept = min_sim, linetype = "dashed", color = "grey30")
  p
}

#' Plot bin count, chaining and singletons across binning thresholds
#'
#' Sweeps `min_sim` over the complete-linkage tree [bin_proteins()] already
#' computed -- no reclustering. The lower panel brackets the two ways binning
#' stops being informative: proteins collapsing into one bin as the threshold
#' drops is chaining, and proteins falling out as singletons as it rises is no
#' aggregation at all. A usable threshold sits in the valley between the two
#' curves, on a plateau of the bin count.
#'
#' @param bins Result of [bin_proteins()], or its `tree` ([stats::hclust]).
#' @param min_sim Optional threshold in use, marked with a vertical line.
#' @return A ggplot object. Needs the `ggplot2` package.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   Z <- matrix(rnorm(200), 50, 4, dimnames = list(paste0("p", 1:50), NULL))
#'   Z <- Z / sqrt(rowSums(Z^2))
#'   b <- bin_proteins(repdist_matrix(Z), min_sim = 0.7)
#'   plot_bin_profile(b, min_sim = 0.7)
#' }
#' @export
plot_bin_profile <- function(bins, min_sim = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_bin_profile() needs the ggplot2 package.", call. = FALSE)
  tree <- if (inherits(bins, "hclust")) bins else bins$tree
  if (!inherits(tree, "hclust"))
    stop("`bins` must be a bin_proteins() result or an hclust tree.",
         call. = FALSE)

  grid <- seq(0.01, 1, by = 0.01)
  cut <- stats::cutree(tree, h = 1 - grid)
  n <- nrow(cut)
  share <- function(f) 100 * apply(cut, 2, function(b) sum(f(tabulate(b)))) / n
  df <- data.frame(
    min_sim = rep(grid, 3L),
    value = c(apply(cut, 2, max), share(max), share(function(k) k == 1L)),
    curve = rep(c("bins", "largest bin", "singletons"), each = length(grid)),
    metric = rep(c("number of bins", "% of proteins", "% of proteins"),
                 each = length(grid)))
  df$metric <- factor(df$metric, unique(df$metric))

  p <- ggplot2::ggplot(df, ggplot2::aes(min_sim, value, color = curve)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::facet_wrap(~metric, ncol = 1L, scales = "free_y") +
    ggplot2::scale_color_discrete(
      name = NULL, breaks = c("largest bin", "singletons")) +
    ggplot2::labs(x = "minimum within-bin similarity (min_sim)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")
  if (!is.null(min_sim))
    p <- p + ggplot2::geom_vline(
      xintercept = min_sim, linetype = "dashed", color = "grey30")
  p
}

#' Plot the similarity matrix inside one bin
#'
#' Heatmap of every pairwise similarity among one bin's members -- for TM-Vec
#' embeddings, predicted TM-scores -- ordered by the clustering tree.
#' [bin_proteins()] guarantees the whole block sits at or above `min_sim`; this
#' is where you see how far above, and whether a bin is one tight group or two
#' sub-groups that only just met the threshold. The fill scale starts at that
#' threshold, so bins are comparable across calls.
#'
#' @details
#' The whole square is drawn, with the diagonal running from the upper left to
#' the lower right. The matrix is symmetric, so each pair appears twice; that
#' redundancy is what lets a row be read straight across, and keeps a
#' sub-group that straddles the diagonal square rather than triangular.
#'
#' @param bins Result of [bin_proteins()].
#' @param bin Bin name, e.g. `"bin1"`.
#' @param D The protein `dist` that produced `bins`.
#' @param label_max Above this many members, accessions and cell values are
#'   dropped and only the block structure is drawn.
#' @param name_max Accessions longer than this are elided in the middle, so a
#'   long naming scheme cannot crowd out the plot.
#' @return A ggplot object. Needs the `ggplot2` package.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   Z <- matrix(rnorm(200), 50, 4, dimnames = list(paste0("p", 1:50), NULL))
#'   Z <- Z / sqrt(rowSums(Z^2))
#'   D <- repdist_matrix(Z)
#'   plot_bin_similarity(bin_proteins(D, min_sim = 0.5), "bin1", D)
#' }
#' @export
plot_bin_similarity <- function(bins, bin, D, label_max = 20L,
                                name_max = 25L) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("plot_bin_similarity() needs the ggplot2 package.", call. = FALSE)
  if (!is.list(bins$clusters) || !inherits(bins$tree, "hclust"))
    stop("`bins` must be a bin_proteins() result.", call. = FALSE)
  if (!is.character(bin) || length(bin) != 1L ||
      is.na(match(bin, names(bins$clusters))))
    stop("`bin` must name one of the ", length(bins$clusters),
         " bins in `bins$clusters`, e.g. \"bin1\".", call. = FALSE)
  if (!inherits(D, "dist")) D <- stats::as.dist(D)
  if (identical(attr(D, "method"), "euclidean"))
    stop("`D` must be on the 1 - similarity scale, not euclidean; ",
         "use repdist_matrix(Z).", call. = FALSE)

  members <- bins$clusters[[bin]]
  n <- length(members)
  if (n < 2L)
    stop("`", bin, "` holds one protein; there is no pair to plot.",
         call. = FALSE)
  members <- members[order(match(members, bins$tree$labels[bins$tree$order]))]

  S <- 1 - as.matrix(usedist::dist_subset(D, members))
  # A discrete y axis starts at the bottom, so `partner`'s levels are reversed
  # to put the first member at the top: the diagonal then falls from the upper
  # left to the lower right, the usual orientation for a matrix.
  df <- data.frame(
    protein = factor(rep(members, times = n), members),
    partner = factor(rep(members, each = n), rev(members)),
    similarity = as.numeric(S))

  # Accessions are long, and the matrix is symmetric, so one axis carries them,
  # elided in the middle: both ends are informative in an assembly-derived name.
  elide <- function(x) ifelse(nchar(x) <= name_max, x, paste0(
    substr(x, 1L, ceiling(name_max / 2) - 1L), "\u2026",
    substring(x, nchar(x) - floor(name_max / 2) + 2L)))

  p <- ggplot2::ggplot(df, ggplot2::aes(protein, partner, fill = similarity)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(limits = c(min(bins$min_sim, S), 1)) +
    ggplot2::scale_y_discrete(labels = elide) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, title = sprintf(
      "%s: %d proteins, lowest pairwise similarity %.3f", bin, n, min(S))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  if (n > label_max)
    return(p + ggplot2::theme(axis.text.y = ggplot2::element_blank()))
  # viridis runs dark to light, so no one text colour reads at both ends
  p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", similarity),
                     color = similarity > mean(range(similarity))),
        size = 2.6, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "white", `TRUE` = "black"))
}
