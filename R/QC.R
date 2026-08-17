utils::globalVariables(c("curve", "pct"))

#' Proportion of proteins, and of abundance, longer than each length
#'
#' Two survival curves on a shared length axis: the fraction of the protein
#' catalog, and the fraction of total abundance, longer than each x. Reading
#' both at a candidate cutoff -- an embedding model's training range, say --
#' shows how much catalog and abundance lie beyond it.
#'
#' `counts` is converted to per-sample relative abundance first, so the
#' abundance curve reflects length rather than sequencing depth. Rarefied input
#' is not required.
#'
#' @param counts Matrix of counts or relative abundances, samples in rows,
#'   proteins in columns. Every row must have a positive total.
#' @param lengths Named numeric vector of protein lengths, or a named character
#'   vector of sequences (lengths via `nchar()`); names must cover
#'   `colnames(counts)`.
#' @param max_length Optional length to mark with a vertical reference line.
#' @return A ggplot object. Needs the `ggplot2` package.
#' @examples
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   set.seed(1)
#'   lengths <- stats::setNames(round(exp(rnorm(100, log(300), 0.6))) + 20,
#'                              paste0("p", 1:100))
#'   counts <- matrix(rpois(400, 5) + 1L, 4, 100,
#'                    dimnames = list(paste0("s", 1:4), names(lengths)))
#'   length_abundance_profile(counts, lengths, max_length = 1000)
#' }
#' @export
length_abundance_profile <- function(counts, lengths, max_length = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("length_abundance_profile() needs the ggplot2 package.", call. = FALSE)
  counts <- as.matrix(counts)
  if (is.character(lengths)) lengths <- nchar(lengths)
  stopifnot(all(colnames(counts) %in% names(lengths)), all(rowSums(counts) > 0))

  len <- lengths[colnames(counts)]
  share <- colSums(counts / rowSums(counts))
  share <- share / sum(share)
  grid <- sort(unique(len))
  pct_protein <- vapply(grid, function(x) mean(len > x), numeric(1))
  pct_abund <- vapply(grid, function(x) sum(share[len > x]), numeric(1))

  df <- data.frame(
    length = rep(grid, 2),
    pct = 100 * c(pct_protein, pct_abund),
    curve = rep(c("proteins", "abundance"), each = length(grid)))

  xmax <- max(grid[pct_protein > 0.01 | pct_abund > 0.01])
  xmax <- min(max(grid), xmax * 1.1)
  p <- ggplot2::ggplot(df, ggplot2::aes(length, pct, color = curve)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::coord_cartesian(xlim = c(min(grid), xmax)) +
    ggplot2::scale_color_manual(
      name = NULL, values = c(proteins = "steelblue", abundance = "firebrick")) +
    ggplot2::labs(x = "protein length (aa)", y = "% longer than this length") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "top")

  if (!is.null(max_length))
    p <- p + ggplot2::geom_vline(
      xintercept = max_length, linetype = "dashed", color = "grey30")
  p
}
