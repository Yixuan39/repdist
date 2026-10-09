utils::globalVariables(c(
    "similarity", "protein", "partner", "tsne1", "tsne2",
    "x_from", "y_from", "x_to", "y_to"
))

#' Plot protein embeddings with t-SNE
#'
#' Two-dimensional t-SNE of the protein embeddings on cosine distance, one
#' point per protein. For TM-Vec embeddings, cosine distance is one minus the
#' predicted TM-score. Read it for neighbourhoods: proteins that sit together
#' are close in embedding space, but t-SNE does not preserve the distances
#' between clusters or their sizes. The layout is stochastic; call
#' [set.seed()] first to reproduce it.
#'
#' @param x Proteins with embeddings, from [embed_proteins()] or
#'   [read_embeddings()].
#' @param perplexity t-SNE perplexity. Defaults to 30, or the largest value
#'   t-SNE accepts when there are fewer than 94 proteins.
#' @param color Optional labels to colour by, such as a functional class: a
#'   vector named by protein. The 12 most frequent labels get colours, the
#'   rest are drawn as `"other"`, and proteins without a label are grey.
#' @return A ggplot object whose `data` holds the t-SNE coordinates (and the
#'   `label` used for colour). Needs the `ggplot2` and `Rtsne` packages.
#' @examples
#' d <- readRDS(system.file("extdata", "diet_metaproteome.rds",
#'     package = "repdist"
#' ))
#' x <- d$proteins[1:100]
#' if (requireNamespace("ggplot2", quietly = TRUE) &&
#'     requireNamespace("Rtsne", quietly = TRUE)) {
#'     set.seed(1)
#'     plot_tsne(x)
#' }
#' @export
plot_tsne <- function(x, perplexity = NULL, color = NULL) {
    Z <- .embeddings(x)
    n <- nrow(Z)
    if (!is.null(color) && is.null(names(color))) {
        stop("`color` must be named by protein.", call. = FALSE)
    }
    if (n < 4L) {
        stop("t-SNE needs at least 4 proteins.", call. = FALSE)
    }
    if (is.null(perplexity)) perplexity <- min(30, floor((n - 1) / 3))
    # on unit-norm rows, 1 - cosine == squared Euclidean distance / 2
    fit <- Rtsne::Rtsne(stats::dist(Z)^2 / 2, is_distance = TRUE,
        perplexity = perplexity)
    df <- data.frame(protein = rownames(Z), tsne1 = fit$Y[, 1],
        tsne2 = fit$Y[, 2])
    size <- if (n > 2000L) 0.4 else 1.5 # a full catalog would be one blob
    p <- ggplot2::ggplot(df, ggplot2::aes(tsne1, tsne2)) +
        ggplot2::labs(x = "t-SNE 1", y = "t-SNE 2",
            subtitle = sprintf("cosine distance, perplexity %g", perplexity)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text = ggplot2::element_blank())
    if (is.null(color)) {
        return(p + ggplot2::geom_point(alpha = 0.6, size = size))
    }
    label <- as.character(color[df$protein])
    top <- utils::head(names(sort(table(label), decreasing = TRUE)), 12L)
    label[!is.na(label) & !label %in% top] <- "other"
    df$label <- factor(label, c(top, "other"))
    # coloured points on top of the grey ones
    p$data <- df[order(!is.na(df$label), df$label != "other"), ]
    p + ggplot2::geom_point(ggplot2::aes(colour = label), alpha = 0.7,
        size = size) +
        ggplot2::scale_colour_manual(values = c(stats::setNames(
            grDevices::hcl.colors(length(top), "Dark 3"), top),
            other = "grey55"), na.value = "grey88", drop = FALSE,
            labels = function(l) ifelse(is.na(l), "unlabelled", l)) +
        ggplot2::guides(colour = ggplot2::guide_legend(
            override.aes = list(size = 2.5, alpha = 1))) +
        ggplot2::labs(colour = NULL)
}

#' Plot the similarity matrix inside one bin
#'
#' Heatmap of every pairwise similarity among one bin's members, the
#' representative first and the rest by decreasing similarity to it. Colour
#' is the similarity before any gate, from the embeddings; a cross marks each
#' pair that is not an edge of `S`, the graph MCL clustered. An MCL bin joins
#' members through shared neighbours, so a pair inside it need not be an edge,
#' and a pair above `min_sim` is not one if it failed the coverage gate.
#'
#' @param bins Result of [bin_proteins()].
#' @param bin Bin name, e.g. `"bin1"`.
#' @param S The similarity matrix given to [bin_proteins()], gates applied.
#' @param x The proteins `S` was computed from; the colours come from their
#'   embeddings, so no ungated matrix is needed.
#' @param label_max Above this many members, accessions and cell values are
#'   dropped.
#' @param name_max Accessions longer than this are elided in the middle.
#' @return A ggplot object whose `data` holds each pair's `similarity` and
#'   whether it is an `edge`. Needs the `ggplot2` package.
#' @examples
#' d <- readRDS(system.file("extdata", "diet_metaproteome.rds",
#'     package = "repdist"
#' ))
#' x <- d$proteins[1:100]
#' S <- protein_similarity(x)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'     plot_bin_similarity(bin_proteins(S), "bin1", S, x)
#' }
#' @export
plot_bin_similarity <- function(bins, bin, S, x, label_max = 20L,
    name_max = 25L) {
    S <- .check_sim(S)
    if (!is.list(bins$clusters) || is.null(bins$representatives)) {
        stop("`bins` must be a bin_proteins() result.", call. = FALSE)
    }
    if (length(bin) != 1L || !bin %in% names(bins$clusters)) {
        stop("`bin` must name one bin, e.g. \"bin1\".", call. = FALSE)
    }
    # an ungated matrix here would make every pair an edge
    n_edges <- sum(S[upper.tri(S)] > 0)
    if (nrow(S) != bins$graph[["proteins"]] ||
        n_edges != bins$graph[["edges"]]) {
        stop("`S` has ", n_edges, " edges but `bins` was built from ",
            bins$graph[["edges"]], "; pass the matrix given to bin_proteins().",
            call. = FALSE)
    }
    members <- bins$clusters[[bin]]
    if (length(members) < 2L) {
        stop("`", bin, "` holds one protein; no pair to plot.", call. = FALSE)
    }
    Z <- .embeddings(x)
    if (!all(members %in% rownames(Z))) {
        stop("`x` must hold every member of `", bin, "`.", call. = FALSE)
    }
    # as protein_similarity(x, min_sim = 0, min_coverage = 0) gives it
    R <- pmax(pmin(tcrossprod(Z[members, , drop = FALSE]), 1), 0)
    medoid <- bins$representatives[[bin]]
    members <- members[order(members != medoid, -R[medoid, members])]
    R <- R[members, members]
    E <- S[members, members] > 0
    diag(E) <- FALSE
    n <- length(members)
    # y levels reversed so the diagonal runs from the upper left
    df <- data.frame(
        protein = factor(rep(members, times = n), members),
        partner = factor(rep(members, each = n), rev(members)),
        similarity = as.numeric(R),
        edge = as.vector(E)
    )
    # each pair that is not an edge is crossed corner to corner
    self <- as.character(df$protein) == as.character(df$partner)
    off <- df[!df$edge & !self, ]
    at_x <- as.numeric(off$protein)
    at_y <- as.numeric(off$partner)
    h <- 0.42
    cross <- data.frame(x_from = at_x - h, x_to = at_x + h,
        y_from = c(at_y - h, at_y + h), y_to = c(at_y + h, at_y - h))
    half <- name_max %/% 2
    elide <- function(x) ifelse(nchar(x) <= name_max, x, paste0(
        substr(x, 1L, half), "\u2026", substring(x, nchar(x) - half + 2L)))
    p <- ggplot2::ggplot(df,
        ggplot2::aes(protein, partner, fill = similarity)) +
        ggplot2::geom_tile() +
        ggplot2::geom_segment(data = cross, ggplot2::aes(x = x_from,
            y = y_from, xend = x_to, yend = y_to), inherit.aes = FALSE,
            linewidth = if (n > label_max) 0.3 else 0.6) +
        ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
        ggplot2::scale_y_discrete(labels = elide) +
        ggplot2::coord_fixed() +
        ggplot2::labs(x = NULL, y = NULL,
            title = sprintf("%s: %d proteins", bin, n),
            subtitle = sprintf("%.0f%% of pairs are edges",
                100 * mean(E[upper.tri(E)]))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())
    if (n > label_max) {
        return(p + ggplot2::theme(axis.text.y = ggplot2::element_blank()))
    }
    # viridis runs dark to light, so no one text colour reads at both ends
    p + ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", similarity),
        color = similarity > 0.5), size = 2.6, show.legend = FALSE) +
        ggplot2::scale_color_manual(
            values = c(`FALSE` = "white", `TRUE` = "black"))
}
