unit <- function(W) W / sqrt(rowSums(W^2))

set.seed(3)
Z <- unit(matrix(rnorm(90), 30, 3, dimnames = list(paste0("p", 1:30), NULL)))

test_that("plot_bin_similarity shows one bin's pairwise similarities", {
    skip_if_not_installed("ggplot2")
    W <- unit(rbind(
        p1 = c(1, 0), p2 = c(0.99, 0.14), p3 = c(0.97, 0.24), p4 = c(0, 1)
    ))
    x <- as_proteins(W)
    G <- protein_similarity(x, min_sim = 0.9)
    b <- bin_proteins(G)
    p <- plot_bin_similarity(b, "bin1", G, x)
    expect_s3_class(p, "ggplot")
    members <- b$clusters$bin1
    S <- G[members, members]
    # the whole square, every cell a similarity
    expect_equal(nrow(p$data), length(members)^2)
    expect_setequal(as.character(p$data$protein), members)
    expect_equal(sort(p$data$similarity), sort(as.numeric(S)))
    expect_s3_class(p$layers[[1]]$geom, "GeomTile")
    expect_equal(p$scales$get_scales("fill")$limits, c(0, 1))
    expect_identical(p$scales$get_scales("colour")$palette(2),
        c(`FALSE` = "white", `TRUE` = "black"))
    # the representative comes first, and y is reversed against x so the
    # diagonal falls from the upper left
    expect_identical(levels(p$data$protein)[1], b$representatives[["bin1"]])
    expect_identical(levels(p$data$partner), rev(levels(p$data$protein)))
    expect_match(p$labels$title, sprintf("bin1: %d proteins", length(members)))

    expect_error(plot_bin_similarity(b, 1L, G, x), "must name one")
    expect_error(plot_bin_similarity(b, "bin99", G, x), "must name one")
    expect_error(
        plot_bin_similarity(list(clusters = 1), "bin1", G, x), "bin_proteins"
    )
    singleton <- names(which(lengths(b$clusters) == 1L))[[1]]
    expect_error(plot_bin_similarity(b, singleton, G, x), "no pair")
    expect_error(plot_bin_similarity(b, "bin1", G, x[4]), "every member")
})

test_that("an MCL bin can hold pairs that are not edges", {
    skip_if_not_installed("ggplot2")
    # a chain: neighbours clear the floor, the two ends do not
    ang <- c(0, 0.5, 1, 1.5)
    W <- cbind(cos(ang), sin(ang))
    rownames(W) <- paste0("c", 1:4)
    x <- as_proteins(W)
    G <- protein_similarity(x, min_sim = 0.8)
    b <- bin_proteins(G, inflation = 1.4)
    expect_length(b$clusters, 1L)
    p <- plot_bin_similarity(b, "bin1", G, x)
    # every cell carries its ungated similarity, the 3 edges are outlined
    S_all <- protein_similarity(x, min_sim = 0, min_coverage = 0)
    expect_equal(p$data$similarity, as.numeric(S_all[
        levels(p$data$protein), levels(p$data$protein)]))
    expect_gt(min(p$data$similarity), 0)
    expect_equal(sum(p$data$edge), 6L)
    expect_equal(nrow(p$layers[[2]]$data), 6L)
    expect_match(p$labels$subtitle, "50% of pairs are edges \\(outlined\\)")
    # the ungated matrix is not the graph MCL clustered
    expect_error(plot_bin_similarity(b, "bin1", S_all, x),
        "6 edges but `bins` was built from 3")
})

test_that("plot_tsne embeds every protein on cosine distance", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("Rtsne")
    x <- as_proteins(Z)
    set.seed(1)
    p <- plot_tsne(x)
    expect_s3_class(p, "ggplot")
    expect_equal(p$data$protein, rownames(Z))
    expect_equal(dim(p$data[c("tsne1", "tsne2")]), c(30L, 2L))
    # 30 proteins is too few for perplexity 30; the default backs off to the
    # largest Rtsne accepts
    expect_match(p$labels$subtitle, "perplexity 9")
    # it is Rtsne on cosine distance, which on unit-norm rows is half the
    # squared Euclidean distance. (Layouts match only for the identical
    # distance object: float noise at 1e-16 reshuffles a t-SNE.)
    Zn <- .embeddings(x)
    expect_equal(as.matrix(stats::dist(Zn)^2 / 2), 1 - tcrossprod(Zn),
        ignore_attr = TRUE)
    set.seed(1)
    fit <- Rtsne::Rtsne(stats::dist(Zn)^2 / 2, is_distance = TRUE,
        perplexity = 9)
    expect_equal(as.matrix(p$data[c("tsne1", "tsne2")]), fit$Y,
        ignore_attr = TRUE, tolerance = 1e-6)
    expect_error(plot_tsne(as_proteins(Z[1:3, ])), "at least 4")
    expect_error(plot_tsne(Z), "AAStringSet")
})

test_that("plot_tsne colours the 12 most frequent labels", {
    skip_if_not_installed("ggplot2")
    skip_if_not_installed("Rtsne")
    labels <- c(rep(LETTERS[1:13], each = 2), "N", NA)
    names(labels) <- rownames(Z)[1:28] # p29 and p30 get no label
    set.seed(1)
    p <- plot_tsne(as_proteins(Z), color = labels)
    lab <- p$data$label[match(rownames(Z), p$data$protein)]
    expect_identical(levels(lab), c(LETTERS[1:12], "other"))
    # M and N fall outside the top 12; NA and unnamed proteins stay unlabelled
    expect_equal(as.character(lab[25:28]), c("other", "other", "other", NA))
    expect_true(all(is.na(lab[29:30])))
    expect_error(plot_tsne(as_proteins(Z), color = unname(labels)), "named")
})
