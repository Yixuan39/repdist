set.seed(0)

base <- rbind(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1))
Z <- do.call(rbind, lapply(1:3, function(k) {
    matrix(base[k, ], 5, 3, byrow = TRUE) + matrix(rnorm(15, 0, 0.03), 5, 3)
}))
Z <- Z / sqrt(rowSums(Z^2))
rownames(Z) <- paste0("p", 1:15)
family <- rep(1:3, each = 5)
S <- protein_similarity(as_proteins(Z), min_sim = 0.5)

unit <- function(W) W / sqrt(rowSums(W^2))

# Plain dense MCL with the same settings and no pruning: the reference
# bin_proteins() must agree with on small graphs.
dense_mcl <- function(S, inflation) {
    A <- S
    diag(A) <- 0
    linked <- rowSums(A > 0) > 0
    A <- A[linked, linked, drop = FALSE]
    diag(A) <- apply(A, 2, max)
    M <- sweep(A, 2, colSums(A), "/")
    for (k in 1:200) {
        E <- (M %*% M)^inflation
        E <- sweep(E, 2, colSums(E), "/")
        if (max(abs(E - M)) < 1e-12) break
        M <- E
    }
    to <- apply(M, 2, which.max)
    membership <- repdist:::.components(seq_along(to), to, length(to))
    out <- stats::setNames(rep(NA_integer_, nrow(S)), rownames(S))
    out[linked] <- membership
    out[!linked] <- max(membership, 0L) + seq_len(sum(!linked))
    out
}

same_partition <- function(a, b) {
    identical(outer(a, a, "=="), outer(b, b, "=="))
}

test_that("bin_proteins returns MCL clusters, membership, medoids, graph", {
    out <- bin_proteins(S)
    expect_named(out, c(
        "clusters", "membership", "representatives", "inflation", "graph"
    ))
    expect_equal(out$inflation, 2)
    expect_named(out$graph, c(
        "proteins", "edges", "isolated", "components", "largest_component"
    ))
    expect_false(is.unsorted(rev(lengths(out$clusters))))
    expect_equal(unname(lengths(out$clusters)), c(5L, 5L, 5L))
    for (fam in 1:3) {
        expect_length(unique(out$membership[family == fam]), 1L)
    }
})

test_that("every protein appears in exactly one bin", {
    set.seed(5)
    W <- unit(matrix(rnorm(80 * 6), 80, 6))
    rownames(W) <- paste0("q", 1:80)
    for (t in c(0.3, 0.6, 0.9)) {
        out <- bin_proteins(protein_similarity(as_proteins(W), min_sim = t))
        members <- unlist(out$clusters, use.names = FALSE)
        expect_setequal(members, rownames(W))
        expect_false(anyDuplicated(members) > 0)
        expect_identical(names(out$membership), rownames(W))
        expect_identical(
            out$membership,
            stats::setNames(
                rep(names(out$clusters), lengths(out$clusters)), members
            )[rownames(W)]
        )
    }
})

test_that("sparse MCL agrees with a dense reference", {
    set.seed(8)
    centres <- unit(matrix(rnorm(4 * 8), 4, 8))
    W <- unit(centres[rep(1:4, each = 12), ] + matrix(rnorm(48 * 8, 0, 0.35),
        48, 8))
    rownames(W) <- paste0("q", 1:48)
    for (t in c(0.3, 0.5)) {
        G <- protein_similarity(as_proteins(W), min_sim = t)
        for (inflation in c(1.4, 2, 4)) {
            out <- bin_proteins(G, inflation = inflation)
            expect_true(same_partition(
                out$membership, dense_mcl(G, inflation)[rownames(W)]
            ))
        }
    }
})

test_that("isolated proteins stay separate singleton bins", {
    W <- unit(rbind(
        a1 = c(1, 0, 0), a2 = c(0.99, 0.14, 0), b = c(0, 1, 0), c = c(0, 0, 1)
    ))
    out <- bin_proteins(protein_similarity(as_proteins(W), min_sim = 0.9))
    expect_setequal(
        lapply(out$clusters, sort),
        list(c("a1", "a2"), "b", "c")
    )
    expect_equal(out$graph[["edges"]], 1)
    expect_equal(out$graph[["isolated"]], 2)
    expect_equal(out$graph[["components"]], 3)
})

test_that("MCL never merges disconnected components", {
    # two tight groups per component, the components far apart; a low
    # inflation merges within a component but can never cross between them
    set.seed(3)
    centres <- unit(rbind(
        c(1, 0.3, 0, 0), c(1, -0.3, 0, 0), c(0, 0, 1, 0.3), c(0, 0, 1, -0.3)
    ))
    W <- unit(centres[rep(1:4, each = 6), ] + matrix(rnorm(96, 0, 0.02), 24))
    rownames(W) <- paste0("q", 1:24)
    side <- rep(1:2, each = 12)
    G <- protein_similarity(as_proteins(W), min_sim = 0.5)
    for (inflation in c(1.1, 1.4, 2, 6)) {
        out <- bin_proteins(G, inflation = inflation)
        expect_equal(out$graph[["components"]], 2)
        for (members in out$clusters) {
            expect_length(unique(side[match(members, rownames(W))]), 1L)
        }
    }
})

test_that("a higher inflation splits what a lower one merges", {
    # two groups of eight joined by every cross pair at similarity ~0.6
    set.seed(1)
    W <- unit(rbind(
        matrix(rep(c(1, 0, 0), 8), 8, byrow = TRUE),
        matrix(rep(c(0.6, 0.8, 0), 8), 8, byrow = TRUE)
    ) + matrix(rnorm(48, 0, 0.02), 16))
    rownames(W) <- paste0("q", 1:16)
    G <- protein_similarity(as_proteins(W), min_sim = 0.5)
    count <- function(inflation) length(bin_proteins(G, inflation)$clusters)
    expect_equal(count(1.2), 1)
    expect_equal(count(6), 2)
})

test_that("clustering is deterministic", {
    set.seed(6)
    W <- unit(matrix(rnorm(120 * 5), 120, 5))
    rownames(W) <- paste0("q", 1:120)
    G <- protein_similarity(as_proteins(W), min_sim = 0.4)
    expect_identical(bin_proteins(G), bin_proteins(G))
})

test_that("representatives are members and within-bin medoids", {
    set.seed(7)
    W <- unit(matrix(rnorm(60 * 4), 60, 4))
    rownames(W) <- paste0("q", 1:60)
    G <- protein_similarity(as_proteins(W), min_sim = 0.6)
    out <- bin_proteins(G)
    expect_named(out$representatives, names(out$clusters))
    for (bin in names(out$clusters)) {
        members <- out$clusters[[bin]]
        expect_true(out$representatives[[bin]] %in% members)
        expect_equal(
            out$representatives[[bin]],
            members[which.max(colSums(G[members, members, drop = FALSE]))]
        )
    }
})

test_that("bin_proteins agrees with the mcl program", {
    skip_if(!nzchar(Sys.which("mcl")), "mcl is not installed")
    set.seed(9)
    centres <- unit(matrix(rnorm(6 * 10), 6, 10))
    W <- unit(centres[rep(1:6, each = 20), ] + matrix(rnorm(120 * 10, 0, 0.4),
        120, 10))
    rownames(W) <- paste0("q", 1:120)
    G <- protein_similarity(as_proteins(W), min_sim = 0.5)
    e <- which(upper.tri(G) & G > 0, arr.ind = TRUE)
    abc <- tempfile(fileext = ".abc")
    out <- tempfile()
    utils::write.table(
        data.frame(rownames(W)[e[, 1]], rownames(W)[e[, 2]], G[e]), abc,
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    for (inflation in c(1.4, 2, 4)) {
        system2("mcl", c(abc, "--abc", "-I", inflation, "-o", out),
            stdout = FALSE, stderr = FALSE
        )
        cl <- strsplit(readLines(out), "\t")
        ref <- stats::setNames(rep(seq_along(cl), lengths(cl)), unlist(cl))
        ours <- bin_proteins(G, inflation)$membership[names(ref)]
        expect_true(same_partition(ours, ref))
    }
})

test_that("an MCL failure is an error, not a message", {
    # the MCL package returns this string instead of signalling
    testthat::local_mocked_bindings(
        mcl = function(...) "An Error occurred at iteration 100",
        .package = "MCL"
    )
    expect_error(bin_proteins(S), "did not converge")
})

test_that("protein_similarity is cosine similarity gated at min_sim", {
    set.seed(14)
    W <- matrix(rnorm(50 * 4), 50, 4, dimnames = list(paste0("p", 1:50), NULL))
    W <- W / sqrt(rowSums(W^2))
    cosine <- tcrossprod(W)
    off <- upper.tri(cosine)
    for (t in c(0, 0.3, 0.7)) {
        G <- protein_similarity(as_proteins(W), min_sim = t)
        expect_true(isSymmetric(G))
        expect_identical(dimnames(G), list(rownames(W), rownames(W)))
        expect_equal(unname(diag(G)), rep(1, 50))
        expect_equal(G[off], ifelse(cosine[off] >= t, cosine[off], 0))
    }
    # rows are normalised, with similarity clipped at 1
    rounded <- rbind(a = c(1 + 1e-7, 0), b = c(3, 0), c = c(-1, 0))
    expect_equal(unname(protein_similarity(as_proteins(rounded))),
        rbind(c(1, 1, 0), c(1, 1, 0), c(0, 0, 1)))
    expect_error(protein_similarity(W), "AAStringSet")
    expect_error(protein_similarity(as_proteins(unname(W))), "names")
    expect_error(protein_similarity(as_proteins(rbind(W[1:3, ], dead = 0))),
        "is.finite")
    for (bad in list(-0.1, 1.5, NA, c(0.5, 0.7), "0.5")) {
        expect_error(protein_similarity(as_proteins(W), min_sim = bad), "min_sim")
        expect_error(protein_similarity(as_proteins(W), min_coverage = bad),
            "min_coverage")
        expect_error(protein_similarity(as_proteins(W), n_cores = bad),
            "n_cores")
    }
})

test_that("coverage is the share of the less-covered sequence aligned", {
    aa <- Biostrings::AAStringSet(c(long = "MKTAYIAKQRQISFVKSHFSRQ",
        short = "AYIAKQRQIS"))
    ij <- rbind(c(1L, 2L), c(1L, 1L))
    expect_equal(repdist:::.coverage(1:2, ij, aa, "BLOSUM62"), c(10 / 22, 1))
})

test_that("the coverage gate needs both sequences half covered", {
    W <- rbind(a = c(1, 0), b = c(0.95, 0.05), c = c(0.9, 0.1), d = c(0, 1))
    seqs <- c(a = "MKTAYIAKQRQISFVKSHFSRQ", b = "MKTAYIAKQRQISFVKSHFSRA",
        c = "AYIAKQRQIS", d = "WWWWWW")
    x <- as_proteins(W, seqs)
    plain <- protein_similarity(x, min_coverage = 0)
    gated <- protein_similarity(x)
    ab <- c("a", "b")
    expect_equal(gated[ab, ab], plain[ab, ab])
    # c is wholly covered, but covers less than half of a and b
    expect_equal(gated[ab, "c"], c(a = 0, b = 0))
    expect_gt(protein_similarity(x, min_coverage = 0.4)["a", "c"], 0)
    expect_true(isSymmetric(gated))
    # only pairs that pass min_sim are aligned
    aligned <- 0L
    testthat::local_mocked_bindings(.coverage = function(k, ...) {
        aligned <<- aligned + length(k)
        rep(1, length(k))
    })
    expect_equal(protein_similarity(x), plain)
    expect_equal(aligned, sum(plain[upper.tri(plain)] > 0))
})

test_that("parallel and serial coverage gates agree", {
    skip_on_os("windows")
    set.seed(2)
    W <- matrix(rnorm(12 * 3), 12, 3, dimnames = list(paste0("p", 1:12), NULL))
    seqs <- vapply(1:12, function(i) {
        paste(sample(Biostrings::AA_STANDARD, 20 + 5 * i, TRUE), collapse = "")
    }, "")
    names(seqs) <- rownames(W)
    x <- as_proteins(abs(W), seqs)
    expect_identical(
        protein_similarity(x, 0.5, 0.3, n_cores = 2),
        protein_similarity(x, 0.5, 0.3)
    )
})

test_that("bin_glom sums members and conserves totals", {
    out <- bin_proteins(S)
    set.seed(4)
    counts <- matrix(rpois(4 * 15, 20), 4, 15,
        dimnames = list(paste0("s", 1:4), sample(rownames(Z)))
    )
    agg <- bin_glom(counts, out)
    expect_equal(rowSums(agg), rowSums(counts))
    expect_identical(colnames(agg), names(out$clusters))
    for (bin in colnames(agg)) {
        expect_equal(
            agg[, bin], rowSums(counts[, out$clusters[[bin]], drop = FALSE])
        )
    }
    sparse <- bin_glom(Matrix::Matrix(counts, sparse = TRUE), out)
    expect_s4_class(sparse, "sparseMatrix")
    expect_equal(as.matrix(sparse), agg)
    # a subset of proteins keeps only the bins it touches
    part <- bin_glom(counts[, out$clusters$bin2], out)
    expect_identical(colnames(part), "bin2")
    expect_error(bin_glom(cbind(counts, extra = 1), out), "in no bin")
    expect_error(bin_glom(counts, list(clusters = 1)), "bin_proteins")
})

test_that("invalid inflation and similarity matrices are rejected", {
    for (bad in list(1, 0.5, c(1.5, 2), Inf, NA, "2")) {
        expect_error(bin_proteins(S, inflation = bad))
    }
    bad <- S
    bad[1, 2] <- 0.3
    expect_error(bin_proteins(bad), "symmetric")
    expect_error(bin_proteins(S[, -1]), "symmetric")
    bad <- S
    bad[1, 2] <- bad[2, 1] <- NA
    expect_error(bin_proteins(bad), "symmetric")
    expect_error(bin_proteins(2 * S), "\\[0, 1\\]")
    expect_error(bin_proteins(-S), "\\[0, 1\\]")
    expect_error(bin_proteins(unname(S)), "names")
})
