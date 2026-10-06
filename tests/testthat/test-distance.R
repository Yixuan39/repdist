set.seed(0)

unit <- function(W) W / sqrt(rowSums(W^2))
# three tight clusters of unit-norm embeddings
Z <- unit(rbind(
    matrix(rnorm(10 * 4, 0, 0.05), 10, 4) + rep(c(1, 0, 0, 0), each = 10),
    matrix(rnorm(10 * 4, 0, 0.05), 10, 4) + rep(c(0, 1, 0, 0), each = 10),
    matrix(rnorm(10 * 4, 0, 0.05), 10, 4) + rep(c(0, 0, 1, 1), each = 10)
))
rownames(Z) <- paste0("p", 1:30)
x <- as_proteins(Z)
counts <- matrix(
    rpois(6 * 30, 20) + 1L, 6, 30,
    dimnames = list(paste0("s", 1:6), rownames(Z))
)

test_that("sample_mmd takes the proteins, not a matrix or dist", {
    expect_error(sample_mmd(counts, Z), "AAStringSet")
    expect_error(sample_mmd(counts, stats::dist(Z)), "AAStringSet")
})

test_that("sample_mmd matches a dense reference", {
    G <- as.matrix(stats::dist(Z))
    sigma <- stats::median(G[upper.tri(G)])
    P <- counts / rowSums(counts)
    PKP <- P %*% exp(-G^2 / (2 * sigma^2)) %*% t(P)
    M2 <- outer(diag(PKP), diag(PKP), "+") - 2 * PKP
    D <- sample_mmd(counts, x)
    expect_equal(attr(D, "sigma"), sigma)
    expect_equal(as.numeric(D), as.numeric(stats::as.dist(sqrt(pmax(M2, 0)))))
})

test_that("embeddings are normalised to unit rows", {
    expect_equal(sample_mmd(counts, as_proteins(3 * Z)),
        sample_mmd(counts, x))
})

test_that("sample_mmd preserves sparse counts", {
    sparse_counts <- counts
    sparse_counts[(row(sparse_counts) + col(sparse_counts)) %% 4 != 0] <- 0L
    sparse <- Matrix::Matrix(sparse_counts, sparse = TRUE)
    for (weighted in c(TRUE, FALSE)) {
        expect_equal(sample_mmd(sparse, x, weighted = weighted),
            sample_mmd(sparse_counts, x, weighted = weighted))
    }
})

test_that("embedding rows are aligned to the retained proteins", {
    cc <- counts[1:3, 1:7]
    cc[, 7] <- 0L
    for (weighted in c(TRUE, FALSE)) {
        expect_equal(sample_mmd(cc, x[30:1], weighted = weighted),
            sample_mmd(cc[, 1:6], x[1:6], weighted = weighted))
    }
})

test_that("sample_mmd returns a dist that adonis2 accepts", {
    skip_if_not_installed("vegan")
    D <- sample_mmd(counts, x)
    expect_s3_class(D, "dist")
    expect_equal(attr(D, "Size"), nrow(counts))
    expect_equal(labels(D), rownames(counts))

    md <- data.frame(group = rep(c("A", "B"), each = 3))
    a <- vegan::adonis2(D ~ group, data = md, permutations = 99)
    expect_true(is.finite(a$F[1]))
})

test_that("RBF-MMD is symmetric, non-negative, and zero on the diagonal", {
    for (weighted in c(TRUE, FALSE)) {
        D <- as.matrix(sample_mmd(counts, x, weighted = weighted))
        expect_true(isSymmetric(D))
        expect_true(all(D >= 0))
        expect_equal(diag(D), setNames(rep(0, nrow(D)), rownames(D)))
    }
})

test_that("identical composition gives zero distance", {
    cc <- rbind(a = c(10, 0, 0), b = c(10, 0, 0), d = c(0, 10, 0))
    colnames(cc) <- rownames(Z)[c(1, 11, 21)]
    D <- as.matrix(sample_mmd(cc, x))
    expect_equal(D["a", "b"], 0)
    expect_gt(D["a", "d"], 0)
})

test_that("weighted = FALSE uses support rather than abundance", {
    xb <- as_proteins(rbind(p1 = c(1, 0), p2 = c(0, 1), p3 = c(-1, 0)))
    cc <- rbind(a = c(98, 1, 1), b = c(1, 1, 98))
    colnames(cc) <- names(xb)

    expect_gt(as.numeric(sample_mmd(cc, xb)), 0.1)
    expect_equal(as.numeric(sample_mmd(cc, xb, weighted = FALSE)), 0)
})

test_that("weighted = FALSE is exactly sign(abundance)", {
    cc <- rbind(a = c(90, 5, 5, 0), b = c(1, 1, 1, 97))
    colnames(cc) <- rownames(Z)[c(1, 11, 21, 2)]

    unweighted <- sample_mmd(cc, x, weighted = FALSE)
    expect_equal(unweighted, sample_mmd(sign(cc), x))
    expect_false(isTRUE(all.equal(
        as.numeric(unweighted), as.numeric(sample_mmd(cc, x))
    )))
})

test_that("RBF-MMD sees mass rearranged at a fixed weighted mean", {
    xb <- as_proteins(rbind(
        p1 = c(1, 0), p2 = c(-1, 0), p3 = c(0, 1), p4 = c(0, -1)
    ))
    cc <- rbind(a = c(50, 50, 0, 0), b = c(0, 0, 50, 50))
    colnames(cc) <- names(xb)
    # both samples have mean embedding zero
    expect_gt(as.numeric(sample_mmd(cc, xb)), 0)
})

test_that("missing or unnamed proteins fail loudly", {
    expect_error(sample_mmd(counts, x[1:5]))
    expect_error(sample_mmd(counts, unname(x)), "names")
    expect_error(sample_mmd(counts, as_proteins(rbind(Z, Z[1, ]))),
        "names")
})

test_that("sample_mmd rejects malformed counts", {
    bad_names <- counts
    colnames(bad_names) <- NULL
    expect_error(sample_mmd(bad_names, x))

    negative <- counts
    negative[1, 1] <- -1
    expect_error(sample_mmd(negative, x))

    nonfinite <- counts
    nonfinite[1, 1] <- NA
    expect_error(sample_mmd(nonfinite, x))

    empty_row <- counts
    empty_row[1, ] <- 0
    expect_error(sample_mmd(empty_row, x))
    expect_error(sample_mmd(counts, x, weighted = NA))
})

test_that("globally zero-count proteins are dropped before the bandwidth", {
    extra <- as_proteins(rbind(Z, junk = c(0, 0, 0, 1)))
    cc <- cbind(counts, junk = 0L)
    expect_equal(sample_mmd(cc, extra), sample_mmd(counts, x))
})

test_that("one retained protein gives zero sample distance", {
    cc <- rbind(a = 10, b = 20)
    colnames(cc) <- "p1"
    expect_equal(as.numeric(sample_mmd(cc, x)), 0)
})

test_that("a degenerate bandwidth is an error, not a NaN kernel", {
    same <- Z[rep(1, nrow(Z)), ] # every protein identical
    rownames(same) <- rownames(Z)
    expect_error(sample_mmd(counts, as_proteins(same)),
        "identical embeddings")
    near <- same
    near[, 2] <- near[, 2] + seq_len(nrow(near)) * 1e-6
    expect_true(all(is.finite(sample_mmd(counts, as_proteins(near)))))
    expect_error(sample_mmd(counts, x, sigma = 0), "positive")
    expect_error(sample_mmd(counts, x, sigma = -1), "positive")
    for (bad in list("1", TRUE, NaN, Inf, c(1, 2), NA)) {
        expect_error(sample_mmd(counts, x, sigma = bad), "sigma")
    }
})

test_that("sample and protein IDs are valid before alignment", {
    for (bad in list(c("x", "x"), c(NA, "y"), c(" ", "y"))) {
        cc <- counts[1:2, , drop = FALSE]
        rownames(cc) <- bad
        expect_error(sample_mmd(cc, x), "names")
        cc <- counts[, 1:2, drop = FALSE]
        colnames(cc) <- bad
        expect_error(sample_mmd(cc, x), "names")
    }
})

test_that("sample_mmd reports its bandwidth and holds it when passed", {
    D <- sample_mmd(counts, x)
    expect_equal(sample_mmd(counts, x, sigma = attr(D, "sigma")), D)

    # The default moves with the retained proteins; an explicit sigma does not.
    cc <- counts
    cc[, 11:30] <- 0 # six samples on one cluster
    cc7 <- rbind(cc, s7 = c(rep(0, 10), rep(5, 20))) # a seventh adds the rest
    d6 <- as.matrix(sample_mmd(cc, x))
    expect_false(isTRUE(all.equal(
        d6, as.matrix(sample_mmd(cc7, x))[1:6, 1:6]
    )))
    sigma6 <- attr(sample_mmd(cc, x), "sigma")
    fixed <- sample_mmd(cc7, x, sigma = sigma6)
    expect_equal(as.matrix(fixed)[1:6, 1:6], d6)
})
