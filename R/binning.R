#' Protein similarity matrix from embeddings
#'
#' Cosine similarity between proteins -- for TM-Vec embeddings, the predicted
#' TM-score -- with every pair below `min_sim` set to 0: the weighted graph
#' [bin_proteins()] clusters.
#'
#' Pairs that pass `min_sim` are then aligned (local, `substitutionMatrix`,
#' gap opening 11, extension 1, via [pwalign::pairwiseAlignment()]) and set to
#' 0 unless the alignment covers at least `min_coverage` of both sequences.
#' Pairs below `min_sim` are never aligned.
#'
#' @param x Proteins with embeddings, from
#'   [embed_proteins()] or
#'   [read_embeddings()].
#' @param min_sim Similarity floor in `[0, 1]` (default 0.7).
#' @param min_coverage Coverage floor in `[0, 1]` (default 0.5); 0 skips
#'   alignment.
#' @param substitutionMatrix Substitution matrix for the alignments, by name
#'   or as a matrix; see [pwalign::pairwiseAlignment()].
#' @param n_cores Number of cores for the alignments (default 1). Above 1 the
#'   alignments run in forked processes, which Windows does not support: there
#'   they run on one core with a warning.
#' @return Symmetric numeric matrix of similarities, 1 on the diagonal and 0
#'   for every gated pair.
#' @usage protein_similarity(
#'     x,
#'     min_sim = 0.7,
#'     min_coverage = 0.5,
#'     substitutionMatrix = "BLOSUM62",
#'     n_cores = 1L
#' )
#' @examples
#' x <- Biostrings::AAStringSet(c(
#'     p1 = "MKTAYIAKQRQISFVKSHFSRQ", p2 = "AYIAKQRQIS", p3 = "MKVLAA"
#' ))
#' S4Vectors::mcols(x)$embedding <- rbind(c(1, 0), c(0.9, 0.1), c(0, 1))
#' protein_similarity(x, min_coverage = 0)
#' protein_similarity(x)
#' @export
protein_similarity <- function(x, min_sim = 0.7, min_coverage = 0.5,
    substitutionMatrix = "BLOSUM62", n_cores = 1L) {
    Z <- .embeddings(x)
    stopifnot(
        is.numeric(min_sim), length(min_sim) == 1L,
        min_sim >= 0, min_sim <= 1,
        is.numeric(min_coverage), length(min_coverage) == 1L,
        min_coverage >= 0, min_coverage <= 1,
        is.numeric(n_cores), length(n_cores) == 1L,
        n_cores >= 1, n_cores == round(n_cores)
    )
    S <- pmin(tcrossprod(Z), 1)
    S[S < min_sim] <- 0
    diag(S) <- 1
    if (min_coverage > 0) {
        S4Vectors::mcols(x) <- NULL # pwalign drops them with a warning
        ij <- which(upper.tri(S) & S > 0, arr.ind = TRUE)
        BPPARAM <- if (n_cores > 1) {
            BiocParallel::MulticoreParam(n_cores)
        } else {
            BiocParallel::SerialParam()
        }
        cov <- BiocParallel::bpvec(seq_len(nrow(ij)), .coverage, ij = ij,
            seqs = x, substitutionMatrix = substitutionMatrix,
            BPPARAM = BPPARAM)
        low <- ij[unlist(cov) < min_coverage, , drop = FALSE]
        S[rbind(low, low[, 2:1])] <- 0
    }
    S
}

# Alignment coverage of the protein pairs in rows `rows` of `ij`, a two-column
# matrix of indices into `seqs`: for each pair, the share of the less-covered
# sequence that their local alignment spans. Top level so workers are not sent
# the similarity matrix; no argument may be named `i`, which bpvec() uses.
.coverage <- function(rows, ij, seqs, substitutionMatrix) {
    p <- seqs[ij[rows, 1]]
    s <- seqs[ij[rows, 2]]
    a <- pwalign::pairwiseAlignment(p, s,
        substitutionMatrix = substitutionMatrix, type = "local",
        gapOpening = 11, gapExtension = 1)
    pmin(Biostrings::width(pwalign::pattern(a)) / Biostrings::width(p),
        Biostrings::width(pwalign::subject(a)) / Biostrings::width(s))
}

# A protein_similarity() result, or any symmetric named similarity matrix in
# [0, 1]: nonzero off-diagonal entries are the graph's edges.
.check_sim <- function(S) {
    S <- as.matrix(S)
    if (!is.numeric(S) || anyNA(S) || !isSymmetric(S) || any(S < 0 | S > 1)) {
        stop("`S` must be a symmetric similarity matrix in [0, 1], as ",
            "returned by protein_similarity().", call. = FALSE)
    }
    .check_ids(rownames(S), "Protein")
    S
}

# Connected components of the graph on nodes 1..n with edges (i, j), which
# bin_proteins() clusters one at a time and reports in `graph`. Hooks each
# root to its smallest neighbouring root, then pointer-jumps to the new roots,
# until no edge joins two roots. Returns an integer component per node,
# numbered by first appearance.
.components <- function(i, j, n) {
    parent <- seq_len(n)
    repeat {
        ri <- parent[i]
        rj <- parent[j]
        cross <- ri != rj
        if (!any(cross)) break
        lo <- pmin(ri[cross], rj[cross])
        hi <- pmax(ri[cross], rj[cross])
        o <- order(hi, lo)
        first <- !duplicated(hi[o])
        parent[hi[o][first]] <- lo[o][first]
        repeat {
            up <- parent[parent]
            if (identical(up, parent)) break
            parent <- up
        }
    }
    match(parent, unique(parent))
}

#' Cluster proteins into structural bins
#'
#' Partitions the similarity graph from
#' [protein_similarity()]
#' by Markov clustering (MCL, via the `MCL` package). Every nonzero pair is an
#' edge weighted by its similarity, so the edge floor is the matrix's
#' `min_sim`. Bins are named `bin1`, `bin2`, ... from largest to
#' smallest; proteins with no edge are singleton bins.
#'
#' The floor applies to edges, not to all pairs: a bin can hold pairs joined
#' only through shared neighbours (see
#' [plot_bin_similarity()]).
#' Collapse a counts table to bins with [bin_glom()].
#'
#' @param S Protein similarity matrix from
#'   [protein_similarity()].
#' @param inflation MCL inflation, greater than 1 (default 2). Higher values
#'   give more, finer bins.
#' @return A list:
#'   \describe{
#'     \item{`clusters`}{Named list of member accessions per bin.}
#'     \item{`membership`}{Bin of each protein, in the order of `S`.}
#'     \item{`representatives`}{Medoid of each bin: the member with the
#'       highest total similarity to the others.}
#'     \item{`inflation`}{The inflation used.}
#'     \item{`graph`}{Counts of `proteins`, `edges`, `isolated`
#'       proteins, connected `components` and the `largest_component`
#'       size.}
#'   }
#' @references van Dongen S (2000). Graph clustering by flow simulation.
#'   PhD thesis, University of Utrecht.
#' @examples
#' d <- readRDS(system.file("extdata", "diet_metaproteome.rds",
#'     package = "repdist"
#' ))
#' x <- d$proteins[1:100]
#' b <- bin_proteins(protein_similarity(x))
#' head(lengths(b$clusters))
#' @export
bin_proteins <- function(S, inflation = 2) {
    S <- .check_sim(S)
    stopifnot(is.numeric(inflation), length(inflation) == 1L,
        is.finite(inflation), inflation > 1)
    labels <- rownames(S)
    n <- length(labels)
    ij <- which(upper.tri(S) & S > 0, arr.ind = TRUE)
    comp <- .components(ij[, 1], ij[, 2], n)
    size <- tabulate(comp)

    # MCL one connected component at a time: flow never crosses components,
    # MCL::mcl() holds dense matrices, and a protein with no edge would be an
    # all-zero column. As in the mcl program, each node's self-loop takes the
    # weight of its heaviest edge; with that, the clusters match mcl's.
    mcl_id <- integer(n)
    for (idx in split(seq_len(n), comp)) {
        cl <- 1L
        if (length(idx) > 1L) {
            A <- S[idx, idx]
            diag(A) <- 0
            diag(A) <- apply(A, 2, max)
            # convergence is exact equality, which low inflations take more
            # than the package's default 100 iterations to reach
            res <- MCL::mcl(A, addLoops = FALSE, inflation = inflation,
                allow1 = TRUE, max.iter = 1000L)
            # the package returns a message, not an error, when it fails
            if (is.character(res)) {
                stop("MCL did not converge: ", res, call. = FALSE)
            }
            cl <- match(res$Cluster, unique(res$Cluster))
        }
        mcl_id[idx] <- max(mcl_id) + cl
    }
    clusters <- unname(split(labels, mcl_id))
    clusters <- clusters[order(-lengths(clusters))]
    names(clusters) <- paste0("bin", seq_along(clusters))
    membership <- stats::setNames(
        rep(names(clusters), lengths(clusters)), unlist(clusters)
    )[labels]
    list(
        clusters = clusters,
        membership = membership,
        representatives = vapply(clusters, function(m) {
            m[[which.max(colSums(S[m, m, drop = FALSE]))]]
        }, character(1)),
        inflation = inflation,
        graph = c(
            proteins = n, edges = nrow(ij), isolated = sum(size[comp] == 1L),
            components = length(size), largest_component = max(size)
        )
    )
}

#' Sum protein counts within bins
#'
#' Collapses a sample-by-protein table to sample-by-bin by summing each
#' bin's members, as `tax_glom()` in phyloseq does for taxa, so totals per
#' sample are unchanged. Use the raw counts for
#' count-based differential abundance.
#'
#' @param counts Numeric matrix, or sparse `Matrix`, with samples in rows and
#'   proteins in columns named as in `bins`.
#' @param bins Result of [bin_proteins()].
#' @return Samples-by-bins matrix (sparse if `counts` is), with the bins that
#'   hold at least one column of `counts`, in `bins` order.
#' @examples
#' x <- Biostrings::AAStringSet(c(p1 = "MKVLAA", p2 = "MKVLAG", p3 = "WWPPGG"))
#' S4Vectors::mcols(x)$embedding <- rbind(c(1, 0), c(0.95, 0.05), c(0, 1))
#' b <- bin_proteins(protein_similarity(x, min_coverage = 0))
#' counts <- rbind(s1 = c(p1 = 5, p2 = 3, p3 = 1), s2 = c(0, 2, 7))
#' bin_glom(counts, b)
#' @export
bin_glom <- function(counts, bins) {
    if (!is.character(bins$membership) || is.null(names(bins$membership))) {
        stop("`bins` must be a bin_proteins() result.", call. = FALSE)
    }
    if (!inherits(counts, "sparseMatrix")) counts <- as.matrix(counts)
    ids <- colnames(counts)
    .check_ids(ids, "Protein")
    missing <- setdiff(ids, names(bins$membership))
    if (length(missing)) {
        stop(length(missing), " protein(s) in `counts` are in no bin, e.g. ",
            missing[[1]], ".",
            call. = FALSE
        )
    }
    bin <- droplevels(factor(bins$membership[ids], names(bins$clusters)))
    A <- Matrix::sparseMatrix(
        i = seq_along(ids), j = as.integer(bin), x = 1,
        dims = c(length(ids), nlevels(bin)),
        dimnames = list(ids, levels(bin))
    )
    out <- counts %*% A
    if (inherits(counts, "sparseMatrix")) out else as.matrix(out)
}
