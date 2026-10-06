# Simulation of a coherent functional perturbation carried by non-overlapping
# accessions, shared by benchmark_simulation.R and test-simulation.R.
# Every protein is present in every sample; only abundance changes.

# Which proteins each sample upregulates. Group A: every sample takes `n_up`
# private proteins, all from the `target` EC3. Group B: every sample takes
# `n_up` private proteins from the other EC3s, dealt round-robin over a
# shuffled EC3 order, so a sample's proteins come from distinct EC3s and no
# EC3 is favoured across the group. No accession is upregulated twice.
design_upregulation <- function(members, target, n_per_group = 10L,
                                n_up = 5L) {
    k <- n_per_group * n_up
    ids_a <- members$protein[members$ec3 == target]
    stopifnot(length(ids_a) >= k)
    blocks <- rep(seq_len(n_per_group), each = n_up)
    up_a <- unname(split(sample(ids_a, k), blocks))

    other <- members[members$ec3 != target, ]
    ec3s <- sample(unique(other$ec3))
    stopifnot(length(ec3s) >= n_up)
    dealt <- rep(ec3s, length.out = k)
    pools <- lapply(split(other$protein, other$ec3), sample)
    taken <- ave(seq_along(dealt), dealt, FUN = seq_along)
    up_b <- unname(split(
        vapply(seq_along(dealt), function(i) pools[[dealt[i]]][taken[i]],
            character(1)),
        blocks
    ))

    up <- c(up_a, up_b)
    names(up) <- sprintf("s%02d", seq_along(up))
    list(
        group = factor(rep(c("A", "B"), each = n_per_group)),
        up = up,
        target = target
    )
}

# Negative-binomial counts. Protein j has baseline relative abundance pi_j
# (lognormal, normalised) and baseline mean mu_j = depth * pi_j; sample i has
# size factor s_i (lognormal); an upregulated protein's mean is multiplied by
# `fold` (31 by default). With 50 target proteins each raised in exactly one
# of ten A samples, the size-factor-normalised group mean target aggregate
# is exactly 1 + (31 - 1) / 10 = 4 times baseline, for any baseline abundances.
# Individual samples vary with the baseline mass of their five proteins.
# Dispersion follows a DESeq2-style trend, phi_j = asymptote +
# extra / mu_j, with lognormal scatter. The upregulated proteins of the k-th
# sample of group B take the baseline abundances of those of the k-th sample
# of group A, so both groups move the same amount of abundance and differ only
# in which proteins carry it. One raw count matrix is returned, samples in
# rows, with every parameter needed to inspect the ground truth.
simulate_counts <- function(proteins, design, fold = 31, depth = 5e5,
                            sd_abundance = 1, sd_size = 0.25,
                            asymptote = 0.1, extra = 1, sd_dispersion = 0.5) {
    up <- design$up
    n <- length(proteins)
    m <- length(up)
    pi <- stats::setNames(exp(stats::rnorm(n, 0, sd_abundance)), proteins)
    a <- which(design$group == "A")
    b <- which(design$group == "B")
    for (k in seq_along(a)) {
        pi[up[[b[k]]]] <- sample(pi[up[[a[k]]]])
    }
    pi <- pi / sum(pi)
    mu0 <- depth * pi
    phi <- (asymptote + extra / mu0) * exp(stats::rnorm(n, 0, sd_dispersion))
    size_factor <- exp(stats::rnorm(m, 0, sd_size))
    fc <- matrix(1, m, n, dimnames = list(names(up), proteins))
    for (i in seq_len(m)) fc[i, up[[i]]] <- fold
    mu <- outer(size_factor, mu0) * fc
    counts <- matrix(
        stats::rnbinom(m * n, mu = mu, size = rep(1 / phi, each = m)),
        m, n,
        dimnames = dimnames(fc)
    )
    list(
        counts = counts,
        truth = list(
            baseline = pi, dispersion = stats::setNames(phi, proteins),
            size_factor = stats::setNames(size_factor, names(up)),
            fold = fold, up = up
        )
    )
}

# Scores against the EC3 labels. The clustering never sees them; they are
# used here only, to ask how the representation-based bins relate to EC3.
score_bins <- function(bins, ec3, target) {
    membership <- bins$membership[names(ec3)]
    size <- table(membership)
    multi <- membership %in% names(size)[size > 1L]
    plurality <- tapply(ec3[multi], membership[multi],
        function(x) max(table(x)))
    in_target <- table(membership[ec3 == target])
    c(
        bins = length(size),
        in_multi_bins = mean(multi),
        ec3_purity = sum(plurality) / sum(multi),
        target_bins = length(in_target),
        target_completeness = max(in_target) / sum(ec3 == target)
    )
}

# DESeq2 Wald test of group A against B. Outlier replacement and Cook's
# filtering are off at both levels: a single upregulated sample is the
# design, not an outlier to remove.
deseq_ab <- function(counts, group) {
    dds <- DESeq2::DESeqDataSetFromMatrix(
        t(round(as.matrix(counts))), data.frame(group = group), ~group
    )
    dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf, quiet = TRUE)
    res <- DESeq2::results(dds, contrast = c("group", "A", "B"),
        cooksCutoff = FALSE)
    res <- as.data.frame(res)
    res$padj[is.na(res$padj)] <- 1
    res
}

# Differential abundance at the protein and bin levels against the design:
# `up_a` are the target proteins upregulated in group A.
score_da <- function(
    protein_res, bin_res, bins, ec3, target, up_a, alpha = 0.05) {
    sig <- function(r, sign) r$padj < alpha & sign * r$log2FoldChange > 0
    target_bin <- names(which.max(table(bins$membership[names(ec3)][
        ec3 == target])))
    rank_p <- rank(bin_res$pvalue, ties.method = "min")
    names(rank_p) <- rownames(bin_res)
    bin_up <- rownames(bin_res)[sig(bin_res, 1)]
    plurality <- vapply(bins$clusters[bin_up], function(x) {
        names(which.max(table(ec3[x])))
    }, character(1))
    c(
        protein_sig = sum(protein_res$padj < alpha),
        protein_sig_target = sum(sig(protein_res[up_a, ], 1)),
        protein_recall = mean(sig(protein_res[up_a, ], 1)),
        bin_sig = sum(bin_res$padj < alpha),
        bin_sig_up_a = length(bin_up),
        bin_sig_up_a_target = sum(plurality == target),
        bin_sig_up_b = sum(sig(bin_res, -1)),
        bin_recall = mean(bins$membership[up_a] %in% bin_up),
        target_bin_padj = bin_res[target_bin, "padj"],
        target_bin_rank = rank_p[[target_bin]],
        target_bin_lfc = bin_res[target_bin, "log2FoldChange"]
    )
}

# Group separation on one distance: PERMANOVA R2 and p; PERMDISP (betadisper
# permutation p, and the mean distance to each group's median); and the mean
# and median A-A, B-B and A-B distances. A location shift puts A-B above both
# within-group values; A-A far below B-B ~ A-B is a dispersion difference.
group_separation <- function(D, group, permutations = 999) {
    set.seed(1)
    a <- vegan::adonis2(D ~ group, data.frame(group = group),
        permutations = permutations)
    bd <- vegan::betadisper(D, group)
    set.seed(1)
    pd <- vegan::permutest(bd, permutations = permutations)
    M <- as.matrix(D)
    in_a <- group == levels(group)[1]
    aa <- M[in_a, in_a][upper.tri(diag(sum(in_a)))]
    bb <- M[!in_a, !in_a][upper.tri(diag(sum(!in_a)))]
    ab <- as.vector(M[in_a, !in_a])
    c(
        R2 = a$R2[[1]], p = a$`Pr(>F)`[[1]], permdisp_p = pd$tab[1, "Pr(>F)"],
        disp_A = mean(bd$distances[in_a]), disp_B = mean(bd$distances[!in_a]),
        AA_mean = mean(aa), BB_mean = mean(bb), AB_mean = mean(ab),
        AA_median = stats::median(aa), BB_median = stats::median(bb),
        AB_median = stats::median(ab)
    )
}

# Validation only: preserve a failed fit as missing scores, never as zero
# discoveries or a different dispersion model. `template` names the scores.
record_fit <- function(expr, template) {
    tryCatch(
        list(scores = force(expr), error = NA_character_),
        error = function(e) list(scores = template * NA_real_,
            error = conditionMessage(e))
    )
}
