# Benchmark the EC3 simulation (simulate_ec3.R) over every EC3 in the
# universe as the target, and a sensitivity grid over the MCL parameters.
# Writes validation/simulation_benchmark.rds for maintainer validation.
#
# Run from the source root after installing repdist; needs DESeq2 and vegan.
#
# The clustering never sees the EC3 labels, and the illustrative settings
# (min_sim = 0.5, inflation = 2, no coverage gate) were fixed before this ran. The grid is
# scored against the labels to show how the parameters move the result; it
# is not a tuning procedure, and no setting from it is reported as if chosen
# without labels.
library(repdist)
source(file.path("validation", "simulate_ec3.R"))

u <- readRDS(file.path("validation", "microbial_ec_universe.rds"))
members <- u$members
x <- Biostrings::AAStringSet(u$sequences)
S4Vectors::mcols(x)$embedding <- u$embeddings
ec3 <- stats::setNames(members$ec3, members$protein)
targets <- sort(unique(members$ec3))
grid <- expand.grid(min_sim = c(0.4, 0.5, 0.6, 0.7),
    inflation = c(1.4, 2, 4))
grid_bins <- lapply(seq_len(nrow(grid)), function(k) {
    S <- protein_similarity(x, grid$min_sim[k], min_coverage = 0)
    bin_proteins(S, grid$inflation[k])
})
reference <- which(grid$min_sim == 0.5 & grid$inflation == 2)
bins <- grid_bins[[reference]]
# EC3 itself as the grouping: the reference a perfect functional grouping
# would reach, available only because the labels are known.
oracle <- list(membership = ec3, clusters = split(names(ec3), ec3))

relative <- function(x) x / rowSums(x)
separation <- function(D, group) {
    group_separation(D, group)
}

per_target <- list()
sensitivity <- list()
for (i in seq_along(targets)) {
    target <- targets[[i]]
    message(i, "/", length(targets), ": ", target)
    set.seed(1000L + i)
    design <- design_upregulation(members, target)
    counts <- simulate_counts(members$protein, design)$counts
    up_a <- unlist(design$up[design$group == "A"])
    bin_cnt <- bin_glom(counts, bins)

    D <- list(
        `Bray-Curtis (protein)` = vegan::vegdist(relative(counts), "bray"),
        `RBF-MMD` = sample_mmd(counts, x),
        `Bray-Curtis (MCL bins)` = vegan::vegdist(relative(bin_cnt), "bray"),
        `Bray-Curtis (EC3 oracle)` = vegan::vegdist(
            relative(bin_glom(counts, oracle)), "bray"
        )
    )
    sep <- do.call(rbind, lapply(names(D), function(n) {
        data.frame(target = target, distance = n,
            t(separation(D[[n]], design$group)))
    }))

    protein_res <- deseq_ab(counts, design$group)
    da <- score_da(protein_res, deseq_ab(bin_cnt, design$group), bins, ec3,
        target, up_a)
    oracle_res <- deseq_ab(bin_glom(counts, oracle), design$group)
    per_target[[i]] <- list(
        separation = sep,
        scores = data.frame(target = target,
            t(c(score_bins(bins, ec3, target), da)),
            oracle_padj = oracle_res[target, "padj"],
            oracle_lfc = oracle_res[target, "log2FoldChange"])
    )

    sensitivity[[i]] <- do.call(rbind, lapply(seq_len(nrow(grid)), function(k) {
        b <- grid_bins[[k]]
        fit <- record_fit({
            res <- deseq_ab(bin_glom(counts, b), design$group)
            score_da(protein_res, res, b, ec3, target, up_a)
        }, da)
        data.frame(target = target, grid[k, ],
            t(c(score_bins(b, ec3, target), fit$scores)), error = fit$error)
    }))
    # Keep completed targets if a later target or the process fails.
    saveRDS(list(per_target = per_target, sensitivity = sensitivity),
        file.path("validation", "simulation_benchmark_partial.rds"))
}

benchmark <- list(
    separation = do.call(rbind, lapply(per_target, `[[`, "separation")),
    scores = do.call(rbind, lapply(per_target, `[[`, "scores")),
    sensitivity = do.call(rbind, sensitivity),
    settings = list(
        n_per_group = 10L, n_up = 5L, fold = 31,
        seeds = 1000L + seq_along(targets),
        min_sim = 0.5, inflation = 2, alpha = 0.05,
        universe = u$provenance[c("uniprot_release", "counts")]
    ),
    session = utils::capture.output(utils::sessionInfo())
)
rownames(benchmark$separation) <- NULL
rownames(benchmark$scores) <- NULL
rownames(benchmark$sensitivity) <- NULL
saveRDS(benchmark, file.path("validation", "simulation_benchmark.rds"),
    compress = "xz"
)
