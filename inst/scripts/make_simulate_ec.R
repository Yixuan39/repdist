# Writes data/simulate_ec.rda: contrast 1 of the website's sequence-disjoint
# EC simulation (simulation.Rmd on the gh-pages branch), on a protein universe
# cut down to fit in the package: every protein of the two enzyme classes plus
# a seeded background of 1000 others. Laid out like the case study's curated
# object (proteins, counts, samples, annotation, provenance). See
# DATA_PROVENANCE.md. Run from the source repository root with repdist
# installed.
#
# Usage: Rscript inst/scripts/make_simulate_ec.R \
#     [simulation_universe.rds simulation_mmseqs2.rds]

library(repdist) # also makes the readRDS() AAStringSet subsettable
src <- commandArgs(TRUE)[1:2]
if (anyNA(src)) {
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    src <- vapply(c("simulation_universe.rds", "simulation_mmseqs2.rds"),
        function(f) BiocFileCache::bfcrpath(bfc, paste0(
            "https://github.com/Yixuan39/repdist/raw/gh-pages/curated/", f)),
        "")
}
u <- readRDS(src[1])
mm <- readRDS(src[2])
stopifnot(identical(mm$provenance$universe_md5, unname(tools::md5sum(src[1]))),
    identical(mm$clusters$protein, u$annotation$protein))
ann <- u$annotation
ann$mmseqs50 <- mm$clusters$mmseqs50
ann$mmseqs30 <- mm$clusters$mmseqs30
mm50 <- setNames(ann$mmseqs50, ann$protein)

# Eligible classes and the seeded pairing, as on the simulation page: an EC3
# needs 20 distinct 50% clusters; contrast 1 is the first shuffled pair.
n_per_group <- 10L
n_enriched <- 2L
k <- n_per_group * n_enriched
n50 <- tapply(ann$mmseqs50, ann$ec3, function(v) length(unique(v)))
set.seed(20261008)
ec <- sample(sort(names(n50)[n50 >= k]))[1:2]

# Enriched proteins of contrast 1: one per 50% cluster, one cluster per EC4
# before a second from any EC4; B never reuses a cluster A used.
draw_enriched <- function(ec, exclude = character()) {
    cand <- unique(ann[ann$ec3 == ec & !ann$mmseqs50 %in% exclude,
        c("mmseqs50", "ec4")])
    cand <- cand[sample.int(nrow(cand)), ]
    cand <- cand[!duplicated(cand$mmseqs50), ]
    cand <- cand[order(ave(seq_len(nrow(cand)), cand$ec4, FUN = seq_along)), ]
    vapply(seq_len(k), function(r) {
        p <- ann$protein[ann$mmseqs50 == cand$mmseqs50[r] &
            ann$ec4 == cand$ec4[r]]
        p[sample.int(length(p), 1L)]
    }, "")
}
set.seed(1)
a <- draw_enriched(ec[1])
b <- draw_enriched(ec[2], exclude = mm50[a])
stopifnot(!anyDuplicated(mm50[a]), !anyDuplicated(mm50[b]))
samples <- sprintf("s%02d", seq_len(2 * n_per_group))
up <- setNames(unname(split(c(a, b), rep(samples, each = n_enriched))),
    samples)

# The cut-down universe, in universe order
n_background <- 1000L
background <- sample(ann$protein[!ann$ec3 %in% ec], n_background)
keep <- ann$ec3 %in% ec | ann$protein %in% background
ann <- ann[keep, ]
ids <- ann$protein

# Negative-binomial counts, the simulation page's model: sample i puts rho of
# its expected abundance on its enriched proteins over a log-normal baseline;
# the k-th enriched protein of B takes the baseline of the k-th of A.
rho <- 0.05
depth <- 2e4
pi <- setNames(rlnorm(length(ids)), ids)
pi[b] <- pi[a]
pi <- pi / sum(pi)
s <- rlnorm(length(samples), 0, 0.25)
e <- rlnorm(length(ids), 0, 0.5)
p <- matrix((1 - rho) * pi, length(samples), length(ids), byrow = TRUE,
    dimnames = list(samples, ids))
for (j in samples) p[j, up[[j]]] <- p[j, up[[j]]] + rho / n_enriched
mu <- s * depth * p
phi <- (0.1 + 1 / mu) * rep(e, each = length(samples))
counts <- matrix(rnbinom(length(mu), mu = mu, size = 1 / phi),
    length(samples), dimnames = dimnames(p))
storage.mode(counts) <- "integer"

# Rebuilt rather than subset, which would keep all 4920 sequences in the pool;
# embeddings rounded to 4 decimals to fit the package
proteins <- Biostrings::AAStringSet(as.character(u$proteins[ids]))
S4Vectors::mcols(proteins)$embedding <- round(
    S4Vectors::mcols(u$proteins)$embedding[keep, ], 4)

group <- factor(rep(c("A", "B"), each = n_per_group))
simulate_ec <- list(
    proteins = proteins,
    counts = counts,
    samples = data.frame(sample = samples, group = group,
        enriched_ec3 = ec[as.integer(group)], row.names = samples),
    annotation = data.frame(protein = ids, name = ann$name,
        organism = ann$organism, length = ann$length, ec3 = ann$ec3,
        ec4 = ann$ec4, mmseqs50 = ann$mmseqs50, mmseqs30 = ann$mmseqs30,
        enriched_in = rep(samples, each = n_enriched)[match(ids, unlist(up))],
        row.names = ids),
    provenance = list(
        source = u$provenance$source,
        uniprot_release = u$provenance$uniprot_release,
        license = "CC BY 4.0",
        universe = "curated/simulation_universe.rds (gh-pages branch)",
        universe_md5 = unname(tools::md5sum(src[1])),
        universe_steps = u$provenance$steps,
        model = u$provenance$model,
        mmseqs = mm$provenance[c("mmseqs_version", "cluster_50",
            "cluster_30")],
        design = list(contrast = 1L, ec3 = ec, rho = rho, depth = depth,
            n_per_group = n_per_group, n_enriched = n_enriched,
            n_background = n_background, seeds = c(pairing = 20261008,
                draw = 1)),
        script_md5 = unname(tools::md5sum(sub("^--file=", "",
            grep("^--file=", commandArgs(), value = TRUE)))),
        repdist = as.character(utils::packageVersion("repdist")),
        built = format(Sys.Date()),
        session = utils::capture.output(utils::sessionInfo())
    )
)
dir.create("data", showWarnings = FALSE)
save(simulate_ec, file = file.path("data", "simulate_ec.rda"),
    compress = "xz")
