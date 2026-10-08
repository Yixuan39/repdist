# Build curated/simulation_mmseqs2.rds: MMseqs2 sequence clusters and
# all-vs-all sequence hits for the proteins of curated/simulation_universe.rds.
# The simulation uses the clusters as the conventional sequence-family
# baseline and to keep enriched proteins in different sequence families.
#
# Run from the root of the gh-pages branch after simulation_universe.R:
#   Rscript data-raw/simulation_mmseqs2.R
# Needs the `mmseqs` binary (MMseqs2) on the PATH.

universe <- file.path("curated", "simulation_universe.rds")
out <- file.path("curated", "simulation_mmseqs2.rds")
x <- readRDS(universe)$proteins

tmp <- tempfile("mmseqs")
dir.create(tmp)
fasta <- file.path(tmp, "universe.fasta")
Biostrings::writeXStringSet(x, fasta)
mmseqs <- function(...) {
    args <- c(...)
    stopifnot(system2("mmseqs", args, stdout = FALSE, stderr = FALSE) == 0L)
    paste("mmseqs", paste(args, collapse = " "))
}

# 1. Clusters ---------------------------------------------------------------------
# Cascaded clustering (greedy set cover) at 50% identity, the usual protein
# family level (as UniRef50), and at 30% as a permissive sensitivity level.
# Both require 80% coverage of query and target (--cov-mode 0) and use the
# most sensitive prefilter (-s 7.5), so the sequence baseline merges as much
# as these thresholds allow.
cluster <- function(id) {
    prefix <- file.path(tmp, paste0("c", id))
    cmd <- mmseqs("easy-cluster", fasta, prefix, file.path(tmp, "work"),
        "--min-seq-id", id, "-c", 0.8, "--cov-mode", 0, "--cluster-mode", 0,
        "-s", 7.5)
    m <- utils::read.delim(paste0(prefix, "_cluster.tsv"), header = FALSE,
        col.names = c("representative", "member"))
    list(membership = m$representative[match(names(x), m$member)],
        command = gsub(tmp, "<tmp>", cmd, fixed = TRUE))
}
c50 <- cluster(0.5)
c30 <- cluster(0.3)
clusters <- data.frame(protein = names(x), mmseqs50 = c50$membership,
    mmseqs30 = c30$membership)
stopifnot(!anyNA(clusters))

# 2. All-vs-all hits ----------------------------------------------------------------
# Local alignments with E-value <= 1e-3; the best hit per unordered pair.
# Pairs absent from the table have no detectable sequence similarity at this
# sensitivity.
m8 <- file.path(tmp, "hits.m8")
search <- mmseqs("easy-search", fasta, fasta, m8, file.path(tmp, "work"),
    "-s", 7.5, "-e", "1e-3", "--max-seqs", 10000,
    "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue")
h <- utils::read.delim(m8, header = FALSE, col.names = c("query", "target",
    "fident", "alnlen", "qcov", "tcov", "evalue"))
h <- h[h$query != h$target, ]
h$p1 <- pmin(h$query, h$target)
h$p2 <- pmax(h$query, h$target)
h <- h[order(-h$fident), ]
hits <- h[!duplicated(h[c("p1", "p2")]),
    c("p1", "p2", "fident", "alnlen", "qcov", "tcov", "evalue")]
rownames(hits) <- NULL

saveRDS(list(
    clusters = clusters,
    hits = hits,
    provenance = list(
        universe = universe,
        universe_md5 = unname(tools::md5sum(universe)),
        mmseqs_version = system2("mmseqs", "version", stdout = TRUE),
        cluster_50 = c50$command,
        cluster_30 = c30$command,
        search = gsub(tmp, "<tmp>", search, fixed = TRUE),
        n_clusters = c(mmseqs50 = length(unique(clusters$mmseqs50)),
            mmseqs30 = length(unique(clusters$mmseqs30))),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
unlink(tmp, recursive = TRUE)
