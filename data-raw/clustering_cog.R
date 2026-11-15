# Build curated/clustering_cog.rds: 2,000 NCBI COG 2024 proteins, 20 from
# each of 100 random COGs, with TM-Vec embeddings, MMseqs2 30% clusters and
# all-vs-all MMseqs2 hits. clustering.Rmd scores repdist bins, MMseqs2
# clusters and a TribeMCL-like graph against the COG labels.
#
# Run from the root of the gh-pages branch with repdist installed:
#   Rscript data-raw/clustering_cog.R
# Needs the `mmseqs` binary (MMseqs2) on the PATH, the curl package and the
# default `tmvec` model. The COG release (about 2.2 GB) is downloaded once
# into data/cog/, which is never committed.
library(repdist)

out <- file.path("curated", "clustering_cog.rds")
n_cog <- 100L
per_cog <- 20L
length_range <- c(50L, 1000L)
cog_dir <- file.path("data", "cog")
dir.create(cog_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Source: NCBI COG 2024 ---------------------------------------------------------
base <- "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/"
files <- c("checksums.md5", "cog-24.def.tab", "cog-24.cog.csv",
    "COGorg24.faa.gz")
for (f in files) {
    p <- file.path(cog_dir, f)
    if (!file.exists(p)) curl::curl_download(paste0(base, f), p)
}
md5 <- utils::read.table(file.path(cog_dir, "checksums.md5"),
    col.names = c("md5", "file"))
md5 <- md5[md5$file %in% files, ]
stopifnot(identical(unname(tools::md5sum(file.path(cog_dir, md5$file))),
    md5$md5))

# 2. One-label pool ----------------------------------------------------------------
# The rule of data_simulation.Rmd: exactly one COG, every assignment row in
# membership class 0 (full-length label), 50-1000 residues.
a <- utils::read.csv(file.path(cog_dir, "cog-24.cog.csv"), header = FALSE,
    colClasses = c("NULL", "NULL", "character", "integer", "NULL", "NULL",
        "character", "NULL", "integer", rep("NULL", 4)))
names(a) <- c("protein", "length", "cog", "class")
pc <- unique(a[c("protein", "cog")])
multi <- unique(pc$protein[duplicated(pc$protein)])
pool <- unique(a[!a$protein %in% c(multi, a$protein[a$class != 0L]),
    c("protein", "length", "cog")])
pool <- pool[pool$length >= length_range[1] &
    pool$length <= length_range[2], ]
rm(a, pc)
steps <- c(one_label_pool = nrow(pool))

# 3. Draw: 100 random COGs, 20 random proteins each ---------------------------------
# 150 COGs with at least 40 pool proteins, 40 candidates each (slack for the
# sequence filters); the first 100 COGs that keep 20 are used.
size <- table(pool$cog)
set.seed(20261114L)
cogs <- sample(names(size)[size >= 2L * per_cog], 1.5 * n_cog)
cand <- do.call(rbind, lapply(cogs, function(g) {
    m <- pool[pool$cog == g, ]
    m[sample.int(nrow(m), 2L * per_cog), ]
}))
ids <- tempfile(fileext = ".txt")
fa <- tempfile(fileext = ".fa")
writeLines(cand$protein, ids)
stopifnot(system(sprintf(paste("gzip -dc %s | awk 'NR==FNR {k[$1]; next}",
    "/^>/ {p = (substr($1, 2) in k)} p' %s - > %s"),
    shQuote(file.path(cog_dir, "COGorg24.faa.gz")), shQuote(ids),
    shQuote(fa))) == 0L)
s <- Biostrings::readAAStringSet(fa)
s <- s[!duplicated(sub(" .*", "", names(s)))]
hdr <- names(s)
names(s) <- sub(" .*", "", hdr)
cand$sequence <- as.character(s)[cand$protein]
cand$organism <- sub(".*\\[(.*)\\]$", "\\1",
    hdr[match(cand$protein, names(s))])
# Standard residues only; identical sequences keep their first draw.
cand <- cand[!is.na(cand$sequence) &
    grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cand$sequence), ]
cand <- cand[!duplicated(cand$sequence), ]
kept <- names(which(table(cand$cog) >= per_cog))
cogs <- head(cogs[cogs %in% kept], n_cog)
stopifnot(length(cogs) == n_cog)
annotation <- do.call(rbind, lapply(cogs, function(g) {
    head(cand[cand$cog == g, ], per_cog)
}))
def <- utils::read.delim(file.path(cog_dir, "cog-24.def.tab"), header = FALSE,
    quote = "", colClasses = c("character", "character", "character",
        rep("NULL", 4)), col.names = c("cog", "category", "name",
        rep("x", 4)))
annotation <- cbind(annotation, def[match(annotation$cog, def$cog),
    c("category", "name")])
rownames(annotation) <- NULL
seqs <- setNames(annotation$sequence, annotation$protein)
annotation$sequence <- NULL

# 4. Embed with the default model ----------------------------------------------------
# Sorted by length so each batch pads little; the order is restored after.
o <- order(nchar(seqs))
e <- S4Vectors::mcols(embed_proteins(seqs[o], model = "tmvec"))$embedding
x <- Biostrings::AAStringSet(seqs)
S4Vectors::mcols(x)$embedding <- round(e[match(names(x), names(seqs)[o]), ],
    6)
stopifnot(!anyNA(S4Vectors::mcols(x)$embedding))

# 5. MMseqs2 -------------------------------------------------------------------------
# 30% clusters with the settings of the simulation (80% coverage of both,
# greedy set cover, -s 7.5), and every hit of the all-vs-all search at
# E <= 0.01, the cut-off of protein_similarity(), in both directions.
tmp <- tempfile("mmseqs")
dir.create(tmp)
fasta <- file.path(tmp, "proteins.fasta")
Biostrings::writeXStringSet(x, fasta)
mmseqs <- function(...) {
    args <- c(...)
    stopifnot(system2("mmseqs", args, stdout = FALSE, stderr = FALSE) == 0L)
    gsub(tmp, "<tmp>", paste("mmseqs", paste(args, collapse = " ")),
        fixed = TRUE)
}
cluster_cmd <- mmseqs("easy-cluster", fasta, file.path(tmp, "c"),
    file.path(tmp, "work"), "--min-seq-id", 0.3, "-c", 0.8, "--cov-mode", 0,
    "--cluster-mode", 0, "-s", 7.5)
m <- utils::read.delim(file.path(tmp, "c_cluster.tsv"), header = FALSE,
    col.names = c("representative", "member"))
annotation$mmseqs30 <- m$representative[match(annotation$protein, m$member)]
stopifnot(!anyNA(annotation$mmseqs30))
m8 <- file.path(tmp, "hits.m8")
search_cmd <- mmseqs("easy-search", fasta, fasta, m8, file.path(tmp, "work"),
    "-s", 7.5, "-e", 0.01, "--max-seqs", 10000,
    "--format-output", "query,target,fident,qcov,tcov,evalue")
h <- utils::read.delim(m8, header = FALSE, col.names = c("query", "target",
    "fident", "qcov", "tcov", "evalue"))
h <- h[h$query != h$target, ]
h <- h[!duplicated(h[c("query", "target")]), ] # best alignment per direction
rownames(h) <- NULL
unlink(tmp, recursive = TRUE)

saveRDS(list(
    proteins = x,
    annotation = annotation,
    hits = h,
    provenance = list(
        source = paste("NCBI COG 2024,", base),
        accessed = format(Sys.Date()),
        pool_rule = paste("exactly one COG, all rows membership class 0,",
            paste(length_range, collapse = "-"), "residues"),
        draw = sprintf("%d random COGs with >= %d pool proteins, %d random",
            n_cog, 2L * per_cog, per_cog),
        seed = 20261114L,
        steps = steps,
        model = c(key = "tmvec", known_models()[["tmvec"]]),
        repdist = as.character(utils::packageVersion("repdist")),
        mmseqs_version = system2("mmseqs", "version", stdout = TRUE),
        cluster = cluster_cmd,
        search = search_cmd,
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
