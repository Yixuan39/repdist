# Build curated/simulation_kofam.rds: one COG (NCBI COG 2024) is the treatment
# label; every other protein belongs to background COGs. The page enriches
# each condition-A sample for treatment proteins from 30% clusters no other
# sample uses, so the test is remote homology inside one function.
#
# Run from the root of this branch:
#   KOFAM_COG=COGxxxx Rscript data-raw/simulation_kofam.R
#   KOFAM_SMOKE=1 KOFAM_COG=COGxxxx Rscript data-raw/simulation_kofam.R
#
# Needs ftp.ncbi.nih.gov (COG 2024 release, ~2.2 GB once, cached untracked in
# data/cog/), the `mmseqs` binary, and gzip/awk. KOFAM_SMOKE=1 builds tens of
# proteins with random unit-vector embeddings: a pipeline check, not a result.

smoke <- identical(Sys.getenv("KOFAM_SMOKE", "0"), "1")
treatment_cog <- Sys.getenv("KOFAM_COG")
if (!grepl("^COG[0-9]{4}$", treatment_cog)) {
    stop("Set KOFAM_COG to the treatment COG id (see KOFAM_SIM_PLAN.md).", call. = FALSE)
}
out <- file.path("curated", "simulation_kofam.rds")
cog_dir <- file.path("data", "cog")
dir.create(cog_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters, fixed before any result --------------------------------------------
length_range <- c(50L, 1000L)
n_treatment <- if (smoke) 20L else 300L   # seeded draw from the treatment COG
n_per_cog <- if (smoke) 5L else 20L       # proteins per background COG
global_cap <- if (smoke) 60L else 3000L
n_max <- if (smoke) 4L else 20L           # largest per-group n on the page
n_enriched <- 2L                          # enriched proteins per condition-A sample
seed <- 20261112L

# 1. COG 2024 release files, md5-checked --------------------------------------------
base <- "https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/"
files <- c("checksums.md5", "cog-24.def.tab", "cog-24.cog.csv", "COGorg24.faa.gz")
for (f in files) {
    p <- file.path(cog_dir, f)
    if (!file.exists(p)) curl::curl_download(paste0(base, f), p, quiet = FALSE)
}
md5 <- utils::read.table(file.path(cog_dir, "checksums.md5"), col.names = c("md5", "file"))
md5 <- md5[md5$file %in% files, ]
bad <- md5$file[unname(tools::md5sum(file.path(cog_dir, md5$file))) != md5$md5]
if (length(bad)) stop("md5 mismatch, delete and re-run: ", toString(bad), call. = FALSE)

# 2. One-label pool: exactly one COG, every row membership class 0 ------------------
# (class 0: the COG footprint covers most of the protein and most of the
# profile, i.e. a full-length label), 50-1000 aa.
message("2. Reading cog-24.cog.csv")
a <- utils::read.csv(file.path(cog_dir, "cog-24.cog.csv"), header = FALSE,
    colClasses = c("NULL", "NULL", "character", "integer", "NULL", "NULL",
        "character", "NULL", "integer", rep("NULL", 4)))
names(a) <- c("protein", "length", "cog", "class")
steps <- c(assignment_rows = nrow(a), proteins = length(unique(a$protein)))
pc <- unique(a[c("protein", "cog")])
multi <- unique(pc$protein[duplicated(pc$protein)])
pool <- unique(a[!a$protein %in% c(multi, a$protein[a$class != 0L]), c("protein", "length", "cog")])
rm(a, pc)
steps[["one_COG_class0"]] <- nrow(pool)
pool <- pool[pool$length >= length_range[1] & pool$length <= length_range[2], ]
steps[["length_50_1000"]] <- nrow(pool)
if (!treatment_cog %in% pool$cog) stop(treatment_cog, " has no qualifying proteins.", call. = FALSE)

# 3. Candidates, with slack for the sequence filters and the leakage step ------------
set.seed(seed)
shuffle <- function(v) v[sample.int(length(v))]
size <- table(pool$cog)
bg_cogs <- shuffle(setdiff(names(size)[size >= 3L * n_per_cog], treatment_cog))
n_bg_cogs <- ceiling((global_cap - n_treatment) / n_per_cog)
bg_cogs <- bg_cogs[seq_len(min(length(bg_cogs), 2L * n_bg_cogs))]
cand <- rbind(
    pool[pool$cog == treatment_cog, ],
    do.call(rbind, lapply(bg_cogs, function(cg) {
        d <- pool[pool$cog == cg, ]
        d[shuffle(seq_len(nrow(d)))[seq_len(3L * n_per_cog)], ]
    }))
)

message("3. Extracting ", nrow(cand), " candidate sequences from COGorg24.faa.gz")
ids <- tempfile(fileext = ".txt")
fa <- tempfile(fileext = ".fa")
writeLines(cand$protein, ids)
stopifnot(system(sprintf(paste("gzip -dc %s | awk 'NR==FNR {k[$1]; next}",
    "/^>/ {p = (substr($1, 2) in k)} p' %s - > %s"),
    shQuote(file.path(cog_dir, "COGorg24.faa.gz")), shQuote(ids), shQuote(fa))) == 0L)
s <- Biostrings::readAAStringSet(fa)
names(s) <- sub(" .*", "", names(s))
cand$sequence <- as.character(s)[cand$protein]
cand <- cand[!is.na(cand$sequence) & grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", cand$sequence), ]
cand <- cand[order(cand$protein), ]
cand <- cand[!duplicated(cand$sequence), ] # one protein id per unique sequence
steps[["treatment_COG_qualifying"]] <- sum(cand$cog == treatment_cog)

# Seeded draws in the shuffled order: n_treatment from the treatment COG,
# n_per_cog from each background COG that still has that many.
set.seed(seed)
tr <- cand[cand$cog == treatment_cog, ]
tr <- tr[shuffle(seq_len(nrow(tr)))[seq_len(min(nrow(tr), n_treatment))], ]
bg <- do.call(rbind, lapply(bg_cogs, function(cg) {
    d <- cand[cand$cog == cg, ]
    if (nrow(d) >= n_per_cog) d[seq_len(n_per_cog), ]
}))
tr$role <- "treatment"
bg$role <- "background"

# 4. MMseqs2 (settings of data-raw/simulation_mmseqs2.R) ----------------------------
mmseqs <- function(seqs) {
    tmp <- tempfile("mmseqs")
    dir.create(tmp)
    f <- file.path(tmp, "u.fasta")
    Biostrings::writeXStringSet(Biostrings::AAStringSet(seqs), f)
    run <- function(...) {
        args <- c(...)
        stopifnot(system2("mmseqs", args, stdout = FALSE, stderr = FALSE) == 0L)
        gsub(tmp, "<tmp>", paste("mmseqs", paste(args, collapse = " ")), fixed = TRUE)
    }
    c30 <- run("easy-cluster", f, file.path(tmp, "c30"), file.path(tmp, "w"),
        "--min-seq-id", 0.3, "-c", 0.8, "--cov-mode", 0, "--cluster-mode", 0, "-s", 7.5)
    m <- utils::read.delim(file.path(tmp, "c30_cluster.tsv"), header = FALSE)
    search <- run("easy-search", f, f, file.path(tmp, "hits.m8"), file.path(tmp, "w"),
        "-s", 7.5, "-e", "1e-3", "--max-seqs", 10000,
        "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue")
    h <- utils::read.delim(file.path(tmp, "hits.m8"), header = FALSE,
        col.names = c("query", "target", "fident", "alnlen", "qcov", "tcov", "evalue"))
    h <- h[h$query != h$target, ]
    # One row per unordered pair (its lowest-E direction), with flags over both.
    key <- paste(pmin(h$query, h$target), pmax(h$query, h$target))
    full <- h$qcov >= 0.8 & h$tcov >= 0.8
    h$cov80 <- as.logical(ave(full, key, FUN = any))
    h$id30_cov80 <- as.logical(ave(full & h$fident >= 0.3, key, FUN = any))
    o <- order(h$evalue)
    h <- h[o, ][!duplicated(key[o]), ]
    names(h)[1:2] <- c("p1", "p2")
    rownames(h) <- NULL
    unlink(tmp, recursive = TRUE)
    list(c30 = setNames(m$V1, m$V2)[names(seqs)], hits = h,
        provenance = list(cluster_30 = c30, search = search,
            version = system2("mmseqs", "version", stdout = TRUE)))
}

message("4. Leakage: drop background proteins with an E <= 1e-3, 80%-coverage hit to the treatment")
u <- rbind(tr, bg)
mm <- mmseqs(setNames(u$sequence, u$protein))
h <- mm$hits[mm$hits$cov80, ]
tr_ids <- tr$protein
leak <- unique(c(h$p1[h$p2 %in% tr_ids], h$p2[h$p1 %in% tr_ids]))
leak <- intersect(leak, bg$protein)
steps[["background_leakage_removed"]] <- length(leak)
bg <- bg[!bg$protein %in% leak, ]
# Whole background COGs, in the shuffled order, up to the cap.
keep_cogs <- unique(bg$cog)[cumsum(table(factor(bg$cog, unique(bg$cog)))) <=
    global_cap - nrow(tr)]
bg <- bg[bg$cog %in% keep_cogs, ]
u <- rbind(tr, bg)
rownames(u) <- NULL
mm <- mmseqs(setNames(u$sequence, u$protein))
u$mmseqs30 <- unname(mm$c30[u$protein])
steps[["universe"]] <- nrow(u)
steps[["treatment"]] <- nrow(tr)
steps[["background_COGs"]] <- length(keep_cogs)
steps[["treatment_30pct_clusters"]] <- length(unique(u$mmseqs30[u$role == "treatment"]))
if (steps[["treatment_30pct_clusters"]] < n_enriched * n_max) {
    stop(treatment_cog, " spans ", steps[["treatment_30pct_clusters"]], " 30% clusters; ",
        "the design needs ", n_enriched * n_max, ".", call. = FALSE)
}

# 5. Embeddings -------------------------------------------------------------------
x <- Biostrings::AAStringSet(setNames(u$sequence, u$protein))
if (smoke) {
    set.seed(seed)
    Z <- matrix(rnorm(length(x) * 512), length(x), 512, dimnames = list(names(x), NULL))
    Z <- Z / sqrt(rowSums(Z^2))
    embed_note <- "KOFAM_SMOKE=1: random unit vectors, not TM-Vec"
} else {
    library(repdist)
    basilisk::setBasiliskShared(FALSE)
    basilisk::setBasiliskFork(FALSE)
    o <- order(Biostrings::width(x)) # length-sorted batches embed ~3x faster
    Z <- S4Vectors::mcols(embed_proteins(x[o], model = "tmvec"))$embedding
    embed_note <- "embed_proteins(model = 'tmvec')"
}
S4Vectors::mcols(x)$embedding <- round(Z[names(x), ], 6)

# 6. Similarity graph (quadratic alignments, so built here, not at knit) ------------
message("6. protein_similarity() on ", length(x), " proteins")
t0 <- Sys.time()
S <- repdist::protein_similarity(x, n_cores = 1L)
align_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
ij <- which(upper.tri(S) & S > 0, arr.ind = TRUE)
edges <- data.frame(i = ij[, 1], j = ij[, 2], s = S[ij])

def <- utils::read.delim(file.path(cog_dir, "cog-24.def.tab"), header = FALSE,
    quote = "", fileEncoding = "latin1")
u$cog_name <- def$V3[match(u$cog, def$V1)]
u$sequence <- NULL
saveRDS(list(
    proteins = x,
    annotation = u,
    similarity = edges, # upper-triangle edges of protein_similarity(proteins)
    mmseqs = list(hits = mm$hits, provenance = mm$provenance),
    design = list(treatment_cog = treatment_cog, n_enriched = n_enriched, n_max = n_max,
        n_treatment = n_treatment, n_per_cog = n_per_cog, global_cap = global_cap,
        seed = seed),
    steps = steps,
    provenance = list(
        source = paste0(base, " (COG 2024 release)"),
        md5 = setNames(md5$md5, md5$file),
        accessed = format(Sys.Date()),
        length_range = length_range,
        smoke = smoke,
        embedding = embed_note,
        align_minutes = align_minutes,
        repdist = as.character(utils::packageVersion("repdist")),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
print(steps)
message(sprintf("protein_similarity(): %.1f min for %d proteins", align_minutes, length(x)))
