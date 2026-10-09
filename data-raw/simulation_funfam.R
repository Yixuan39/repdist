# Build curated/simulation_funfam.rds: prokaryotic CATH FunFam domain
# sequences, one per 50%-identity MMseqs2 cluster, embedded from their
# sequences with the default repdist model. Read by simulation_funfam.Rmd.
# Which FunFams are simulation targets is decided on the page (enough
# distinct clusters); every FunFam kept here is background.
#
# Run from the root of the gh-pages branch with repdist installed:
#   Rscript data-raw/simulation_funfam.R
# Needs network access to www.cathdb.info and download.cathdb.info, the curl
# and jsonlite packages, the `mmseqs` binary (MMseqs2) on the PATH, and the
# default `tmvec` model. Downloads are kept in data/cath_v4_3/ (not
# committed) and reused on a rerun; cathdb.info rate-limits bursts, so
# requests are spaced.
library(repdist)

out <- file.path("curated", "simulation_funfam.rds")
dl <- file.path("data", "cath_v4_3")
dir.create(dl, recursive = TRUE, showWarnings = FALSE)
cath <- "https://www.cathdb.info/version/v4_3_0"
sf_list_url <- paste0("https://download.cathdb.info/cath/releases/all-releases/",
    "v4_3_0/cath-classification-data/cath-superfamily-list-v4_3_0.txt")

# Every rule is fixed here, before any sequence is embedded.
n_superfamilies <- 300L # most structurally diverse superfamilies (S35 reps)
min_dops <- 70 # CATH's threshold for an informative FunFam alignment
kingdoms <- c("Bacteria", "Archaea")
length_range <- c(50L, 1000L) # the TM-Vec head's training range
max_members <- 100L # cap per FunFam, a seeded draw
seed <- 20261008L

fetch <- function(url, file) {
    if (!file.exists(file)) {
        for (attempt in 1:10) {
            Sys.sleep(3)
            r <- curl::curl_fetch_disk(url, paste0(file, ".part"))
            if (r$status_code == 200L) break
            # rate limited (429) or a transient server error (5xx): wait, retry
            if (r$status_code != 429L && r$status_code < 500L) break
            Sys.sleep(90)
        }
        if (r$status_code != 200L) {
            stop("Download failed (HTTP ", r$status_code, "): ", url,
                call. = FALSE)
        }
        file.rename(paste0(file, ".part"), file)
    }
    file
}

# 1. Superfamilies and one FunFam each ----------------------------------------------
sf <- utils::read.delim(fetch(sf_list_url, file.path(dl, basename(sf_list_url))),
    comment.char = "#", header = FALSE, quote = "",
    col.names = c("superfamily", "s35_reps", "domains", "name"))
sf <- head(sf[order(-sf$s35_reps, sf$superfamily), ], n_superfamilies)

# The largest FunFam (most members) with DOPS >= 70 in each superfamily; one
# per superfamily, so different FunFams are different folds or superfamilies.
families <- do.call(rbind, lapply(sf$superfamily, function(s) {
    f <- fetch(sprintf("%s/api/rest/superfamily/%s/funfam", cath, s),
        file.path(dl, paste0(s, "_funfams.json")))
    j <- jsonlite::fromJSON(f)$data
    j <- data.frame(superfamily = s, funfam = as.integer(j$funfam_number),
        name = j$name, members = as.integer(j$num_members_in_funfam),
        dops = as.numeric(j$seed_dops_score))
    j <- j[!is.na(j$dops) & j$dops >= min_dops, ]
    j[order(-j$members, j$funfam)[1], ]
}))
families <- families[!is.na(families$funfam), ]
families$family <- sprintf("%s/FF/%06d", families$superfamily, families$funfam)

# 2. Domain sequences from the FunFam seed alignments --------------------------------
# Ungapping a Stockholm row gives the domain's residues; lower case marks
# insert columns and is kept.
read_funfam <- function(file) {
    l <- readLines(file)
    gs <- grep("^#=GS \\S+ +DR ORG;", l, value = TRUE)
    org <- data.frame(domain = sub("^#=GS (\\S+) .*", "\\1", gs),
        kingdom = sub(".*DR ORG; ([^;]+);.*", "\\1", gs),
        organism = sub(".*; ([^;]+); *$", "\\1", gs))
    ac <- grep("^#=GS \\S+ +AC ", l, value = TRUE)
    ac <- setNames(sub(".* AC +", "", ac), sub("^#=GS (\\S+) .*", "\\1", ac))
    rows <- l[!grepl("^(#|//)", l) & nzchar(trimws(l))]
    id <- sub(" .*", "", rows)
    aln <- tapply(sub("^\\S+ +", "", rows), factor(id, unique(id)), paste,
        collapse = "")
    seq <- toupper(gsub("[-.]", "", aln))
    data.frame(domain = names(seq), accession = unname(ac[names(seq)]),
        sequence = unname(seq))[, c("domain", "accession", "sequence")] |>
        merge(org, by = "domain")
}
dom <- do.call(rbind, lapply(seq_len(nrow(families)), function(i) {
    f <- families[i, ]
    file <- fetch(sprintf("%s/superfamily/%s/funfam/%d/files/stockholm", cath,
        f$superfamily, f$funfam), file.path(dl, sprintf("%s_FF%06d.sto",
        f$superfamily, f$funfam)))
    cbind(family = f$family, read_funfam(file))
}))
steps <- c(seed_alignment_domains = nrow(dom))
dom <- dom[dom$kingdom %in% kingdoms, ]
steps[["bacteria_archaea"]] <- nrow(dom)
dom$length <- nchar(dom$sequence)
dom <- dom[dom$length >= length_range[1] & dom$length <= length_range[2], ]
steps[["length_50_1000"]] <- nrow(dom)
# 20 standard residues only (BLOSUM62 alignments in protein_similarity())
dom <- dom[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", dom$sequence), ]
steps[["standard_residues"]] <- nrow(dom)
dom <- dom[order(dom$domain), ]
dom <- dom[!duplicated(dom$sequence), ]
steps[["unique_sequences"]] <- nrow(dom)

# 3. One domain per 50%-identity cluster ---------------------------------------------
# Same settings as the Simulation 2 primary clustering. A cluster that spans
# two FunFams is dropped; otherwise its lowest domain id represents it, so no
# two retained domains share a 50% cluster.
tmp <- tempfile("mmseqs")
dir.create(tmp)
fasta <- file.path(tmp, "domains.fasta")
Biostrings::writeXStringSet(Biostrings::AAStringSet(setNames(dom$sequence,
    dom$domain)), fasta)
mmseqs_args <- c("easy-cluster", fasta, file.path(tmp, "c50"),
    file.path(tmp, "work"), "--min-seq-id", 0.5, "-c", 0.8, "--cov-mode", 0,
    "--cluster-mode", 0, "-s", 7.5)
stopifnot(system2("mmseqs", mmseqs_args, stdout = FALSE, stderr = FALSE) == 0L)
m <- utils::read.delim(file.path(tmp, "c50_cluster.tsv"), header = FALSE,
    col.names = c("representative", "member"))
dom$mmseqs50 <- m$representative[match(dom$domain, m$member)]
stopifnot(!anyNA(dom$mmseqs50))
mixed <- tapply(dom$family, dom$mmseqs50, function(v) length(unique(v)) > 1)
dom <- dom[!mixed[dom$mmseqs50], ]
steps[["clusters_within_one_funfam"]] <- nrow(dom)
dom <- dom[!duplicated(dom$mmseqs50), ] # dom is sorted by domain id
steps[["one_per_50pct_cluster"]] <- nrow(dom)
# One domain per UniProt protein: a multi-domain protein would otherwise be
# one accession in two FunFams (or twice in one, for repeats).
dom <- dom[!duplicated(dom$accession), ]
steps[["one_per_accession"]] <- nrow(dom)
unlink(tmp, recursive = TRUE)

# 4. FunFams with any domain left, capped ------------------------------------------
n <- table(dom$family)
families$candidates <- as.integer(n[families$family])
families$candidates[is.na(families$candidates)] <- 0L
families <- families[families$candidates > 0L, ]
set.seed(seed)
dom <- do.call(rbind, lapply(families$family, function(f) {
    d <- dom[dom$family == f, ]
    d[sort(sample.int(nrow(d), min(nrow(d), max_members))), ]
}))
families$kept <- as.integer(table(dom$family)[families$family])
steps[["capped_per_funfam"]] <- nrow(dom)
dom$label <- families$name[match(dom$family, families$family)]
rownames(dom) <- NULL

# 5. Embed every sequence with the default model -------------------------------------
# In chunks of 500, each in a fresh worker process: on Apple GPUs (MPS) one
# long run can exhaust device memory and fall back to the CPU.
embed_chunked <- function(seqs, chunk = 500L) {
    basilisk::setBasiliskShared(FALSE)
    basilisk::setBasiliskFork(FALSE)
    parts <- lapply(split(seq_along(seqs), ceiling(seq_along(seqs) / chunk)),
        function(i) {
            message("embedding ", max(i), " / ", length(seqs))
            S4Vectors::mcols(embed_proteins(seqs[i], model = "tmvec"))$embedding
        })
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- do.call(rbind, unname(parts))
    x
}
x <- embed_chunked(setNames(dom$sequence, dom$domain))
stopifnot(identical(names(x), dom$domain),
    identical(rownames(S4Vectors::mcols(x)$embedding), dom$domain))
S4Vectors::mcols(x)$embedding <- round(S4Vectors::mcols(x)$embedding, 6)
dom$sequence <- NULL

saveRDS(list(
    proteins = x,
    annotation = dom,
    families = families,
    provenance = list(
        source = paste("CATH v4.3 FunFam seed alignments,",
            "https://www.cathdb.info (CC BY 4.0)"),
        alignment_url = paste0(cath,
            "/superfamily/{superfamily}/funfam/{number}/files/stockholm"),
        funfam_list_url = paste0(cath, "/api/rest/superfamily/{superfamily}/funfam"),
        superfamily_list = sf_list_url,
        accessed = format(Sys.Date()),
        rules = list(n_superfamilies = n_superfamilies, min_dops = min_dops,
            kingdoms = kingdoms, length_range = length_range,
            max_members = max_members, seed = seed),
        superfamilies = sf,
        mmseqs_version = system2("mmseqs", "version", stdout = TRUE),
        mmseqs = paste("mmseqs", paste(gsub(tmp, "<tmp>", mmseqs_args,
            fixed = TRUE), collapse = " ")),
        steps = steps,
        model = c(key = "tmvec", known_models()[["tmvec"]]),
        repdist = as.character(utils::packageVersion("repdist")),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
print(steps)
print(families)
