# Build curated/diet_series.rds: every sample of the mouse dietary-protein
# metaproteome series of Blakeley-Ruiz et al. (2025),
# ISME J 19(1), wraf048 (doi:10.1093/ismejo/wraf048; PRIDE PXD041586),
# with protein embeddings computed from the published sequences.
#
# Run from the root of the gh-pages branch with repdist installed:
#   Rscript data-raw/diet_series.R
# The two supplementary tables (CC BY 4.0) are downloaded from Europe PMC into
# data/wraf048/ (not committed) and checked against their MD5 sums.
library(repdist)

out <- file.path("curated", "diet_series.rds")
src <- file.path("data", "wraf048")
MAX_AA <- 1500L # analysis choice, see the real-data page

# 1. Source tables -------------------------------------------------------------
# Data Set 1: PSM counts per protein and sample, with sequences.
# Data Set 6: annotation of microbial proteins (columns 1-19 only; its %NSAF
# abundances are not used).
files <- c(
    table1 = "data_set_1_wraf048.txt", # 028c39efa8ec944484c093bf48c4365b
    table6 = "data_set_6_wraf048.txt" # e571367e516a1139ae002b883a27e1f0
)
md5 <- c(table1 = "028c39efa8ec944484c093bf48c4365b",
    table6 = "e571367e516a1139ae002b883a27e1f0")
dir.create(src, showWarnings = FALSE, recursive = TRUE)
for (k in names(files)) {
    f <- file.path(src, files[[k]])
    if (!file.exists(f)) {
        curl::curl_download(paste0(
            "https://europepmc.org/api/fulltextRepo?pprId=PMC12066410",
            "&type=FILE&mimeType=text/plain&fileName=", files[[k]]), f)
    }
    stopifnot(unname(tools::md5sum(f)) == md5[[k]])
}
read_table <- function(f) {
    utils::read.delim(file.path(src, f), skip = 2, quote = "\"",
        check.names = FALSE, colClasses = "character", comment.char = "",
        fileEncoding = "latin1")
}

# 2. Counts and sequences (Data Set 1, microbial proteins) ---------------------
t1 <- read_table(files[["table1"]])
names(t1)[1:5] <- c("code", "protein", "aa_length", "prokka", "sequence")
t1 <- t1[t1$code == "Microbiome", ]
counts <- vapply(t1[, -(1:5)], function(v) {
    blank <- !nzchar(trimws(v)) # a blank cell is a protein not seen
    n <- suppressWarnings(as.numeric(v))
    stopifnot(!anyNA(n[!blank]), all(n[!blank] >= 0), all(n[!blank] %% 1 == 0))
    n[blank] <- 0
    n
}, numeric(nrow(t1)))
rownames(counts) <- t1$protein
counts <- t(counts) # samples in rows

# 3. Samples ----------------------------------------------------------------------
# Names are {cage}_{mouse}_{diet}, e.g. C2_1_20YST. Sex is fixed by cage
# (C1, C3 male; C2, C4 female). Period follows the published feeding order.
parts <- regmatches(rownames(counts),
    regexec("^(C[1-4])_([1-3])_([0-9A-Z]+)_?$", rownames(counts)))
stopifnot(lengths(parts) == 4L)
samples <- data.frame(
    sample = rownames(counts),
    subject = vapply(parts, function(p) paste(p[2], p[3], sep = "_"), ""),
    cage = vapply(parts, `[`, "", 2),
    token = vapply(parts, `[`, "", 4)
)
feeding_order <- c("T0", "20SOY", "20CAS", "20RIC", "40SOY", "20YST", "40CAS",
    "20PEA", "20EGG", "CTL") # the 20% chicken-bone week yielded no samples
stopifnot(samples$token %in% feeding_order)
# Diet tokens as in the sample names (percentage and protein source), in
# feeding order.
samples$diet <- factor(samples$token, feeding_order)
samples$period <- match(samples$token, feeding_order) - 1L
samples$sex <- factor(ifelse(samples$cage %in% c("C1", "C3"), "male", "female"))
samples$subject <- factor(samples$subject)
samples$token <- NULL
samples <- samples[order(samples$diet, samples$subject), ]
rownames(samples) <- samples$sample
counts <- counts[samples$sample, ]

# 4. Protein universe -----------------------------------------------------------
# Every microbial protein with at least one PSM in any sample, up to MAX_AA.
len <- nchar(t1$sequence)
seen <- colSums(counts) > 0
steps <- c(detected = sum(seen), longer_than_max = sum(seen & len > MAX_AA))
keep <- seen & len <= MAX_AA
long_share <- sum(counts[, seen & len > MAX_AA]) / sum(counts[, seen])
counts <- counts[, keep]
storage.mode(counts) <- "integer"
t1 <- t1[keep, ]
steps[["proteins"]] <- ncol(counts)

# 5. Annotation (Data Set 6) ----------------------------------------------------
# Join on Protein Identifier, else on the species-code identifier, keeping a
# match only where the residue counts agree.
t6 <- read_table(files[["table6"]])[, 1:19]
len <- nchar(t1$sequence)
agrees <- function(r) !is.na(r) & suppressWarnings(as.integer(t6$AA[r])) == len
r_id <- match(t1$protein, t6[["Protein Identifier"]])
r_code <- match(t1$protein, t6[["Protein Identifier with Species Code"]])
row <- ifelse(agrees(r_id), r_id, ifelse(agrees(r_code), r_code, NA))
field <- function(k) {
    v <- trimws(t6[[k]][row])
    v[is.na(v) | !nzchar(v) | v == "Unknown"] <- NA
    v
}
# The first token of an identifier is its species code: a dRep species group
# (AB<n>-<n>), BT for B. thetaiotaomicron, or Unbin/amb/lowQ for unbinned,
# ambiguous and low-quality-bin proteins, which have no species.
species_code <- sub("_.*", "", t1$protein)
code_lineage <- tapply(t6[["Taxa Lineage"]],
    sub("_.*", "", t6[["Protein Identifier with Species Code"]]),
    function(v) unique(c(stats::na.omit(v), NA))[1])
lineage <- unname(code_lineage[species_code])
lineage[species_code == "BT"] <- paste0("d__Bacteria;p__Bacteroidota;",
    "c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;",
    "s__Bacteroides thetaiotaomicron")
# Readable name in the study's own style: the species where GTDB names one,
# else "<genus> sp. <code>" or "<higher taxon> bacterium <code>".
taxon <- vapply(seq_along(lineage), function(i) {
    if (is.na(lineage[i])) return(NA_character_)
    ranks <- sub("^[a-z]__", "", strsplit(lineage[i], ";")[[1]])
    named <- which(nzchar(ranks))
    last <- max(named)
    if (last == 7L) return(ranks[7])
    paste(ranks[last], if (last == 6L) "sp." else "bacterium", species_code[i])
}, character(1))
cazy <- vapply(strsplit(t6$CAZy[row], ";"), function(v) {
    v <- unique(v[!is.na(v) & nzchar(v) & v != "-"])
    if (length(v)) paste(v, collapse = ";") else NA_character_
}, character(1))
ko_mantis <- field("MantisKO")
ko_microbeannotator <- field("MicrobeKO")
annotation <- data.frame(
    protein = t1$protein,
    length = len,
    species_code = species_code,
    taxon = taxon,
    lineage = lineage,
    prokka = t1$prokka,
    table6_match = ifelse(is.na(row), NA,
        ifelse(agrees(r_id), "Protein Identifier", "species-code identifier")),
    broad_function = field("Broad Function"),
    detailed_function = field("Detailed Function"),
    consensus_annotation = field("Consensus annotation"),
    ko_mantis = ko_mantis,
    ko_microbeannotator = ko_microbeannotator,
    ko = ifelse(is.na(ko_mantis), ko_microbeannotator, ko_mantis),
    cazy = cazy,
    row.names = t1$protein
)
steps[["with_table6_row"]] <- sum(!is.na(row))
steps[["with_ko"]] <- sum(!is.na(annotation$ko))

# 6. Embed every sequence with the default model --------------------------------
# In chunks of 500, each in a fresh worker process: on Apple GPUs (MPS) one
# long run can exhaust device memory and fall back to the CPU for most
# sequences. Sequences are embedded shortest first so that a batch pads to a
# similar length; neither the chunks nor the order change an embedding, only
# the run time.
embed_chunked <- function(seqs, chunk = 500L) {
    basilisk::setBasiliskShared(FALSE)
    basilisk::setBasiliskFork(FALSE)
    o <- order(nchar(seqs))
    parts <- lapply(split(o, ceiling(seq_along(o) / chunk)), function(i) {
        message("embedding ", length(i), " more, up to ", max(nchar(seqs[i])),
            " residues")
        S4Vectors::mcols(embed_proteins(seqs[i], model = "tmvec"))$embedding
    })
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- do.call(rbind, unname(parts))[names(seqs), ]
    x
}
x <- embed_chunked(setNames(t1$sequence, t1$protein))
stopifnot(identical(names(x), colnames(counts)),
    identical(rownames(S4Vectors::mcols(x)$embedding), colnames(counts)))
# Six decimals keep the file small; the cosine error is below 1e-5.
S4Vectors::mcols(x)$embedding <- round(S4Vectors::mcols(x)$embedding, 6)

saveRDS(list(
    proteins = x,
    counts = counts,
    samples = samples,
    annotation = annotation,
    provenance = list(
        study = paste("Blakeley-Ruiz JA et al. (2025) Dietary protein source",
            "alters gut microbiota composition and function. ISME J 19(1):",
            "wraf048. doi:10.1093/ismejo/wraf048"),
        source = c(
            counts = "Supplementary Data Set 1 (Extended Data Table 1)",
            annotation = "Supplementary Data Set 6 (Extended Data Table 6)"),
        source_url = "https://europepmc.org/article/MED/40116459",
        source_md5 = md5,
        license = "CC BY 4.0",
        pride = "PXD041586",
        diets = feeding_order,
        max_aa = MAX_AA,
        long_protein_psm_share = long_share,
        steps = steps,
        ko_rule = "MantisKO, falling back to MicrobeKO (MicrobeAnnotator)",
        model = c(key = "tmvec", known_models()[["tmvec"]]),
        repdist = as.character(utils::packageVersion("repdist")),
        built = format(Sys.Date()),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
print(steps)
print(table(samples$diet))
