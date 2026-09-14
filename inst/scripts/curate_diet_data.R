# Reconstruct the diet study intermediate from public supplementary tables.
# See DATA_PROVENANCE.md for downloads, names, license and interpretation.
# Run from the source repository root; no network access or embedding here.
suppressPackageStartupMessages(library(phyloseq))
HERE <- file.path("data", "test_study")

raw <- read.delim(file.path(HERE, "ExtendedDataTable1CompleteMetaproteome.txt"),
    skip = 2, quote = "\"", check.names = FALSE,
    colClasses = "character", comment.char = ""
)
names(raw)[1:5] <- c("code", "protein_id", "aa_length", "description", "sequence")
raw <- raw[raw$code == "Microbiome", ]

counts <- vapply(raw[, -(1:5)], function(x) {
    blank <- !nzchar(trimws(x)) # blank cell = protein not seen
    v <- suppressWarnings(as.numeric(x))
    bad <- is.na(v) & !blank # a parse failure that ISN'T blank is
    if (any(bad)) { # unexpected data, not a documented zero
        stop("Non-numeric, non-blank count cell(s): ",
            paste(unique(x[bad]), collapse = ", "),
            call. = FALSE
        )
    }
    v[blank] <- 0
    v
}, numeric(nrow(raw)))
stopifnot(all(counts >= 0), all(counts == round(counts))) # PSM counts only
storage.mode(counts) <- "integer"
dimnames(counts) <- list(raw$protein_id, names(raw)[-(1:5)])

c(proteins = nrow(counts), samples = ncol(counts), total_psm = sum(counts))

diets <- data.frame(
    token = c(
        "T0", "20SOY", "20CAS", "20RIC", "40SOY", "20YST", "40CAS", "20PEA",
        "20EGG", "CTL"
    ),
    period = 0:9,
    diet = c(NA, "Soy", "Casein", "Rice", "Soy", "Yeast", "Casein", "Pea", "Egg", NA),
    dose = c(NA, "20", "20", "20", "40", "20", "40", "20", "20", "20"),
    phase = c("baseline", rep("treatment", 8), "control")
)

f <- strsplit(colnames(counts), "_")
cage <- vapply(f, `[`, "", 1)
token <- vapply(f, `[`, "", 3)
stopifnot(token %in% diets$token)

md <- data.frame(diets[match(token, diets$token), c("period", "phase", "diet", "dose")],
    sample_id = colnames(counts), cage = cage,
    subject = paste(cage, vapply(f, `[`, "", 2), sep = "_"),
    # sex is not in the name and is not free to vary; it is also the
    # paper's "mouse group", the sexes having come from different rooms
    sex = ifelse(cage %in% c("C1", "C3"), "Male", "Female"),
    row.names = colnames(counts)
)
# control re-feeds take their source from the cage
md$diet[md$phase == "control"] <- ifelse(cage[md$phase == "control"] %in% c("C1", "C2"),
    "Soy", "Casein"
)
# adonis2 on a continuous term is invariant to a constant shift, so weeks-since-T0
# stands in for age. The bone week still elapsed, so CTL sits two weeks after period 8.
md$age_week <- ifelse(md$period == 9, 10L, md$period)

table(md$period, md$diet, useNA = "ifany")

bin <- sub("_[0-9]+$", "", rownames(counts))

t6 <- read.delim(file.path(HERE, "annotated", "ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt"),
    skip = 2, quote = "\"", check.names = FALSE, colClasses = "character",
    comment.char = "", fileEncoding = "latin1"
)
lineage <- setNames(t6[["Taxa Lineage"]], t6[["Protein Identifier"]])[rownames(counts)]

ranks <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
rk <- t(vapply(strsplit(lineage, ";"), function(x) {
    x <- sub("^[a-z]__", "", x)
    x[!nzchar(x)] <- NA
    length(x) <- length(ranks)
    x
}, character(length(ranks))))
dimnames(rk) <- list(rownames(counts), ranks)

tax <- cbind(Bin = bin, rk)
c(
    bins = length(unique(bin)), genera = length(unique(na.omit(tax[, "Genus"]))),
    unmatched = sum(is.na(lineage)), proteins = nrow(tax)
)

ps <- phyloseq(
    otu_table(counts, taxa_are_rows = TRUE),
    sample_data(md),
    tax_table(tax),
    refseq(Biostrings::AAStringSet(setNames(raw$sequence, raw$protein_id)))
)
ps

saveRDS(ps, file.path(HERE, "ps.rds"))
