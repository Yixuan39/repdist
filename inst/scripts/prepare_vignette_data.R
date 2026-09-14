# Build the vignette's example dataset from case study 1 (diet metaproteome).
#
# The full study is 30,561 proteins x 111 samples; embedding all of them takes
# about an hour and the resulting matrix is ~120 MB, neither of which belongs in
# a package build. This keeps the 2000 most abundant proteins after rarefaction,
# which hold 80% of the total mass and reproduce the diet PERMANOVA (Bray-Curtis
# diet R2 0.275 against 0.253 on the full catalog).
#
# First run curate_diet_data.R; see DATA_PROVENANCE.md for source tables.
# Run from the repository root:  Rscript inst/scripts/prepare_vignette_data.R
# Output: inst/extdata/diet_metaproteome.rds

suppressMessages({
    library(phyloseq)
    library(repdist)
})

N_KEEP <- 2000L
HERE <- file.path("data", "test_study")
OUT <- file.path("inst", "extdata", "diet_metaproteome.rds")
# Resume point for a re-run; not shipped.
CACHE <- file.path(HERE, "vignette_embeddings.rds")

set.seed(1)
ps <- readRDS(file.path(HERE, "ps.rds"))

# Rarefy to the shallowest sample (drops none), then drop the pre-intervention
# baseline, for which the diet term is undefined. Rarefying first is deliberate:
# it subsamples raw counts, so it has to see them.
psd <- rarefy_even_depth(ps,
    sample.size = min(sample_sums(ps)), rngseed = 1,
    replace = FALSE, verbose = FALSE
)
psd <- prune_samples(sample_data(psd)$phase != "baseline", psd)
psd <- prune_taxa(taxa_sums(psd) > 0, psd)

psu <- prune_samples(sample_data(ps)$phase != "baseline", ps) # unrarefied twin
psu <- prune_taxa(taxa_sums(psu) > 0, psu)

cnt_all <- t(as(otu_table(psd), "matrix")) # samples x proteins
keep <- colnames(cnt_all)[order(-colSums(cnt_all))][seq_len(N_KEEP)]

counts <- cnt_all[, keep, drop = FALSE]
counts_raw <- t(as(otu_table(psu), "matrix"))[rownames(counts), keep, drop = FALSE]

md <- data.frame(sample_data(psd))[rownames(counts), ]
md$diet <- factor(md$diet)
md$dose <- factor(md$dose)
md <- md[, c("sample_id", "subject", "cage", "diet", "dose", "sex", "age_week", "period")]

seqs <- setNames(as.character(refseq(psd)[keep]), keep)

# Extended Data Table 6 supplies the curated consensus annotation; the lineage
# already lives in the phyloseq tax_table, joined there during data curation.
t6 <- read.delim(
    file.path(
        HERE, "annotated",
        "ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt"
    ),
    skip = 2, quote = "\"", check.names = FALSE,
    colClasses = "character", comment.char = "", fileEncoding = "latin1"
)
tx <- as.data.frame(as(tax_table(psd), "matrix"))[keep, ]
taxa <- data.frame(
    protein = keep,
    genus = tx$Genus,
    species = tx$Species,
    annotation = unname(setNames(
        t6[["Consensus annotation"]],
        t6[["Protein Identifier"]]
    )[keep]),
    length = Biostrings::width(refseq(psd)[keep]),
    row.names = keep, stringsAsFactors = FALSE
)
taxa$annotation[!nzchar(trimws(taxa$annotation))] <- NA_character_

# Embedding is the slow step, so it is cached outside the knitr/build cycle and
# saved in blocks -- an interrupted run resumes from the last one.
spec <- known_models()[["tmvec"]]
fingerprint <- digest::digest(list(seqs, spec, packageVersion("repdist")))
cache <- if (file.exists(CACHE)) readRDS(CACHE) else NULL
if (!is.null(cache) && !identical(cache$fingerprint, fingerprint)) {
    stop("Embedding cache does not match sequences, model or package version.")
}
cached <- cache$embeddings
need <- setdiff(keep, rownames(cached))
if (length(need)) {
    s <- seqs[need]
    s <- s[order(nchar(s))] # length-sorted keeps batches tight
    for (i in split(seq_along(s), ceiling(seq_along(s) / 500))) {
        message(sprintf("embedding %d/%d", max(i), length(s)))
        cached <- rbind(cached, embed_proteins(s[i]))
        saveRDS(list(fingerprint = fingerprint, embeddings = cached), CACHE)
    }
}
Z <- cached[keep, , drop = FALSE]

# Round for storage; callers renormalize before cosine-distance calculations.
# This is lossy compression and can affect pairs close to a bin threshold.
Z <- round(Z, 4)

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
saveRDS(
    list(
        counts = counts, counts_raw = counts_raw, metadata = md,
        embeddings = Z, taxa = taxa, sequences = seqs,
        provenance = list(
            source = "https://doi.org/10.1093/ismejo/wraf048",
            license = "CC BY 4.0",
            counts = "Data Set 1; microbial PSM counts",
            annotation = "Data Set 6; lineage and consensus annotation",
            processing = "Full-catalog rarefaction seed 1, remove baseline, top 2000",
            embedding_model = spec,
            package_version = as.character(packageVersion("repdist")),
            session = capture.output(sessionInfo()),
            rounding_digits = 4L
        )
    ),
    OUT,
    compress = "xz"
)

message(sprintf(
    "wrote %s (%.2f MB): %d samples x %d proteins",
    OUT, file.size(OUT) / 1048576, nrow(counts), ncol(counts)
))
