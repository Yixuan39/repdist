# Writes inst/extdata/metaproteomics.fasta, 100 proteins of the diet study,
# and metaproteomics_{counts,samples}.tsv, their PSM counts and sample
# metadata on the Yeast and Casein diets. See DATA_PROVENANCE.md. Run from the
# source repository root with repdist installed.
#
# Usage: Rscript inst/scripts/make_metaproteomics_example.R [diet_series.rds]

library(repdist) # also makes the readRDS() AAStringSet subsettable
src <- commandArgs(TRUE)[1]
if (is.na(src)) {
    src <- BiocFileCache::bfcrpath(BiocFileCache::BiocFileCache(ask = FALSE),
        paste0("https://github.com/Yixuan39/repdist/raw/gh-pages/",
            "curated/diet_series.rds"))
}
d <- readRDS(src)

# the 12 GH92 proteins seen on the Yeast and Casein diets, then the 88 most
# abundant other proteins on those diets
yc <- colSums(d$counts[d$samples$diet %in% c("20CAS", "20YST"), ])
gh92 <- names(yc)[yc > 0 & grepl("GH92", d$annotation[names(yc), "cazy"])]
rest <- setdiff(names(sort(yc, decreasing = TRUE)), gh92)
ids <- c(gh92, rest[seq_len(100 - length(gh92))])

x <- d$proteins[ids]
S4Vectors::mcols(x) <- NULL
Biostrings::writeXStringSet(x, "inst/extdata/metaproteomics.fasta")

s <- d$samples[d$samples$diet %in% c("20CAS", "20YST"), ]
utils::write.table(d$counts[rownames(s), ids],
    "inst/extdata/metaproteomics_counts.tsv",
    sep = "\t", quote = FALSE, col.names = NA)
utils::write.table(s, "inst/extdata/metaproteomics_samples.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE)
