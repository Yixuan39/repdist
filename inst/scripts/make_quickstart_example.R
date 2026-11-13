# Writes inst/extdata/quickstart_{proteins.fasta,counts.tsv,metadata.tsv}: a
# small trio of typical repdist input files (sequences with no embeddings, a
# counts table, sample metadata), cut from data(simulate_ec) for the README
# and vignette file-reading examples. See DATA_PROVENANCE.md. Run from the
# source repository root with repdist installed.
#
# Usage: Rscript inst/scripts/make_quickstart_example.R

library(repdist)
data(simulate_ec, package = "repdist")
s <- simulate_ec

samples <- c(s$samples$sample[s$samples$group == "A"][1:3],
    s$samples$sample[s$samples$group == "B"][1:3])
enriched <- s$annotation$protein[s$annotation$enriched_in %in% samples]
set.seed(1)
background <- sample(setdiff(s$annotation$protein, enriched), 18L)
ids <- c(enriched, background)

counts <- s$counts[samples, ids]
metadata <- s$samples[samples, c("sample", "group")]
seqs <- s$proteins[ids]
S4Vectors::mcols(seqs) <- NULL # raw sequences only, as if never embedded

Biostrings::writeXStringSet(seqs,
    "inst/extdata/quickstart_proteins.fasta")
utils::write.table(counts, "inst/extdata/quickstart_counts.tsv",
    sep = "\t", quote = FALSE, col.names = NA)
utils::write.table(metadata, "inst/extdata/quickstart_metadata.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE)
