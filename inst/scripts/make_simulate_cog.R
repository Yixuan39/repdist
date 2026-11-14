# Writes data/simulate_cog.rda and inst/extdata/simulate_cog.fasta: the
# one-COG test dataset of the website (data_simulation.Rmd on the gh-pages
# branch) with TM-Vec embeddings. 20 samples, 110 proteins: 100 background
# proteins in every sample, and 10 GrpE (COG0576) proteins from 10 different
# MMseqs2 30% clusters, each spiked into one group-A sample with 5% of its
# abundance. See DATA_PROVENANCE.md. Run from the source repository root with
# repdist installed.
#
# Usage: Rscript inst/scripts/make_simulate_cog.R [path/to/data_simulation]

library(repdist) # also makes the readRDS() AAStringSet subsettable
src <- commandArgs(TRUE)[1]
files <- c("proteins.fasta", "labels.tsv", "samples.tsv", "counts.tsv")
path <- if (is.na(src)) {
    bfc <- BiocFileCache::BiocFileCache(ask = FALSE)
    vapply(files, function(f) BiocFileCache::bfcrpath(bfc, paste0(
        "https://github.com/Yixuan39/repdist/raw/gh-pages/data_simulation/", f)), "")
} else {
    setNames(file.path(src, files), files)
}

proteins <- Biostrings::readAAStringSet(path[["proteins.fasta"]])
ids <- names(proteins)
annotation <- read.delim(path[["labels.tsv"]], na.strings = "")
samples <- read.delim(path[["samples.tsv"]], na.strings = "")
stopifnot(identical(annotation$protein, ids))
rownames(annotation) <- ids
rownames(samples) <- samples$sample
counts <- as.matrix(read.delim(path[["counts.tsv"]], row.names = 1, check.names = FALSE))
storage.mode(counts) <- "integer"
stopifnot(identical(colnames(counts), ids), identical(rownames(counts), samples$sample))

# TM-Vec embeddings, length-sorted batches; rounded to 4 decimals for size
basilisk::setBasiliskShared(FALSE)
basilisk::setBasiliskFork(FALSE)
o <- order(Biostrings::width(proteins))
Z <- S4Vectors::mcols(embed_proteins(proteins[o]))$embedding
S4Vectors::mcols(proteins)$embedding <- round(Z[ids, ], 4)

simulate_cog <- list(
    proteins = proteins,
    counts = counts,
    samples = samples,
    annotation = annotation,
    provenance = list(
        source = "NCBI COG 2024 release, https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/",
        license = paste("NCBI places no restrictions on the use or distribution",
            "of its data: https://www.ncbi.nlm.nih.gov/home/about/policies/"),
        simulation = "data_simulation.Rmd (gh-pages branch)",
        model = "embed_proteins(model = 'tmvec')",
        design = list(treatment_cog = "COG0576", n_per_group = 10L,
            n_background = 100L, rho = 0.05, depth = 2e4,
            seeds = c(proteins = 20261112, counts = 1)),
        script_md5 = unname(tools::md5sum(sub("^--file=", "",
            grep("^--file=", commandArgs(), value = TRUE)))),
        repdist = as.character(utils::packageVersion("repdist")),
        built = format(Sys.Date()),
        session = utils::capture.output(utils::sessionInfo())
    )
)
dir.create("data", showWarnings = FALSE)
save(simulate_cog, file = file.path("data", "simulate_cog.rda"), compress = "xz")
fasta <- proteins
S4Vectors::mcols(fasta) <- NULL
Biostrings::writeXStringSet(fasta, file.path("inst", "extdata", "simulate_cog.fasta"))
