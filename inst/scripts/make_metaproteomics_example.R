# Writes inst/extdata/metaproteomics.{npz,fasta}: 100 proteins of the diet
# study, as a tmvec build-db database plus the FASTA it was built from. See
# DATA_PROVENANCE.md. Run from the source repository root with repdist
# installed; the npz is written with numpy from repdist's own environment.
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
Z <- S4Vectors::mcols(x)$embedding
S4Vectors::mcols(x) <- NULL
Biostrings::writeXStringSet(x, "inst/extdata/metaproteomics.fasta")

basilisk::basiliskRun(env = repdist:::repdist_env, fun = function(ids, Z) {
    np <- reticulate::import("numpy", convert = FALSE)
    np$savez_compressed("inst/extdata/metaproteomics.npz",
        headers = np$array(ids, dtype = "U"),
        # row-major, as tmvec writes it
        embeddings = np$ascontiguousarray(Z, dtype = "float32"),
        input_fasta = np$array("metaproteomics.fasta", dtype = "U"),
        tm_vec_weights = np$array("", dtype = "U"),
        protrans_model_path = np$array("", dtype = "U"))
    NULL
}, ids = ids, Z = unname(Z))
