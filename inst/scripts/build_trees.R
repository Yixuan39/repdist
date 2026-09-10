# Builds the two trees the beta-diversity comparisons need, for one dataset:
# one fitted to sequence composition, one to the TM-Vec structure metric.
# Expensive (tens of minutes, gigabytes), so it is a cached script rather than a
# knitr chunk.
#
#   Rscript inst/scripts/build_trees.R test_study
#   Rscript inst/scripts/build_trees.R case_study2
#
# Both trees use the same complete-linkage fit, so the
# comparison downstream isolates the representation (sequence vs structure)
# rather than confounding it with the tree-fitting algorithm. Complete linkage
# also preserves bin_proteins()'s all-pairs TM-score threshold guarantee.
#
# The sequence tree is alignment-free (3-mer composition), NOT an inferred
# phylogeny. That is deliberate, not a shortcut: these ~19k proteins are the
# functional repertoire of a whole community, not a homologous family, so a
# global MSA of them is ~99% gaps (longest sequence 32,877aa against a median of
# 365) and any tree inferred from it is noise. An alignment-free composition
# distance is the honest sequence-space counterpart to the structure-space one,
# and it is the same functional form -- cosine on a profile vector -- in both.

suppressMessages({
  library(phyloseq)
  library(Matrix)
  library(repdist)
})

dataset <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(dataset %in% c("test_study", "case_study2"))

path  <- here::here("data", dataset)
out   <- file.path(path, "trees_kmer_struct.rds")
K <- 3L  # 3-mers: 20^3 = 8000 features, ~360 nonzero for a median protein

set.seed(1)
cluster_tree <- function(D) ape::as.phylo(fastcluster::hclust(D, method = "complete"))

# ---- the same psd the corresponding Rmd builds --------------------------
MAX_AA <- 3000L
ps <- readRDS(file.path(path, "ps.rds"))
if (dataset == "case_study2") {
  ps <- prune_samples(data.frame(sample_data(ps))$substrate != "PLA", ps)
  ps <- prune_taxa(taxa_sums(ps) > 0, ps)
}
psd <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 1,
                         replace = FALSE, verbose = FALSE)
if (dataset == "test_study")
  psd <- prune_samples(sample_data(psd)$phase != "baseline", psd)
psd  <- prune_taxa(taxa_sums(psd) > 0, psd)
psd  <- prune_taxa(Biostrings::width(refseq(psd)) <= MAX_AA, psd)
prot <- taxa_names(psd)
message(sprintf("%s: %d proteins, %d samples", dataset, length(prot), nsamples(psd)))

# ---- sequence side: 3-mer composition profile ---------------------------
AA   <- strsplit("ACDEFGHIKLMNPQRSTVWY", "")[[1]]
seqs <- as.character(refseq(psd)[prot])

kmer_rows <- function(s, idx) {
  a <- match(strsplit(s, "", fixed = TRUE)[[1]], AA)   # non-standard residues -> NA
  n <- length(a)
  if (n < K) return(NULL)
  j <- (a[1:(n - 2L)] - 1L) * 400L + (a[2:(n - 1L)] - 1L) * 20L + a[3:n]
  j <- j[!is.na(j)]                                     # drop k-mers spanning X/B/Z/U
  if (!length(j)) return(NULL)
  tj <- tabulate(j, nbins = 8000L)
  nz <- which(tj > 0L)
  cbind(idx, nz, tj[nz])
}

message("counting k-mers ...")
trip <- do.call(rbind, lapply(seq_along(seqs), function(i) kmer_rows(seqs[i], i)))
Kmat <- sparseMatrix(i = trip[, 1], j = trip[, 2], x = as.numeric(trip[, 3]),
                     dims = c(length(prot), 8000L), dimnames = list(prot, NULL))
rm(trip); invisible(gc(verbose = FALSE))
message(sprintf("k-mer matrix: %d x %d, %.1f%% dense", nrow(Kmat), ncol(Kmat),
                100 * nnzero(Kmat) / prod(dim(Kmat))))

# cosine distance, the same form repdist_matrix() uses on embeddings
nrm <- sqrt(Matrix::rowSums(Kmat^2))
stopifnot(all(nrm > 0))
Kmat  <- Kmat / nrm
D_seq <- as.dist(1 - as.matrix(Matrix::tcrossprod(Kmat)))
rm(Kmat); invisible(gc(verbose = FALSE))

message("fitting sequence tree ...")
tree_seq <- cluster_tree(D_seq)
fid_seq  <- cor(as.vector(D_seq),
                as.vector(as.dist(ape::cophenetic.phylo(tree_seq)[prot, prot])))
rm(D_seq); invisible(gc(verbose = FALSE))
message(sprintf("sequence tree fidelity: %.4f", fid_seq))

# ---- structure side: TM-Vec cosine --------------------------------------
Z <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))$Z[prot, , drop = FALSE]
D_str <- repdist_matrix(Z)

message("fitting structure tree ...")
tree_str <- cluster_tree(D_str)
fid_str  <- cor(as.vector(D_str),
                as.vector(as.dist(ape::cophenetic.phylo(tree_str)[prot, prot])))
rm(D_str); invisible(gc(verbose = FALSE))
message(sprintf("structure tree fidelity: %.4f", fid_str))

saveRDS(list(tree_seq = tree_seq, tree_str = tree_str,
             fidelity = c(sequence = fid_seq, structure = fid_str),
             k = K, proteins = prot), out)
message("wrote ", out)
