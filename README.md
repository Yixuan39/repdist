# repdist

**Measure, analyze and compare protein representation distances**

`repdist` represents each sample as an abundance-weighted distribution over
protein embeddings and compares samples with maximum mean discrepancy
(MMD), giving a beta-diversity metric that is aware of protein representation
similarity rather than exact sequence identity alone. The same protein
distances also cluster proteins into structural bins, with quality-control
plots for choosing the similarity threshold and reading a bin.

The default `TM-Vec` model predicts structural similarity: cosine similarity
between protein embeddings approximates the TM-score. The corresponding cosine
distance is `1 - similarity`.

![PCoA of 20 simulated samples under four distances: RBF-MMD, structure-tree UniFrac, bin-level Bray-Curtis and protein-level Bray-Curtis](man/figures/README-simulation-pcoa.png)

Two groups of simulated samples differing by *function*, carried on
non-overlapping homologs. MMD, a structure tree and structural bins all separate
them; protein-level Bray-Curtis, which can only ask whether two samples hold the
same accession, sees nothing. Built in `vignette("simulation", package = "repdist")`.

## Install

After acceptance into Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("repdist")
```

Until then, install the development version with Bioconductor dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("basilisk", "BiocFileCache", "Biostrings"))
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Yixuan39/repdist")
```

## Quick start

```r
library(repdist)

seqs <- read_fasta("proteins.fasta")
Z <- embed_proteins(seqs, model = "tmvec")

# Beta diversity between samples
D <- sample_repdist(counts, Z)
vegan::adonis2(D ~ condition, data = meta)

# Structure-level protein clusters, e.g. for biomarker discovery
G <- repdist_matrix(Z)                    # 1 - predicted TM-score
bins <- bin_proteins(G, min_sim = 0.7)
```

`bin_proteins()` cuts a complete-linkage tree, so every pair inside a bin meets
`min_sim`. Structural similarity alone does not guarantee shared function, so
weigh annotation agreement against how much of the catalog the non-singleton
bins hold, and how far families fragment, before transferring labels.

## License

Package code: MIT. Example data: CC BY 4.0; see
[the data provenance and attribution](inst/scripts/DATA_PROVENANCE.md).
