# repdist

**Measure, analyze and compare protein representation distances**

`repdist` represents each sample as an abundance-weighted distribution over
protein embeddings and compares samples with maximum mean discrepancy
(MMD), giving a beta-diversity metric that is aware of protein representation
similarity rather than exact sequence identity alone. The same embeddings also
cluster proteins into structural bins by Markov clustering (MCL).

The default `TM-Vec` model predicts structural similarity: cosine similarity
between protein embeddings approximates the TM-score. The corresponding cosine
distance is `1 - similarity`.

![PCoA of 20 simulated samples under three distances: RBF-MMD, MCL-bin Bray-Curtis and protein-level Bray-Curtis](man/figures/README-simulation-pcoa.png)

Two groups of simulated samples differing by *function*, carried on
non-overlapping homologs. MMD and MCL bins both separate them; protein-level
Bray-Curtis, which can only ask whether two samples hold the same accession,
sees nothing. Built in the
[simulation study](https://yixuan39.github.io/repdist/simulation.html) on the
package website.

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
BiocManager::install(c("basilisk", "BiocFileCache", "BiocParallel",
    "Biostrings", "pwalign", "rhdf5", "S4Vectors", "MCL"))
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Yixuan39/repdist")
```

## Quick start

```r
library(repdist)

x <- embed_proteins("proteins.fasta", model = "tmvec")
# or, from the tmvec command-line tool:
# x <- read_embeddings("db/proteins.npz", seqs = "proteins.fasta")

# Beta diversity between samples
D <- sample_mmd(counts, x)
vegan::adonis2(D ~ condition, data = meta)

# Structure-level protein clusters, e.g. for biomarker discovery
S <- protein_similarity(x, min_sim = 0.7, min_coverage = 0.5)
bins <- bin_proteins(S, inflation = 2)
bin_glom(counts, bins)            # samples x bins
```

`x` is an `AAStringSet` of the sequences with the embeddings in
`mcols(x)$embedding`. `protein_similarity()` sets pairs below `min_sim` to 0,
and pairs whose alignment covers less than `min_coverage` of either sequence.
`bin_proteins()` runs MCL on the remaining edges. The floor applies to edges,
not all pairs: a bin can hold pairs below it. Structural similarity alone does
not guarantee shared function, so weigh annotation agreement against how much of the catalog the
non-singleton bins hold, and how far families fragment, before transferring
labels.

## License

Package code: MIT. Example data: CC BY 4.0; see
[the data provenance and attribution](inst/scripts/DATA_PROVENANCE.md).
