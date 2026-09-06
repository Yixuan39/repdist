# repdist

**Measure, analyze and compare protein representation distances**

`repdist` represents each sample as an abundance-weighted distribution over
protein embeddings and compares samples with maximum mean discrepancy
(MMD), giving a beta-diversity metric that is aware of protein representation
similarity rather than exact sequence identity alone.

The package use `TM-Vec` by default, which the protein and sample level embedding
distances represent structural differences. The cosine distance between a pair of 
protein embeddings is defined as an approximation of the TM-score.

![PCoA of 20 simulated samples under four distances: RBF-MMD, structure-tree UniFrac, bin-level Bray-Curtis and protein-level Bray-Curtis](man/figures/README-simulation-pcoa.png)

Two groups of simulated samples differing by *function*, carried on
non-overlapping homologs. MMD, a structure tree and structural bins all separate
them; protein-level Bray-Curtis, which can only ask whether two samples hold the
same accession, sees nothing. Built in `vignette("simulation", package = "repdist")`.

## Install

```r
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

For annotation, compare graph and density clusters on the same distance:

```r
graph <- bin_proteins(G, method = "mcl", min_sim = 0.5, inflation = 2)
dense <- bin_proteins(G, method = "dbscan", min_sim = 0.7, min_pts = 5)
adaptive <- bin_proteins(G, method = "hdbscan", min_pts = 5)
```

MCL needs the [MCL executable](https://micans.org/mcl/) on `PATH` (or set
`mcl_bin`); density methods need `install.packages("dbscan")`. MCL uses
`min_sim` to filter graph edges and inflation to control granularity; DBSCAN
uses radius `1 - min_sim` and defaults to excluding border points (DBSCAN*).
HDBSCAN selects clusters across density levels. Only hclust guarantees an
all-pairs similarity floor. Structural similarity alone does not guarantee
shared function, so compare annotation agreement with coverage and fragmentation
before transferring labels. Density noise remains in separate singleton bins
and is listed in `$noise`; alternative methods return `$tree = NULL`.
`plot_bin_similarity()` supports all methods; `plot_bin_profile()` needs hclust.

Run `Rscript analysis/annotation_clustering.R` for the local EC benchmark grid,
including purity, pair precision/recall, coverage, and protein-level membership
across settings. The input data are local research files, not shipped with the
package; see [the comparison notes](analysis/annotation_clustering.md).

## License

MIT
