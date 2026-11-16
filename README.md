# repdist

[![R-CMD-check](https://github.com/Yixuan39/repdist/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Yixuan39/repdist/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/Yixuan39/repdist/branch/main/graph/badge.svg)](https://codecov.io/gh/Yixuan39/repdist)

**Measure, analyze and compare protein representation distances**

`repdist` represents each sample as an abundance-weighted distribution over
protein embeddings and compares samples with maximum mean discrepancy
(MMD), giving a beta-diversity metric that is aware of protein representation
similarity rather than exact sequence identity alone. The same embeddings also
cluster proteins into structural bins by Markov clustering (MCL).

The default `TM-Vec` model predicts structural similarity: cosine similarity
between protein embeddings approximates the TM-score.

![PCoA of 20 simulated samples under MMD, Bray-Curtis on proteins and Bray-Curtis on bins, each group's outermost samples joined](man/figures/README-simulation-pcoa.png)

Two groups of ten simulated samples share 100 background proteins of the NCBI
COG database. Each group-A sample also puts 5% of its abundance on a GrpE
protein (COG0576) of its own, and the ten GrpE proteins fall in ten different
30%-identity MMseqs2 clusters, so no protein or sequence family is enriched
in more than one sample. All three distances use relative abundances;
Bray-Curtis on bins first sums the proteins of each structural bin
(`bin_glom()`). PERMANOVA, 999 permutations:

| Distance | R² | p |
| --- | --- | --- |
| MMD | 0.150 | 0.001 |
| Bray-Curtis (protein) | 0.051 | 0.482 |
| Bray-Curtis (bins) | 0.165 | 0.001 |

Protein-level Bray-Curtis misses the enrichment; MMD and Bray-Curtis on bins
both detect it. This is the
[vignette](https://yixuan39.github.io/repdist/repdist.html)'s dataset,
`data(simulate_cog)`; the
[simulated data page](https://yixuan39.github.io/repdist/data_simulation.html)
shows how it was built.

## What you need

| Input | Typical file | Read with |
| --- | --- | --- |
| Protein sequences | FASTA, one record per protein, the protein ID as the first word of the header | `embed_proteins()` |
| Protein counts | tab-separated, samples in rows, sample IDs in the first column, protein IDs as column names (e.g. PSM or spectral counts) | `as.matrix(read.delim(f, row.names = 1, check.names = FALSE))` |
| Sample metadata | tab-separated, one row per sample, in the count table's row order | `read.delim()` |

`inst/extdata` ships one of each, 100 proteins of a mouse gut metaproteome in
23 samples on two diets; the end of the [Quick start](#quick-start) reads
them.

Python, conda and git do not need to be installed: the first
`embed_proteins()` call downloads a standalone Python and sets up a private,
version-pinned environment (up to a few GB, once), on Windows, macOS and
Linux.

**New to repdist? Start with the
[vignette](https://yixuan39.github.io/repdist/repdist.html)**
(`vignette("repdist")`), a worked example on simulated data with a known
answer: embeddings, beta diversity, protein bins and differential abundance.

## Simulated example data

```r
data(simulate_cog, package = "repdist")
dim(simulate_cog$counts) # 20 samples by 110 proteins
```

`simulate_cog` holds the protein sequences with precomputed TM-Vec
embeddings, simulated counts, sample groups, and COG annotations marking the
GrpE proteins and the sample each is spiked into. Loading it needs no model
download. [DATA_PROVENANCE.md](inst/scripts/DATA_PROVENANCE.md) explains how
it was built.

## Install

After acceptance into Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("repdist", dependencies = TRUE)
```

Until then, install the development version with Bioconductor dependencies.
You need R >= 4.5; Rtools is not needed on Windows. On Linux, R compiles some
dependencies from source and needs these system libraries first:

```sh
# Debian/Ubuntu
sudo apt-get install -y curl build-essential gfortran cmake libcurl4-openssl-dev libssl-dev libpng-dev libuv1-dev libbz2-dev libicu-dev
# Fedora/RHEL
sudo dnf install -y curl gcc-c++ gcc-gfortran cmake libcurl-devel openssl-devel libpng-devel libuv-devel bzip2-devel libicu-devel
```

Then, in R:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("basilisk", "BiocFileCache", "BiocParallel",
    "Biostrings", "pwalign", "S4Vectors", "MCL"))
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
# dependencies = TRUE also installs vegan and ggplot2 for the vignette, and Rtsne for plot_tsne()
remotes::install_github("Yixuan39/repdist", dependencies = TRUE)
```

## Quick start

This runs as shown, on the bundled `simulate_cog` data, with no model
download and no file of your own:

```r
library(repdist)

data(simulate_cog, package = "repdist")
proteins <- simulate_cog$proteins # AAStringSet, embeddings in mcols()$embedding
counts <- simulate_cog$counts     # 20 samples x 110 proteins
meta <- simulate_cog$samples      # one row per sample; `group` is "A" or "B"

# Beta diversity between samples
D <- sample_mmd(counts, proteins)
vegan::adonis2(D ~ group, data = meta)
plot_pcoa(counts, proteins, meta$group, drop = 0) # the same distance, plotted

# Structure-level protein bins, e.g. for biomarker discovery
S <- protein_similarity(proteins, max_evalue = 0.01)
bins <- bin_proteins(S, inflation = 2)
bin_counts <- bin_glom(counts, bins) # samples x bins
vegan::adonis2(vegan::vegdist(bin_counts / rowSums(bin_counts)) ~ group,
    data = meta) # Bray-Curtis on bins
```

`proteins` is an `AAStringSet` of the sequences with the embeddings in
`mcols(proteins)$embedding`; `simulate_cog` ships with its embeddings already
computed, so loading it needs no model download.

The bins need only `proteins`, not `counts` or any sample, and MMD needs no
bins, so each can be used without the other. `protein_similarity()` aligns
every pair of proteins with a positive embedding similarity, so its time
grows with the square of the number of proteins. A pair becomes an *edge*,
weighted by its predicted TM-score (0 to 1; above about 0.5 usually the same
fold), only if its alignment E-value is at most `max_evalue`; the others get
0. `bin_proteins()` runs MCL on the edges. A bin can hold pairs that are not
edges, joined only through shared neighbours. Structural similarity alone
does not guarantee shared function, so weigh annotation agreement against how
much of the catalog the non-singleton bins hold, and how far families
fragment, before transferring labels. `plot_pcoa()` with its default
`drop = 2` leaves out the most abundant protein regions and is for display
only; test on `sample_mmd()` itself.

The same from files, with the examples in `inst/extdata`; swap in your own
paths. `embed_proteins()` downloads the model on first use. The 100 proteins
were chosen for a GH92 family and the most abundant others, so the example
checks that the workflow runs; it is not a result. Eleven of the twelve mice
are sampled on both diets, so the test permutes within mouse:

```r
ex <- function(f) system.file("extdata", f, package = "repdist")
proteins <- embed_proteins(ex("metaproteomics.fasta"))
counts <- as.matrix(read.delim(ex("metaproteomics_counts.tsv"),
    row.names = 1, check.names = FALSE))
meta <- read.delim(ex("metaproteomics_samples.tsv"))
D <- sample_mmd(counts, proteins)
vegan::adonis2(D ~ diet, data = meta,
    permutations = permute::how(blocks = meta$subject, nperm = 999))
```

## License

Package code: MIT. Example data: the metaproteomics files are CC BY 4.0, and
`simulate_cog` comes from NCBI's COG database, which NCBI distributes without
restrictions; see
[the data provenance and attribution](inst/scripts/DATA_PROVENANCE.md).
