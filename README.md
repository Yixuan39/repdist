# repdist

**Measure, analyze and compare protein representation distances**

`repdist` represents each sample as an abundance-weighted distribution over
protein embeddings and compares samples with maximum mean discrepancy
(MMD), giving a beta-diversity metric that is aware of protein representation
similarity rather than exact sequence identity alone.

The package use `TM-Vec` by default, which the protein and sample level embedding
distances represent structural differences. The cosine distance between a pair of 
protein embeddings is defined as an approximation of the TM-score.

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
Z <- embed_proteins(seqs, model = "tmvec-swissmodel-large")

# Beta diversity between samples
D <- sample_repdist(counts, Z)
vegan::adonis2(D ~ condition, data = meta)

# Structure-level protein clusters, e.g. for biomarker discovery
binned <- repdist_bin(counts, Z, cluster_threshold = 0.5, seqs = seqs)
```

## License

MIT
