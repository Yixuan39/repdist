# repdist (devel)

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

# Impute missing labels from protein distances
# labels is a named character vector: c(protein_id = "known label", ...)
predictions <- impute_labels(G, labels)  # 1NN + spline logistic score
rf_predictions <- impute_labels(G, labels, scoreMethod = "rf")
attr(predictions, "models")             # KNN references and fitted scorer
attr(predictions, "diagnostics")        # calibration size and support ranges
```

`impute_labels()` returns a label for every unannotated protein, together with
`supporting_protein`, `confidence`, supporting similarity and KNN vote counts.
Names align labels to proteins; omitted IDs, `NA` and blanks are unannotated.
KNN uses a fixed `k = 1` by default; the scoring method cannot change its labels.
RF scoring requires the optional `ranger` package. Both scorers use held-out
predictions with entire TM95 single-linkage groups excluded from references.
Supply named `groups` for another grouping protocol. Final label propagation
uses all original annotations; predictions never become training references.

`confidence` estimates agreement with original annotations. It is `NA`, with
a warning, when there are too few independent groups or correct/incorrect
calibration outcomes. `score_extrapolation` flags predictions outside the
calibration feature ranges; it is not a complete overlap test. Full coverage
can include incorrect assignments, and biological correctness still requires
independent validation. The `seed` argument preserves R's random-number state.

Use `biasAdjust = TRUE` to fit the scorer with capped selection weights.
This adjusts observed feature differences between annotated and unannotated
calibration records; it cannot resolve hidden annotation bias. Changing the
target catalog may change these weights and scores.

Run `Rscript inst/scripts/annotation_imputation.R` for the nested case study2
comparison of spline/RF scoring with and without selection adjustment;
see [the performance report](inst/scripts/annotation_imputation.md).

`bin_proteins()` cuts a complete-linkage tree, so every pair inside a bin meets
`min_sim`. Structural similarity alone does not guarantee shared function, so
compare annotation agreement with coverage and fragmentation before
transferring labels.

## License

MIT
