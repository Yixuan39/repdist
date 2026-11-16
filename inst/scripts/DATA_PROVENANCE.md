# Example data provenance

## metaproteomics.fasta, metaproteomics_counts.tsv, metaproteomics_samples.tsv

100 microbial proteins from the mouse gut metaproteomes of Blakeley-Ruiz et
al. (2025), *Dietary protein source alters gut microbiota composition and
function*, ISME Journal 19, wraf048, https://doi.org/10.1093/ismejo/wraf048.
The article and its supplementary material are distributed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/); raw proteomics accession
PXD041586.

`metaproteomics.fasta` holds their sequences. It is cut from the curated
object of the package website's case study, `curated/diet_series.rds` on the
repository's `gh-pages` branch, built by `data-raw/diet_series.R` there from
the paper's supplementary Data Set 1 (PSM counts and sequences) and Data
Set 6 (annotation). The 100 proteins are the 12
GH92-annotated proteins seen on the 20% yeast and 20% casein diets, then the 88
most abundant other proteins on those diets.
`metaproteomics_counts.tsv` holds their PSM counts in the 23 samples of
those two diets (samples in rows, sample IDs in the first column), and
`metaproteomics_samples.tsv` the samples' metadata, one row per sample.
`inst/scripts/make_metaproteomics_example.R` writes all three files.

## data/simulate_cog.rda and simulate_cog.fasta

The vignette's dataset (`data(simulate_cog)`): 20 simulated samples over 110
bacterial and archaeal proteins of the NCBI COG 2024 release
(https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/; Galperin et al. 2025,
https://doi.org/10.1093/nar/gkae983). NCBI places no restrictions on the use
or distribution of its data (https://www.ncbi.nlm.nih.gov/home/about/policies/).

The proteins, labels and counts come from `data_simulation.Rmd` on the
repository's `gh-pages` branch: proteins with exactly one COG, assigned over
their full length (membership class 0), 50-1000 residues; 10 GrpE (COG0576)
proteins from 10 different MMseqs2 30%-identity clusters, none sharing a 30%
alignment over 70% of both proteins; 100 random background proteins from
other COGs with no MMseqs2 hit (E <= 1e-3, 80% coverage) to a GrpE protein.
Each group-A sample puts 5% of its abundance on its own GrpE protein, with
negative-binomial counts. `inst/scripts/make_simulate_cog.R` reads those
files, adds TM-Vec embeddings from `embed_proteins()` (rounded to 4
decimals), and writes `data/simulate_cog.rda` and
`inst/extdata/simulate_cog.fasta` (the sequences alone, for the vignette's
embedding step). Design, seeds and session are in `simulate_cog$provenance`.
