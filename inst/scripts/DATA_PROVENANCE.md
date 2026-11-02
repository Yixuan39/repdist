# Example data provenance

## metaproteomics.npz and metaproteomics.fasta

100 microbial proteins from the mouse gut metaproteomes of Blakeley-Ruiz et
al. (2025), *Dietary protein source alters gut microbiota composition and
function*, ISME Journal 19, wraf048, https://doi.org/10.1093/ismejo/wraf048.
The article and its supplementary material are distributed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/); raw proteomics accession
PXD041586.

`metaproteomics.npz` is a database in the layout `tmvec build-db` writes
(tmvec 1.0.5): protein IDs in `headers`, one float32 TM-Vec embedding per row
in `embeddings`. It stores no sequences; `metaproteomics.fasta` holds them,
which is what `read_embeddings()` takes as `seqs`.

Both files are cut from the curated object of the package website's case
study, `curated/diet_series.rds` on the repository's `gh-pages` branch, built
by `data-raw/diet_series.R` there from the paper's supplementary Data Set 1
(PSM counts and sequences) and Data Set 6 (annotation), with embeddings from
`embed_proteins()` at the default `"tmvec"` model. The 100 proteins are the 12
GH92-annotated proteins seen on the 20% yeast and 20% casein diets, then the 88
most abundant other proteins on those diets.
`inst/scripts/make_metaproteomics_example.R` writes both files.

## data/simulate_ec.rda

The vignette's dataset (`data(simulate_ec)`): 20 simulated samples over 1,160
bacterial and archaeal enzymes from UniProtKB/Swiss-Prot
(https://www.uniprot.org), distributed under CC BY 4.0
(https://creativecommons.org/licenses/by/4.0/).

The proteins, their TM-Vec embeddings and MMseqs2 clusters come from the
package website's EC simulation universe, `curated/simulation_universe.rds`
and `curated/simulation_mmseqs2.rds` on the repository's `gh-pages` branch,
built by `data-raw/simulation_universe.R` and `data-raw/simulation_mmseqs2.R`
there: reviewed Bacteria and Archaea entries whose single complete EC number
has experimental evidence, 50-1000 residues, one protein per UniRef50 cluster,
embedded with `embed_proteins()` at the default `"tmvec"` model.
`inst/scripts/make_simulate_ec.R` takes contrast 1 of the website's EC
simulation (EC 2.7.7 vs EC 2.3.2, 5% of abundance on two enriched enzymes per
sample), keeps every enzyme of the two classes plus 1,000 seeded background
enzymes, and simulates negative-binomial counts on that universe. Embeddings
are rounded to 4 decimals. Seeds, design, source checksum and commands are in
`simulate_ec$provenance`. The website's simulated data page describes each
step.
