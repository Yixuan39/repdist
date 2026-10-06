# Example data provenance

The RDS file in `extdata` is a fixed, rounded analysis example. Its numeric
payload is unchanged from the original run; metadata were added. The original
embedding run did not record immutable backbone/head revisions. Those missing
historical facts cannot be reconstructed from the coordinates. New inference
pins revisions in `extdata/models.json`.

## Diet metaproteome

Source: Blakeley-Ruiz et al. (2025), *Dietary protein source alters gut
microbiota composition and function*, ISME Journal 19, wraf048.
https://doi.org/10.1093/ismejo/wraf048

Use the paper's supplementary Data Set 1 (complete metaproteome, counts and
sequences) and Data Set 6 (annotated microbial metaproteomes). The article and
supplementary material are distributed under CC BY 4.0; attribution is retained
here and in the data help page. License: https://creativecommons.org/licenses/by/4.0/
Raw proteomics accession: PXD041586. Raw mass spectra need not be reprocessed to
reconstruct this example from the published tables.

Download/extract the tab-separated tables and place them under these names:

* `data/test_study/ExtendedDataTable1CompleteMetaproteome.txt`
* `data/test_study/annotated/ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt`

From the source repository root, run `Rscript inst/scripts/curate_diet_data.R`
(requires phyloseq and Biostrings). This reads the header after two title lines,
keeps `Microbiome` rows, converts blank abundance cells to zero, and rejects
other non-numeric counts. Table 1 supplies PSM spectral counts and sequences;
Table 6 supplies taxonomy and consensus annotation, **not** its %NSAF values.
Protein identifiers join the tables. Missing taxonomy stays missing. Sample
names encode cage, mouse and dietary treatment. The script explicitly records
feeding order, dose, sex by cage and the terminal control diet; `age_week` is
elapsed weeks since baseline, not the animal's absolute age.

Then run `Rscript inst/scripts/prepare_vignette_data.R` (also requires repdist).
Rarefy the complete catalog to its smallest sample depth with seed 1, without
replacement; remove baseline samples and zero-abundance proteins; retain the
2000 proteins with greatest summed rarefied counts. Subsetting changes sample
totals, so the final 99-by-2000 table is not equal-depth. Retain the corresponding
unrarefied counts for differential abundance. Embed whole sequences using the
current pinned `tmvec` preset, preserve row order, round to four decimal places,
and compress with xz. Cache identity includes sequences, model and package
version. This is a new reproducible run and need not reproduce the historical
coordinates bit for bit. The shipped historical data used repdist 0.99.0 and
the default Figshare head (file 49181530; SHA-256 recorded in the RDS).

The dataset stores embeddings rounded to four decimals. Renormalize rows
before cosine distances; rounding can move pairs across a sufficiently close
clustering threshold. The data help page documents dimensions, units and
identifier correspondence. These preparation scripts are deliberately outside
vignette execution because inference requires multi-GB model downloads.

