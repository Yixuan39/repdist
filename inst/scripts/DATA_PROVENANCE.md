# Example data provenance

The two RDS files in `extdata` are fixed, rounded analysis examples. Their
numeric payloads are unchanged from the original runs; metadata and the
simulation sequences were added. Original embedding runs did not record immutable backbone/head
revisions. Those missing historical facts cannot be reconstructed from the
coordinates. New inference pins revisions in `extdata/models.json`.

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

## FunFam simulation universe

Source: CATH v4.3 FunFam Stockholm alignments, accessed 2026-08-23.
CATH database content is CC BY 4.0 (https://www.cathdb.info/wiki/doku/?id=start).
Attribution: CATH, Orengo and colleagues, University College London.

`funfam_universe.rds` retains 150 domains from each of these families:

| Superfamily | FunFam | Label |
| --- | --- | --- |
| 3.40.50.720 | 1 | GAPDH |
| 3.40.50.300 | 1 | FtsH protease |
| 1.10.10.10 | 1 | LysR regulator |
| 3.30.420.10 | 2 | RuvC resolvase |
| 2.40.50.140 | 1 | Ribosomal S12 |
| 3.20.20.70 | 1 | PdxS synthase |
| 1.20.1250.20 | 1 | MFS transporter |
| 3.40.640.10 | 1 | SHMT |
| 3.50.50.60 | 1 | Dihydrolipoyl DH |
| 3.40.190.10 | 2 | ABC binding protein |

Public alignment URL template:
`https://www.cathdb.info/version/v4_3_0/superfamily/{superfamily}/funfam/{number}/files/stockholm`

To recover the exact selected set, download each alignment, concatenate
Stockholm sequence fragments by record identifier, remove alignment gaps (`-`
and `.`), and select the identifiers in `members$protein` for that family.
The shipped `sequences` supplies the exact ungapped strings and order as a
reconstruction manifest. Historical selection metadata (seed 20260822, 150 per
family and length filter 60--510) is retained in `provenance`; reproducing the
fixed set uses its recorded identifiers rather than depending on sampling order
in a potentially rearranged server response.

The historical embeddings used ProtT5 XL UniRef50 and the `scikit-bio/tmvec`
head, not the current default Figshare head or `scikit-bio/tmvec-swissmodel`.
`prepare_simulation_data.R` can recompute the exact sequence set with the same
head identity pinned to commit 405ecd7651b378cdf3e3cb544d5477c73ac384bd and the
registry's pinned ProtT5 backbone. This revision was verified against
the current registry; it is not asserted to be the unrecorded historical revision.
New results record the complete model specification and R session.

Both datasets store embeddings rounded to four decimals. Renormalize rows
before cosine distances; rounding can move pairs across a sufficiently close
clustering threshold. The RDS data help pages document dimensions, units and
identifier correspondence. These preparation scripts are deliberately outside
vignette execution because inference requires multi-GB model downloads.

