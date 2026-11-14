# Curated inputs for the site

The case study reads one curated object in `curated/`, rebuilt by the script
here, run from the branch root with the current `repdist` installed; it
computes protein embeddings with `embed_proteins(model = "tmvec")` (see
`known_models()`). The page computes every result from it at render time.
The vignette's and reference's dataset, `simulate_cog`, is built by
`data_simulation.Rmd` at the branch root instead.

| Object | Script | Used by |
| --- | --- | --- |
| `curated/diet_series.rds` | `diet_series.R` | `case_study.Rmd` |

```sh
Rscript data-raw/diet_series.R
```

Embeddings are stored rounded to six decimal places (cosine error below
1e-5). Each object carries a `provenance` list with the source, the filtering
counts, the model specification and the R session that built it.

## `diet_series.rds`

Source: Blakeley-Ruiz JA *et al.* (2025) Dietary protein source alters gut
microbiota composition and function. *ISME J* 19(1): wraf048.
doi:10.1093/ismejo/wraf048. Supplementary Data Set 1 (Extended Data Table 1:
PSM counts and sequences) and Data Set 6 (Extended Data Table 6: annotation),
CC BY 4.0, downloaded from Europe PMC (PMC12066410) into `data/wraf048/` and
checked against their MD5 sums. Raw spectra: PRIDE PXD041586.

1. Keep `Microbiome` proteins of Data Set 1; blank cells are zero counts.
2. Parse cage, mouse and diet from the sample names and keep every sample
   (12 mice, ten diet periods, 111 samples). Diet is the token of the sample
   name, ordered by feeding period. Sex is fixed by cage (C1, C3 male; C2, C4
   female).
3. Keep proteins with at least one PSM in any sample and at most 1500
   residues (`MAX_AA`).
4. Join Data Set 6 on `Protein Identifier`, falling back to `Protein
   Identifier with Species Code`, keeping a match only where the residue
   counts agree. The species code is the first token of the protein
   identifier; its GTDB lineage comes from Data Set 6 (`BT` is
   *B. thetaiotaomicron*). KO is `MantisKO`, falling back to `MicrobeKO`.
5. Embed with `tmvec`.

Fields: `proteins` (AAStringSet with embeddings), `counts` (raw PSM counts,
samples × proteins), `samples` (sample, subject, cage, diet, period, sex),
`annotation` (species code, taxon, lineage, PROKKA description, curated broad
and detailed function, consensus annotation, KOs, CAZy), `provenance`.

The case study selects its Yeast–Casein pairs from this object; another diet
pair needs no rebuild.
