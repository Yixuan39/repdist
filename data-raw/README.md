# Curated inputs for the site

The analysis pages read four curated objects in `curated/`. Each is rebuilt
by one script here, run from the branch root with the current `repdist`
installed. The three source scripts compute protein embeddings from the
sequences with `embed_proteins(model = "tmvec")` (see `known_models()`); the
MMseqs2 script clusters the simulation universe's sequences. The pages then
compute every result from these objects at render time. Nothing else is
cached.

| Object | Script | Used by |
| --- | --- | --- |
| `curated/simulation_universe.rds` | `simulation_universe.R` | `reference.Rmd`, `simulation_ec.Rmd`, `parameters.Rmd` |
| `curated/simulation_mmseqs2.rds` | `simulation_mmseqs2.R` | `simulation_ec.Rmd` |
| `curated/simulation_funfam.rds` | `simulation_funfam.R` | `simulation_funfam.Rmd` |
| `curated/diet_series.rds` | `diet_series.R` | `case_study.Rmd` |

```sh
Rscript data-raw/simulation_universe.R
Rscript data-raw/simulation_mmseqs2.R
Rscript data-raw/simulation_funfam.R
Rscript data-raw/diet_series.R
Rscript -e 'rmarkdown::render_site()'
```

Embeddings are stored rounded to six decimal places (cosine error below
1e-5). Each object carries a `provenance` list with the source, the filtering
counts, the model specification and the R session that built it.

## `simulation_universe.rds`

Source: UniProtKB/Swiss-Prot through the REST API
(<https://rest.uniprot.org>), CC BY 4.0. The release is recorded in
`provenance$uniprot_release`.

1. Query `(reviewed:true) AND (fragment:false) AND (taxonomy_id:2 OR
   taxonomy_id:2157) AND (cc_catalytic_activity_exp:*)`.
2. Keep entries whose EC numbers (all catalytic-activity reactions and the
   protein name, including multi-domain parts) are one single EC number,
   complete to four levels, carried by at least one reaction with
   experimental evidence (ECO:0000269).
3. Keep 50-1000 residues.
4. Map to UniRef50 (UniProt ID mapping) and keep the lowest accession per
   cluster.
5. Embed with `tmvec`.

Fields: `proteins` (AAStringSet with `mcols(proteins)$embedding`),
`annotation` (accession, EC4, EC3, name, organism, taxon id, superkingdom,
length, UniRef50 cluster), `provenance`.

## `simulation_mmseqs2.rds`

Source: the sequences of `simulation_universe.rds`, with MMseqs2 (version
recorded in `provenance$mmseqs_version`) on the `PATH`. Rerun it whenever
`simulation_universe.rds` changes; `provenance$universe_md5` records the file
it was computed from, and `simulation_ec.Rmd` stops if they differ.

1. `mmseqs easy-cluster` at 50% identity (primary) and at 30% (sensitivity),
   both with 80% coverage of query and target (`-c 0.8 --cov-mode 0`),
   greedy set cover (`--cluster-mode 0`) and `-s 7.5`.
2. `mmseqs easy-search` all against all (`-s 7.5 -e 1e-3 --max-seqs
   10000`), keeping the best hit per unordered pair.

Fields: `clusters` (accession, `mmseqs50` and `mmseqs30` cluster
representatives), `hits` (p1, p2, fident, alnlen, qcov, tcov, evalue),
`provenance` (the exact commands, cluster counts, universe MD5, session).

## `simulation_funfam.rds`

Source: CATH v4.3 [FunFam](https://www.cathdb.info) seed alignments (CC BY
4.0), through `https://www.cathdb.info/version/v4_3_0`; downloads are kept in
`data/cath_v4_3/`. Needs MMseqs2 on the `PATH`.

1. Take the 300 superfamilies with the most S35 representatives
   (`cath-superfamily-list-v4_3_0.txt`).
2. In each, take the FunFam with the most members among those whose seed
   alignment has DOPS >= 70 (`api/rest/superfamily/<sf>/funfam`).
3. Read its Stockholm seed alignment; an ungapped row is the domain
   sequence. Keep Bacteria and Archaea, 50-1000 residues, the 20 standard
   residues, unique sequences.
4. Cluster all domains with `mmseqs easy-cluster --min-seq-id 0.5 -c 0.8
   --cov-mode 0 --cluster-mode 0 -s 7.5`; drop clusters spanning two
   FunFams and keep the lowest domain id per cluster, then one domain per
   UniProt accession.
5. Keep at most 100 domains per FunFam (seeded draw); embed with `tmvec`.

Fields: `proteins` (AAStringSet with embeddings), `annotation` (domain,
family, UniProt accession, kingdom, organism, length, 50% cluster, FunFam
name), `families` (superfamily, FunFam number, name, members, DOPS, domains
before and after the cap), `provenance`.

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
