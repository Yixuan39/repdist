# Simulation 4: implementation notes (2026-11-13)

Branch `claude/affectionate-goldwasser-2c8db4`, cut from `sim-kofam`
(`732ebd2`, itself on `origin/gh-pages` `e2f4059`). Design: `KOFAM_SIM_PLAN.md`.
Nothing pushed; no PR.

## Changed

| File | Change |
| --- | --- |
| `data-raw/simulation_kofam.R` | rewritten: one COG 2024 label; ko-rn, BRITE, "subunit" and cross-KO filters deleted |
| `simulation_kofam.Rmd` | rewritten: bins vs labels, DA at 4 unit levels, MMD vs Bray-Curtis |
| `_site.yml`, `index.Rmd`, `simulation.Rmd` | menu and hub text |
| `references.bib` | `galperin2025cog` |
| `KOFAM_SIM_PLAN.md`, `data-raw/README.md`, this file | rewritten |

Not committed: `curated/simulation_kofam.rds` (smoke), `simulation_kofam.html`,
`data/` (COG, UniProt, KEGG caches). `kofam-official-sources.md` (KEGG licence
research) is kept unchanged for the record.

## Open: treatment COG

`KOFAM_COG` has no default until Yixuan picks the label. The smoke build used
COG0576 (GrpE) as a stand-in. n = 20 with 2 enriched proteins per sample
needs >= 40 30% clusters among the 300 drawn treatment proteins:

| COG | Function | 30% clusters in a 300 draw | mean aa |
| --- | --- | --- | --- |
| COG1589 | FtsQ | 90 | 270 |
| COG1826 | TatA/B/E (three subunits) | 63 | 90 |
| COG1381 | RecO | 41 | 249 |
| COG0576 | GrpE | 31 | 196 |
| COG0711 | ATP synthase b | 31 | 169 |
| COG1758, COG0712, COG0037, COG0360, COG0319 | omega, ATP synthase delta, TilS/TtcA, S6, YbeY | 16-24 | 74-304 |

Superfamily COGs (TetR, MFS, LysR...) have 700-2,400 clusters but are broad
labels.

## Smoke (run)

`KOFAM_SMOKE=1 KOFAM_COG=COG0576 Rscript data-raw/simulation_kofam.R`: 55 s.

| Step | n |
| --- | --- |
| COG 2024 assignment rows / proteins | 6,266,336 / 5,872,258 |
| one COG, all class 0 / 50-1000 aa | 3,936,892 / 3,906,828 |
| GrpE qualifying | 1,936 |
| universe (20 treatment in 9 30% clusters + 8 COGs x 5) | 60 |
| leakage removed | 0 |

`rmarkdown::render("simulation_kofam.Rmd")`: 21 s. Random embeddings, n = 4:
numbers are not results (protein-level Type I is inflated by DESeq2's
gene-wise fallback on 60 proteins).

## Full build: not run (estimate)

| Stage | Measured | Full estimate |
| --- | --- | --- |
| COG 2024 download (2.2 GB, cached) | ~2 min | 0 (cached) |
| pool + sequence extraction + MMseqs2 | 55 s (smoke) | ~2 min |
| TM-Vec embedding | 2,533 proteins 7.6 min (MPS) | ~9 min |
| `protein_similarity(n_cores = 1)` | 2,956 pairs/s (GrpE, 196 aa) | 4.5M pairs: 25 min at 196 aa, ~60 min at 300 aa |
| knit: 400 replicates | DESeq2 3000 x 40 3.9 s; `sample_mmd` 3.7 s; `adonis2` 0.3 s | ~1.5 h |

```sh
KOFAM_COG=COGxxxx Rscript data-raw/simulation_kofam.R   # ~40-70 min
Rscript -e 'rmarkdown::render("simulation_kofam.Rmd")'  # ~1.5 h
```

## Step-1 headroom check

Scripts and the full table are in the side-session notes folder
(`headroom-step1.md`, `headroom-table.md`); per-COG 30% counts for COG 2024
were run on the BRC cluster (`~/project_data/cog24`) plus locally.
