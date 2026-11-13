# Simulation 4 plan: one function across remote homologs (2026-11-13)

Replaces the ko-rn design (several KOs of one KEGG reaction, dropped by
Yixuan on 2026-11-13). The file names keep `kofam`; the label is now a COG.

## Question

Do `repdist` bins recover one function whose proteins share < 30% identity,
and does that help DA and beta diversity, against proteins, MMseqs2 30%
clusters and the label itself (oracle)?

## Label source

| Option | Licence | Pool for one specific COG | Used |
| --- | --- | --- | --- |
| KO (KEGG `link/ko` on Swiss-Prot genes) | no stated redistribution terms | 100-400 | no |
| COG via UniProt `xref_eggnog` (Swiss-Prot) | CC BY 4.0 | 100-500 | no |
| NCBI COG 2024 release | NCBI | 900-2,200, class-0 (full-length) flag | **yes** |

## Headroom check (step 1)

Top labels by MMseqs2 30% clusters, <= 150 proteins each, package defaults,
binned per label (recall only; purity needs the background, below).

| Set | same 30% cluster | same bin |
| --- | --- | --- |
| Swiss-Prot COG top 10 | 0.01-0.19 | 0.58-1.00 |
| Swiss-Prot KO top 10 | 0.03-0.27 | 0.91-1.00 |
| COG 2024 top 10 (superfamilies) | 0.01-0.07 | 0.89-1.00 |
| COG 2024, COGs matching the top KOs | 0.04-0.31 | 0.71-1.00 |

E <= 0.01 alignments alone join 0.39-1.00 of same-label pairs directly; bins
follow the gated graph, so the recall gain comes from the alignment gate.

## Universe (`data-raw/simulation_kofam.R`)

| Item | Full | Smoke |
| --- | --- | --- |
| Pool | exactly one COG, all rows class 0, 50-1000 aa, standard residues, unique sequence | same |
| Treatment | `KOFAM_COG`, seeded draw of 300 | 20 |
| Background | random COGs (>= 60 members), 20 each, whole COGs to the cap | 5 each |
| Leakage | drop background with an MMseqs2 hit to treatment (E <= 1e-3, 80% cov both) | same |
| Cap | 3,000 proteins | 60 |
| Embeddings | TM-Vec, length-sorted | random unit vectors |
| Graph | `protein_similarity(x, n_cores = 1)` at build | same |
| Seed | `20261112L` | same |

The build stops if the treatment spans fewer than `2 x max(n)` 30% clusters.

## Design (`simulation_kofam.Rmd`)

| Item | Value |
| --- | --- |
| Count model | NB, D = 2e4, s ~ LN(0, 0.25), log-normal baseline, phi = (0.1 + 1/mu) LN(0, 0.5) |
| Enrichment | each A sample: rho over 2 treatment proteins from 30% clusters no other sample uses; B: rho = 0 |
| Grid | rho {0, .01, .02, .05} x n {10, 20} x 50 replicates (smoke: n = 4, 2 replicates) |
| DA units | protein, MMseqs2 30%, bin (`bin_glom`), COG (oracle); DESeq2 poscounts Wald, padj <= 0.05 |
| True positive | unit holds an enriched protein and > 50% of its abundance is treatment |
| Beta | PERMANOVA (999) on MMD, Bray-Curtis (protein), Bray-Curtis (MMseqs2 30%) |
| Clustering | pairwise recall/purity of bins and MMseqs2 30% vs COG (treatment; all COGs) |
| Regime | share of same-COG pairs with no >= 30% identity / 80% coverage hit |
