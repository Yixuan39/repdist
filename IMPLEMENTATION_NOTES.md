# KOfam simulation: implementation notes

Branch `sim-kofam`, worktree
`/Users/yixuanyang/Documents/Projects/repdist-worktrees/sim-kofam`, based on
`origin/gh-pages` (`e2f4059`). Follows `KOFAM_SIM_PLAN.md`. Official KEGG
and UniProt URL/licence research is in `kofam-official-sources.md`.

The dirty `gh-pages` checkout at
`/Users/yixuanyang/Documents/Projects/repdist` and the `review-ben-2026-11`
worktree were not modified. Nothing was pushed, and no PR was opened.

## Files on this branch

- `data-raw/simulation_kofam.R` (new): builds `curated/simulation_kofam.rds`.
- `simulation_kofam.Rmd` (new): four method variants and the planned
  metrics. Knit-time recomputation from the RDS.
- `_site.yml`, `index.Rmd`, `simulation.Rmd`: FunFam and EC removed from
  navigation and hub links only; those files and their curated data are
  still in the worktree.
- `data-raw/README.md`, `.gitignore` (`/data/kegg/` is ignored; already
  covered by `/data`), `references.bib` (`kanehisa2000kegg`).
- `KOFAM_SIM_PLAN.md`, `kofam-official-sources.md`, this file.

Not committed: `curated/simulation_kofam.rds` (a `KOFAM_SMOKE=1` object
with random embeddings and KEGG-derived labels), `simulation_kofam.html`
(a smoke knit, not a published result), `data/kegg/` (raw KEGG/UniProt
cache), and `CLAUDE_PROMPT.md`.

## How grouping works (and one mapping change)

Primary grouping is still KEGG ko-rn: reactions with two or more KOs are
candidates; BRITE `ko00001` is only used to prefer Metabolism-annotated
candidates and to annotate; `ko01000` is never fetched or used as a
grouping key. No treatment KO or reaction ID is hard-coded.

The documented endpoint `https://rest.kegg.jp/link/rn/ko` returned HTTP
400 when checked live. `https://rest.kegg.jp/link/reaction/ko` returns the
same two-column table and is what the build uses; the script tries the
`rn` alias first and records the working path in
`provenance$ko_rn_endpoint` (smoke run: `link/reaction/ko`).

Mapping route is Swiss-Prot first, then KEGG (plan §2.1), not KO-first.
An earlier KO-first draft (`link/genes/ko` then `conv/uniprot`) mostly
hit unreviewed TrEMBL genes and found no qualifying proteins. The script
now downloads reviewed Bacteria/Archaea Swiss-Prot entries with a KEGG
cross-reference and maps their genes with `link/ko/<genes>`. Grouping is
still applied to the ko-rn table; only the protein lookup changed.

## Defaults interpreted (plan §11)

Adopted as proposed unless noted:

- length 50–1000, ≥ 3 proteins/KO, ≤ 20 proteins/KO (10 in smoke),
  universe cap 3000 (60 in smoke), seed `20261112L`.
- MMseqs2 settings copied from `data-raw/simulation_mmseqs2.R`
  (`-s 7.5`, 80% coverage, 50%/30% clusters, E ≤ 1e-3 hits). Cross-KO
  hits or a shared 30% cluster exclude a treatment group; leaking
  background proteins are dropped.
- Homology-only weights default to 1. Bit-score sensitivity was not
  implemented; alignments are not cached in the RDS.
- `k_target = 1`: only the first surviving candidate is accepted;
  every probed reaction is recorded in `candidates`.
- Count model ports Simulations 1–2: depth `D = 2e4`,
  `phi = (0.1 + 1/mu) * lognormal(0, 0.5)`. Grid
  `rho in {0, 0.01, 0.02, 0.05}`, `n` per group 10 and 20 (smoke: `n = 4`,
  `n_rep = 2`). Baseline composition is uniform. Condition B has
  `rho = 0`.
- Step D (subunit filter) is a case-insensitive text search for
  "subunit" in the official `get/rn:` and `get/ko:` entries. Cheap and
  recorded, not a semantic review.
- DA on the page follows the case study (DESeq2, `poscounts`, Wald,
  adjusted p ≤ 0.05) rather than the FunFam/EC limma/DEqMS helper. A
  unit is a true positive only if a majority of its abundance comes
  from treatment-reaction proteins.

## Smoke build (what was actually run)

`KOFAM_SMOKE=1 Rscript data-raw/simulation_kofam.R` completed and wrote
`curated/simulation_kofam.rds` (kept local, not committed).

Recorded in that object, not hard-coded in the scripts:

- 15,000 Swiss-Prot entries were fetched (smoke cap), 8,787 unique
  standard-residue sequences with a single KO.
- 12,416 well-formed ko-rn rows; 2,206 candidate reactions with ≥ 2 KOs.
- After the subunit / protein-count / cross-KO homology filters, the
  first accepted treatment was reaction `R02251` (KOs `K00635` and
  `K18851`, 16 proteins). 47 earlier candidates in the ranked list were
  rejected (subunit language, too few qualifying proteins in the smoke
  sample, or cross-KO homology). This is a smoke outcome, not a chosen
  published contrast.
- Final universe: 60 proteins, 46 KOs (16 treatment, 44 background,
  one protein each after the cap). Leakage removed: 0.
- Embeddings: random unit vectors (512-d), clearly flagged in
  `provenance$embedding`. Seed `20261112`. Accessed 2026-11-13
  (America/New_York).

## Knit

A first knit failed in two ways, both in `simulation_kofam.Rmd`, after
the build had already succeeded:

1. DESeq2's parametric and `fitType = "mean"` dispersion fits both fail
   on this small smoke universe. The page now falls back to DESeq2's
   documented gene-wise-as-final-dispersions path, then `nbinomWaldTest`.
2. `aricode::ARI` / `NMI` coerce labels with `as.numeric()`, so character
   labels like `"bin1"` or `"treatment"` became NA. Labels are now mapped
   to integer codes first.

After those fixes, `rmarkdown::render("simulation_kofam.Rmd")` succeeded
locally (about 30 seconds, including `protein_similarity()` on 60
proteins). The HTML is a smoke render and was not committed. The page
banners that every number on it is computed but is not a result.

## Full build: not started

A non-smoke run was not started. Without `KOFAM_SMOKE=1`, the script
issues one UniProt probe and stops. The query (reviewed, non-fragment,
Bacteria/Archaea, length 50–1000, any KEGG cross-reference) matched
334,791 entries (live header `x-total-results` on 2026-11-12/13, of
355,745 reviewed Bacteria/Archaea). That is a large sequence download,
over the plan's 30-minute / large-download rule. To proceed anyway,
remove or edit that `stop()` deliberately, then:

```sh
Rscript data-raw/simulation_kofam.R
```

Quadratic `pwalign` in `protein_similarity()` (knit-time, not in the
build script) is then bounded by the 3000-protein cap.

## Review residual issues (not all fixed)

- No separate label-permutation Type I experiment; Type I is estimated
  from the `rho = 0` grid only. MMD uses 999 permutations inside
  `adonis2`, which is the test, not a pipeline-level permutation.
- Pairwise E-values / bit scores are not stored, so a bit-score
  homology-only sensitivity cannot be added without re-alignment.
- The "subunit" text heuristic can be both too strict and too loose.
- Ranking still prefers reactions with the most KOs, so a smoke sample
  of 15,000 proteins had to probe many candidates before one passed.
  A full Swiss-Prot table may change which reaction is accepted.
- `protein_similarity()` has no `db_size` argument; this simulation
  uses one fixed universe, so all variants see the same database size.
- The distances / MCL-parameters pages still describe the old FunFam
  universe; they remain in the menu, as allowed by the plan.
- KEGG-derived labels in a public curated RDS / page still need
  Yixuan's written scope confirmation (plan §11 item 1).

## Claude session

Implementation was done with `claude -p` in this worktree. That session
hit a usage limit mid-knit (reset 3:10am America/New_York) and exited
without a commit or a successful knit. The DESeq2 and ARI/NMI fixes
above, the successful smoke knit, this notes update, and the local
commit were finished in review. Codex completed the official-source
research; the plan document was assembled from that research after two
`codex exec` write sessions hung mid-reasoning.
