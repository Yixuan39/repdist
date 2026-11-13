# KOfam simulation plan (Korect / repdist)

> 中文摘要：用一个基于 KEGG Orthology 的模拟，替换网站上现有的 FunFam 和 EC 模拟。treatment 组来自 **KEGG KO→REACTION（ko-rn）映射**：同一个反应编号对应多个 KO，这些 KO 被视为可互相替代的同反应酶，不要求结构相近。BRITE（建议用 ko00001）只做辅助过滤和功能类别标注；**不要用 ko01000（Enzymes）做分组**，因为它按 EC 组织。蛋白质序列和元数据优先用 UniProtKB/Swiss-Prot（CC BY 4.0）；KO 标签经 Swiss-Prot 的 `xref_kegg` 再查 KEGG REST 得到。所有 KEGG 映射必须在构建时从官方 REST/文件查得，不得硬编码 KO 或反应编号。对照（背景 KO）的效应量在设计上就是 0，并另做全局零效应和标签置换以估计 Type I 错误。丰度用带过离散的负二项，扫描样本量和效应量。方法对比为：两者结合（`max_evalue=0.01` + 余弦权重）、只用 embedding（`max_evalue=Inf`）、只用同源性（同一批 E≤0.01 的边，权重为 1，或把 bit score 缩放到 (0,1]）、以及按 KO 计数的朴素差异丰度基线。评估包括簇层面的 AUPRC / ARI / NMI，差异丰度的功效和经验 FDR，以及 MMD 的 Type I 错误和功效。E 值门槛本身这次不改。全对比对是平方级，蛋白质子集最多几千条；脚本必须带 DRY_RUN/smoke 模式。KEGG 不是开放许可，公开网站上能放哪些派生标签需要 Yixuan 确认。

This document is a build and evaluation plan, not a result. No KO identifiers, reaction identifiers, accessions, or wall-clock times are chosen from memory. Numbers that appear below are **proposed caps, seeds, or defaults to be confirmed**, and are labelled as such.

Codex researched official KEGG/UniProt URLs and licence texts on 2026-11-12; those findings are in `kofam-official-sources.md` beside this file. A later `codex exec` write session hung after that research, so this plan was assembled from those findings and the already-decided design. Implementation should follow this file, not invent a different grouping scheme.

---

## 1. Goal and relationship to the retired FunFam / EC simulations

The website currently presents three protein-level simulations (committed `origin/gh-pages` at `e2f4059`):

- FunFam (`simulation_funfam`): a positive control on CATH domain sequences. Retired because Korect’s unit is the **full protein**, while FunFam’s unit is a **domain**.
- EC (`simulation_ec`): sequence-disjoint enzyme-class contrasts on Swiss-Prot. Retired because a four-level EC number can group structurally different, convergent enzymes, which does not match the method’s goal, and the existing page was not a strong demonstration.
- Distances / MCL parameters pages that reuse those universes.

They are replaced by **one** KOfam / KO simulation:

- The functional unit is a **KEGG Orthology (KO) group of complete proteins**.
- A treatment set is several KOs that catalyse the **same reaction** and are mutually interchangeable alternatives. Structural similarity is **not** required.
- KOs that merely sit in the same module or pathway, but catalyse **different reaction steps**, are not a valid treatment group.

The new page should make this contrast explicit: the question is whether TM-Vec cosine graphs, with or without the homology E-value gate, can group interchangeable enzymes of one reaction and support differential abundance / MMD, relative to a naive “count by KO” baseline. It is not a claim that the method transcends sequence homology in general; a later E-value-threshold benchmark is out of scope here.

---

## 2. Data sources, URLs, licences, and provenance

### 2.1 UniProtKB/Swiss-Prot (sequences and protein metadata)

- Licence: [CC BY 4.0](https://www.uniprot.org/help/license) (also machine-readable at `https://rest.uniprot.org/help/license`). Attribute UniProt, link the licence, and record transformations.
- API style should follow the existing `data-raw/simulation_universe.R`: UniProt REST (`https://rest.uniprot.org`), paginated queries, recorded release.
- Official return fields (`https://rest.uniprot.org/configure/uniprotkb/result-fields`) include `xref_kegg` under genome annotation and **do not** currently list a direct `xref_ko` field. Official ID mapping (`https://rest.uniprot.org/configure/idmapping/fields`) maps to KEGG **gene** identifiers, not KO. Do not assume a Swiss-Prot `KO` column exists.
- Default mapping route (already-annotated KOs; no KOfam HMM redistribution):
  1. Download reviewed Swiss-Prot records with sequences and `xref_kegg`.
  2. Map KEGG gene IDs to KO with KEGG REST `https://rest.kegg.jp/link/ko/<KEGG-gene-identifiers>` (batched, ≤ 3 calls/s).
  3. Fallback where needed, as documented in the API manual: `https://rest.kegg.jp/conv/genes/uniprot:<accession>`.
- Record unmapped and conflicting mappings. Proposed default, to be confirmed: exclude proteins with no KO, and exclude proteins that map to **more than one KO** (ambiguous labels).
- Provenance fields (build-time): UniProt release, query string, retrieval timestamp, response hashes, accession and sequence versions, filter counts at each step, R/`repdist` session info.

### 2.2 KEGG (KO–reaction links and BRITE)

Confirmed official endpoints (validate live at build time; require HTTP 200 and a nonempty, well-formed payload; record URL, timestamp, hash, and `https://rest.kegg.jp/info/kegg`, `info/ko`, `info/reaction`):

| Role | Official endpoint |
| --- | --- |
| **Primary grouping**: KO → reaction (ko-rn) | `https://rest.kegg.jp/link/rn/ko` (documented as `/link/rn/ko` in the [KEGG announcement](https://www.kegg.jp/kegg/docs/announce.html) and the [API manual](https://www.kegg.jp/kegg/rest/keggapi.html); `reaction` alias `rn`; tab-delimited) |
| Auxiliary BRITE hierarchy (KO tree) | `https://rest.kegg.jp/get/br:ko00001/json`, or the official download `https://www.kegg.jp/kegg-bin/download_htext?filedir=&format=json&htext=ko00001` ([BRITE catalog](https://www.kegg.jp/kegg/brite.html), [live ko00001](https://www.kegg.jp/kegg-bin/get_htext?htext=ko00001)) |
| Explicitly **not** a grouping key | BRITE `ko01000` (Enzymes), which descends from EC classes to KO leaves: [official hierarchy](https://www.kegg.jp/kegg-bin/get_htext?filedir=&hier=4&htext=ko01000) |

Rate limit: academic REST access is limited to **at most 3 calls/second** ([KEGG API policy](https://www.kegg.jp/kegg/rest/)).

Licence (research findings, not legal advice):

- KEGG is copyright Kanehisa Laboratories and is **not** an open-licensed database ([copyright page](https://www.kegg.jp/kegg/legal.html)).
- `rest.kegg.jp` is limited to academic use by academic users. The API policy does not say that ordinary academic REST access requires a signed agreement.
- Academic users “providing services” are requested to obtain an academic service-provider licence, bundled with the academic FTP subscription ([academic subscription](https://www.pathway.jp/en/academic.html)). The current end-user option does not allow setting up a local copy. The service-provider option covers research groups developing software or web services that utilize KEGG.
- The current academic service-provider agreement (SP-202507, [PDF](https://www.pathway.jp/doc/agreement_sp.pdf)) grants specified research/publication/service rights, limits redistribution of FTP data, and allows user retrieval only in small quantities. Buying a subscription must not be assumed to grant unlimited public rehosting.
- **No explicit blanket permission** was found for public redistribution of REST-retrieved raw files or derived KO/reaction/BRITE mapping tables. Converting a table to RDS, or obtaining gene IDs via UniProt, does **not** make KEGG-derived relations open-licensed.

Conservative, licence-safe workflow (proposed default until Yixuan gets written scope confirmation):

1. Download KEGG mappings at **build time** into a local **untracked** cache (e.g. `data/kegg/`, listed in `.gitignore`). Never commit raw KEGG files or complete ko-rn / BRITE tables.
2. Preserve hashes, URLs, and version metadata in `provenance`.
3. A public curated RDS may contain UniProt-sourced sequences/metadata (CC BY 4.0), embeddings, and numerical analysis objects. Inclusion of KEGG-derived labels or reverse-decodable encodings on the **public** gh-pages site needs written confirmation (see §11).
4. To preserve knit-time recomputation, an authorized **local/private** curated RDS can render a public page that shows only aggregate plots and counts. Public rebuilds would re-download mappings under the visitor’s own academic use.
5. Do not redistribute KOfam HMM profiles; they are unnecessary because labels come from already-annotated Swiss-Prot ↔ KEGG mappings.

### 2.3 Software

- Installed `repdist` from local `origin/main` at `d423023` (`protein_similarity: E-value gate (pwalign, E<=0.01) with cosine weights, no TM threshold`). There is no remaining local branch named `evalue-links`; that algorithm is already on main. Signature: `protein_similarity(x, max_evalue = 0.01, n_cores = 1L)`. There is **no** `db_size` argument.
- Default embeddings: `embed_proteins(model = "tmvec")` (ProtT5-XL + `tmvec1-large` head), same as the existing site scripts.
- MMseqs2 on `PATH` (version recorded in provenance), following `data-raw/simulation_mmseqs2.R`.
- R packages already used on the site: `Biostrings`, `pwalign` (via the package), `MCL`, `limma` / `DEqMS` if the DA section follows the existing vignette, `vegan` for PERMANOVA if MMD distances are tested that way.

---

## 3. Reproducible KO-group selection algorithm

All steps run against downloaded data. No KO or reaction ID is hard-coded as the chosen treatment. Illustrative names, if mentioned on the page, must be marked as examples and must not drive selection.

**Step A — Load ko-rn.** Fetch `https://rest.kegg.jp/link/rn/ko`. Parse two-column KO / reaction identifiers. Drop malformed rows. Group KOs by reaction ID.

**Step B — Candidate treatment groups.** A reaction with **two or more** distinct KOs is a candidate group: those KOs are interchangeable alternatives for that reaction. Reactions with a single KO are not treatments (they may remain background). Same-module / same-pathway membership is ignored as a grouping key.

**Step C — BRITE is auxiliary only.** Attach `ko00001` paths to each KO for filtering and annotation. Preserve **all** paths when a KO appears in more than one place; do not silently collapse categories. Proposed default filters, labelled as defaults to confirm:

- Prefer candidate groups whose KOs are annotated under Metabolism (or another explicitly recorded top-level category) in `ko00001`, so the page can report functional context.
- **Never** group, rank, or exclude candidates by `ko01000` (EC-organised Enzymes). `ko01000` may be stored only as an optional annotation, if at all.

**Step D — Reject complexes and unclear cases.** A bare multi-KO reaction link can include **subunits of one complex** rather than independently interchangeable full-protein enzymes. At build time, inspect official KO / reaction entries for each surviving candidate (rate-limited REST `get`), and reject required multisubunit components and unclear cases with a recorded reason. This is a **filter after** reaction grouping, not a different grouping scheme.

**Step E — Restrict to the protein universe.** Keep a candidate group only if, after the Swiss-Prot filters of §5, **at least two KOs** of that reaction still have at least a proposed minimum number of proteins each (proposed default, to confirm: **≥ 3 proteins per KO**, otherwise the group is too thin for a contrast).

**Step F — Cross-KO homology and leakage filters** (details in §4). Drop a candidate if its member KOs are sequence-homologous to each other, or if treatment proteins leak into the control set.

**Step G — Tie-breaks and which group(s) to simulate.** After all filters, rank remaining candidate groups by a fully specified key and keep the top *k*. Proposed default, to confirm:

- Rank by (i) number of distinct KOs (descending), then (ii) total qualifying proteins (descending), then (iii) reaction ID (lexicographic, for stability).
- Proposed *k*: **one primary treatment group** for the representative figures, plus **all other groups that pass the caps** for a sensitivity/robustness table, *or* the top 5, whichever is smaller. The exact *k* is a cap, not a result.
- Record the full ranked table (reaction ID, KO list, counts, filter reasons for exclusions) in the private build log / RDS provenance, even for groups not simulated.

**Step H — Seeds.** All random draws (protein caps per KO, if any; simulation of counts; permutations) use a single recorded seed. Proposed default seed, to confirm: `20261112L` (date of this plan). Additional stream seeds must be derived from it, not drawn ad hoc.

---

## 4. Homology / leakage filters

Purpose: the scientific claim “these treatment KOs were not already grouped by sequence homology” is only valid if we **check** it, and control proteins must not be hidden homologs of the treatment.

Proposed defaults, to be confirmed (chosen to stay close to the existing `simulation_mmseqs2.R` settings, not because they have been tuned for KO):

1. Cluster or search the simulation universe with MMseqs2 `easy-search` all-vs-all, sensitivity `-s 7.5`, E-value `1e-3`, and 80% coverage of query and target (`-c 0.8 --cov-mode 0`), keeping the best hit per unordered pair. Also compute `easy-cluster` at 50% and 30% identity with the same coverage, as in the existing script, for reporting.
2. **Within a treatment group:** if any pair of proteins from **different KOs** of that reaction has an MMseqs2 hit under the filter above, or falls in the same 30% cluster, the group is excluded (the alternatives are not sequence-disjoint). A stricter optional report at 50% identity should be tabulated but is not the exclusion rule unless Yixuan prefers it.
3. **Treatment vs background leakage:** background proteins that hit any treatment protein under the same filter are **removed from the background** (not silently kept). If this removal empties a background KO, that is fine; do not promote new treatments.
4. Optionally, as a secondary check on the shortlisted groups only, run `pwalign` (same BLOSUM62, gap 11/1, Karlin–Altschul E as the package) on cross-KO pairs. Proposed default exclusion: any cross-KO pair with E ≤ 0.01. This is a consistency check with the method’s own gate; it is applied only to shortlisted candidates so the cost stays bounded.

All thresholds above are defaults to confirm, not empirical results. The build must print how many candidate groups were dropped at each filter.

---

## 5. Universe construction and embeddings

Follow the existing site conventions unless a change is recorded:

- Taxonomy: reviewed Swiss-Prot, **Bacteria or Archaea** (`taxonomy_id:2` or `taxonomy_id:2157`), non-fragment, as in `simulation_universe.R`. (Open question if Yixuan wants a broader taxon set.)
- Length: **50–1000** residues (TM-Vec head training range used by the existing scripts). Longer proteins are excluded, not padded.
- Alphabet: the 20 standard residues; drop sequences with non-standard residues.
- Unique sequences: keep one accession per unique sequence (lowest accession as tie-break).
- Proposed optional cap, to confirm: at most **N_per_KO** proteins per KO, drawn with the recorded seed after uniqueness (proposed `N_per_KO = 20`). This keeps rare and common KOs from being dominated by a few huge families, and helps the all-vs-all budget.
- Global cap: the final universe must have **at most a few thousand proteins**. Proposed hard cap, to confirm: **≤ 3000** proteins. If the rules yield more, downsample background KOs first (evenly / by recorded rank), never by post-hoc significance. If even the treatment proteins alone exceed the cap, stop and report the counts; do not silently thin the treatment.
- Embeddings: `embed_proteins(model = "tmvec")`, stored rounded to six decimal places as on the existing site. Skip embedding in `DRY_RUN` (see §10).
- Do not reuse Atlas / meta-PUF TM-Vec vectors in this website simulation unless a separate, checksummed identity check is added; the site scripts should remain self-contained.

The curated object (`curated/simulation_kofam.rds`, see licence note in §2 and §11) should contain at least:

- `proteins`: `AAStringSet` with embeddings in `mcols`.
- `annotation`: accession, length, organism, taxon, UniProt name, KO, reaction ID(s), BRITE path(s), treatment/background flag, filter provenance.
- `ko_rn`, `candidates`, `filters`: tabulated selection, or a pointer to a sidecar if KEGG-derived tables must stay private.
- `mmseqs`: cluster IDs and cross-hits used in §4.
- `alignments` (recommended): for pairs actually aligned by `protein_similarity()`, store E-value and bit score. The package itself **discards** these (`scoreOnly = TRUE` inside `.evalue()` and returns only the cosine matrix). Caching them at build time avoids repeating quadratic `pwalign` at knit, and is required if homology-only bit-score weights are ever used.
- `provenance`.

---

## 6. Simulation design

**Treatment.** Proteins belonging to the selected same-reaction KO group are the enriched set in the treatment condition. All of those KOs share the same designed effect (they are alternatives for one reaction). Do not require them to share a fold.

**Background / control.** Every other protein in the universe has **effect size 0 by design**. Do not post-select “non-significant” KOs as controls.

**Count model.** Reuse the existing site’s negative-binomial model (FunFam / EC pages): each sample has a composition

\[
p = (1-\rho)\,\pi + \rho\, q
\]

where \(\pi\) is a baseline composition (proposed default: uniform over proteins, or a recorded Dirichlet draw — pick one and fix it; do not fit it to real data after seeing outcomes), and \(q\) is uniform over the enrichment set of that sample’s condition. Counts are drawn from a negative binomial with mean \(\text{depth} \times p\) and a recorded overdispersion (the existing pages implement this as `nb_counts()`; port that helper rather than inventing a new parameterisation). Proposed depth and dispersion should be copied from the existing scripts when those values are already fixed there; otherwise mark them as proposed defaults and stop rather than inventing “realistic” numbers.

**Grids (proposed, to confirm):**

- Effect size \(\rho \in \{0, 0.01, 0.02, 0.05\}\). \(\rho=0\) is the global null.
- Sample size per condition \(n \in \{10, 20\} (planning caps; larger *n* only if the smoke/pilot budget allows).
- Replicates of the whole contrast: proposed `n_rep = 50` for the full run, `n_rep = 2` for smoke, as in Simulation 3’s “50 replicates” convention where applicable.

**Two conditions.** A two-group design: group A enriched for the treatment-reaction proteins, group B not enriched (or enriched for a different, non-homologous background set with \(\rho=0\) in the primary design). The primary design is **asymmetric enrichment vs null background**, which matches a differential-abundance setting. A two-treatment swap is optional and must be specified before looking at results if used.

**Nulls for Type I error.**

1. Global null: \(\rho = 0\), same pipeline, estimate false-positive rate for DA and for MMD.
2. Label permutation: permute condition labels on \(\rho>0\) data (or permute protein treatment labels before DA), estimate empirical Type I / FDR behaviour. Record whether permutation is of sample labels (for MMD / global DA) or of protein-level treatment membership (for clustering metrics).

**Seeds.** See §3.H. Each replicate has a derived seed.

---

## 7. Method variants (no package changes)

All graph methods should start from **one** `protein_similarity()` call on the universe, or from cached alignments, so differences are not confounded by re-alignment.

The current API on `d423023`:

```r
S <- protein_similarity(x, max_evalue = 0.01, n_cores = n_cores)
S_embed <- protein_similarity(x, max_evalue = Inf)
```

Implementation facts that the benchmark must respect:

- Cosine is clipped to `[0, 1]`; pairs with cosine `> 0` are aligned when `max_evalue` is finite. That is effectively all-vs-all among non-orthogonal embeddings, hence quadratic.
- `.evalue()` uses BLOSUM62, local alignment, gap 11/1, \(\lambda=0.267\), \(K=0.041\), database size `sum(width(seqs))`, query length = shorter protein, no length correction, no `db_size` argument. **Do not redesign this gate** in this work.
- The return value is only the cosine matrix with E > `max_evalue` edges zeroed. Bit scores are not returned.
- `bin_proteins(S)` requires a symmetric similarity in `[0, 1]`. Raw bit scores **cannot** be passed in.

Variants:

1. **Combined (default method):** `max_evalue = 0.01`, edge weights = TM-Vec cosine; `bin_proteins(S)`.
2. **Embedding only:** `max_evalue = Inf`; same cosine weights; `bin_proteins(S_embed)`.
3. **Homology only:** take the combined graph’s support (edges with E ≤ 0.01) and replace weights by **1** (unweighted). This is the proposed default because it fits `bin_proteins()` without rescaling. An optional sensitivity uses bit scores **min–max scaled to (0, 1]** from the cached alignments; do not feed raw bits to `bin_proteins()`. If alignments were not cached and a full re-alignment would exceed the 30-minute rule, skip the bit-score sensitivity and report that.
4. **Naive KO baseline:** ignore embeddings. Aggregate counts by the protein’s KO annotation (`rowsum` / equivalent), then run the **same** differential-abundance test as the bin-level analysis. This is the “just count KOs” comparator. It is expected to be strong when labels are complete and treatment KOs are the true unit; the interesting case is whether MCL bins recover the reaction group without being told the KO labels, and how MMD behaves on the embedding graph.

Do not add a `db_size` parameter or change `protein_similarity()` in this branch. A known limitation is that E-values depend on the input set’s total length; this simulation uses one fixed universe, so all variants see the same database size. Subset-vs-full-Atlas comparability is out of scope here.

---

## 8. Metrics

Define labels **before** computing scores.

**Clustering labels (primary):** one positive class = “member of the selected treatment reaction group” (all interchangeable KOs pooled); all other proteins are background. Report:

- **AUPRC** of a protein-pair or protein-level score against this label (specify one and stick to it; proposed: pairwise edges / bin co-membership as the predictor of “same treatment-reaction membership”, plus singleton handling recorded in the script).
- **ARI** and **NMI** between MCL bins and the three-or-more-way partition {KO_1, KO_2, …, background}, **and** between MCL bins and the two-way partition {treatment-reaction, background}. Primary headline = two-way treatment-reaction vs background; the per-KO partition is a secondary diagnostic (the method is not required to split interchangeable KOs).

**Differential abundance:**

- Bin-level (or KO-level for the naive baseline) test following the existing site/vignette practice (DEqMS / limma on depth-normalised log counts, or the same helper the current pages use — port, do not invent a new test).
- A bin (or KO) is a true positive if it is significantly called and is enriched for treatment-reaction proteins above a recorded purity rule (proposed default: a called bin is TP if a majority of its abundance in the contrast comes from treatment proteins; record the rule).
- Report **power** (fraction of replicates with at least one true-positive call, and/or recall of treatment proteins) and **empirical FDR** (among calls, the fraction that fail the purity rule). On \(\rho=0\) and on permutations, report the Type I rate.

**MMD:**

- `sample_mmd(counts, proteins)` as on the existing distances page, plus a permutation or \(\rho=0\) Type I estimate, and power at each \(\rho\) and *n* (e.g. PERMANOVA on the MMD distance, matching the existing page). Record PERMDISP if the existing page does.

All metrics must be computed from the curated object at **knit time** when that is cheap (clustering scores, DA, MMD on stored counts and stored graphs). Quadratic `pwalign` and TM-Vec embedding belong in the **build** script, not in the Rmd.

---

## 9. Website files, navigation, knit-time vs build-time

Claude will create a **new git worktree and branch** `sim-kofam` from `origin/gh-pages` (`e2f4059`), in a directory **outside** the dirty checkout, for example:

```sh
git -C /Users/yixuanyang/Documents/Projects/repdist worktree add -b sim-kofam \
  /Users/yixuanyang/Documents/Projects/repdist-worktrees/sim-kofam origin/gh-pages
```

**Do not** touch `/Users/yixuanyang/Documents/Projects/repdist` (dirty `gh-pages` working tree with unrelated local edits). **Do not** push, open a PR, modify `main`, or modify the existing `gh-pages` worktree. Do not submit Hazel/HPC jobs.

Files to add or edit in the new worktree only:

| File | Action |
| --- | --- |
| `data-raw/simulation_kofam.R` | **Add.** Builds `curated/simulation_kofam.rds` (+ untracked `data/kegg/` cache). |
| `simulation_kofam.Rmd` | **Add.** Page; reads the RDS; recomputes stats and figures at knit. |
| `_site.yml` | Edit the Simulations menu: add the KOfam page; **remove** FunFam and EC (and the distances page if it only exists as a FunFam/EC companion — if it can stand alone, keep it; if not, drop it from the menu but do not delete the file). Keep “MCL parameters” if it still makes sense, or note that it still points at the old universe. |
| `index.Rmd` | Replace FunFam/EC blurbs with a short KOfam simulation link. |
| `simulation.Rmd` | Hub: point to the new page; stop promoting FunFam/EC as current benchmarks. |
| `data-raw/README.md` | Document the new object, licence caveats, and build command. |
| `.gitignore` | Ensure `data/kegg/` (and any other raw KEGG cache) is ignored. |

**Keep** the old `simulation_funfam.*`, `simulation_ec.*`, and related `data-raw` / `curated` files on disk; only remove them from navigation and hub links.

Build commands (from the worktree root, with the `d423023` package installed):

```sh
Rscript data-raw/simulation_kofam.R          # full build (may stop; see §10)
KOFAM_SMOKE=1 Rscript data-raw/simulation_kofam.R
# After a successful curated RDS exists:
Rscript -e 'rmarkdown::render("simulation_kofam.Rmd")'
```

Do not run `render_site()` over the whole site unless it is fast and does not overwrite unrelated HTML unnecessarily.

---

## 10. Compute budget, smoke mode, and what is deferred

- All-vs-all `pwalign` in `protein_similarity()` is \(O(n^2)\) pairwise local alignments. Keep \(n \le 3000\) (proposed cap). A formula, not a measured time: cost ≈ \(n(n-1)/2\) alignments of typical length ≤ 1000. **Do not invent wall-clock estimates.** A cheap pilot (the smoke mode) may report alignments/second on a handful of pairs; only then may a full-run estimate be written.
- **DRY_RUN / smoke mode** (`KOFAM_SMOKE=1` or `--smoke`): same rules, but cap to tens of proteins (proposed: ≤ 30 proteins, ≤ 2 KOs if possible, `n_rep = 2`, skip TM-Vec if weights are unavailable and use random unit vectors **only** in smoke, clearly labelled as not for publication). Target: seconds to a few minutes, no large downloads beyond a tiny UniProt/KEGG probe (still rate-limited). This mode is what a reviewer should run.
- If a full run would download a large database or is projected (from the pilot or from payload size) to take **> 30 minutes**, **stop** and write an estimate plus the exact command to resume. Do not start Hazel jobs. Do not embed the full universe in an unsupervised overnight run as part of the first implementation pass.
- Deferred: E-value threshold sweep; `db_size` parameter; package changes; pushing / PRs; deleting old FunFam/EC files; any Atlas/meta-PUF embedding reuse; bit-score homology-only if alignments cannot be cached in budget.

---

## 11. Open questions for Yixuan

Short list only:

1. **KEGG publication scope.** Confirm academic eligibility and ask Pathway/Kanehisa (official form on [the academic subscription page](https://www.pathway.jp/en/academic.html); do not contact them automatically) whether the public build script, a curated RDS that includes selected KO/reaction labels, and a static results page require a service-provider agreement or written redistribution permission, and which derived fields may be public. Until then, follow the conservative workflow in §2.2 (raw KEGG untracked; public page prefers aggregate results).
2. **Confirm proposed defaults:** homology filters in §4 (MMseqs2 E `1e-3`, 80% coverage, exclude treatment groups with any cross-KO 30% cluster or hit; leaky background proteins dropped); protein caps in §5 (50–1000 aa, prokaryotes only, ≤ 20 proteins/KO, ≤ 3000 total); grids in §6 (\(\rho\), *n*, `n_rep`); ranking key and *k* in §3.G; seed `20261112L`.
3. **Homology-only weights:** accept the default of unweighted edges (1), with optional scaled bit scores as a sensitivity?
4. **MCL parameters page / distances page:** leave them in the menu pointing at the old universe, or later retarget them to the KOfam object (out of scope for the first patch unless trivial)?

The E-value threshold of 0.01 is **not** an open question for this simulation; it stays as implemented, to be benchmarked later.

---

## 12. Implementation checklist for Claude

Working rules:

- [ ] Create a new worktree + branch `sim-kofam` from `origin/gh-pages` **outside** `/Users/yixuanyang/Documents/Projects/repdist`. Do not dirty that checkout, do not checkout over it, do not push, do not open a PR, do not edit `main`.
- [ ] Copy this plan into the worktree (e.g. `data-raw/KOFAM_SIM_PLAN.md` or the root) so the branch is self-contained.
- [ ] Implement `data-raw/simulation_kofam.R` with the algorithm in §§3–6, untracked KEGG cache, provenance, `KOFAM_SMOKE` mode, and a hard stop with an estimate if the full run would exceed 30 minutes or a large download.
- [ ] Do not hard-code KO or reaction IDs. Look up ko-rn and BRITE from the official URLs. Do not use `ko01000` as a grouping key.
- [ ] Implement homology / leakage checks with MMseqs2; record versions and commands.
- [ ] Implement `simulation_kofam.Rmd` with the four method variants and the metrics in §§7–8. Port count-model and DA/MMD helpers from the committed FunFam/EC pages rather than inventing new ones. If `curated/simulation_kofam.rds` does not yet exist (full build stopped), the Rmd should knit a clearly labelled placeholder / smoke result, or skip chunks behind a flag, rather than fail on missing fake numbers.
- [ ] Update `_site.yml`, `index.Rmd`, `simulation.Rmd`, `data-raw/README.md`, and `.gitignore`. Remove FunFam and EC from the **navigation and hub links only**; do not delete the old files.
- [ ] Local commit(s) on `sim-kofam` only, with a clear message. No `git push`. No force operations on other branches.
- [ ] Run the smoke mode. If it passes and a small knit is possible, knit `simulation_kofam.Rmd`. Do not start the full quadratic embedding+alignment unless it is clearly under 30 minutes.
- [ ] Report: files changed, commit hash, smoke output, anything stopped for budget, and any default that had to be interpreted.

Package install note: the worktree is a website branch and may not contain `R/`. Install `repdist` from `/Users/yixuanyang/Documents/Projects/repdist/.claude/worktrees/main` (commit `d423023` or its successor on that worktree), without committing unrelated review-ben edits.
