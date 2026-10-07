# Simulation data and the EC3 benchmark (not part of the package)

Experimental work on the `gh-pages` branch, outside the package; it did not
choose any MCL setting. Run everything from the branch root with the current
`repdist` installed:

* `prepare_microbial_ec_universe.R` builds `microbial_ec_universe.rds`
  (needs network access and the TM-Vec model).
* `simulate_ec3.R` holds the design and count simulation.
* `benchmark_simulation.R` writes `simulation_benchmark.rds`, which
  `benchmark_report.Rmd` renders.
* `testthat::test_file("validation/test-simulation.R")` checks the design.
* `case_study1.R` runs the full diet metaproteome and writes
  `case_study1/` (not committed), which `case_study1_report.Rmd` renders.

Render a report with `rmarkdown::render("validation/benchmark_report.Rmd")`;
the home page links both HTML files.

## Microbial experimental-EC universe

Source: UniProtKB REST API (<https://rest.uniprot.org>), release 2026_03,
accessed 2026-09-25. UniProt content is CC BY 4.0
(<https://www.uniprot.org/help/license>). Built by
`prepare_microbial_ec_universe.R`; every rule below is also recorded in the
RDS `provenance`.

1. **Candidates.** `(cc_catalytic_activity_exp:*) AND (fragment:false) AND
   (taxonomy_id:2 OR taxonomy_id:2157)`: Bacteria and Archaea, non-fragment,
   with at least one catalytic-activity comment carrying experimental
   evidence. 18,779 entries.
2. **Experimental EC.** An EC number counts as experimental only when it sits
   on a CATALYTIC ACTIVITY reaction whose own evidence list includes
   ECO:0000269 (experimental evidence used in manual assertion). EC numbers in
   the protein name, and reactions with only curator-inferred (ECO:0000305),
   sequence-similarity or automatic evidence, never qualify. Review status is
   recorded but is not a filter; every entry returned happens to be Swiss-Prot,
   because only curators assign ECO:0000269. 6,739 entries have an experimental
   reaction EC.
3. **One EC3.** The experimental EC numbers, reduced to their first three
   levels, must give exactly one EC3, and every EC number on the entry —
   experimental or not, in reactions or in the protein name, including
   multi-domain parts — must agree with it on every level it specifies
   (`3.1.-.-` agrees with 3.1.3; `2.7.-.-` does not). 6,432 entries.
4. **Length.** 50 to 1000 residues, the range TM-Vec's head was trained on.
   A choice for this universe, not a package limit. 6,313 entries.
5. **Redundancy.** Mapped to UniRef50 with the UniProt ID-mapping service and
   reduced to one protein per UniRef50 cluster (lowest accession), so no two
   proteins are near-duplicates that sequence identity alone would group.
   5,324 proteins over 229 EC3s.
6. **Eligibility.** A simulation target needs 10 samples x 5 private
   upregulated proteins = 50 proteins. Every EC3 with at least 50 proteins
   after step 5 is kept, and 50 of each are drawn with seed 20260924. The rule
   is a count, fixed before embedding, so no EC3 is in or out because of how
   TM-Vec treats it. 28 EC3s, 1,400 proteins.
7. **Embedding.** `embed_proteins(model = "tmvec")` on the full sequences,
   rounded to four decimals.

## Simulated effect size and example target

Ten A samples each upregulate five private target proteins, covering the 50
target proteins exactly once. A 31-fold per-protein mean increase gives a
size-factor-normalised group-mean target aggregate of
`1 + (31 - 1) / 10 = 4` times baseline for any baseline abundances. Individual
sample effects vary. The effect size is defined by this ground-truth algebra,
never calibrated against a distance or test result; test-simulation.R checks it.

Before simulating counts, the worked target is the EC3 whose median within-EC3
predicted TM-score is closest to the 75th percentile across the 28 EC3s, ties
in sorted EC3 order. This selects EC 2.1.1 (median similarity 0.430790) without
using p-values, sample distances or clustering outcomes. EC 2.7.1 remains a
less coherent comparison. Neither this example selection nor any EC label
determines the graph edge floor.

## Simulation benchmark

`simulation_benchmark.rds` is the output of `benchmark_simulation.R`:
the EC3 simulation rerun with every EC3 of the
universe as the target (seed 1000 + target index), plus a grid of
`min_sim` and inflation values scored against EC3. It records its
settings and R session.

The universe stores embeddings rounded to four decimals. Renormalize rows
before cosine distances; rounding can move pairs across a sufficiently close
edge floor.

Run `Rscript validation/benchmark_simulation.R` after installing the current
package. Its `simulation_benchmark.rds` records all 28 target EC3s,
five distance methods, PERMANOVA R2/p, dispersion and within/between-group
distances, protein/bin differential abundance, and a graph-floor/inflation
sensitivity grid. This detail is kept outside the vignettes.

The benchmark reference setting (edge floor 0.5, inflation 2) is an explicit
analysis setting, not a package default or a recommendation. Sensitivity is
scored against labels for validation only; never select the floor using EC
agreement or downstream detection.

The 31-fold private-protein effect gives exactly fourfold group-mean EC3
ground truth before independent sample size factors and count noise. It is
fixed by algebra, with a runnable check in test-simulation.R. PERMANOVA may
reflect dispersion as well as group location in this concentrated-versus-
diffuse design; inspect the stored PERMDISP and pairwise-distance summaries.

## FunFam universe (`funfam_universe.rds`, the simulation page)

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
The stored `proteins` supply the exact ungapped strings and order as a
reconstruction manifest. Historical selection metadata (seed 20260822, 150 per
family and length filter 60--510) is retained in `provenance`; reproducing the
fixed set uses its recorded identifiers rather than depending on sampling order
in a potentially rearranged server response.

The historical embeddings used ProtT5 XL UniRef50 and the `scikit-bio/tmvec`
head, not the current default Figshare head or `scikit-bio/tmvec-swissmodel`.
`prepare_funfam_universe.R` can recompute the exact sequence set with the same
head identity pinned to commit 405ecd7651b378cdf3e3cb544d5477c73ac384bd and the
registry's pinned ProtT5 backbone. This revision was verified against
the current registry; it is not asserted to be the unrecorded historical revision.
New results record the complete model specification and R session.

Embeddings are rounded to four decimals; renormalize rows before cosine
distances.
