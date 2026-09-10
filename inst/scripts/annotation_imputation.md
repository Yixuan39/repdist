# Case study 2: label propagation and reliability scoring

Updated 2026-09-09. The rewritten `impute_labels()` separates KNN label
propagation from a shared binary model of whether the propagated label agrees
with the annotation. Selection-weighted random forest had the lowest pooled
Brier score in this evaluation, 0.13668, compared with 0.14208 for unweighted
spline logistic regression. The advantage is modest and depends on the metric;
it does not establish calibration on genuinely unannotated proteins.

## API and model design

```r
spline <- impute_labels(D, labels, k = 1L, scoreMethod = "spline")
rf <- impute_labels(D, labels, k = 1L, scoreMethod = "rf",
                    biasAdjust = TRUE)
attr(rf, "models")
attr(rf, "calibration")
attr(rf, "diagnostics")
```

`D` is a named cosine-distance `dist` object or symmetric matrix; `labels`
is a named character vector or factor containing original annotations.
Each returned row is an unannotated protein. The two model blocks are `knn`
and `score`; the latter also records the auxiliary annotation-selection model.
The scoring method never changes the propagated label or rejects a prediction.
Original annotations are the only propagation references. The former SVM
label classifier and automatic KNN parameter search have been removed.

The default `k = 1` transfers the nearest reference label. Exact distance ties
are included; vote ties use nearest support, then protein ID. General `k`
remains available, but this experiment fixes it at 1. The output preserves the
supporting protein, similarity and neighbour counts, and adds the following
scoring evidence:

| Feature | Definition |
|---|---|
| Similarity | `1 - D` to the supporting reference |
| Margin | Supporting similarity minus the nearest alternative-label similarity |
| Label support | `log1p()` of the candidate label's reference count |
| Local similarity | Similarity to the nearest catalog protein, excluding self |
| Agreement | Candidate-label fraction among the top five references, including distance ties |
| Alternative flag | Whether any reference has a different label |

Spline logistic regression uses a three-degree-of-freedom natural cubic spline
for supporting similarity and linear terms for the remaining features.
`quasibinomial()` provides the same fitted logit mean while allowing fractional
weights. Natural splines extrapolate linearly beyond their boundary knots on
the linear-predictor scale; this is an extrapolation assumption, not additional
validation. See the [R natural-spline documentation](https://stat.ethz.ch/R-manual/R-devel/library/splines/html/ns.html).

RF uses `ranger` probability forests with 500 trees, `min.node.size = 20`,
default feature subsampling, and one thread. The node-size setting controls
when nodes can split, not a guaranteed minimum terminal-node size. Forest
probabilities average tree probability estimates. Selection weights enter
through bootstrap sampling via `case.weights`. No post-hoc Platt or isotonic
calibrator, tuning grid, or outcome-class balancing was added; the outer
evaluation directly tests these probability estimates. See the
[ranger documentation](https://imbs-hl.github.io/ranger/reference/ranger.html).

The function uses existing `dist` and `data.frame` representations, validates
inputs, preserves R's random-number state, and records fit failures instead
of inventing confidence values. Insufficient scoring data produce `NA` scores
with a warning while retaining the propagated labels. These choices follow
the [Bioconductor R-code guidance](https://contributions.bioconductor.org/r-code.html)
on reusing representations and dependencies, small internal helpers, argument
validation, and reproducible behavior. This design statement is not a claim
of Bioconductor submission approval.

## Data and independent evaluation

The cached case-study-2 catalog contains 5,245 proteins: 2,033 with consensus
KO annotations and 3,212 without them. There are 701 distinct annotated KO
labels, including 440 singleton labels. Consensus KO requires agreement
between annotation pipelines; it is not an experimentally curated gold
standard. The EC benchmark is not included in this report.

Five outer folds hold out complete single-linkage components at TM-Vec cosine
similarity at least 0.95. Components are constructed from the entire catalog,
so unannotated bridge proteins cannot connect different folds. There are
3,586 components, of which 1,171 contain annotated proteins; the largest has
99 proteins. Component assignment uses annotation counts, protein IDs and
seed 20260909, without label identities. TM95 here is a threshold on the
representation's similarity, not 95% sequence identity or a measured TM-score.

Each outer fit receives only the 1,626 or 1,627 remaining annotations. Its
five inner grouped folds generate correctness-training records using references
outside the corresponding fold. Annotated and unannotated inner queries use
the same reference pool. No outer-test annotation enters inner propagation,
correctness fitting, or selection fitting. Test distances and annotation
availability remain observable, as in the intended transductive catalog task.

All four scoring variants use the same outer folds, features, `k = 1`, inner
fold seed and weight cap of 10. Baselines are the inner-OOF mean correctness
and raw supporting similarity clipped to [0, 1]. No model or threshold is
tuned on the outer test outcomes. Every annotated protein is tested once.

## Performance

All methods propagate exactly the same labels: 965 of 2,033 are correct,
giving 47.47% exact-label accuracy and 100% propagation coverage. Among 1,261
test proteins whose true label occurs in the outer references, accuracy is
76.53%. The remaining 772 test proteins, 37.97%, have unseen true labels and
cannot be assigned their correct label by closed-reference propagation.
Outer-fold accuracy ranges from 39.56% to 59.21%.

| Score model | Brier | Log loss | Selection-weighted Brier | Selection-weighted log loss |
|---|---:|---:|---:|---:|
| Constant baseline | 0.25336 | 0.70002 | 0.22506 | 0.64287 |
| Raw similarity baseline | 0.30786 | 0.85414 | 0.36606 | 0.97662 |
| Spline logistic | 0.14208 | 0.44172 | 0.12556 | 0.40076 |
| Random forest | 0.14443 | 0.44941 | 0.12472 | 0.40607 |
| Selection-weighted spline | 0.14158 | 0.44804 | 0.12439 | 0.40041 |
| Selection-weighted RF | 0.13668 | 0.43395 | 0.12187 | 0.40261 |

Lower loss is better. Losses are pooled over all held-out proteins; log-loss
probabilities are clipped to [1e-6, 1 - 1e-6] for evaluation only. Selection-
weighted evaluation uses a common spline selection model from the unweighted
spline run, applied to each test protein's actual outer-prediction features.
All models receive the same evaluation weights, with pooled effective sample
size 1,415.2. These weights describe a mixture of masked annotated proteins
and genuinely unannotated proteins, so the weighted columns are a covariate-
shift sensitivity analysis, not measured accuracy on the unannotated target.

Selection-weighted RF improves ordinary Brier score over unweighted spline
in four of five outer folds, by 0.00540 overall. It is a reasonable candidate
for this case study. However, weighted spline has the lowest weighted log
loss, and differences between fitted scorers are much smaller than their
improvement over the baselines. This is one grouped partition and seed, with
overlapping training sets; no significance or general superiority claim is
made. The API retains unweighted spline as its default and exposes RF and
selection weighting explicitly.

![Outer held-out score reliability and selection-weighted sensitivity analysis](../../results/annotation_imputation_scoring/calibration.png)

The curves use fixed probability bins of width 0.1. They show useful score
separation but imperfect calibration, including underconfidence in several
middle-score bins. Bin counts and values are available in `calibration.csv`.
For example, among weighted-RF predictions with score at least 0.9, 213 of
219 labels are correct: 97.26% accuracy at 10.77% coverage. Unweighted spline
at the same fixed threshold selects 299 proteins with 94.65% accuracy.
These are descriptive held-out results, not guaranteed deployment precision.

Mean outer-fit runtime is 4.71 seconds for spline, 5.04 for RF, 4.70 for
weighted spline, and 5.11 for weighted RF. Runtime includes feature extraction,
inner calibration, selection fitting, scorer fitting and prediction; it
excludes initial distance loading and outer grouping. These are single-run
wall-clock timings on the recorded local environment.

## Annotation bias and full-catalog predictions

The similarity shift is substantial. In the full-fit calibration records,
median nearest-any-protein similarity is 0.95964 for annotated proteins and
0.87878 for unannotated proteins. With matched grouped reference exclusion,
median supporting similarity is 0.85703 versus 0.71296. The final unannotated
predictions, which can use every annotation, have median supporting similarity
0.72967. These quantities distinguish biological neighborhood density from
the availability of usable annotated references.

The auxiliary spline model estimates annotation probability `e(X)` with
grouped cross-fitting. Scorer-training weights are

```text
min((1 - e(X)) / e(X) * nAnnotated / nUnannotated, 10)
```

The prior ratio uses each selection model's training population. In the
full fit, annotated training-weight effective sample size is 1,105.3 out
of 2,033; no annotated weight reaches the cap. Weighting therefore reduces
the effective information even though the raw dataset has two thousand
training records. It assumes feature overlap and equal conditional
correctness between annotated and unannotated populations. Unobserved
novel functions or annotation-pipeline errors can violate that assumption.
Unannotated proteins are never coded as incorrect training examples.

All four final fits use all 2,033 original annotations and output the same
3,212 candidate labels. Their scores differ:

| Full-fit score model | Median score | Score >= 0.8 | Score >= 0.9 |
|---|---:|---:|---:|
| Spline logistic | 0.13785 | 505 | 232 |
| Random forest | 0.17272 | 426 | 256 |
| Selection-weighted spline | 0.12803 | 436 | 246 |
| Selection-weighted RF | 0.17810 | 363 | 186 |

All models flag 293 unannotated proteins, 9.12%, as outside at least one
source-calibration feature range. This includes the 282 proteins whose
supporting similarity is at least 0.95; grouped calibration has no such
annotated examples. Of the 186 weighted-RF predictions scoring at least 0.9,
44 carry the extrapolation flag. The flag checks marginal ranges only and
does not establish joint-distribution overlap. Even unflagged predictions
on actual unknowns lack observed correctness outcomes. Full-reference
predictions also have larger reference pools than inner calibration records.
Independent annotations are needed before treating these scores as verified
biological probabilities.

## Reproduction and artifacts

Run from the repository root with the cached research input and the installed
R dependencies, including `ranger`:

```sh
Rscript inst/scripts/annotation_imputation.R
```

The script writes to `results/annotation_imputation_scoring/`. Main outputs
are `test_summary.csv`, `fold_metrics.csv`, `test_predictions.csv`,
`calibration.csv`, `calibration.png`, `similarity_strata.csv`, `thresholds.csv`,
`selection_diagnostics.csv`, `unannotated_predictions.csv`, and `models.rds`.
Full calibration records are saved separately for each scorer. `splits.csv`,
`cohort.csv`, `runtimes.csv`, `diagnostics.rds`, input and code MD5 fingerprints,
and `sessionInfo.txt` preserve the evaluation context.

The completed benchmark checks that all references are original training
annotations, every supporting label matches its recorded reference, no outer
test label appears in calibration outcomes, every test protein receives a
finite score in [0, 1], and all four scorers return identical propagated
labels. Numerical smoke checks verified the metric calculations. Input and
source fingerprints were identical before and after the completed run. A
verification rerun reproduced byte-identical summary, fold-metric, calibration,
threshold, similarity-stratum and unannotated-prediction CSV files.
Previous benchmark result directories remain untouched; their older API and
split-specific numbers do not describe this evaluation.

## Package validation

The final source package was built with both vignettes evaluated and checked
on R 4.5.2, macOS arm64, Bioconductor 3.22 and BiocCheck 1.46.3:

```sh
R CMD build --no-manual /path/to/repdist
R CMD check --no-manual repdist_0.99.0.tar.gz
Rscript -e 'BiocCheck::BiocCheck("repdist_0.99.0.tar.gz", `no-check-bioc-help`=TRUE)'
```

`R CMD check` completed with **0 errors, 0 warnings and 1 note**. All 638 test
expectations passed, with no failures, warnings or skips; examples and vignette
rebuilding also passed. The note reports the package's existing non-staged
installation setting. Rounded embeddings in the two cached vignette inputs
are now explicitly normalized to satisfy the current unit-row contract.

`BiocCheck` completed with **0 errors, 0 warnings and 6 notes**. These suggest
an additional biocView, maintainer ORCID, a funder role if applicable, shorter
functions, shorter lines and four-space indentation. The function-length note
includes the new 75-line `impute_labels()` orchestration function and two
existing functions; package-wide style and author metadata were not rewritten
for this task. Maintainer support-site and mailing-list account checks were
explicitly skipped. These checks therefore document local package validation,
not Bioconductor submission approval or cross-platform certification.

Imputation tests cover input validation, KNN ties, protein-order invariance,
group exclusion including unannotated bridges, score-method label invariance,
selection weights, RNG preservation, insufficient calibration, failed fits,
and non-estimable predictions. Check and test logs are retained in
`results/annotation_imputation_scoring/checks/`.

The repository's required `graphify update .` maintenance command also
completed. Its installed parser does not support the R source files, so this
index update is not evidence of R dependency coverage.
