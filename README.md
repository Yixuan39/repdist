# repdist analyses

The `gh-pages` branch: rendered re-analyses of published metaproteomes with
[repdist](https://github.com/Yixuan39/repdist), served at
<https://yixuan39.github.io/repdist/>.

| Report | Source |
| --- | --- |
| [Data curation](https://yixuan39.github.io/repdist/data_curation.html) | `inst/scripts/data_curation.Rmd` |
| [Case study — mouse faecal metaproteomes](https://yixuan39.github.io/repdist/case_study.html) | `inst/scripts/case_study.Rmd` |
| [Case study 2 — marine plastic substrates](https://yixuan39.github.io/repdist/case_study2.html) | `inst/scripts/case_study2.Rmd` |
| [Case study 2 — annotation-oriented clustering](https://yixuan39.github.io/repdist/case_study2_annotation_report.html) | `inst/scripts/case_study2_annotation_report.Rmd` |

Supporting scripts: `build_trees.R` (the k-mer and TM-Vec trees the UniFrac
comparisons need), `case_study2_annotation.R` (writes
`results/case_study2_annotation_clustering/` for the annotation report) and
`case_study_bin_ordination.R`.

The sources expect to run from a repdist checkout with the study inputs under
`data/`; those inputs are not redistributed here. The package itself, its
vignettes and its tests are on `main`.
