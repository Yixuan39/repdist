# repdist site

The `gh-pages` branch holds the [repdist](https://github.com/Yixuan39/repdist)
documentation site, served at <https://yixuan39.github.io/repdist/>. The
package itself, its vignette and its tests are on `main`.

| Page | Source | Reads |
| --- | --- | --- |
| Home | `index.Rmd` | |
| Reference | `reference.Rmd` | `curated/simulation_universe.rds` |
| Simulation | `simulation.Rmd` | `curated/simulation_funfam.rds`, `curated/simulation_universe.rds`, `curated/simulation_mmseqs2.rds` |
| MCL parameters | `parameters.Rmd` | `curated/simulation_universe.rds` |
| Case study | `case_study.Rmd` | `curated/diet_yeast.rds` |

```
data-raw/   scripts that build curated/ from the public sources (see its README)
curated/    curated input objects read by the pages
supplement/ full result tables written by the pages at render time
data/       downloaded source tables (not committed)
```

Rebuild from scratch, from the branch root, with the current `repdist`
installed:

```sh
Rscript data-raw/simulation_universe.R
Rscript data-raw/simulation_mmseqs2.R
Rscript data-raw/simulation_funfam.R
Rscript data-raw/diet_yeast.R
Rscript -e 'rmarkdown::render_site()'
```

The pages compute every result at render time from the curated objects; no
knitr cache is used.
