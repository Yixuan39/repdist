# repdist site

The `gh-pages` branch holds the [repdist](https://github.com/Yixuan39/repdist)
documentation site, served at <https://yixuan39.github.io/repdist/>. The
package itself, its vignette and its tests are on `main`.

Navbar: Install | Reference | Case study | Simulations (Hidden coherence,
Sequence-disjoint functions, MCL parameters) | GitHub.

| Page | Source | Reads |
| --- | --- | --- |
| Home and Install | `index.Rmd` | |
| Reference | `reference.Rmd` | `curated/simulation_universe.rds` |
| Case study | `case_study.Rmd` | `curated/diet_series.rds` |
| Simulation 1: Hidden coherence | `simulation_funfam.Rmd` | `curated/simulation_funfam.rds` |
| Simulation 2: Sequence-disjoint functions | `simulation_ec.Rmd` | `curated/simulation_universe.rds`, `curated/simulation_mmseqs2.rds` |
| MCL parameters | `parameters.Rmd` | `curated/simulation_universe.rds` |
| Old simulation link | `simulation.Rmd` | nothing; links to the three pages above |

Each page renders on its own; no page reads another page's output.

```
data-raw/   scripts that build curated/ from the public sources (see its README)
curated/    curated input objects read by the pages
data/       downloaded source tables (not committed)
```

Rebuild from scratch, from the branch root, with the current `repdist`
installed:

```sh
Rscript data-raw/simulation_universe.R
Rscript data-raw/simulation_mmseqs2.R
Rscript data-raw/simulation_funfam.R
Rscript data-raw/diet_series.R
Rscript -e 'rmarkdown::render_site()'
```

The pages compute every result at render time from the curated objects; no
knitr cache is used.
