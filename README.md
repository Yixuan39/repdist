# repdist site

The `gh-pages` branch holds the [repdist](https://github.com/Yixuan39/repdist)
documentation site, served at <https://yixuan39.github.io/repdist/>. The
package itself, its vignette and its tests are on `main`.

Navbar: Get started | Install | Reference | Case study | Simulated data |
GitHub.

| Page | Source | Reads |
| --- | --- | --- |
| Home and Install | `index.Rmd` | |
| Get started | `repdist.html`, rendered from `main`'s `vignettes/repdist.Rmd` | the package's `simulate_cog` |
| Reference | `reference.Rmd` | the package's `simulate_cog` |
| Case study | `case_study.Rmd` | `curated/diet_series.rds` |
| Simulated data | `data_simulation.Rmd` | the NCBI COG 2024 release, downloaded into `data/cog/` |

Each page renders on its own; no page reads another page's output.
`data_simulation.Rmd` writes `data_simulation/`, which the package's
`inst/scripts/make_simulate_cog.R` turns into `simulate_cog`.

```
data-raw/   script that builds curated/ from the public source (see its README)
curated/    curated input object read by the case study
data/       downloaded source tables (not committed)
```

Rebuild from scratch, from the branch root, with the current `repdist`
installed:

```sh
Rscript data-raw/diet_series.R
Rscript -e 'rmarkdown::render_site()'
git show main:vignettes/repdist.Rmd > repdist.Rmd
Rscript -e 'rmarkdown::render("repdist.Rmd")' && rm repdist.Rmd
```

`repdist.Rmd` is not committed: the vignette's source stays on `main`, and
only its rendered page lives here. `data_simulation.Rmd` downloads about
2.2 GB once.

The pages compute every result at render time from these inputs; no knitr
cache is used.
