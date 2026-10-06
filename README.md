# repdist analyses

The `gh-pages` branch: the [repdist](https://github.com/Yixuan39/repdist) site,
served at <https://yixuan39.github.io/repdist/>.

An `rmarkdown` website (`_site.yml`, Bootstrap `cosmo`) rendered in place, so
all pages share one theme and navbar:

| Page | Source |
| --- | --- |
| [Home](https://yixuan39.github.io/repdist/) | `index.Rmd` |
| [Simulation](https://yixuan39.github.io/repdist/simulation.html) | `simulation.Rmd` |
| [Real-data validation](https://yixuan39.github.io/repdist/case_study.html) | `case_study.Rmd` |
| [Data curation](https://yixuan39.github.io/repdist/data_curation.html) | `data_curation.Rmd` |

Rebuild with:

```r
rmarkdown::render_site()
```

That writes the page HTMLs, `site_libs/` and `*_files/` into the branch
root. `README.md`, `figures/` and `data/` are excluded in `_site.yml`.

The simulation page reads `validation/funfam_universe.rds`; `validation/`
also holds the EC3 benchmark and case-study scripts, which are not rendered.
The real-data pages use the mouse faecal metaproteome of Blakeley-Ruiz
*et al.*, *ISME J* **19**(1) wraf048 and read `data/test_study/` relative to
the working directory. Those inputs are not redistributed here, so a rebuild needs them
put in place first; `data_curation.Rmd` must run before `case_study.Rmd`,
which reads the phyloseq object it writes. The package itself, its vignette
and its tests are on `main`.
