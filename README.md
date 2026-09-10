# repdist analyses

The `gh-pages` branch: the [repdist](https://github.com/Yixuan39/repdist) site,
served at <https://yixuan39.github.io/repdist/>.

An `rmarkdown` website (`_site.yml`, Bootstrap `cosmo`) rendered in place, so
all three pages share one theme and navbar:

| Page | Source |
| --- | --- |
| [Home](https://yixuan39.github.io/repdist/) | `index.Rmd` |
| [Data curation](https://yixuan39.github.io/repdist/data_curation.html) | `data_curation.Rmd` |
| [Case study](https://yixuan39.github.io/repdist/case_study.html) | `case_study.Rmd` |

Rebuild with:

```r
rmarkdown::render_site()
```

That writes `index.html`, the two report HTMLs, `site_libs/` and
`*_files/` into the branch root. `README.md`, `figures/` and `data/` are
excluded in `_site.yml`.

Both reports use the mouse faecal metaproteome of Blakeley-Ruiz *et al.*,
*ISME J* **19**(1) wraf048 and read `data/test_study/` relative to the working
directory. Those inputs are not redistributed here, so a rebuild needs them
put in place first; `data_curation.Rmd` must run before `case_study.Rmd`,
which reads the phyloseq object it writes. The package itself, its vignettes
and its tests are on `main`.
