# repdist analyses

The `gh-pages` branch: the [repdist](https://github.com/Yixuan39/repdist) site,
served at <https://yixuan39.github.io/repdist/>.

The landing page is an `rmarkdown` website (`_site.yml` + `index.Rmd`, Bootstrap
`cosmo`), rendered in place. To rebuild it after editing `index.Rmd`:

```r
rmarkdown::render_site()
```

That writes `index.html` and `site_libs/` into the branch root. It renders
`index.Rmd` only — `README.md` and `inst/` are excluded in `_site.yml`, and the
report HTMLs are committed renders, not rebuilt by the site.

| Report | Source |
| --- | --- |
| [Data curation](https://yixuan39.github.io/repdist/data_curation.html) | `inst/scripts/data_curation.Rmd` |
| [Case study](https://yixuan39.github.io/repdist/case_study.html) | `inst/scripts/case_study.Rmd` |

Both use the mouse faecal metaproteome of Blakeley-Ruiz *et al.*, *ISME J*
**19**(1) wraf048, under `data/test_study/`.

The reports expect to run from a repdist checkout with those inputs in place;
they are not redistributed here, so the committed HTML cannot be regenerated
from this branch alone. The package itself, its vignettes and its tests are on
`main`.
