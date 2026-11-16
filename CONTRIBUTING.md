# Contributing to repdist

Standard R package development: `devtools::load_all()`, `devtools::test()`,
`devtools::check()`. `BiocCheck::BiocCheck()` runs the Bioconductor-specific
checks.

## Tests

`devtools::test()` runs everything in `tests/testthat/` with the fabricated
`AAStringSet`/embedding fixtures in `helper-repdist.R` -- no network, no
model download. The one exception is `test-embed-model.R`, which calls
`embed_proteins()` against a real registered model through `basilisk`. It is
skipped unless the `REPDIST_CHECK_MODEL` environment variable names a model
(see the comment at the top of that file), since that is the only test that
downloads model weights.

## Continuous integration

`.github/workflows/R-CMD-check.yaml` runs `R CMD check` on every push and pull
request to `main`. A passing push to `main` also uploads coverage to Codecov.
The site's `repdist.html` is rendered by hand on the `gh-pages` branch.

## Vignette

`tools::buildVignettes()` or `devtools::build_vignettes()`.
`vignettes/repdist.Rmd` runs on `data(simulate_cog)` (110 proteins) and
builds in seconds, `protein_similarity()` included. `inst/scripts` has the
scripts that build the shipped data.

The README figure is the vignette's `pcoa` chunk; after changing it, from
the package root with the package loaded:

```r
source(knitr::purl("vignettes/repdist.Rmd", output = tempfile()))
ggplot2::ggsave("man/figures/README-simulation-pcoa.png", width = 9,
    height = 3.4, dpi = 150)
```
