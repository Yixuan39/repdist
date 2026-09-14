# Recompute the fixed FunFam sequence set with the historical head identity.
# See DATA_PROVENANCE.md: the original run did not record model revisions.
# This produces new embeddings, not a promise of bitwise historical recovery.
# Run from the source root after installing repdist; requires cached/downloaded
# ProtT5 (several GB) and scikit-bio/tmvec weights. No private intermediates.
library(repdist)
u <- readRDS(file.path("inst", "extdata", "funfam_universe.rds"))
stopifnot(identical(names(u$sequences), u$members$protein))
spec <- known_models()[["tmvec-300"]]
spec$description <- "Historical simulation head: scikit-bio/tmvec"
spec$head_repo <- "scikit-bio/tmvec"
spec$head_revision <- "405ecd7651b378cdf3e3cb544d5477c73ac384bd"
# The historical head is intentionally not added as a new public preset.
proc <- basilisk::basiliskStart(get("repdist_env", asNamespace("repdist")))
Z <- tryCatch(
    basilisk::basiliskRun(
        proc, get(".repdist_embed_worker", asNamespace("repdist")),
        as.list(unname(u$sequences)),
        as.character(jsonlite::toJSON(spec, auto_unbox = TRUE)),
        spec$head_repo, "cpu", 16L
    ),
    finally = basilisk::basiliskStop(proc)
)
stopifnot(
    is.matrix(Z), nrow(Z) == length(u$sequences), all(is.finite(Z)),
    max(abs(rowSums(Z^2) - 1)) < 1e-6
)
dimnames(Z) <- list(names(u$sequences), paste0("d", seq_len(ncol(Z))))
u$embeddings <- round(Z, 4)
u$provenance$embedding_model <- spec
u$provenance$model_revision <- spec$head_revision
u$provenance$repdist <- as.character(packageVersion("repdist"))
u$provenance$session <- capture.output(sessionInfo())
saveRDS(u, file.path("inst", "extdata", "funfam_universe.rds"), compress = "xz")
