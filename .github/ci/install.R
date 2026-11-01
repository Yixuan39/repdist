# The README's development-version install, verbatim apart from the ref.
t0 <- Sys.time()
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("basilisk", "BiocFileCache", "BiocParallel",
    "Biostrings", "pwalign", "rhdf5", "S4Vectors", "MCL"),
    ask = FALSE, update = FALSE)
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("Yixuan39/repdist", ref = Sys.getenv("REF"),
    dependencies = TRUE, upgrade = "never")
stopifnot(requireNamespace("repdist"))
message("install took ", format(round(Sys.time() - t0, 1)))
