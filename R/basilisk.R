# The Python side of embed_proteins() lives in a basilisk-managed virtual
# environment rather than in whatever Python reticulate happens to find. Every
# version is pinned, so the stack Bioconductor's build machines install is the
# stack a user gets, and it cannot be perturbed by the rest of their system.
#
# Bumping any pin here changes the embeddings, so treat it as a version bump of
# the ground metric, not as routine maintenance: the coordinates cached in
# inst/extdata were produced by the stack below. Model weights are pinned
# separately, to commit revisions in inst/extdata/models.json.
#
# basilisk is also why DESCRIPTION sets `StagedInstall: no`: its environments
# cannot carry a staged install path. R CMD check's NOTE about it is expected.
#
# This file is sourced standalone by configureBasiliskEnv() from ./configure, so
# it must not assume the rest of the package is loaded. There is deliberately no
# configure.win: R runs it through sh, which Windows has only with Rtools or
# git, so it made a GitHub install fail on a bare Windows machine.
# `repdist_env` carries no leading dot on purpose: configureBasiliskEnv()
# collects environments with ls(), which skips dot-prefixed names, and would
# silently build nothing.
#
# One environment, not one per model: every entry in inst/extdata/models.json
# runs on these pins. Split only when a model needs pins that genuinely cannot
# coexist with these -- configure builds every environment declared here, so an
# extra one costs another multi-GB torch on the Bioconductor builders.
# 3.11.9 is the last 3.11 release with a python.org Windows installer, which
# basilisk's pyenv route needs there; 3.11.10+ are source-only.
repdist_python <- "3.11.9"

repdist_env <- basilisk::BasiliskEnvironment(
    envname = "repdist-embed",
    pkgname = "repdist",
    packages = paste0("python=", repdist_python),
    pip = c(
        "torch==2.5.1",
        "transformers==4.46.3",
        "safetensors==0.4.5",
        "huggingface_hub==0.26.5",
        "numpy==1.26.4",
        # Rostlab/prot_t5_xl_uniref50 ships only a SentencePiece
        # `spiece.model`, so transformers has to convert it into a fast
        # tokenizer at load time and needs both of these to do it. Without
        # protobuf the conversion silently falls back to a TikToken extractor
        # and dies on a file that is not one.
        "sentencepiece==0.2.0",
        "protobuf==5.29.1"
    )
)

# basilisk >= 1.22 builds the environment on a Python from pyenv, which needs
# git on every OS and compiles Python from source on macOS and Linux. A uv
# standalone build needs neither (reticulate fetches uv itself), so hand one to
# basilisk through its documented BASILISK_CUSTOM_PYTHON_* override before the
# environment is first built. An override already set (e.g. by an admin) wins;
# if uv fails, basilisk falls back to pyenv.
.standalone_python <- function() {
    var <- paste0("BASILISK_CUSTOM_PYTHON_",
        gsub(".", "_", repdist_python, fixed = TRUE))
    built <- file.path(basilisk::getExternalDir(), "repdist",
        utils::packageVersion("repdist"), "repdist-embed")
    if (nzchar(Sys.getenv(var)) || basilisk::useSystemDir() ||
        file.exists(built)) {
        return(invisible())
    }
    py <- tryCatch(
        reticulate::uv_run_tool("python",
            c("-c", shQuote("import sys; print(sys._base_executable)")),
            python_version = repdist_python, stdout = TRUE),
        error = identity, warning = identity
    )
    why <- NULL
    if (inherits(py, "condition")) {
        why <- conditionMessage(py)
        py <- character()
    }
    py <- py[length(py)]
    if (length(py) && file.exists(py)) {
        do.call(Sys.setenv, stats::setNames(list(py), var))
    } else {
        # Otherwise the user sees only pyenv's "requires git", which hides the
        # real cause (on Linux, usually no curl or wget for uv's installer).
        message("Could not get a standalone Python through uv",
            if (length(why)) c(" (", why, ")"),
            "; falling back to basilisk's pyenv, which needs git.")
    }
    invisible()
}
