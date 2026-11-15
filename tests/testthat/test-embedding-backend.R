test_that("R dispatch preserves sequences and validates backend output", {
    s <- c(p1 = "ACDE", p2 = "KLMN")
    result <- diag(2)
    testthat::local_mocked_bindings(.standalone_python = function() NULL)
    testthat::local_mocked_bindings(
        basiliskStart = function(...) NULL,
        basiliskStop = function(...) NULL,
        basiliskRun = function(proc, fun, seqs, spec, weights, device, batch) {
            expect_identical(seqs, as.list(unname(s)))
            expect_identical(device, "cpu")
            expect_identical(batch, 2L)
            expect_match(
                jsonlite::fromJSON(spec)$backbone_revision,
                "^[0-9a-f]{40}$"
            )
            result
        }, .package = "basilisk"
    )
    x <- embed_proteins(s, "esm2-8m", "cpu", 2, progress = FALSE)
    expect_s4_class(x, "AAStringSet")
    expect_identical(as.character(x), s)
    expect_identical(rownames(S4Vectors::mcols(x)$embedding), names(s))
    expect_equal(unname(S4Vectors::mcols(x)$embedding), diag(2))
    for (bad in list(matrix(0, 2, 2), matrix(NA_real_, 2, 2), diag(3))) {
        result <- bad
        expect_error(embed_proteins(s, "esm2-8m", "cpu", 2, progress = FALSE),
            "backend")
    }
})

test_that("chunked embedding keeps input order", {
    s <- stats::setNames(strrep("A", 1:250), paste0("p", 1:250))
    calls <- 0L
    testthat::local_mocked_bindings(.standalone_python = function() NULL)
    testthat::local_mocked_bindings(
        basiliskStart = function(...) NULL,
        basiliskStop = function(...) NULL,
        basiliskRun = function(proc, fun, seqs, ...) {
            calls <<- calls + 1L
            n <- nchar(unlist(seqs))
            cbind(cos(n), sin(n))
        }, .package = "basilisk"
    )
    x <- embed_proteins(s, "esm2-8m", "cpu", 1, progress = FALSE)
    expect_gt(calls, 1L)
    expect_equal(unname(S4Vectors::mcols(x)$embedding),
        cbind(cos(1:250), sin(1:250)))
})

test_that("all remote model artifacts have immutable revisions", {
    for (model in known_models()) {
        expect_match(model$backbone_revision, "^[0-9a-f]{40}$")
        if (!is.null(model$head_repo)) {
            expect_match(model$head_revision, "^[0-9a-f]{40}$")
        }
    }
})

test_that("a uv Python is handed to basilisk only before the env is built", {
    var <- paste0("BASILISK_CUSTOM_PYTHON_",
        gsub(".", "_", repdist_python, fixed = TRUE))
    keep <- c(var, "BASILISK_EXTERNAL_DIR", "BASILISK_USE_SYSTEM_DIR")
    old <- Sys.getenv(keep, unset = NA)
    on.exit({
        Sys.unsetenv(keep[is.na(old)])
        if (any(!is.na(old))) do.call(Sys.setenv, as.list(old[!is.na(old)]))
    })
    Sys.unsetenv(c(var, "BASILISK_USE_SYSTEM_DIR"))
    Sys.setenv(BASILISK_EXTERNAL_DIR = tempfile())
    py <- tempfile()
    file.create(py)
    calls <- 0
    testthat::local_mocked_bindings(uv_run_tool = function(...) {
        calls <<- calls + 1
        c("Downloading cpython", py)
    }, .package = "reticulate")

    .standalone_python()
    expect_identical(Sys.getenv(var), py)
    .standalone_python() # an override already set wins
    expect_identical(calls, 1)

    Sys.unsetenv(var)
    dir.create(file.path(basilisk::getExternalDir(), "repdist",
        utils::packageVersion("repdist"), "repdist-embed"), recursive = TRUE)
    .standalone_python() # a built env needs no Python
    expect_identical(calls, 1)
    expect_identical(Sys.getenv(var), "")
})
