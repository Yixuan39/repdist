.write_npz <- function(path, embeddings, ids = NULL, headers = NULL,
                       accession = NULL, singular = FALSE) {
  proc <- basilisk::basiliskStart(repdist:::.repdist_io_env)
  on.exit(basilisk::basiliskStop(proc), add = TRUE)
  basilisk::basiliskRun(proc, function(path, embeddings, ids, headers,
                                      accession, singular) {
    np <- reticulate::import("numpy", convert = FALSE)
    key <- if (singular) "embedding" else "embeddings"
    values <- list(embeddings)
    names(values) <- key
    if (!is.null(ids)) {
      values$ids <- np$asarray(ids, dtype = "str")
    } else if (!is.null(headers)) {
      values$headers <- np$asarray(headers, dtype = "str")
    } else if (!is.null(accession)) {
      values$accession <- np$asarray(accession, dtype = "str")
    }
    do.call(np$savez, c(list(path), values))
  }, path, embeddings, ids, headers, accession, singular)
}

test_that("read_embeddings imports Python and official TM-Vec NPZ layouts", {
  Z <- matrix(c(1, 0, 0, 1, 1, 1), 3, 2, byrow = TRUE)
  ids <- c("p1", "p2", "p3")

  python_file <- tempfile(fileext = ".npz")
  .write_npz(python_file, Z, ids = ids)
  expect_identical(read_embeddings(python_file),
                   structure(Z, dimnames = list(ids, c("d1", "d2"))))

  tmvec_file <- tempfile(fileext = ".npz")
  .write_npz(tmvec_file, Z, headers = ids)
  expect_identical(read_embeddings(tmvec_file),
                   structure(Z, dimnames = list(ids, c("d1", "d2"))))

  puf_file <- tempfile(fileext = ".npz")
  .write_npz(puf_file, Z, accession = ids, singular = TRUE)
  expect_identical(read_embeddings(puf_file),
                   structure(Z, dimnames = list(ids, c("d1", "d2"))))
})

test_that("read_embeddings requires identifiers", {
  path <- tempfile(fileext = ".npz")
  .write_npz(path, matrix(1, 1, 2))
  expect_error(read_embeddings(path), "ids.*headers.*accession")
})
