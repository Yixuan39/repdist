# read_fasta only. The model paths in embed_proteins need network + multi-GB
# downloads, so they are exercised by the analyses rather than here.

fasta_file <- function(lines) {
  p <- tempfile(fileext = ".fasta")
  writeLines(lines, p)
  p
}

test_that("read_fasta joins multi-line records and keeps header order", {
  p <- fasta_file(c(">p1 some description", "ACDE", "FGHI",
                    ">p2", "KLMN",
                    ">p3 another", "PQRS", "TVWY"))
  s <- read_fasta(p)
  expect_identical(names(s), c("p1", "p2", "p3"))
  expect_identical(unname(s), c("ACDEFGHI", "KLMN", "PQRSTVWY"))
})

test_that("record order survives past ten records", {
  # split() on a bare integer group sorts "10" before "2", which would pair
  # sequences with the wrong headers -- silently, and only for >9 records
  n <- 12
  p <- fasta_file(as.vector(rbind(paste0(">p", seq_len(n)),
                                  strrep("ACDEFGHIKL", seq_len(n)))))
  s <- read_fasta(p)
  expect_identical(names(s), paste0("p", seq_len(n)))
  expect_identical(nchar(s, type = "chars"),
                   stats::setNames(10L * seq_len(n), paste0("p", seq_len(n))))
})

test_that("read_fasta uppercases, removes terminal stops, and ignores blank lines", {
  p <- fasta_file(c(">p1", "acde", "", "fghi*", ">p2", "klmn"))
  s <- read_fasta(p)
  expect_identical(unname(s), c("ACDEFGHI", "KLMN"))
})

test_that("read_fasta performs embedding-compatible sequence QC", {
  expect_error(read_fasta(fasta_file(c(">p1", "AC DE"))), "invalid one-letter")
  expect_error(read_fasta(fasta_file(c(">p1", "ACJDE"))), "model-incompatible")
  expect_error(read_fasta(fasta_file(c(">p1", "AC*DE"))), "model-incompatible")
  expect_error(read_fasta(fasta_file(c(">p1", ">p2", "ACDE"))), "empty")
  expect_error(read_fasta(fasta_file(c(">p1", "ACDE", ">p1", "FGHI"))), "unique")
})

test_that("read_fasta errors on a file with no header", {
  expect_error(read_fasta(fasta_file(c("ACDE", "FGHI"))), "expected")
})

test_that("embed_proteins accepts a FASTA path in place of a named vector", {
  p <- fasta_file(c(">p1", "ACDE", ">p2", "KLMN"))
  # stop before any model is fetched: a duplicate name must be rejected after
  # the file is read, which proves the path was expanded into sequences
  p_dup <- fasta_file(c(">p1", "ACDE", ">p1", "KLMN"))
  expect_error(embed_proteins(p_dup), "unique")
  expect_identical(names(read_fasta(p)), c("p1", "p2"))
})

test_that("every registered model names a head embed_proteins can build", {
  # a typo'd repo key falls through to the generic path rather than erroring,
  # so a malformed registry would only surface as a wrong embedding later
  reg <- known_models()
  expect_gt(length(reg), 0)
  heads <- vapply(reg, function(m) m$head, character(1))
  expect_true(all(heads %in% c("generic", "tmvec1", "tmvec1-large")))
  expect_true(all(nzchar(vapply(reg, function(m) m$backbone, character(1)))))
  expect_true("tmvec-swissmodel-large" %in% names(reg))  # the default model
  expect_identical(reg$esm2$backbone, "facebook/esm2_t33_650M_UR50D")
  expect_identical(reg$`esm2-small`$backbone, "facebook/esm2_t6_8M_UR50D")

  # a repeated top-level key in models.json parses as two same-named entries,
  # the second silently shadowing the first
  expect_false(anyDuplicated(names(reg)) > 0)

  # a tmvec1-large head has no HuggingFace release, so its entry has to carry
  # where to fetch it, what it should hash to, and -- since a bare checkpoint
  # ships no config.json -- the architecture to rebuild it with. Otherwise the
  # omission only surfaces mid-embedding, after a multi-hundred-MB download.
  for (m in reg[heads == "tmvec1-large"]) {
    expect_match(m$url, "^https://")
    expect_match(m$sha256, "^[0-9a-f]{64}$")
    expect_setequal(names(m$config),
                    c("d_model", "nhead", "num_layers", "dim_feedforward",
                      "out_dim", "dropout", "activation"))
  }
  # a tmvec1 head must NOT carry one: it takes config.json from its own repo,
  # and a stale copy here would silently win over the published architecture
  for (m in reg[heads == "tmvec1"]) expect_null(m$config)
})
