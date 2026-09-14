# Reference values from the FunFHMMer implementation at
# UCLOrengoGroup/eMMA-FunFHMMER: its shipped scorecons and GroupSim outputs
# for the example alignments in t/example_data. These pin the port to the
# published tool's own numbers.

joint <- c(
  "a537db49c8089da45f7c7a64b27f2727/124-233|cluster1" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTYQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "dcc33ac588339221129a304e85ab0769/124-235|cluster1" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "597d8a166408cb27f323762efcb3dd03/124-233|cluster1" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "c54fe0642f3603eaafc4c25bbd03250b/124-233|cluster1" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "a1900ac5253db25e6b59fd30f7559ddd/124-233|cluster1" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "66c9ec71c53d711731af571200f5a5c2/124-235|cluster2" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "296260caf498312cc31f3486e5b73fe7/124-235|cluster2" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "096c52f420044b6ff34a849271c35598/124-235|cluster2" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN")
cluster1 <- c(
  "a537db49c8089da45f7c7a64b27f2727/124-233" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTYQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "dcc33ac588339221129a304e85ab0769/124-235" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "597d8a166408cb27f323762efcb3dd03/124-233" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "c54fe0642f3603eaafc4c25bbd03250b/124-233" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "a1900ac5253db25e6b59fd30f7559ddd/124-233" =
    "LVKQEFRTKVEETAKQKAEEVLLDILLPFPGENKHGSTHQITG---SSQFTEEEDRKTHFLETREFMRKKLKAGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN")
cluster2 <- c(
  "66c9ec71c53d711731af571200f5a5c2/124-235" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "296260caf498312cc31f3486e5b73fe7/124-235" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN",
  "096c52f420044b6ff34a849271c35598/124-235" =
    "LVKQEFRTKVEETAKQKAEEALLDILLPFPGENKHGSG-QITGFATSSTLADEEDRKTHFLETREFMRKKLKTGKLDDQEVELDLPNPSVSQVPMLQVFGAGNLDDLDNQLQN")

# scorecons' per-column scores for the 8-sequence alignment, as its binary
# prints them.
reference_scorecons <- c(1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.657, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.619, 0.174, 1.000, 1.000, 1.000, 1.000, 0.210, 0.210, 0.210, 1.000, 1.000, 0.581, 0.695, 0.695, 0.771, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 0.695, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000)

chars <- function(x) do.call(rbind, strsplit(unname(x), ""))

# The reference summarises GroupSim as the share of columns in each score
# band plus the quartiles. repdist's own rule does not need the bands, so this
# lives here, where it exists only to compare against the published numbers.
bands <- function(scores) {
  L <- length(scores)
  share <- function(ok) 100 * sum(ok, na.rm = TRUE) / L
  filled <- ifelse(is.na(scores), 1.1, scores)
  c(le3 = share(scores <= 0.3),
    le4 = share(scores > 0.3 & scores <= 0.4),
    mid = share(scores > 0.4 & scores < 0.7),
    le8 = share(scores >= 0.7 & scores < 0.8),
    high = share(scores >= 0.8 & scores <= 1),
    none = 100 * sum(is.na(scores)) / L,
    total = L,
    q1 = stats::quantile(filled, 0.25, names = FALSE),
    median = stats::quantile(filled, 0.5, names = FALSE),
    q3 = stats::quantile(filled, 0.75, names = FALSE))
}

test_that("scorecons reproduces the reference per-column scores", {
  scores <- round(repdist:::.repdist_scorecons(chars(joint)), 3)
  expect_equal(sum(scores == reference_scorecons), 109L)
  # the four that differ are gapped columns, where the binary's gap term is
  # not published
  expect_lt(max(abs(scores - reference_scorecons)), 0.005)
})

test_that("DOPS reproduces the reference scores exactly", {
  expect_equal(round(repdist:::.repdist_dops(chars(cluster1)), 3), 10.456)
  expect_equal(round(repdist:::.repdist_dops(chars(cluster2)), 3), 0)
  expect_equal(round(repdist:::.repdist_dops(chars(joint)), 3), 10.456)
})

test_that("GroupSim reproduces the reference band counts and quartiles", {
  group <- sub(".*[|]", "", names(joint))
  b <- bands(repdist:::.repdist_groupsim(chars(joint), group))

  # the reference scored 112 of the 113 columns, so compare counts, not
  # percentages: 101 at or below 0.3, 7 between 0.4 and 0.7, 4 unscorable
  counts <- round(b[1:6] / 100 * b[["total"]])
  expect_equal(unname(counts), c(101, 1, 7, 0, 0, 4))
  expect_lt(max(abs(unname(b[c("q1", "median", "q3")]) -
                     c(0.1728445, 0.2113905, 0.24123))), 0.002)
})

test_that("a column split cleanly between the groups reads as separable", {
  # every column identical within each group and different between them: the
  # signal FunFHMMer refuses to merge across
  split <- c(a1 = "AAAAAAAAAA", a2 = "AAAAAAAAAA", a3 = "AAAAAAAAAA",
             b1 = "KKKKKKKKKK", b2 = "KKKKKKKKKK", b3 = "KKKKKKKKKK")
  group <- rep(c("a", "b"), each = 3)
  scores <- repdist:::.repdist_groupsim(chars(split), group)
  expect_true(all(scores > 0.9))
  expect_true(repdist:::.repdist_gs_separable(scores, c(100, 100)))

  # one group, one residue everywhere: nothing distinguishes them
  same <- c(a1 = "AAAAAAAAAA", a2 = "AAAAAAAAAA",
            b1 = "AAAAAAAAAA", b2 = "AAAAAAAAAA")
  scores <- repdist:::.repdist_groupsim(chars(same), rep(c("a", "b"), each = 2))
  expect_true(all(scores < 0.5))
  expect_false(repdist:::.repdist_gs_separable(scores, c(0, 0)))
})

test_that("funfhmmer with seqs runs the GroupSim test end to end", {
  skip_if_not_installed("DECIPHER")
  seqs <- gsub("-", "", joint)
  names(seqs) <- sub("[|].*", "", names(joint))

  set.seed(3)
  Z <- matrix(rnorm(length(seqs) * 8), length(seqs), 8)
  Z <- Z / sqrt(rowSums(Z^2))
  rownames(Z) <- names(seqs)
  D <- repdist_matrix(Z)

  out <- funfhmmer(D, seqs = seqs)
  expect_equal(out$test, "groupsim")
  expect_setequal(unlist(out$clusters), names(seqs))
  expect_equal(sum(lengths(out$clusters)), length(seqs))
  expect_named(out$representatives, names(out$clusters))
})

test_that("seqs that do not cover every protein are rejected", {
  skip_if_not_installed("DECIPHER")
  seqs <- gsub("-", "", joint)
  names(seqs) <- sub("[|].*", "", names(joint))
  set.seed(3)
  Z <- matrix(rnorm(length(seqs) * 8), length(seqs), 8)
  Z <- Z / sqrt(rowSums(Z^2))
  rownames(Z) <- names(seqs)
  expect_error(funfhmmer(repdist_matrix(Z), seqs = seqs[-1]), "no sequence")
})

test_that("the reference cluster pair merges, as it does upstream", {
  # cluster1 and cluster2 are two starting clusters of one superfamily, which
  # the reference pipeline merges; repdist's rule has to agree
  group <- sub(".*[|]", "", names(joint))
  expect_false(repdist:::.repdist_gs_separable(
    repdist:::.repdist_groupsim(chars(joint), group),
    c(repdist:::.repdist_dops(chars(cluster1)),
      repdist:::.repdist_dops(chars(cluster2)))))
})

# Two families sharing a backbone but fixing different residues at a quarter
# of the columns: the specificity-determining signal FunFHMMer looks for.
# `rate` makes each column mutate at its own pace, which is what gives the
# alignments the positional diversity DOPS measures.
planted_families <- function(seed, n = 12, L = 140, marked = 4) {
  set.seed(seed)
  aa <- repdist:::.repdist_aa
  backbone <- sample(aa, L, replace = TRUE)
  rate <- runif(L, 0, 0.9)
  sdp <- seq(marked, L, by = marked)
  family <- function(marker, tag) {
    s <- vapply(seq_len(n), function(i) {
      x <- backbone
      hit <- runif(L) < rate
      x[hit] <- sample(aa, sum(hit), replace = TRUE)
      x[sdp] <- marker
      paste(x, collapse = "")
    }, character(1))
    stats::setNames(s, paste0(tag, seq_len(n)))
  }
  c(family("K", "a"), family("D", "b"))
}

test_that("informative families with their own marker residues stay apart", {
  skip_if_not_installed("DECIPHER")
  seqs <- Biostrings::AAStringSet(planted_families(5))
  half <- length(seqs) / 2
  A1 <- DECIPHER::AlignSeqs(seqs[seq_len(half)], verbose = FALSE)
  A2 <- DECIPHER::AlignSeqs(seqs[-seq_len(half)], verbose = FALSE)
  dops <- c(repdist:::.repdist_dops(as.matrix(A1)),
            repdist:::.repdist_dops(as.matrix(A2)))
  expect_gte(max(dops), 70)            # the alignments are informative

  scores <- repdist:::.repdist_groupsim(
    as.matrix(DECIPHER::AlignSeqs(seqs, verbose = FALSE)),
    rep(c("a", "b"), each = half))
  scored <- scores[!is.na(scores)]
  expect_gt(100 * mean(scored >= 0.7), 20)   # enough marker columns
  expect_true(repdist:::.repdist_gs_separable(scores, dops))
})

test_that("funfhmmer splits those families and keeps each one whole", {
  skip_if_not_installed("DECIPHER")
  seqs <- Biostrings::AAStringSet(planted_families(5))
  aligned <- DECIPHER::AlignSeqs(seqs, verbose = FALSE)
  D <- stats::as.dist(DECIPHER::DistanceMatrix(
    aligned, verbose = FALSE, includeTerminalGaps = TRUE))

  out <- funfhmmer(D, seqs = seqs, auto_merge = 0)
  membership <- stats::setNames(
    rep(names(out$clusters), lengths(out$clusters)), unlist(out$clusters))
  a <- unique(membership[grep("^a", names(seqs), value = TRUE)])
  b <- unique(membership[grep("^b", names(seqs), value = TRUE)])
  expect_length(a, 1L)
  expect_length(b, 1L)
  expect_false(identical(a, b))
})

test_that("an uninformative alignment cannot block a merge", {
  # same marker columns, but every sequence a near-copy: DOPS stays under 70
  # and FunFHMMer declines to read a functional signal off it
  skip_if_not_installed("DECIPHER")
  aa <- repdist:::.repdist_aa
  set.seed(9)
  backbone <- sample(aa, 140, replace = TRUE)
  one <- function(marker, tag, n = 8) {
    s <- vapply(seq_len(n), function(i) {
      x <- backbone
      x[sample(140, 3)] <- sample(aa, 3, replace = TRUE)
      x[seq(4, 140, by = 4)] <- marker
      paste(x, collapse = "")
    }, character(1))
    stats::setNames(s, paste0(tag, seq_len(n)))
  }
  seqs <- Biostrings::AAStringSet(c(one("K", "a"), one("D", "b")))
  A1 <- DECIPHER::AlignSeqs(seqs[1:8], verbose = FALSE)
  A2 <- DECIPHER::AlignSeqs(seqs[9:16], verbose = FALSE)
  dops <- c(repdist:::.repdist_dops(as.matrix(A1)),
            repdist:::.repdist_dops(as.matrix(A2)))
  expect_lt(max(dops), 70)

  scores <- repdist:::.repdist_groupsim(
    as.matrix(DECIPHER::AlignSeqs(seqs, verbose = FALSE)),
    rep(c("a", "b"), each = 8))
  expect_false(repdist:::.repdist_gs_separable(scores, dops))
})
