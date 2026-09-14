# Scoring behind FunFHMMer's merge test, ported from the reference
# implementation at UCLOrengoGroup/eMMA-FunFHMMER: scorecons conservation and
# its DOPS summary (Valdar 2002), and GroupSim specificity-determining
# position scores (Capra & Singh 2008), which that pipeline shells out to
# `scorecons` and `group_sim_sdp.py` for.

.repdist_aa <- c("A", "R", "N", "D", "C", "Q", "E", "G", "H", "I",
                 "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")

# PET91 (Jones, Taylor & Thornton 1991) as scorecons uses it: rescaled so
# every maximum is on the diagonal at 10. Rows and columns in .repdist_aa
# order. Normalised to [0, 1] on first use, which puts identity at 1.
.repdist_pet91 <- local({
  rows <- c(
    "10,-1,0,-1,-1,-1,-1,1,-2,0,-1,-1,-1,-3,1,1,2,-4,-3,1",
    "-1,10,0,-1,-1,2,0,0,2,-3,-3,4,-2,-4,-1,-1,-1,0,-2,-3",
    "0,0,10,2,-1,0,1,0,1,-2,-3,1,-2,-3,-1,1,1,-4,-1,-2",
    "-1,-1,2,10,-3,0,4,1,0,-3,-4,0,-3,-5,-2,0,-1,-5,-2,-3",
    "-1,-1,-1,-3,10,-3,-4,-1,0,-2,-3,-3,-2,0,-2,1,-1,1,2,-2",
    "-1,2,0,0,-3,10,2,-1,3,-3,-2,2,-2,-4,0,-1,-1,-3,-1,-3",
    "-1,0,1,4,-4,2,10,1,0,-3,-4,1,-3,-5,-2,-1,-1,-5,-4,-2",
    "1,0,0,1,-1,-1,1,10,-2,-3,-4,-1,-3,-5,-1,1,0,-2,-4,-2",
    "-2,2,1,0,0,3,0,-2,10,-3,-2,1,-2,0,0,-1,-1,-3,4,-3",
    "0,-3,-2,-3,-2,-3,-3,-3,-3,10,2,-3,3,0,-2,-1,1,-4,-2,4",
    "-1,-3,-3,-4,-3,-2,-4,-4,-2,2,10,-3,3,2,0,-2,-1,-2,-1,2",
    "-1,4,1,0,-3,2,1,-1,1,-3,-3,10,-2,-5,-2,-1,-1,-3,-3,-3",
    "-1,-2,-2,-3,-2,-2,-3,-3,-2,3,3,-2,10,0,-2,-1,0,-3,-3,2",
    "-3,-4,-3,-5,0,-4,-5,-5,0,0,2,-5,0,10,-2,-2,-2,-1,5,0",
    "1,-1,-1,-2,-2,0,-2,-1,0,-2,0,-2,-2,-2,10,1,1,-5,-3,-1",
    "1,-1,1,0,1,-1,-1,1,-1,-1,-2,-1,-1,-2,1,10,1,-3,-1,-1",
    "2,-1,1,-1,-1,-1,-1,0,-1,1,-1,-1,0,-2,1,1,10,-4,-3,0",
    "-4,0,-4,-5,1,-3,-5,-2,-3,-4,-2,-3,-3,-1,-5,-3,-4,10,0,-4",
    "-3,-2,-1,-2,2,-1,-4,-4,4,-2,-1,-3,-3,5,-3,-1,-3,0,10,-3",
    "1,-3,-2,-3,-2,-3,-2,-2,-3,4,2,-3,2,0,-1,-1,0,-4,-3,10")
  M <- matrix(as.numeric(unlist(strsplit(rows, ","))), 20, 20, byrow = TRUE,
              dimnames = list(.repdist_aa, .repdist_aa))
  (M - min(M)) / (max(M) - min(M))
})

# BLOSUM62 background, the distribution group_sim_sdp.py scores columns
# against. .repdist_aa order.
.repdist_bg <- c(0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074,
                 0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
                 0.051, 0.013, 0.032, 0.073)

# Residue counts per alignment column, one row per standard amino acid.
# Everything else -- gaps, X, B, Z, U -- is simply absent from the counts.
.repdist_counts <- function(A) {
  vapply(.repdist_aa, function(a) colSums(A == a), numeric(ncol(A)))
}

# scorecons, verified against the reference binary's shipped test alignment:
# 109 of its 113 columns to the printed three decimals, the four gapped ones
# within 0.004 because the binary's exact gap term is not published. Each
# column is the mean normalised PET91 score over its pairs, with gapped pairs
# contributing nothing rather than being dropped from the denominator.
.repdist_scorecons <- function(A) {
  n <- nrow(A)
  if (n < 2L) return(rep(1, ncol(A)))
  C <- t(.repdist_counts(A))
  (colSums(C * (.repdist_pet91 %*% C)) - colSums(C)) / (n * (n - 1))
}

# DOPS: 0 when every column has the same conservation score, 100 when no two
# share one. Entropy over the distinct scores, normalised by log(n columns).
# All-gap columns are dropped, which is what makes the reference binary score
# its own test clusters 10.456 and 0.000.
.repdist_dops <- function(A) {
  A <- A[, colSums(A != "-") > 0L, drop = FALSE]
  scores <- round(.repdist_scorecons(A), 3)
  L <- length(scores)
  if (L < 2L) return(0)
  p <- table(scores) / L
  100 * -sum(p * log(p)) / log(L)
}

# Jensen-Shannon divergence of a column against the background, in bits and
# damped by the column's gap fraction -- group_sim_sdp.py's conservation term.
.repdist_jsd <- function(C, gaps) {
  p <- (C + 0.001) / rep(colSums(C) + 0.02, each = nrow(C))
  q <- matrix(.repdist_bg, nrow(p), ncol(p))
  r <- (p + q) / 2
  kl <- function(x) colSums(ifelse(x > 0, x * log2(x / r), 0))
  (1 - gaps) * (kl(p) + kl(q)) / 2
}

# GroupSim: how much more alike a column's residues are within each group than
# across them, so a high score marks a position that separates the two groups.
# Ported from group_sim_sdp.py with the identity matrix and FunFHMMer's own
# `-c 0.3 -g 0.5`; gaps count as a 21st symbol, as they do there.
.repdist_groupsim <- function(A, group, column_gap = 0.3, group_gap = 0.5,
                              window = 3L, lambda = 0.7) {
  L <- ncol(A)
  split_A <- lapply(sort(unique(group)), function(g) A[group == g, , drop = FALSE])
  gap_fraction <- function(M) colSums(M == "-") / nrow(M)

  # Identity matrix over 21 symbols: a pair scores 1 only when it matches, so
  # every mean collapses to a count of agreeing pairs over valid pairs.
  counts <- lapply(split_A, function(M)
    rbind(t(.repdist_counts(M)), colSums(M == "-")))
  sizes <- lapply(counts, colSums)
  within <- Reduce(`+`, Map(function(C, m) {
    pairs <- m * (m - 1) / 2
    ifelse(pairs > 0, (colSums(C * (C - 1)) / 2) / pmax(pairs, 1), 0)
  }, counts, sizes)) / length(counts)
  between <- colSums(counts[[1L]] * counts[[2L]]) / (sizes[[1L]] * sizes[[2L]])

  sdp <- (1 - gap_fraction(A)) * (within - between)
  sdp[gap_fraction(A) > column_gap] <- NA_real_
  for (M in split_A) sdp[gap_fraction(M) > group_gap] <- NA_real_
  if (all(is.na(sdp))) return(sdp)

  # ConsWin: pull each score towards how conserved its neighbourhood is.
  cons <- ifelse(is.na(sdp), NA_real_, .repdist_jsd(t(.repdist_counts(A)),
                                                    gap_fraction(A)))
  neighbours <- vapply(seq_len(L), function(i) {
    j <- setdiff(seq(max(i - window, 1L), min(i + window, L)), i)
    j <- j[!is.na(cons[j]) & cons[j] >= 0]
    if (length(j)) mean(cons[j]) else 0
  }, numeric(1))
  neighbours[is.na(sdp)] <- NA_real_
  span <- diff(range(neighbours, na.rm = TRUE))
  neighbours <- if (span == 0) ifelse(is.na(neighbours), NA_real_, 1) else
    (neighbours - min(neighbours, na.rm = TRUE)) / span
  ifelse(is.na(sdp), NA_real_, (1 - lambda) * neighbours + lambda * sdp)
}

# The merge decision. FunFHMMer's own rule is a cascade of hardcoded score
# bands, percentage cutoffs and quartile tests with no published derivation,
# two of whose branches cannot fire; this replaces it with the one question
# that cascade is circling -- does the joint alignment carry enough
# specificity-determining signal to call these two functions apart -- asked
# with four parameters that can each be defended and swept.
#
# A pair that cannot be aligned well enough to score is treated as apart:
# unalignability is itself evidence of difference, and is what the reference's
# unscorable-column branches amount to.
.repdist_gs_separable <- function(scores, dops, min_dops = 70,
                                  sdp_score = 0.7, sdp_fraction = 0.2,
                                  max_unscorable = 0.3) {
  if (mean(is.na(scores)) > max_unscorable) return(TRUE)
  scored <- scores[!is.na(scores)]
  max(dops) >= min_dops && mean(scored >= sdp_score) >= sdp_fraction
}
