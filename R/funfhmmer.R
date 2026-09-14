# --- distance merge test (this package's own, unpublished) -----------------

# A live cluster is its member indices plus the sum of the distances within
# it, carried along so the merge test never rescans a cluster it has already
# paid for -- without it the traversal is cubic in the number of proteins.
.funfhmmer_spread <- function(cl) {
  k <- length(cl$i)
  if (k < 2L) 0 else cl$w / (k * (k - 1) / 2)
}

.funfhmmer_distance_test <- function(D, separation) {
  between <- function(a, b)
    sum(usedist::dist_get(D, rep(a$i, each = length(b$i)),
                          rep(b$i, times = length(a$i))))
  list(
    leaf = function(i) list(i = i, w = 0),
    join = function(clusters) {
      i <- unlist(lapply(clusters, `[[`, "i"))
      list(i = i, w = if (length(clusters) == 2L)
             clusters[[1L]]$w + clusters[[2L]]$w +
               between(clusters[[1L]], clusters[[2L]])
           else sum(usedist::dist_subset(D, i)))
    },
    apart = function(a, b, merged) {
      spread <- max(.funfhmmer_spread(a), .funfhmmer_spread(b))
      spread > 0 && (merged$w - a$w - b$w) /
        (length(a$i) * length(b$i)) > separation * spread
    },
    # A cluster with no internal spread is this test's uninformative
    # alignment: it can neither block a merge nor vouch for a three-way one.
    informative = function(cl) .funfhmmer_spread(cl) > 0)
}

# --- GroupSim merge test (the published one) -------------------------------

.funfhmmer_seqs <- function(seqs, labels) {
  if (!requireNamespace("DECIPHER", quietly = TRUE))
    stop("`seqs` needs the DECIPHER package for the alignments.", call. = FALSE)
  seqs <- Biostrings::AAStringSet(seqs)
  missing <- setdiff(labels, names(seqs))
  if (length(missing))
    stop(length(missing), " of the proteins in `D` have no sequence in ",
         "`seqs`, starting with ", missing[[1L]], ".", call. = FALSE)
  seqs[labels]
}

# Every merge realigns its members from scratch, as the reference does by
# concatenating the two clusters and handing them back to MAFFT. Profile
# alignment would be cheaper but freezes each side's columns, and the whole
# point of the test is where the columns fall once the two are read together.
# The size tiers stand in for the reference's own MAFFT tiers, which drop
# accuracy as clusters grow because this runs at every node.
.funfhmmer_align <- function(seqs) {
  if (length(seqs) < 2L) return(seqs)
  accuracy <- if (length(seqs) <= 500L) list(iterations = 2L, refinements = 2L)
              else if (length(seqs) <= 2000L) list(iterations = 1L, refinements = 0L)
              else list(iterations = 0L, refinements = 0L)
  do.call(DECIPHER::AlignSeqs,
          c(list(seqs, verbose = FALSE, processors = NULL), accuracy))
}

.funfhmmer_groupsim_test <- function(seqs, rule) {
  cluster <- function(i) {
    aln <- .funfhmmer_align(seqs[i])
    list(i = i, aln = aln, dops = .repdist_dops(as.matrix(aln)))
  }
  list(
    leaf = function(i) cluster(i),
    join = function(clusters) cluster(unlist(lapply(clusters, `[[`, "i"))),
    apart = function(a, b, merged) {
      group <- ifelse(names(merged$aln) %in% names(seqs)[a$i], "a", "b")
      do.call(.repdist_gs_separable,
              c(list(.repdist_groupsim(as.matrix(merged$aln), group),
                     c(a$dops, b$dops)), rule))
    },
    informative = function(cl) cl$dops >= rule$min_dops)
}

# --- the traversal ---------------------------------------------------------

#' Cut a tree into functional families with FunFHMMer
#'
#' Walks a protein tree bottom-up in merge order, the way FunFHMMer traverses
#' a GeMMA trace (Das et al. 2015), and refuses a merge once the two clusters
#' meeting at a node are distinguishable. Every node below a cut stays a
#' family of its own, so each branch is cut where its own members stop
#' resembling each other rather than at one height for the whole tree -- the
#' difference from [bin_proteins()], which cuts everywhere at `1 - min_sim`.
#'
#' @section Merge tests:
#' With `seqs`, the published test runs: the two clusters' members are
#' realigned together from scratch with DECIPHER, as the reference hands the
#' concatenated pair back to MAFFT at every node, GroupSim scores every column
#' of that alignment for how much better it separates the two clusters than it
#' varies within them, and the merge is blocked when enough of those columns
#' score high enough -- specificity-determining positions, the signal that two
#' functions are about to be pooled.
#'
#' The scoring is ported from the reference implementation and pinned to its
#' own published numbers. The decision it feeds is not: FunFHMMer settles a
#' merge with a cascade of hardcoded score bands, percentage cutoffs and
#' quartile tests that appear in no paper and two of whose branches cannot
#' fire. In its place, a merge is blocked when `sdp_fraction` of the scorable
#' columns reach `sdp_score` and at least one side's alignment is informative
#' (DOPS `min_dops`), or when more than `max_unscorable` of the columns cannot
#' be scored at all -- two clusters that will not align are already apart.
#' Sweep the four rather than trusting them: they are anchors from the
#' literature, not fitted values.
#'
#' Without `seqs` a distance test runs instead: a merge is blocked when the
#' mean distance between the two clusters exceeds `separation` times the
#' larger of their mean within-cluster distances, and a cluster with no
#' internal spread cannot block anything, which is the role the DOPS gate
#' plays above. It needs no alignments and runs in seconds where the published
#' test takes minutes.
#'
#' This test is this package's own and has no published counterpart. It is in
#' particular not what CATH-eMMA does: eMMA replaces the tree building --
#' MMseqs2 for CD-HIT, embedding or Foldseek distances for the all-against-all
#' HMM comparison -- and leaves FunFHMMer itself alone but for deleting the
#' E-value rail, so the published pipeline still realigns and scores every
#' node with GroupSim. Cite it as a substitution of this package, not as a
#' published method.
#'
#' Either way the test only adjudicates a band: `auto_merge` takes the lowest
#' merges outright, `min_sim` refuses the highest outright, and only what
#' falls between is scored. So the result refines `bin_proteins(D, min_sim)` --
#' same rail, with families split further wherever the test finds a reason.
#'
#' @param D Protein `dist` from [repdist_matrix()], or a square distance
#'   matrix on the same `1 - similarity` scale.
#' @param min_sim Similarity floor in `(0, 1]`, on the same scale as
#'   [bin_proteins()]'s. No merge is taken above `1 - min_sim`, whatever the
#'   merge test says, so every family keeps the all-pairs similarity guarantee
#'   and the result is a refinement of `bin_proteins(D, min_sim)`.
#'
#'   FunFHMMer always rails -- it refuses any merge whose E-value exceeds a
#'   fixed threshold before scoring anything, so the test only adjudicates a
#'   band -- and you should set this. It defaults to `NULL`, no rail, because
#'   an HHsearch E-value has no scale-free translation and the right floor
#'   depends entirely on what `D` measures: 0.5 is the conventional fold cut
#'   for TM-score, but is far too strict for a sequence-identity distance over
#'   a diverse family, where it fragments families from the inside. Read it
#'   off your own within-family distances rather than taking a number from
#'   here.
#' @param seqs Named character vector or `AAStringSet` of unaligned protein
#'   sequences covering every label in `D`, which switches on the published
#'   GroupSim merge test and needs the DECIPHER package. Every node is
#'   realigned from scratch, which is where nearly all the time goes: about
#'   20 seconds for 160 proteins against hundredths of a second for the
#'   distance test, and hours rather than minutes once there are thousands.
#' @param separation Positive multiplier on within-cluster spread (default 2),
#'   used by the distance test only. Higher values merge more, giving fewer
#'   and larger families. A uniformly spread cluster produces a ratio near 1
#'   at every one of its own merges, so 1 is the knife edge at which nothing
#'   grows past a pair; useful values sit above it. Tune it against known
#'   families rather than trusting the default.
#' @param min_dops Diversity of positions score, 0 to 100, at which a
#'   cluster's alignment is informative enough to read a functional signal off
#'   (default 70, CATH's own bar for an informative FunFam). Below it on both
#'   sides the merge is taken. GroupSim only.
#' @param sdp_score Column score, 0 to 1, at which a column counts as
#'   specificity-determining (default 0.7, the conventional cut for a highly
#'   conserved position and the band the reference reports against). GroupSim
#'   only.
#' @param sdp_fraction Share of the scorable columns that must reach
#'   `sdp_score` before a merge is blocked (default 0.2, the one live
#'   threshold in the reference's cascade). GroupSim only.
#' @param max_unscorable Share of columns GroupSim may leave unscored --
#'   too gappy in one group or overall -- before the pair is held apart on
#'   that alone (default 0.3). GroupSim only.
#' @param auto_merge Fraction of the merges, lowest first, taken without
#'   testing (default 0.1), as FunFHMMer merges the low-E-value head of a
#'   trace outright. 0 tests every merge, 1 tests none.
#' @param tree A [stats::hclust] tree over the same proteins, e.g. the one
#'   [bin_proteins()] returns. Built with complete linkage when `NULL`.
#' @return A list shaped like [bin_proteins()]'s, with families named
#'   `funfam1`, `funfam2`, ... from largest to smallest:
#'   \describe{
#'     \item{`clusters`}{Named list, one character vector of member
#'       accessions per family.}
#'     \item{`similarity`}{Family-by-family similarity between the
#'       representatives.}
#'     \item{`representatives`}{Family-named medoid accessions.}
#'     \item{`tree`}{The tree that was traversed.}
#'     \item{`separation`, `auto_merge`}{The settings used.}
#'     \item{`test`}{Which merge test ran, `"groupsim"` or `"distance"`.}
#'   }
#' @references Das S, Lee D, Sillitoe I, Dawson NL, Lees JG, Orengo CA (2015).
#'   Functional classification of CATH superfamilies. *Bioinformatics*
#'   31:3460-7. Capra JA, Singh M (2008). Characterization and prediction of
#'   residues determining protein functional specificity. *Bioinformatics*
#'   24:1473-80. Valdar WSJ (2002). Scoring residue conservation. *Proteins*
#'   48:227-41. Bordin N, Scholes H, Rauer C, Roca-Martinez J, Sillitoe I,
#'   Orengo C (2024). Clustering protein functional families at large scale
#'   with hierarchical approaches. *Protein Science* 33:e5140.
#' @seealso [bin_proteins()] for a single-threshold cut of the same tree.
#' @examples
#' Z <- rbind(p1 = c(1, 0), p2 = c(0.95, 0.05), p3 = c(0, 1))
#' Z <- Z / sqrt(rowSums(Z^2))
#' funfhmmer(repdist_matrix(Z))
#' @export
funfhmmer <- function(D, seqs = NULL, min_sim = NULL, separation = 2,
                      auto_merge = 0.1,
                      min_dops = 70, sdp_score = 0.7, sdp_fraction = 0.2,
                      max_unscorable = 0.3, tree = NULL) {
  D <- .repdist_check_dist(D)
  labels <- attr(D, "Labels")
  n <- attr(D, "Size")
  stopifnot(is.numeric(separation), length(separation) == 1L,
            is.finite(separation), separation > 0)
  stopifnot(is.numeric(auto_merge), length(auto_merge) == 1L,
            is.finite(auto_merge), auto_merge >= 0, auto_merge <= 1)
  if (!is.null(min_sim))
    stopifnot(is.numeric(min_sim), length(min_sim) == 1L,
              is.finite(min_sim), min_sim > 0, min_sim <= 1)
  if (is.null(tree)) tree <- fastcluster::hclust(D, method = "complete")
  if (!inherits(tree, "hclust") || length(tree$order) != n)
    stop("`tree` must be an hclust over the same ", n, " proteins as `D`.",
         call. = FALSE)

  rule <- list(min_dops = min_dops, sdp_score = sdp_score,
               sdp_fraction = sdp_fraction, max_unscorable = max_unscorable)
  stopifnot(vapply(rule, function(x)
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 0, logical(1)),
    min_dops <= 100, sdp_score <= 1, sdp_fraction <= 1, max_unscorable <= 1)
  test <- if (is.null(seqs)) .funfhmmer_distance_test(D, separation)
          else .funfhmmer_groupsim_test(.funfhmmer_seqs(seqs, labels), rule)

  # hclust heights are already ascending, so the lowest merges are the first
  # rows and this count is the whole of FunFHMMer's auto-merge phase; the
  # original also needs a flag because GeMMA traces are not sorted.
  auto_merged <- auto_merge * (n - 1L)
  # ...and the rail at the other end: FunFHMMer refuses any merge above a
  # fixed significance threshold before it scores anything, so the test only
  # ever adjudicates a band. Heights ascend, so every node past this one is
  # refused too and nothing above it is ever aligned.
  ceiling <- if (is.null(min_sim)) Inf else 1 - min_sim

  # live[[k]]: the clusters node k still holds. One element means the node
  # merged cleanly and has an alignment in FunFHMMer's terms; more than one
  # means a cut happened below it and the parts are candidate families.
  live <- vector("list", n - 1L)
  blocked <- logical(n - 1L)
  members <- function(i) if (i < 0L) list(test$leaf(-i)) else live[[i]]

  for (k in seq_len(n - 1L)) {
    a <- tree$merge[k, 1L]
    b <- tree$merge[k, 2L]
    A <- members(a)
    B <- members(b)
    # Every node is consumed exactly once, so its alignment can go now.
    if (a > 0L) live[a] <- list(NULL)
    if (b > 0L) live[b] <- list(NULL)

    if ((a > 0L && blocked[[a]]) || (b > 0L && blocked[[b]]) ||
        (length(A) > 1L && length(B) > 1L)) {
      # Neither side is a single cluster any more, so there is nothing left to
      # test: FunFHMMer never reconsiders such a node, nor any node above it.
      live[[k]] <- c(A, B)
      blocked[[k]] <- TRUE
      next
    }
    if (tree$height[[k]] > ceiling) {
      live[[k]] <- c(A, B)
      next
    }
    if (length(A) == 1L && length(B) == 1L) {
      merged <- test$join(list(A[[1L]], B[[1L]]))
      # The auto-merge head is taken on the merged alignment alone, as the
      # reference copies GeMMA's own merge-node alignment without scoring it.
      live[[k]] <- if (k <= auto_merged || !test$apart(A[[1L]], B[[1L]], merged))
                     list(merged) else c(A, B)
      next
    }
    # One child was cut. Test the intact cluster against each part of the cut
    # one, as FunFHMMer does before deciding whether to take one of them, all
    # of them, or none -- generalised past the two parts it ever looked at.
    intact <- if (length(A) == 1L) A[[1L]] else B[[1L]]
    parts <- if (length(A) == 1L) B else A
    with_intact <- lapply(parts, function(part) test$join(list(part, intact)))
    ok <- !vapply(seq_along(parts), function(j)
      test$apart(parts[[j]], intact, with_intact[[j]]), logical(1))
    # Taking more than one part at once needs the intact cluster to be
    # informative in its own right; without that the reference takes none.
    if (sum(ok) > 1L && !test$informative(intact)) ok[] <- FALSE
    live[[k]] <- if (!any(ok)) c(list(intact), parts)
                 else if (sum(ok) == 1L)
                   c(list(with_intact[[which(ok)]]), parts[!ok])
                 else c(list(test$join(c(list(intact), parts[ok]))), parts[!ok])
  }

  clusters <- lapply(live[[n - 1L]], function(cl) labels[sort(cl$i)])
  clusters <- clusters[order(-lengths(clusters))]
  names(clusters) <- paste0("funfam", seq_along(clusters))
  reps <- vapply(clusters, .repdist_medoid, character(1), D = D)
  similarity <- 1 - as.matrix(usedist::dist_subset(D, unname(reps)))
  dimnames(similarity) <- list(names(reps), names(reps))
  list(clusters = clusters,
       similarity = similarity,
       representatives = reps,
       tree = tree,
       separation = separation,
       auto_merge = auto_merge,
       test = if (is.null(seqs)) "distance" else "groupsim")
}
