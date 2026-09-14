#' Structure-aware protein distances, clustering and beta diversity
#'
#' Embed protein sequences with [read_fasta()] and
#' [embed_proteins()], or load the documented
#' [diet_metaproteome] and [funfam_universe]
#' examples. Everything else reads those embeddings.
#'
#' [repdist_matrix()] builds the
#' protein-by-protein distance, and
#' [bin_proteins()] cuts it into
#' complete-linkage structural bins;
#' [plot_similarity_profile()],
#' [plot_bin_profile()] and
#' [plot_bin_similarity()] are the
#' quality-control plots for choosing a threshold and reading a bin.
#'
#' [sample_repdist()] is the beta-diversity
#' side: it represents each sample as an abundance-weighted distribution over
#' the same embeddings and compares samples by Gaussian-kernel maximum mean
#' discrepancy. Counts have samples in rows and proteins in columns;
#' embedding rows are matched by protein identifier.
#'
#' See `vignette("repdist")` and `vignette("simulation")`.
#'
#' @references Gretton et al. (2012). A Kernel Two-Sample Test. JMLR 13,
#' 723-773. <https://jmlr.org/papers/v13/gretton12a.html>.
#'
#' Hamamsy et al. (2024). Protein remote homology detection and structural
#' alignment using deep learning. \doi{10.1038/s41587-023-01917-2}.
#'
#' Elnaggar et al. (2022). ProtTrans: Toward Understanding the Language of
#' Life Through Self-Supervised Learning. \doi{10.1109/TPAMI.2021.3095381}.
#'
#' Lin et al. (2023). Evolutionary-scale prediction of atomic-level protein
#' structure with a language model. \doi{10.1126/science.ade2574}.
#'
#' @author Yixuan Yang \email{yixuanyang309@icloud.com}
#'   (\href{https://orcid.org/0009-0003-5064-6512}{ORCID})
#'
#' @seealso Source and issue tracker:
#'   \url{https://github.com/Yixuan39/repdist}
#'
#' @importFrom utils globalVariables
"_PACKAGE"

#' Dietary metaproteome example
#'
#' A subset of the mouse dietary-protein study of Blakeley-Ruiz et al. (2025).
#' Stored as `extdata/diet_metaproteome.rds`; load with `readRDS()`.
#'
#' @format A list with `counts` (99 samples by 2000 proteins; spectral counts
#'   after full-catalog rarefaction and protein selection), `counts_raw` (same
#'   dimensions, original integer counts), `metadata` (99 rows: sample_id,
#'   subject, cage, diet, dose in percent, sex, age_week measured since
#'   baseline,
#'   period), `embeddings` (2000 by 512, TM-Vec coordinates rounded to 4
#'   decimals),
#'   `taxa` (protein, genus, species, annotation, length in residues),
#'   `sequences`
#'   (named amino-acid strings), and `provenance` (source and processing).
#'   Count column names, embedding row names, taxa row names and sequence names
#'   identify the same proteins; metadata rows match count rows. Missing
#'   taxonomy
#'   is `NA`. The final subset has unequal row totals. Renormalize rounded
#'   embeddings before using [repdist_matrix()].
#' @source Blakeley-Ruiz et al. (2025), Data Sets 1 and 6.
#'   \doi{10.1093/ismejo/wraf048}. CC BY 4.0. See installed
#'   `scripts/DATA_PROVENANCE.md` and `scripts/curate_diet_data.R`.
#' @examples
#' d <- readRDS(system.file("extdata", "diet_metaproteome.rds",
#'     package = "repdist"
#' ))
#' dim(d$counts)
#' @docType data
#' @name diet_metaproteome
NULL

#' CATH FunFam protein universe
#'
#' Ten CATH v4.3 functional families with 150 domains each, used for
#' simulation. Stored as `extdata/funfam_universe.rds`; load with
#' `readRDS()`.
#'
#' @format A list with `members` (1500 rows: protein identifier, family, label),
#'   `embeddings` (1500 by 512, historical scikit-bio/tmvec coordinates rounded
#'   to 4 decimals), `sequences` (named domain amino-acid strings), and
#'   `provenance` (source, selection parameters and model identity).
#'   Member identifiers match embedding row names and sequence names in order.
#'   Renormalize rounded rows before computing cosine distances.
#' @details The historical head differs from the current default `tmvec`.
#'   Original model revisions were not recorded; these fixed coordinates are
#'   an example, not a bitwise inference reference. The installed
#'   `scripts/DATA_PROVENANCE.md` describes recovery of the exact
#'   sequence set and recomputation using a pinned version of the historical
#'   head.
#' @source CATH v4.3 FunFam Stockholm alignments, accessed 2026-08-23.
#'   <https://www.cathdb.info/>. CATH database content: CC BY 4.0.
#' @examples
#' u <- readRDS(system.file("extdata", "funfam_universe.rds",
#'     package = "repdist"
#' ))
#' table(u$members$label)
#' @docType data
#' @name funfam_universe
NULL
