#' Structure-aware protein distances, clustering and beta diversity
#'
#' Embed protein sequences with [read_fasta()] and
#' [embed_proteins()],
#' read embeddings made with the tmvec tool with
#' [read_embeddings()], or load the documented
#' [diet_metaproteome] example. Each gives an `AAStringSet` with the embeddings in
#' `mcols()`.
#'
#' [protein_similarity()] builds the
#' protein-by-protein similarity graph, gated on similarity and alignment
#' coverage, and [bin_proteins()] groups proteins into
#' structural bins by Markov clustering of it;
#' [bin_glom()] sums a count table to those bins.
#' [plot_tsne()] maps the embeddings and [plot_bin_similarity()] shows the
#' pairwise similarities inside one bin.
#'
#' [sample_mmd()]
#' is the beta-diversity side: it represents each sample
#' as an abundance-weighted distribution over the same embeddings and compares
#' samples by Gaussian-kernel maximum mean discrepancy. Counts have samples in
#' rows and proteins in columns, matched to the proteins by name.
#'
#' See `vignette("repdist")`.
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

# Loads Biostrings with the package. Without it, an AAStringSet that readRDS()
# brings back, such as the example data's, cannot be subset with `[`.
#' @importFrom Biostrings AAStringSet
NULL

#' Dietary metaproteome example
#'
#' A subset of the mouse dietary-protein study of Blakeley-Ruiz et al. (2025).
#' Stored as `extdata/diet_metaproteome.rds`; load with `readRDS()`.
#'
#' @format A list with `counts` (99 samples by 2000 proteins; spectral counts
#'   after full-catalog rarefaction and protein selection), `counts_raw` (same
#'   dimensions, original integer counts), `metadata` (99 rows: sample_id,
#'   subject, cage, diet, dose in percent, sex, age_week measured since
#'   baseline, period), `proteins` (2000 sequences as an `AAStringSet`,
#'   with TM-Vec embeddings rounded to 4 decimals in
#'   `mcols(proteins)$embedding`),
#'   `taxa` (protein, genus, species, annotation, length in residues), and
#'   `provenance` (source and processing). Count column names, protein names
#'   and taxa row names identify the same proteins; metadata rows match count
#'   rows. Missing taxonomy is `NA`. The final subset has unequal row
#'   totals.
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
