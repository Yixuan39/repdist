#' Structure-aware protein distances, clustering and beta diversity
#'
#' Embed protein sequences with [read_fasta()] and
#' [embed_proteins()],
#' or read embeddings made with the tmvec tool with
#' [read_embeddings()]. Each gives an `AAStringSet` with the embeddings in
#' `mcols()`.
#'
#' [protein_similarity()] builds the
#' protein-by-protein similarity graph, gated on alignment E-value, and
#' [bin_proteins()] groups proteins into structural bins by Markov clustering
#' of it; [bin_glom()] sums a count table to those bins.
#' [plot_tsne()] maps the embeddings and [plot_bin_similarity()] shows the
#' pairwise similarities inside one bin.
#'
#' [sample_mmd()]
#' is the beta-diversity side: it represents each sample
#' as an abundance-weighted distribution over the same embeddings and compares
#' samples by Gaussian-kernel maximum mean discrepancy. Counts have samples in
#' rows and proteins in columns, matched to the proteins by name.
#'
#' See `vignette("repdist")` for beta diversity, protein bins and
#' differential abundance.
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
# brings back, such as the vignette's data, cannot be subset with `[`.
#' @importFrom Biostrings AAStringSet
NULL

#' Simulated one-COG metaproteome
#'
#' Twenty simulated samples, ten in group A and ten in group B, over 110
#' bacterial and archaeal proteins of the NCBI COG 2024 release, each with
#' exactly one COG assigned over its full length. 100 background proteins,
#' drawn at random from other COGs, are in every sample. Each group-A sample
#' also holds one GrpE protein (COG0576, the nucleotide exchange factor of
#' HSP-70) of its own, and the ten GrpE proteins come from ten different
#' 30%-identity MMseqs2 clusters: group A is enriched for one function, while
#' no sequence family is enriched in more than one sample. Each group-A
#' sample puts 5% of its abundance on its GrpE protein. Load with
#' `data(simulate_cog, package = "repdist")`; the sequences alone are in
#' `system.file("extdata", "simulate_cog.fasta", package = "repdist")`.
#'
#' @format A list of five elements:
#' - `proteins`: `AAStringSet` of the sequences, with TM-Vec embeddings
#'   rounded to 4 decimals in `mcols(proteins)$embedding` (110 x 512).
#' - `counts`: integer matrix, 20 samples x 110 proteins.
#' - `samples`: data frame with `sample`, `group` (A or B) and
#'   `grpe_protein`, the GrpE protein of a group-A sample.
#' - `annotation`: data frame, one row per protein: `protein` (NCBI protein
#'   id), `role` (GrpE or background), `enriched_in` (the sample a GrpE
#'   protein is spiked into), `cog`, `cog_name`, `gene`, `category` and
#'   `category_name` (COG functional category), `length`, `organism`, and
#'   `mmseqs30` (its MMseqs2 30%-identity cluster representative).
#' - `provenance`: source, license, design, seeds and session.
#'
#' @source NCBI COG 2024, <https://ftp.ncbi.nih.gov/pub/COG/COG2024/data/>.
#'   Built by `inst/scripts/make_simulate_cog.R`; see
#'   `inst/scripts/DATA_PROVENANCE.md`.
#'
#' @references Galperin MY, Vera Alvarez R, Karamycheva S, Makarova KS, Wolf
#' YI, Landsman D, Koonin EV (2025). COG database update 2024. Nucleic Acids
#' Research 53(D1), D356-D363. \doi{10.1093/nar/gkae983}.
#'
#' @examples
#' data(simulate_cog, package = "repdist")
#' dim(simulate_cog$counts)
#' table(simulate_cog$annotation$role)
"simulate_cog"
