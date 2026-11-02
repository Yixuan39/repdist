#' Structure-aware protein distances, clustering and beta diversity
#'
#' Embed protein sequences with [read_fasta()] and
#' [embed_proteins()],
#' or read embeddings made with the tmvec tool with
#' [read_embeddings()]. Each gives an `AAStringSet` with the embeddings in
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
# brings back, such as the vignette's data, cannot be subset with `[`.
#' @importFrom Biostrings AAStringSet
NULL

#' Simulated EC metaproteome
#'
#' Twenty simulated samples over 1,160 bacterial and archaeal Swiss-Prot
#' enzymes with experimentally supported EC numbers: every enzyme of EC 2.7.7
#' and EC 2.3.2, and 1,000 other enzymes drawn at random as background. In
#' group A each sample puts 5% of its expected abundance on two EC 2.7.7
#' enzymes, in group B on two EC 2.3.2 enzymes. No two enriched enzymes share a
#' 50%-identity MMseqs2 cluster, so the groups share no enriched protein or
#' sequence family. Contrast 1 of the website's EC simulation, on a smaller
#' protein universe. Load with `data(simulate_ec, package = "repdist")`.
#'
#' @format A list of five elements:
#' - `proteins`: `AAStringSet` of the sequences, with TM-Vec
#'   embeddings rounded to 4 decimals in `mcols(proteins)$embedding`
#'   (1160 x 512).
#' - `counts`: integer matrix, 20 samples x 1160 proteins.
#' - `samples`: data frame with `sample`, `group` (A or B) and
#'   `enriched_ec3`.
#' - `annotation`: data frame, one row per protein: `protein`
#'   (UniProt accession), `name`, `organism`, `length`, `ec3`,
#'   `ec4`, `mmseqs50` and `mmseqs30` (MMseqs2 cluster representatives
#'   at 50% and 30% identity), and `enriched_in`, the sample a protein
#'   is enriched in (`NA` for baseline proteins).
#' - `provenance`: source, UniProt release, license, universe filtering
#'   counts, embedding model, MMseqs2 commands, design and seeds.
#'
#' @source UniProtKB/Swiss-Prot, <https://www.uniprot.org>, CC BY 4.0. Built
#'   by `inst/scripts/make_simulate_ec.R`; see
#'   `inst/scripts/DATA_PROVENANCE.md`.
#'
#' @references The UniProt Consortium (2025). UniProt: the Universal Protein
#' Knowledgebase in 2025. \doi{10.1093/nar/gkae1010}.
#'
#' @examples
#' data(simulate_ec, package = "repdist")
#' dim(simulate_ec$counts)
#' table(simulate_ec$samples$group, simulate_ec$samples$enriched_ec3)
"simulate_ec"
