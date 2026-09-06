# From the repository root:
# Rscript analysis/case_study2_annotation.R [PUP.csv sequences.fasta annotations.xlsx]
# Rebuild the input object and compare annotation candidates, without refitting DA.
suppressPackageStartupMessages(library(phyloseq))
source("R/qc.R")
source("R/distances.R")
source("R/binning.R")
source("analysis/annotation_clustering.R")

files <- commandArgs(trailingOnly = TRUE)
if (!length(files)) files <- file.path(path.expand("~/Downloads"), c(
  "ProteomeMeanPUP.csv", "ProteomOnlyFasta.fasta",
  "ProteomeMeanPUPNSAFPlastic_DepolyAB_TIGERFAM_FinalAnot_Master_MantisSorted_Current.xlsx"))
stopifnot(length(files) == 3L, all(file.exists(files)))
out <- "results/case_study2_annotation_clustering"
path <- "data/case_study2"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Blank abundance cells denote absence. Reject unexpected text, negative or
# infinite values before using the technical-replicate means.
raw <- read.csv(files[1], check.names = FALSE, colClasses = "character")
ids <- raw[[1]]
stopifnot(!anyNA(ids), all(nzchar(ids)), !anyDuplicated(ids), names(raw)[2] == "AA")
mean_pup <- vapply(raw[, -(1:2)], function(x) {
  blank <- is.na(x) | !nzchar(trimws(x))
  value <- suppressWarnings(as.numeric(x))
  stopifnot(all(blank | is.finite(value)), all(blank | value >= 0))
  value[blank] <- 0
  value
}, numeric(nrow(raw)))
rownames(mean_pup) <- ids
seqs <- Biostrings::readAAStringSet(files[2])
stopifnot(!anyDuplicated(names(seqs)), all(names(seqs) %in% ids))
keep <- ids[ids %in% names(seqs)] # Preserve CSV order for deterministic rarefaction.
fasta_matched <- length(keep)
stopifnot(all(as.numeric(raw[match(keep, ids), "AA"]) ==
                Biostrings::width(seqs[keep])))
write.csv(data.frame(protein = setdiff(ids, keep)),
           file.path(out, "abundance_without_sequence.csv"), row.names = FALSE)

sheet <- as.data.frame(readxl::read_excel(files[3], skip = 1,
                                         col_types = "text", .name_repair = "minimal"))
stopifnot(names(sheet)[25] == "Protein ID", names(sheet)[36] == "Protein ID")
names(sheet) <- make.unique(names(sheet))
workbook_rows <- nrow(sheet)
key <- trimws(sheet[[36]])
stopifnot(!anyNA(key), !anyDuplicated(key), all(key %in% keep),
            all(as.numeric(sheet[[37]]) == Biostrings::width(seqs[key])))
# The second ID column is unique and agrees with sequence lengths. Do not join
# on the first column: it contains two copied IDs and a free-text note.
bad <- trimws(sheet[[25]]) != key
write.csv(data.frame(excel_row = which(bad) + 2L,
                     first_id = sheet[[25]][bad], protein = key[bad],
                     genome = sheet[[21]][bad], MantisName = sheet[[7]][bad]),
           file.path(out, "annotation_id_audit.csv"), row.names = FALSE)
sheet <- sheet[match(keep, key), , drop = FALSE]
rownames(sheet) <- keep
tax <- as.matrix(sheet[, 1:24])
tax <- cbind(tax, GenomeID = sheet[["Bin ID"]],
              Depolymerase = sheet[["Potential plastic depolymerase"]])

# Substrate mapping is written in the individual-sample header of the workbook.
substrates <- c(A = "PHBV", B = "PLA", C = "HDPE", D = "PLA/PCL",
                 F = "Filter paper", G = "PHBH", I = "PHBH Green", S = "Sediment")
samples <- colnames(mean_pup)
stopifnot(identical(samples, c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2",
                               "C3", "D1", "D2", "D3", "F2", "F3", "F7", "G1",
                               "G2", "G3", "I1", "I2", "I3", "S1", "S2", "S3")))
md <- data.frame(substrate = factor(unname(substrates[substr(samples, 1, 1)])),
                  replicate = substring(samples, 2), row.names = samples)
mean_pup <- mean_pup[keep, , drop = FALSE]
# The legacy report's PLA totals (1832, 766, 92) identify round(), rather than
# truncation, of mean PUP. Preserve the original means separately: these rounded
# values reproduce a legacy preprocessing convention, not raw spectral counts.
counts <- round(mean_pup)
ps <- phyloseq(otu_table(counts, taxa_are_rows = TRUE), sample_data(md),
                tax_table(tax), refseq(seqs[keep]))
saveRDS(ps, file.path(path, "ps.rds"))
saveRDS(mean_pup, file.path(path, "mean_pup.rds"))

ps <- prune_samples(sample_data(ps)$substrate != "PLA", ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
set.seed(1)
depth <- min(sample_sums(ps))
psd <- rarefy_even_depth(ps, sample.size = depth, rngseed = 1,
                         replace = FALSE, verbose = FALSE)
before_length <- ntaxa(psd)
long <- Biostrings::width(refseq(psd)) > 3000L
long_share <- sum(taxa_sums(psd)[long]) / sum(taxa_sums(psd))
psd <- prune_taxa(!long, psd)
cached <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))
missing <- setdiff(taxa_names(psd), rownames(cached$Z))
missing_abundance <- sum(taxa_sums(psd)[missing]) / sum(taxa_sums(psd))
write.csv(data.frame(protein = missing, length = Biostrings::width(refseq(psd)[missing]),
                      rarefied_abundance = unname(taxa_sums(psd)[missing])),
           file.path(out, "proteins_without_embeddings.csv"), row.names = FALSE)
psd <- prune_taxa(taxa_names(psd) %in% rownames(cached$Z), psd)
keep <- taxa_names(psd)
saveRDS(psd, file.path(path, "annotation_cohort.rds"))
stopifnot(all(keep %in% rownames(cached$Z)),
            identical(unname(as.character(refseq(psd)[keep])), unname(cached$seq[keep])))
Z <- cached$Z[keep, , drop = FALSE]
rm(cached, seqs, counts, mean_pup)
invisible(gc())

clean <- function(x) {
  x <- trimws(tolower(x))
  x[is.na(x) | x %in% c("", "-", ".", "na", "unknown", "uncharacterized protein",
                         "hypothetical protein", "predicted protein")] <- NA_character_
  x
}
single_ko <- function(x) vapply(regmatches(x, gregexpr("K[0-9]{5}", x)), function(k) {
  k <- unique(k)
  if (length(k) == 1L) k else NA_character_
}, character(1))
stopifnot(identical(single_ko(c("ko:K00001", "K00001,K00002", "-", NA)),
                    c("K00001", NA_character_, NA_character_, NA_character_)))
ann <- sheet[keep, , drop = FALSE]
ko <- cbind(MantisKO = single_ko(ann$MantisKO),
              MicrobeKO = single_ko(ann$MicrobeKO), EggNOG_KO = single_ko(ann[["KO Term"]]))
# Require at least two sources and no disagreement among sources with a single
# KO. This is agreement between annotation pipelines, not a curated gold standard.
consensus <- apply(ko, 1, function(x) {
  x <- x[!is.na(x)]
  if (length(x) >= 2L && length(unique(x)) == 1L) x[1] else NA_character_
})
labels <- data.frame(protein = keep, consensus_KO = consensus, ko,
                      curated_name = clean(ann$Name), MantisName = clean(ann$MantisName),
                      MicrobeName = clean(ann$MicrobeName),
                      broad_function = clean(ann[["Broad Function"]]),
                      detailed_function = clean(ann[["Detailed Function"]]),
                      stringsAsFactors = FALSE, row.names = keep)
labels$KO_source_conflict <- apply(ko, 1, function(x) length(unique(x[!is.na(x)])) > 1L)
labels$genome <- ann[["Bin ID"]]
labels$annotation_notes <- ann[["Annotation notes"]]
write.csv(labels, file.path(out, "annotations.csv"), row.names = FALSE)
write.csv(labels[labels$KO_source_conflict, ],
           file.path(out, "annotation_source_conflicts.csv"), row.names = FALSE)
provenance <- data.frame(csv_proteins = nrow(raw), fasta_matched = fasta_matched,
                          workbook_rows = workbook_rows, id_disagreements = sum(bad),
                          samples = nsamples(psd), depth = depth,
                          before_length = before_length, removed_long = sum(long),
                          removed_long_abundance = long_share, proteins = length(keep),
                          without_embeddings = length(missing),
                          without_embeddings_abundance = missing_abundance,
                          legacy_report_proteins = 5250L,
                          consensus_KO = sum(!is.na(labels$consensus_KO)),
                          KO_source_conflicts = sum(labels$KO_source_conflict))
write.csv(provenance, file.path(out, "provenance.csv"), row.names = FALSE)
write.csv(data.frame(file = normalizePath(files), md5 = unname(tools::md5sum(files))),
           file.path(out, "inputs.csv"), row.names = FALSE)
message("Analysis cohort: ", length(keep), " proteins; depth ", depth,
          "; consensus KO available for ", sum(!is.na(labels$consensus_KO)))

grid <- rbind(
  data.frame(method = "hclust", min_sim = c(.5, .7, .9, .95), min_pts = NA, inflation = NA),
  data.frame(method = "mcl", expand.grid(min_sim = c(.5, .7, .9),
               min_pts = NA, inflation = c(2, 4))),
  data.frame(method = "dbscan", expand.grid(min_sim = c(.7, .9),
               min_pts = c(3, 5), inflation = NA)),
  data.frame(method = "hdbscan", min_sim = NA, min_pts = c(3, 5, 10), inflation = NA))
grid$border_points <- ifelse(grid$method == "dbscan", FALSE, NA)
grid <- rbind(grid, data.frame(method = "dbscan", min_sim = .9,
                               min_pts = c(3, 5), inflation = NA, border_points = TRUE))
# Append the TM80 check so the existing run IDs and exported examples stay stable.
grid <- rbind(grid, data.frame(method = "hclust", min_sim = .8,
                               min_pts = NA, inflation = NA, border_points = NA))
message("Computing the shared cosine distance ...")
G <- repdist_matrix(Z)
result <- compare_annotation_clusters(G, setNames(labels$consensus_KO, keep), grid)
write.csv(result$summary, file.path(out, "summary.csv"), row.names = FALSE)
write.csv(result$membership, file.path(out, "membership.csv"), row.names = FALSE)
saveRDS(list(G = G, result = result, labels = labels), file.path(out, "comparison.rds"))

references <- c("consensus_KO", "MantisKO", "MicrobeKO", "EggNOG_KO", "curated_name")
scores <- do.call(rbind, lapply(references, function(ref) {
  scored <- lapply(seq_len(nrow(grid)), function(i) {
    bin <- setNames(result$membership[[paste0("run", i)]], keep)
    noise <- keep[result$membership[[paste0("noise", i)]]]
    annotation_cluster_scores(bin, setNames(labels[[ref]], keep), noise)
  })
  cbind(reference = ref, result$summary[, c("run", names(grid))],
        as.data.frame(do.call(rbind, scored)))
}))
write.csv(scores, file.path(out, "scores_by_annotation.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(result$summary, row.names = FALSE, digits = 3)
