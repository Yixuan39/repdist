# Structural clustering only, for the *complete* case_study2 protein universe.
# No differential abundance, no rarefaction, no substrate-based filtering (so,
# unlike case_study2.Rmd / export_case_study2_bins.R, PLA is NOT dropped here)
# -- every protein in ps.rds gets embedded and clustered. Numbered steps below
# each print a checkpoint so the run can be sanity-checked stage by stage.
# Output: one flat bin-list spreadsheet, both thresholds side by side.
library(phyloseq)
library(writexl)
devtools::load_all(".", quiet = TRUE)

path <- file.path("data", "case_study2")

## Step 1 -- the full curated protein universe: no PLA drop, no rarefaction
cat("== Step 1: load ps.rds ==\n")
ps <- readRDS(file.path(path, "ps.rds"))
cat(sprintf("  %d proteins x %d samples\n", ntaxa(ps), nsamples(ps)))

## Step 2 -- TM-Vec1 embeddings: reuse the cache, embed only what's missing
cat("== Step 2: TM-Vec1 embeddings ==\n")
full_seq   <- setNames(as.character(refseq(ps)[taxa_names(ps)]), taxa_names(ps))
model_tag  <- paste0(eval(formals(embed_proteins)$model), "-full-length-v1")
cache_file <- file.path(path, "embeddings_tmvec1_large.rds")
cached <- if (file.exists(cache_file)) readRDS(cache_file) else
  list(Z = NULL, seq = character(0), model = character(0))

reusable <- intersect(taxa_names(ps), rownames(cached$Z))
reusable <- reusable[cached$seq[reusable] == full_seq[reusable] &
                      cached$model[reusable] == model_tag]
need <- setdiff(taxa_names(ps), reusable)
cat(sprintf("  %d proteins reusable from cache, %d need embedding\n", length(reusable), length(need)))

Zall <- cached$Z
if (length(need)) {
  s <- full_seq[need]
  s <- s[order(nchar(s))]
  for (i in split(seq_along(s), ceiling(seq_along(s) / 2000))) {
    new_Z <- embed_proteins(s[i])
    Zall  <- rbind(Zall[setdiff(rownames(Zall), rownames(new_Z)), , drop = FALSE], new_Z)
    cached$Z <- Zall
    cached$seq[rownames(new_Z)]   <- full_seq[rownames(new_Z)]
    cached$model[rownames(new_Z)] <- model_tag
    saveRDS(cached, cache_file)
    cat(sprintf("  embedded %d more (checkpoint saved)\n", length(i)))
  }
}
Z <- Zall[taxa_names(ps), ]
stopifnot(identical(rownames(Z), taxa_names(ps)), all(is.finite(Z)))
cat(sprintf("  %d proteins embedded, none missing\n", nrow(Z)))

## Step 3 -- unrarefied counts, all 24 samples -- not analyzed, only used to
## rank bins by total abundance (bin1 = highest-total-abundance bin, ...),
## the same convention repdist_bin() always uses.
cat("== Step 3: counts (unrarefied, all samples, ranking only) ==\n")
cnt <- t(as(otu_table(ps), "matrix"))
cat(sprintf("  %d samples x %d proteins, %d total counts\n", nrow(cnt), ncol(cnt), sum(cnt)))

## Step 4 -- structural clustering at TM >= 0.5 and TM >= 0.7
cat("== Step 4: structural clustering (TM >= 0.5, TM >= 0.7) ==\n")
roll_both <- repdist_bin(cnt, Z, cluster_threshold = c(0.5, 0.7))$bins
for (thresh in names(roll_both)) {
  roll <- roll_both[[thresh]]
  n_singleton <- sum(table(roll$bin) == 1)
  cat(sprintf("  TM >= %s: %d bins for %d proteins (%d singleton bins)\n",
              thresh, ncol(roll$counts), length(roll$bin), n_singleton))
}

## Step 5 -- flat protein -> bin table, both thresholds, sorted by accession
cat("== Step 5: build bin list ==\n")
bin05 <- roll_both[["0.5"]]$bin
bin07 <- roll_both[["0.7"]]$bin
stopifnot(setequal(names(bin05), taxa_names(ps)), setequal(names(bin07), taxa_names(ps)))

bin_list <- data.frame(
  `Protein accession` = names(bin05),
  `Bin number 05` = unname(bin05),
  `Bin number 07` = unname(bin07[names(bin05)]),
  check.names = FALSE, row.names = NULL
)
bin_list <- bin_list[order(bin_list$`Protein accession`), ]
stopifnot(nrow(bin_list) == ntaxa(ps), !anyNA(bin_list))
cat(sprintf("  %d proteins, all bin-assigned at both thresholds, no NAs\n", nrow(bin_list)))

## Step 6 -- write
cat("== Step 6: write spreadsheet ==\n")
out <- file.path(path, "case_study2_bin_list_full.xlsx")
write_xlsx(list(`Protein - Bin list` = bin_list), out)
cat(sprintf("  wrote %s\n", out))
