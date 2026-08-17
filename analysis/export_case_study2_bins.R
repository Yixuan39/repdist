# Standalone export: structural bins for the whole case_study2 protein
# universe (30,997 proteins, PLA already dropped -- same set case_study2.Rmd
# clusters), at TM-score >= 0.5 and >= 0.7. Same repdist_bin()
# + LinDA + Wald-test pipeline as case_study2.Rmd's "Differential abundance"
# section, just run once here and exported to spreadsheets instead of shown
# inline. The pairwise TM-Vec distance matrix repdist_bin()
# builds internally is not saved to disk this time -- only the bins it produces.
library(phyloseq)
library(MicrobiomeStat)
library(writexl)
devtools::load_all(".", quiet = TRUE)

path <- file.path("data", "case_study2")

ps <- readRDS(file.path(path, "ps.rds"))
ps <- prune_samples(data.frame(sample_data(ps))$substrate != "PLA", ps)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)
sample_data(ps)$substrate <- droplevels(sample_data(ps)$substrate)
md <- data.frame(sample_data(ps))

full_seq <- setNames(as.character(refseq(ps)[taxa_names(ps)]), taxa_names(ps))
cached    <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))
Z         <- cached$Z[taxa_names(ps), , drop = FALSE]
cnt_da    <- t(as(otu_table(ps), "matrix"))
tax_all   <- as.data.frame(as(tax_table(ps), "matrix"))

cat(sprintf("[%s] clustering %d proteins...\n", Sys.time(), ncol(cnt_da)))
roll_both <- repdist_bin(cnt_da, Z, cluster_threshold = c(0.5, 0.7))$bins

# Per-bin strongest substrate contrast (max |log2FC|) + joint Wald significance,
# for every bin -- not just ones that reject. NA where LinDA's own prevalence
# filter dropped the bin before fitting (present in <3 samples).
bin_stats <- function(roll, md) {
  fit <- linda(t(roll$counts), md, formula = "~ substrate",
               feature.dat.type = "count", prev.filter = 3 / nrow(roll$counts), verbose = FALSE)

  substrate_vars <- grep("^substrate", fit$variables, value = TRUE)
  volc <- do.call(rbind, lapply(substrate_vars, function(v) {
    o <- fit$output[[v]]
    data.frame(bin = rownames(o), substrate = sub("^substrate", "", v),
               log2FoldChange = o$log2FoldChange, pvalue = o$pvalue,
               stringsAsFactors = FALSE)
  }))
  best <- do.call(rbind, lapply(split(volc, volc$bin), function(d) d[which.max(abs(d$log2FoldChange)), ]))
  rownames(best) <- best$bin

  k <- nlevels(md$substrate) - 1
  L <- matrix(0, nrow = k, ncol = length(fit$variables) + 1)
  for (i in 1:k) L[i, i + 1] <- 1
  wald <- linda.wald.test(fit, L, model = "LM", alpha = 0.05)
  rownames(wald) <- rownames(fit$output[[1]])

  all_bins <- unique(roll$bin)
  data.frame(
    Bin = all_bins,
    `N members` = as.integer(table(roll$bin)[all_bins]),
    `Strongest substrate (vs Sediment)` = best[all_bins, "substrate"],
    Log2FC = round(best[all_bins, "log2FoldChange"], 3),
    `p-value` = best[all_bins, "pvalue"],
    `Significant (joint Wald, padj<=0.05)` = wald[all_bins, "reject"],
    check.names = FALSE, row.names = NULL
  )
}

tax_cols <- function(members) {
  t <- tax_all[members, ]
  data.frame(
    Protein = members,
    `Genome ID` = t$GenomeID,
    `Shared genomes` = t$SharedGenomes,
    `Mantis annotation` = t$MantisName,
    `Microbe annotation` = t$MicrobeName,
    PlasticDB = t$PlasticDB,
    `Possible depolymerase` = t$Depolymerase,
    Description = t$Description,
    `Length (aa)` = nchar(full_seq[members]),
    check.names = FALSE, row.names = NULL
  )
}

protein_bin_list <- list()

for (thresh in c("0.5", "0.7")) {
  tag  <- if (thresh == "0.5") "tm05" else "tm07"
  roll <- roll_both[[thresh]]

  cat(sprintf("[%s] LinDA + Wald on %d bins (%s)...\n", Sys.time(), ncol(roll$counts), tag))
  idx <- bin_stats(roll, md)
  idx <- idx[order(-abs(idx$Log2FC), na.last = TRUE), ]

  # Individual member sheets only for bins with >=2 members -- a sheet for a
  # single protein repeats what "Protein - Bin list" already says.
  # ponytail: skips a sheet-per-singleton-bin (thousands of 1-row sheets add
  # nothing); add if a per-bin tab for singletons specifically is wanted.
  multi_bins <- names(table(roll$bin))[table(roll$bin) >= 2]
  bin_sheets <- setNames(
    lapply(multi_bins, function(b) tax_cols(names(roll$bin)[roll$bin == b])),
    multi_bins)

  full_pbl <- data.frame(`Protein accession` = names(roll$bin), `Bin number` = unname(roll$bin),
                         check.names = FALSE, row.names = NULL)
  protein_bin_list[[tag]] <- setNames(full_pbl, c("Protein accession", paste("Bin number", substr(tag, 3, 4))))

  out <- file.path(path, sprintf("case_study2_bins_%s.xlsx", tag))
  write_xlsx(c(list(Index = idx, `Protein - Bin list` = full_pbl), bin_sheets), out)
  cat(sprintf("[%s] wrote %s: %d bins total, %d with >=2 members, %d proteins\n",
              Sys.time(), out, ncol(roll$counts), length(multi_bins), ncol(cnt_da)))
}

pbl_merged <- merge(protein_bin_list[["tm05"]], protein_bin_list[["tm07"]], by = "Protein accession")
pbl_merged <- pbl_merged[order(pbl_merged$`Protein accession`), ]

pbl_out <- file.path(path, "case_study2_protein_bin_list.xlsx")
write_xlsx(list(`Protein - Bin list` = pbl_merged), pbl_out)
cat(sprintf("[%s] wrote %s\n", Sys.time(), pbl_out))
