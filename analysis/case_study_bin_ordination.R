# Standalone: case_study 1 protein universe (Blakeley-Ruiz mouse diet
# metaproteome), bin_proteins()'d at the same-fold TM-score threshold, then
# three ordinations (PCoA, UMAP, t-SNE) run directly on the structural
# distance between bin medoids. No sample counts, no rarefaction, no
# reticulate: binning depends only on the embeddings, which are read from the
# cache written by case_study.Rmd's Embedding section rather than recomputed.
library(phyloseq)
library(ggplot2)
library(here)
library(repdist)

path <- here("data", "test_study")

ps  <- readRDS(file.path(path, "ps.rds"))
psu <- prune_samples(sample_data(ps)$phase != "baseline", ps)   # same protein universe
psu <- prune_taxa(taxa_sums(psu) > 0, psu)                      # as case_study.Rmd's cnt_da,
proteins <- taxa_names(psu)                                     # unrarefied so no protein is
                                                                 # dropped by subsampling RNG

Z <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))$Z[proteins, ]

b <- bin_proteins(repdist_matrix(Z), min_sim = 0.5)
cat(sprintf("bin_proteins: %d proteins -> %d structural bins\n",
            length(proteins), length(b$clusters)))
membership <- setNames(rep(names(b$clusters), lengths(b$clusters)),
                       unlist(b$clusters))

# same functional-annotation call case_study.Rmd uses (its `is_dark()`, negated):
# Table 6's "Consensus annotation" is NA for a protein Mantis/eggNOG couldn't call
t6 <- read.delim(file.path(path, "annotated", "ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt"),
                 skip = 2, quote = "\"", check.names = FALSE, colClasses = "character",
                 comment.char = "", fileEncoding = "latin1")
annotated <- !is.na(setNames(t6[["Consensus annotation"]], t6[["Protein Identifier"]])[proteins])
pct_annotated <- vapply(split(annotated, membership[proteins]), mean, numeric(1)) * 100

Dbin <- as.dist(1 - b$similarity)
pct_annotated <- pct_annotated[labels(Dbin)]

pc <- ape::pcoa(Dbin)
ve <- round(100 * pc$values$Relative_eig[1:2], 1)
um <- uwot::umap(Dbin, n_neighbors = min(15, attr(Dbin, "Size") - 1))
ts <- Rtsne::Rtsne(as.matrix(Dbin), is_distance = TRUE,
                   perplexity = min(30, floor((attr(Dbin, "Size") - 1) / 3)))

df <- rbind(
  data.frame(method = sprintf("PCoA [%s%%, %s%%]", ve[1], ve[2]), x = pc$vectors[, 1], y = pc$vectors[, 2]),
  data.frame(method = "UMAP", x = um[, 1], y = um[, 2]),
  data.frame(method = "t-SNE", x = ts$Y[, 1], y = ts$Y[, 2])
)
df$pct_annotated <- rep(pct_annotated, 3)
df$method <- factor(df$method, levels = unique(df$method))

p <- ggplot(df, aes(x, y, color = pct_annotated)) +
  geom_point(size = 0.8, alpha = 0.7) +
  facet_wrap(~method, scales = "free") +
  scale_color_viridis_c(name = "% annotated\n(Table 6)", limits = c(0, 100)) +
  labs(x = NULL, y = NULL,
      title = "Structural bins only -- ordination of medoid distances, colored by % functionally annotated") +
  theme_bw()

out <- file.path(path, "case_study_bin_ordination.png")
ggsave(out, p, width = 12, height = 4.2, dpi = 150)
cat("saved:", out, "\n")
