# Standalone: case_study 1 protein universe (Blakeley-Ruiz mouse diet
# metaproteome), repdist_bin()'d at the default TM-score threshold, then
# three ordinations (PCoA, UMAP, t-SNE) run directly on bin_dist -- the
# structural distance between bin medoids. No sample counts, no rarefaction,
# no reticulate: bin_dist never depends on abundance (repdist_bin()'s
# clustering and medoid choice are embeddings-only; counts only rank which
# bin is "bin1"), and embeddings are read from the cache written by
# case_study.Rmd's Embedding section rather than recomputed, so there is no
# Python call to make here.
library(phyloseq)
library(ggplot2)
library(here)
devtools::load_all(".", quiet = TRUE)

path <- here("data", "test_study")

ps  <- readRDS(file.path(path, "ps.rds"))
psu <- prune_samples(sample_data(ps)$phase != "baseline", ps)   # same protein universe
psu <- prune_taxa(taxa_sums(psu) > 0, psu)                      # as case_study.Rmd's cnt_da,
proteins <- taxa_names(psu)                                     # unrarefied so no protein is
                                                                 # dropped by subsampling RNG

Z <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))$Z[proteins, ]

# repdist_bin() takes a counts matrix, but bin_dist doesn't use it: clustering
# and medoid choice come from embeddings alone. A single row of 1s makes that
# explicit -- bin totals become protein counts, so "bin1" just means "largest bin".
presence <- matrix(1L, 1, length(proteins), dimnames = list("all", proteins))
b <- repdist_bin(presence, Z, cluster_threshold = 0.5)$bins
cat(sprintf("repdist_bin: %d proteins -> %d structural bins\n", length(proteins), nrow(b$reps)))

# same functional-annotation call case_study.Rmd uses (its `is_dark()`, negated):
# Table 6's "Consensus annotation" is NA for a protein Mantis/eggNOG couldn't call
t6 <- read.delim(file.path(path, "annotated", "ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt"),
                 skip = 2, quote = "\"", check.names = FALSE, colClasses = "character",
                 comment.char = "", fileEncoding = "latin1")
annotated <- !is.na(setNames(t6[["Consensus annotation"]], t6[["Protein Identifier"]])[proteins])
pct_annotated <- vapply(split(annotated, b$bin[proteins]), mean, numeric(1)) * 100

Dbin <- b$bin_dist
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
      title = "Structural bins only -- ordination of bin_dist, colored by % functionally annotated") +
  theme_bw()

out <- file.path(path, "case_study_bin_ordination.png")
ggsave(out, p, width = 12, height = 4.2, dpi = 150)
cat("saved:", out, "\n")
