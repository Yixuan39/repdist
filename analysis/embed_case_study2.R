# Standalone embedding job for case_study2 -- run once in the background.
# Same caching contract as case_study.Rmd's embed chunk: a cache list with
# Z / seq / model, keyed by exact sequence + model tag, saved incrementally
# so an interrupted run resumes rather than restarting.
devtools::load_all(".", quiet = TRUE)

path <- file.path("data", "case_study2")
ps <- readRDS(file.path(path, "ps.rds"))

model_tag <- paste0(eval(formals(embed_proteins)$model), "-full-length-v1")
cache_file <- file.path(path, "embeddings_tmvec1_large.rds")
cached <- if (file.exists(cache_file)) readRDS(cache_file) else
  list(Z = NULL, seq = character(0), model = character(0))

full_seq <- setNames(as.character(phyloseq::refseq(ps)), phyloseq::taxa_names(ps))

reusable <- intersect(phyloseq::taxa_names(ps), rownames(cached$Z))
reusable <- reusable[cached$seq[reusable] == full_seq[reusable] &
                      cached$model[reusable] == model_tag]
need <- setdiff(phyloseq::taxa_names(ps), reusable)
Zall <- cached$Z

cat(sprintf("[%s] total=%d reusable=%d need=%d\n", Sys.time(),
            length(full_seq), length(reusable), length(need)))

if (length(need)) {
  s <- full_seq[need]
  s <- s[order(nchar(s))]
  blocks <- split(seq_along(s), ceiling(seq_along(s) / 2000))
  for (b in seq_along(blocks)) {
    i <- blocks[[b]]
    t0 <- Sys.time()
    new_Z <- embed_proteins(s[i])
    Zall <- rbind(Zall[setdiff(rownames(Zall), rownames(new_Z)), , drop = FALSE], new_Z)
    cached$Z <- Zall
    cached$seq[rownames(new_Z)] <- full_seq[rownames(new_Z)]
    cached$model[rownames(new_Z)] <- model_tag
    saveRDS(cached, cache_file)
    cat(sprintf("[%s] block %d/%d done (%d proteins, %.1f sec) -- checkpoint saved\n",
                Sys.time(), b, length(blocks), length(i),
                as.numeric(Sys.time() - t0, units = "secs")))
  }
}
cat(sprintf("[%s] embedding complete: %d proteins embedded\n", Sys.time(), nrow(Zall)))
