# repdist 0.99.0

* First Bioconductor submission: protein embeddings from a frozen protein
  language model (`embed_proteins()`) or the tmvec command-line tool
  (`read_embeddings()`), abundance-weighted MMD sample distances under a
  Gaussian kernel (`sample_mmd()`), structural bins by Markov clustering of
  a graph gated on alignment E-value (`protein_similarity()`,
  `bin_proteins()`, `bin_glom()`) and plots (`plot_tsne()`,
  `plot_bin_similarity()`).
* One vignette, `vignette("repdist")`, on `data(simulate_cog)`: 20 simulated
  samples over 110 NCBI COG proteins with precomputed TM-Vec embeddings, in
  which group A is enriched for GrpE proteins from ten different
  30%-identity sequence clusters. It covers embeddings, beta diversity,
  protein bins and differential abundance. `inst/extdata` ships the
  dataset's FASTA and an example of each input file (FASTA, counts, sample
  metadata) with matching `tmvec` embeddings.
