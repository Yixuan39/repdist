# repdist 0.99.0

* First Bioconductor submission: protein embeddings from a frozen protein
  language model (`embed_proteins()`) or the tmvec command-line tool
  (`read_embeddings()`), abundance-weighted RBF-MMD sample distances
  (`sample_mmd()`), structural bins by Markov clustering of a graph
  gated on alignment E-value (`protein_similarity()`, `bin_proteins()`,
  `bin_glom()`) and quality-control plots (`plot_tsne()`,
  `plot_bin_similarity()`, `plot_pcoa()`).
* The package has two vignettes, both on `data(simulate_ec)`, a simulated
  metaproteome with precomputed TM-Vec embeddings: `vignette("repdist")`
  covers beta diversity, and `vignette("clustering")` covers protein bins
  and differential abundance. The clustering vignette ships a precomputed
  `protein_similarity()` result to stay fast.
