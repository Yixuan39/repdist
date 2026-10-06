# repdist 0.99.0

* First Bioconductor submission: protein embeddings from a frozen protein
  language model (`embed_proteins()`) or the tmvec command-line tool
  (`read_embeddings()`), abundance-weighted RBF-MMD sample distances
  (`sample_mmd()`), structural bins by Markov clustering of a graph
  thresholded on similarity and alignment coverage (`protein_similarity()`,
  `bin_proteins()`, `bin_glom()`) and two quality-control plots
  (`plot_tsne()`, `plot_bin_similarity()`).
* The worked metaproteome vignette runs on cached embeddings.
