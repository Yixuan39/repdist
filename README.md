# repdist

**Structure-aware beta diversity for metaproteomics.**

Sample distances that treat each sample as an abundance-weighted distribution over
proteins embedded in a protein language model -- the same functional form as weighted
UniFrac and Bray–Curtis, but with a ground metric defined for any amino-acid sequence,
including proteins with no tree position and no annotation. The default model, TM-Vec1,
is calibrated: cosine similarity between two embeddings is a predicted TM-score, not just
"distance in embedding space."

Two uses follow from the same embeddings: sample-level distances for ordination and
PERMANOVA (`weightedRep()`, returning `dist` objects for `vegan::adonis2()`), and
structural clustering of proteins into bins (`rollupStructural()`) for per-bin
differential-abundance testing (corncob, ALDEx2, ANCOM-BM, ...) when protein-level
resolution is too sparse.

## Install

```r
# install.packages("remotes")
remotes::install_github("Yixuan39/repdist")
```

Requires R ≥ 4.1. `embed_proteins()` needs `torch`, `transformers`, `sentencepiece`,
`safetensors` and `huggingface_hub` importable from whatever Python `reticulate` resolves
to (see `reticulate::py_config()`) -- install them yourself, e.g.

```
pip install torch transformers sentencepiece safetensors huggingface_hub
```

`embed_proteins()` errors with this exact command, naming whatever is actually missing,
rather than trying to provision an environment for you. The model itself is pulled from
HuggingFace on first use and cached under `~/.cache/huggingface`.

## Quick start

```r
library(repdist)

# counts: integer matrix, samples in rows, proteins in columns; rarefy first
#         (e.g. vegan::rrarefy()) if samples weren't collected at a common depth
# seqs:   named character vector, or a path to a FASTA file

Z <- embed_proteins(seqs)                       # TM-Vec1
D <- weightedRep(counts, Z, method = "mmd")

vegan::adonis2(D ~ condition, data = meta,
               permutations = permute::how(blocks = meta$subject, nperm = 9999))
```

## Functions

| function | purpose |
|---|---|
| `weightedRep(counts, embeddings, method)` | sample × sample distance; `method = "mmd"` (default) or `"rbf_mmd"` |
| `embed_proteins(seqs, model)` | one vector per protein; TM-Vec1 by default, or any HuggingFace repo via `model =` |
| `known_models()` | the registered head/backbone pairs `embed_proteins()` knows how to assemble |
| `read_fasta(path)` | FASTA → named character vector |
| `mmd_matrix(P, K)`, `rbf_kernel(C)` | the pieces, for building the kernel yourself |
| `rollupStructural(counts, embeddings, method, cluster_threshold)` | collapse proteins into structural bins (complete linkage at a TM-score cutoff) for per-bin DA testing |

`"mmd"` is the abundance-weighted mean representation,
$\lVert\sum_ip_{si}z_i-\sum_ip_{ti}z_i\rVert$ — one matmul, and blind by construction to
mass rearranged at a fixed mean. `"rbf_mmd"` covers that blind spot at roughly 10× the
cost. On the test study the two agree at Mantel $r=0.99$, so the cheap one is the method
there.

## Analyses

- [Simulation study](analysis/simulation_study.html) -- two groups differing only in
  protein function; Bray–Curtis sits at the exact null, `weightedRep()` recovers it.
- [Data curation](analysis/data_curation.html) -- published table → `phyloseq` object,
  validated against the source paper's PERMANOVA.
- [Case study](analysis/case_study.html) -- the same comparison on 111 real
  metaproteomes, plus an accession-shuffle control.
- [Case study, ESM-2 8M backbone](analysis/case_study_esm2_8m.html) -- same pipeline
  with an alternative, uncalibrated embedding model.
- [Protein embedding t-SNE](analysis/protein_tsne.html) -- visualizing the TM-Vec1
  embedding space.

## Citation

Method unpublished. Cite the test study as Blakeley-Ruiz JA *et al.* (2025) Dietary
protein source alters gut microbiota composition and function. *ISME J* 19(1):wraf048.

## License

MIT
