# The Python side of embed_proteins() lives in a basilisk-managed Conda
# environment rather than in whatever Python reticulate happens to find. Every
# version is pinned, so the stack Bioconductor's build machines install is the
# stack a user gets, and it cannot be perturbed by the rest of their system.
#
# Bumping any pin here changes the embeddings, so treat it as a version bump of
# the ground metric, not as routine maintenance.
#
# This file is sourced standalone by configureBasiliskEnv() from ./configure, so
# it must not assume the rest of the package is loaded. `repdist_env` carries no
# leading dot on purpose: configureBasiliskEnv() collects environments with
# ls(), which skips dot-prefixed names, and would silently build nothing.
#
# One environment, not one per model: every entry in inst/extdata/models.json
# runs on these pins. Split only when a model needs pins that genuinely cannot
# coexist with these -- configure builds every environment declared here, so an
# extra one costs another multi-GB torch on the Bioconductor builders.
repdist_env <- basilisk::BasiliskEnvironment(
  envname = "repdist-embed",
  pkgname = "repdist",
  packages = "python=3.11.11",
  pip = c(
    "torch==2.5.1",
    "transformers==4.46.3",
    "safetensors==0.4.5",
    "huggingface_hub==0.26.5",
    "numpy==1.26.4",
    # Rostlab/prot_t5_xl_uniref50 ships only a SentencePiece `spiece.model`, so
    # transformers has to convert it into a fast tokenizer at load time and
    # needs both of these to do it. Without protobuf the conversion silently
    # falls back to a TikToken extractor and dies on a file that is not one.
    "sentencepiece==0.2.0",
    "protobuf==5.29.1"
  )
)
