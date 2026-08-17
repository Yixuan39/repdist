# The Python side of embed_proteins() lives in a basilisk-managed Conda
# environment rather than in whatever Python reticulate happens to find. Every
# version is pinned, so the stack Bioconductor's build machines install is the
# stack a user gets, and it cannot be perturbed by the rest of their system.
#
# Bumping any pin here changes the embeddings, so treat it as a version bump of
# the ground metric, not as routine maintenance.
.repdist_env <- basilisk::BasiliskEnvironment(
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

# Reading NumPy archives should not install or load the much larger model stack.
.repdist_io_env <- basilisk::BasiliskEnvironment(
  envname = "repdist-io",
  pkgname = "repdist",
  packages = "python=3.11.11",
  pip = "numpy==1.26.4"
)
