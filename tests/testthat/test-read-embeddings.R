# The fixtures were written the way tmvec 1.0.5 writes them:
#   np.savez_compressed("tmvec.npz", headers=("p3", "p1", "p2"),
#       embeddings=np.array([[1, 1, 1, 1], [3, 4, 0, 0], [0, 0, 5, 0]],
#       dtype=np.float32), input_fasta=np.array("/nonexistent/..."), ...)
#   h5py: one group per protein with datasets "seq" and "emb" (residues x d)
npz <- test_path("fixtures", "tmvec.npz")

test_that("a tmvec .npz database is read and matched to its sequences", {
    seqs <- c(p1 = "MKV", p2 = "MAL", p3 = "MKL", extra = "MMM")
    x <- read_embeddings(npz, seqs)
    expect_s4_class(x, "AAStringSet")
    expect_identical(names(x), c("p3", "p1", "p2"))
    expect_identical(as.character(x), seqs[c("p3", "p1", "p2")])
    expect_equal(unname(S4Vectors::mcols(x)$embedding), rbind(
        c(0.5, 0.5, 0.5, 0.5), c(0.6, 0.8, 0, 0), c(0, 0, 1, 0)
    ), tolerance = 1e-7)
    fa <- tempfile(fileext = ".fasta")
    writeLines(paste0(">", names(seqs), "\n", seqs), fa)
    expect_identical(read_embeddings(npz, fa),
        x)
    expect_error(read_embeddings(npz),
        "missing; pass `seqs`")
    expect_error(read_embeddings(npz,
        seqs[1:2]), "no sequence, e.g. p3")
})

test_that("a column-major .npy array is read in its stored order", {
    # np.asfortranarray([[3, 4, 0], [0, 0, 2]]), so fortran_order is True
    x <- read_embeddings(test_path("fixtures", "tmvec_fortran.npz"),
        c(p1 = "MKV", p2 = "MAL"))
    expect_equal(unname(S4Vectors::mcols(x)$embedding),
        rbind(c(0.6, 0.8, 0), c(0, 0, 1)), tolerance = 1e-7)
})

test_that("a tmvec embed HDF5 file is averaged over residues", {
    x <- read_embeddings(test_path("fixtures", "tmvec.h5py"))
    expect_identical(as.character(x), c(p1 = "MK", p2 = "MKV"))
    expect_equal(unname(S4Vectors::mcols(x)$embedding),
        rbind(c(1, 0, 0, 0), c(0, 1, 0, 0)))
})

test_that("the shipped example files run from disk to sample_mmd", {
    ex <- function(f) system.file("extdata", f, package = "repdist")
    x <- read_embeddings(ex("metaproteomics.npz"), ex("metaproteomics.fasta"))
    counts <- as.matrix(read.delim(ex("metaproteomics_counts.tsv"),
        row.names = 1, check.names = FALSE))
    samples <- read.delim(ex("metaproteomics_samples.tsv"))
    expect_identical(rownames(counts), samples$sample)
    D <- sample_mmd(counts, x)
    expect_identical(labels(D), samples$sample)
})
