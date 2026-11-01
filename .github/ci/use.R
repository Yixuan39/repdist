# A new user's first session: embed, then the README quick start.
library(repdist)
timed <- function(label, expr) {
    t0 <- Sys.time()
    out <- force(expr)
    message(label, ": ", format(round(Sys.time() - t0, 1)))
    out
}
seqs <- c(
    ubq = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
    ubq2 = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRG",
    trx = "MSDKIIHLTDDSFDTDVLKADGAILVDFWAEWCGPCKMIAPILDEIADEYQGKLTVAKLNIDQNPGTAPKYGIRGIPTLLLFKNGEVAASKVGALSKGQLKEFLDANLA",
    trx2 = "MVKQIESKTAFQEALDAAGDKLVVVDFSATWCGPCKMIKPFFHSLSEKYSNVIFLEVDVDDCQDVASECEVKCMPTFQFFKKGQKVGEFSGANKEKLEATINELV",
    lyz = "KVFGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCRL",
    gfp = "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

small <- timed("first embed_proteins(esm2-8m), provisions Python",
    embed_proteins(seqs, model = "esm2-8m", device = "cpu"))
timed("second embed_proteins(esm2-8m), env reused",
    embed_proteins(seqs, model = "esm2-8m", device = "cpu"))
proteins <- timed("embed_proteins() default model",
    embed_proteins(seqs, device = "cpu"))
Z <- S4Vectors::mcols(proteins)$embedding
print(dim(Z))
stopifnot(identical(rownames(Z), names(seqs)),
    isTRUE(all.equal(unname(rowSums(Z^2)), rep(1, 6), tolerance = 1e-5)))
print(round(tcrossprod(Z), 2))

set.seed(1)
counts <- matrix(rpois(8 * 6, 20), 8, dimnames = list(paste0("s", 1:8), names(seqs)))
meta <- data.frame(condition = rep(c("a", "b"), 4), row.names = rownames(counts))
D <- sample_mmd(counts, proteins)
print(vegan::adonis2(D ~ condition, data = meta))
S <- protein_similarity(proteins, min_sim = 0.7, min_coverage = 0.5)
bins <- bin_proteins(S, inflation = 2)
print(bins)
print(bin_glom(counts, bins))
message("OK")
