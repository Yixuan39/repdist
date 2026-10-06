# An embed_proteins()-shaped object from an embedding matrix. One shared
# sequence by default, so the coverage gate passes every pair.
as_proteins <- function(Z, seqs = NULL) {
    if (is.null(seqs)) {
        seqs <- rep("MKTAYIAKQRQISFVKSHFSRQ", nrow(Z))
        names(seqs) <- rownames(Z)
    }
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- Z
    x
}
