# An embed_proteins()-shaped object from an embedding matrix. One shared
# sequence by default, so the E-value gate passes every pair.
as_proteins <- function(Z, seqs = NULL) {
    if (is.null(seqs)) {
        seqs <- rep("MKTAYIAKQRQISFVKSHFSRQ", nrow(Z))
        names(seqs) <- rownames(Z)
    }
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- Z
    x
}

# Cosine similarity with pairs below `t` set to 0: a sparse graph for testing
# bin_proteins() without sequences to align.
cosine_graph <- function(Z, t) {
    Z <- Z / sqrt(rowSums(Z^2))
    S <- pmax(pmin(tcrossprod(Z), 1), 0)
    S[S < t] <- 0
    diag(S) <- 1
    S
}
