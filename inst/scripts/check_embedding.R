# Optional real-model check, outside R CMD check because weights need downloading.
# Run after installing the package: Rscript inst/scripts/check_embedding.R
# Set REPDIST_CHECK_MODEL to another registered model to validate that preset.
library(repdist)
model <- Sys.getenv("REPDIST_CHECK_MODEL", "esm2-8m")
s <- c(
    p1 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ",
    p2 = "ACDEFGHIKL", p3 = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
)
z1 <- embed_proteins(s, model = model, device = "cpu", batch_size = 1)
z3 <- embed_proteins(s, model = model, device = "cpu", batch_size = 3)
stopifnot(
    identical(rownames(z1), names(s)), nrow(z1) == 3L,
    all(is.finite(z1)), max(abs(rowSums(z1^2) - 1)) < 1e-6,
    max(abs(z1[1, ] - z1[3, ])) < 1e-6,
    max(abs(z1 - z3)) < 1e-5
)
cat(model, "CPU check passed; max batch difference:", max(abs(z1 - z3)), "\n")
