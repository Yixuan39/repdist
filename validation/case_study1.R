# Case study 1: the full mouse dietary-protein metaproteome (Blakeley-Ruiz
# et al. 2025) through the current repdist API. Writes tables and figures to
# validation/case_study1/; case_study1_report.Rmd renders them.
#
# Run from the source root (or set CS1_ROOT) after installing repdist;
# REPDIST_LIB may point at a private library. Needs phyloseq, vegan, permute,
# DESeq2, Rtsne, ggplot2, and python3 with numpy and h5py for the file-reader
# round trip. The coverage alignment and the dense MCL on the largest
# components take most of the time.
if (nzchar(Sys.getenv("REPDIST_LIB"))) {
    .libPaths(c(Sys.getenv("REPDIST_LIB"), .libPaths()))
}
suppressMessages({
    library(repdist)
    library(phyloseq)
    library(ggplot2)
})
root <- Sys.getenv("CS1_ROOT", ".") # the source root
path <- file.path(root, "data", "test_study")
out <- Sys.getenv("CS1_OUT", file.path(root, "validation", "case_study1"))
dir.create(out, showWarnings = FALSE, recursive = TRUE)
res <- list(started = format(Sys.time()),
    workers = BiocParallel::bpnworkers(BiocParallel::bpparam()))
n_cores <- res$workers # for the alignments
perf <- list()
timed <- function(step, expr) {
    invisible(gc(reset = TRUE))
    t <- system.time(value <- expr)[["elapsed"]]
    # peak R heap as max-used cells x 8 bytes (vector cells dominate)
    perf[[step]] <<- data.frame(step = step, seconds = round(t, 1),
        peak_R_GB = round(sum(gc()[, "max used"]) * 8 / 2^30, 1))
    message(step, ": ", round(t, 1), " s")
    value
}
# Slow steps are cached in out/cache, so a rerun resumes past them.
dir.create(file.path(out, "cache"), showWarnings = FALSE)
cache_step <- function(step, expr) {
    f <- file.path(out, "cache", paste0(gsub("[^A-Za-z0-9]+", "_", step),
        ".rds"))
    if (file.exists(f)) {
        hit <- readRDS(f)
        perf[[step]] <<- hit$perf
        message(step, ": cached")
        return(hit$value)
    }
    value <- timed(step, expr)
    saveRDS(list(value = value, perf = perf[[step]]), f)
    value
}
# MCL failures are results, not a reason to stop
attempt <- function(expr) {
    tryCatch(expr, error = function(e) conditionMessage(e))
}
fig <- function(name, p, width = 7, height = 4) {
    ggsave(file.path(out, paste0(name, ".png")), p, width = width,
        height = height, dpi = 150)
}

# 1. Data: rarefy, drop baseline, 3000 aa ceiling (as on the gh-pages site) --
ps <- readRDS(file.path(path, "ps.rds"))
set.seed(1)
psd <- rarefy_even_depth(ps, sample.size = min(sample_sums(ps)), rngseed = 1,
    replace = FALSE, verbose = FALSE)
psd <- prune_samples(sample_data(psd)$phase != "baseline", psd)
psd <- prune_taxa(taxa_sums(psd) > 0, psd)
n_rarefied <- ntaxa(psd)
too_long <- Biostrings::width(refseq(psd)) > 3000L
long_share <- sum(taxa_sums(psd)[too_long]) / sum(taxa_sums(psd))
psd <- prune_taxa(!too_long, psd)
keep <- taxa_names(psd)
# CS1_SUBSET=n runs a quick smoke test on n random proteins
n_sub <- as.integer(Sys.getenv("CS1_SUBSET", "0"))
if (n_sub > 0) {
    set.seed(3)
    keep <- sort(sample(keep, n_sub))
    psd <- prune_taxa(keep, psd)
}

psu <- prune_samples(sample_data(ps)$phase != "baseline", ps)
psu <- prune_taxa(taxa_names(psu) %in% keep, psu)
counts <- Matrix::Matrix(t(as(otu_table(psd), "matrix")), sparse = TRUE)
counts_raw <- t(as(otu_table(psu), "matrix"))[rownames(counts), keep]
md <- data.frame(sample_data(psd))[rownames(counts), ]
md$diet <- factor(md$diet)
md$dose <- factor(md$dose)
md$subject <- factor(md$subject)
res$data <- c(samples = nrow(counts), proteins_after_rarefaction = n_rarefied,
    removed_over_3000aa = sum(too_long), proteins = length(keep),
    removed_abundance_pct = round(100 * long_share, 3))
res$diets <- table(md$diet)

# Annotations: join on Protein Identifier, fall back to the species-code
# column, and keep a match only where the residue count agrees.
t6 <- read.delim(file.path(path, "annotated",
    "ExtendedDataTable6AnnotatedMicrobialMetaproteomes.txt"), skip = 2,
    quote = "\"", check.names = FALSE, colClasses = "character",
    comment.char = "", fileEncoding = "latin1")
len <- Biostrings::width(refseq(psd)[keep])
r_id <- match(keep, t6[["Protein Identifier"]])
r_code <- match(keep, t6[["Protein Identifier with Species Code"]])
agrees <- function(r) !is.na(r) & suppressWarnings(as.integer(t6$AA[r])) == len
row <- ifelse(agrees(r_id), r_id, ifelse(agrees(r_code), r_code, NA))
annotation <- t6[["Consensus annotation"]][row]
annotation[!nzchar(trimws(annotation)) | tolower(annotation) == "unknown"] <- NA
taxa <- data.frame(row.names = keep,
    genus = as(tax_table(psd), "matrix")[keep, "Genus"],
    annotation = annotation)
# Function from KOfamScan (prokaryote.hal against proteins.fasta, written by
# the run in validation/case_study1/): each protein's top-scoring KO above its
# family threshold, and that KO's KEGG BRITE level-B category (ko00001 from
# the KEGG REST API). A KO in several categories takes the one holding most of
# its pathways, pathway classes before BRITE protein families, and ties go to
# the category KEGG lists first.
ko_file <- file.path(out, "kofamscan.tsv")
ko <- category <- NULL
if (file.exists(ko_file)) {
    hits <- read.delim(ko_file, header = FALSE, comment.char = "#",
        quote = "", colClasses = "character", col.names = c("hit",
            "protein", "ko", "threshold", "score", "evalue", "definition"))
    hits <- hits[hits$hit == "*" & hits$protein %in% keep, ]
    hits <- hits[order(hits$protein, -as.numeric(hits$score)), ]
    hits <- hits[!duplicated(hits$protein), ]
    ko <- stats::setNames(hits$ko, hits$protein)

    tree <- jsonlite::fromJSON(file.path(out, "ko00001.json"),
        simplifyVector = FALSE)
    brite <- list()
    for (a in tree$children) for (b in a$children) for (cc in b$children) {
        if (length(cc$children)) {
            brite[[length(brite) + 1L]] <- data.frame(A = a$name,
                B = b$name, ko = vapply(cc$children,
                    function(d) sub(" .*", "", d$name), ""))
        }
    }
    brite <- do.call(rbind, brite)
    brite$priority <- ifelse(grepl("^091[0-4]", brite$A), 1L,
        ifelse(grepl("^0918", brite$A), 2L, 3L))
    n_b <- aggregate(list(n = brite$ko), brite[c("ko", "B", "priority")],
        length)
    n_b <- n_b[order(n_b$ko, n_b$priority, -n_b$n, n_b$B), ]
    n_b <- n_b[!duplicated(n_b$ko), ]
    category <- stats::setNames(sub("^[0-9]+ ", "", n_b$B[match(ko, n_b$ko)]),
        names(ko))
    # how often the category was a choice among several
    top <- brite[brite$priority == ave(brite$priority, brite$ko, FUN = min), ]
    n_cat <- tapply(top$B, top$ko, function(v) length(unique(v)))
    res$kofam <- c(proteins = length(keep), `with a KO` = length(ko),
        `distinct KOs` = length(unique(ko)),
        `with a BRITE category` = sum(!is.na(category)),
        `KO in more than one category` = sum(n_cat[ko] > 1, na.rm = TRUE))
    res$categories <- utils::head(sort(table(category), decreasing = TRUE), 15)
}
res$annotation <- c(
    `matched by Protein Identifier (old join)` = sum(!is.na(r_id)),
    `matched, length agrees` = sum(agrees(r_id)),
    `rescued by species-code column` = sum(!agrees(r_id) & agrees(r_code)),
    unmatched = sum(is.na(row)),
    `with a consensus annotation` = sum(!is.na(annotation)))

# 2. The protein object, and the file readers at full size -----------------
cached <- readRDS(file.path(path, "embeddings_tmvec1_large.rds"))
x <- refseq(psd)[keep]
stopifnot(identical(as.character(x), cached$seq[keep]))
S4Vectors::mcols(x)$embedding <- cached$Z[keep, ]
Zn <- cached$Z[keep, ] / sqrt(rowSums(cached$Z[keep, ]^2))
rm(cached)

tmp <- file.path(tempdir(), "readers")
dir.create(tmp, showWarnings = FALSE)
Biostrings::writeXStringSet(x, file.path(tmp, "proteins.fasta"))
writeLines(names(x), file.path(tmp, "ids.txt"))
writeLines(as.character(x), file.path(tmp, "seqs.txt"))
writeBin(as.vector(t(Zn)), file.path(tmp, "Z.f32"), size = 4)
writeLines(c(
    "import sys, numpy as np, h5py",
    "d = sys.argv[1]",
    "ids = open(d + '/ids.txt').read().split()",
    "seqs = open(d + '/seqs.txt').read().split()",
    "Z = np.fromfile(d + '/Z.f32', dtype='<f4').reshape(len(ids), -1)",
    "# tmvec build-db layout",
    "np.savez_compressed(d + '/db.npz', headers=tuple(ids), embeddings=Z,",
    "    input_fasta=np.array(d + '/proteins.fasta'),",
    "    tm_vec_weights=np.array(''), protrans_model_path=np.array(''))",
    "# tmvec embed layout, one pooled vector per protein",
    "with h5py.File(d + '/emb.h5py', 'w') as f:",
    "    for h, s, e in zip(ids, seqs, Z):",
    "        g = f.create_group(h)",
    "        g.create_dataset('seq', data=s)",
    "        g.create_dataset('emb', data=e)"
), file.path(tmp, "write.py"))
stopifnot(system2("python3", c(file.path(tmp, "write.py"), tmp)) == 0)
x_npz <- timed("read_embeddings, tmvec .npz",
    read_embeddings(file.path(tmp, "db.npz")))
x_h5 <- timed("read_embeddings, HDF5",
    read_embeddings(file.path(tmp, "emb.h5py")))
check_reader <- function(y) {
    E <- S4Vectors::mcols(y)$embedding
    data.frame(proteins = length(y), same_names = setequal(names(y), names(x)),
        same_sequences = identical(as.character(y)[names(x)], as.character(x)),
        max_abs_embedding_diff = signif(max(abs(E[names(x), ] - Zn)), 2))
}
res$readers <- rbind(`.npz (tmvec build-db)` = check_reader(x_npz),
    `HDF5 (tmvec embed layout)` = check_reader(x_h5))
res$file_MB <- round(file.size(file.path(tmp, c("db.npz", "emb.h5py"))) / 2^20)
rm(x_npz, x_h5)

# 3. Beta diversity ---------------------------------------------------------
D <- list(`Bray-Curtis` = vegan::vegdist(as.matrix(counts), "bray"),
    MMD = cache_step("sample_mmd", sample_mmd(counts, x)))
res$sigma <- attr(D$MMD, "sigma")
perm <- permute::how(blocks = md$subject, nperm = 999)
res$permanova <- do.call(rbind, lapply(names(D), function(n) {
    set.seed(1)
    a <- vegan::adonis2(D[[n]] ~ dose + sex + age_week + diet, data = md,
        by = "terms", permutations = perm)
    data.frame(distance = n, R2 = round(a["diet", "R2"], 3),
        F = round(a["diet", "F"], 2), p = a["diet", "Pr(>F)"])
}))
ord <- do.call(rbind, lapply(names(D), function(n) {
    e <- cmdscale(D[[n]], k = 2, eig = TRUE)
    rel <- 100 * e$eig / sum(e$eig[e$eig > 0])
    data.frame(Axis1 = e$points[, 1], Axis2 = e$points[, 2],
        diet = md[rownames(e$points), "diet"],
        panel = sprintf("%s (PCo1 %.1f%%, PCo2 %.1f%%)", n, rel[1], rel[2]))
}))
fig("pcoa", ggplot(ord, aes(Axis1, Axis2, colour = diet)) +
    geom_point(size = 1.6, alpha = 0.85) +
    stat_ellipse(type = "t", linewidth = 0.3) +
    facet_wrap(~panel, scales = "free") +
    labs(x = "PCo1", y = "PCo2") + theme_bw() +
    theme(legend.position = "bottom"), width = 8, height = 4.2)
rm(D)

tsne <- cache_step("plot_tsne", {
    set.seed(1)
    p <- plot_tsne(x)
    fig("tsne_plain", p, width = 6, height = 5.5)
    p$data
})
if (!is.null(category)) {
    cache_step("plot_tsne, coloured by KEGG category", {
        set.seed(1)
        fig("tsne_ko", plot_tsne(x, color = category), width = 8.5,
            height = 6)
        TRUE
    })
}

# 4. The similarity graph ---------------------------------------------------
# Pairs at each predicted TM-score floor, before any gate
S_all <- timed("protein_similarity, no gates", protein_similarity(x,
    min_sim = 0, min_coverage = 0))
u <- S_all[upper.tri(S_all)]
grid <- c(0.5, 0.6, 0.7, 0.75, 0.8, 0.9)
res$pairs <- data.frame(min_sim = grid,
    pairs = vapply(grid, function(t) sum(u >= t), numeric(1)))
rm(S_all, u)
invisible(gc())

# The package defaults, TM >= 0.7 and coverage 0.5; timed together with the
# sparse copy kept in the cache
S <- as.matrix(cache_step("protein_similarity, defaults",
    Matrix::Matrix(protein_similarity(x, n_cores = n_cores), sparse = TRUE)))
b <- cache_step("bin_proteins, coverage gate", bin_proteins(S, inflation = 2))

summarise_bins <- function(b, setting, t) {
    size <- lengths(b$clusters)
    data.frame(setting = setting, t(b$graph), bins = length(size),
        singletons = sum(size == 1), largest_bin = max(size), mcl_seconds = t)
}
res$bins <- summarise_bins(b, "similarity + coverage",
    perf[["bin_proteins, coverage gate"]]$seconds)
sizes <- lengths(b$clusters)
res$bin_sizes <- table(cut(sizes, c(0, 1, 2, 5, 10, Inf),
    labels = c("1", "2", "3-5", "6-10", ">10")))

# 5. Differential abundance on the bins ------------------------------------
bin_cnt <- bin_glom(counts_raw, b)
dds <- timed("DESeq2 LRT + Wald + shrinkage", {
    dds <- DESeq2::DESeqDataSetFromMatrix(t(bin_cnt), md,
        design = ~ subject + diet)
    dds <- DESeq2::DESeq(dds, test = "LRT", reduced = ~subject,
        sfType = "poscounts", quiet = TRUE)
    dds
})
res_lrt <- DESeq2::results(dds)
lev <- setdiff(levels(md$diet), "Casein")
ddsw <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
shrunk <- lapply(lev, function(l) DESeq2::lfcShrink(ddsw,
    coef = paste0("diet_", l, "_vs_Casein"), type = "normal", quiet = TRUE))
names(shrunk) <- lev
fc <- vapply(shrunk, function(o) o$log2FoldChange, numeric(nrow(res_lrt)))
rownames(fc) <- rownames(res_lrt)
sig <- rownames(res_lrt)[!is.na(res_lrt$padj) & res_lrt$padj <= 0.05]
best <- lev[apply(abs(fc[sig, , drop = FALSE]), 1, which.max)]
top <- data.frame(bin = sig, members = sizes[sig],
    padj = formatC(res_lrt[sig, "padj"], digits = 2, format = "e"),
    diet = best, log2FC = round(fc[cbind(sig, best)], 2))
top <- top[order(-abs(top$log2FC)), ]
res$deseq <- c(bins_tested = nrow(res_lrt), diet_significant = length(sig),
    multi_member_significant = sum(top$members > 1))
res$top <- head(top, 8)
volc <- do.call(rbind, lapply(lev, function(l) {
    o <- shrunk[[l]]
    data.frame(diet = l, log2FoldChange = o$log2FoldChange,
        pvalue = o$pvalue, reject = !is.na(o$padj) & o$padj <= 0.05)
}))
fig("volcano", ggplot(volc, aes(log2FoldChange, -log10(pvalue),
    colour = reject)) +
    geom_point(alpha = 0.6, size = 0.8) + facet_wrap(~diet, nrow = 1) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"),
        guide = "none") +
    labs(x = "log2 fold change vs Casein", y = expression(-log[10](p))) +
    theme_bw(), width = 9, height = 3.2)

# 6. The most responsive annotated bins ------------------------------------
n_ann <- vapply(b$clusters[top$bin],
    function(m) sum(!is.na(taxa[m, "annotation"])), integer(1))
picks <- head(top[top$members > 1 & n_ann >= 2, ], 3)
res$picks <- picks
res$members <- do.call(rbind, lapply(seq_len(nrow(picks)), function(i) {
    m <- b$clusters[[picks$bin[i]]]
    data.frame(bin = picks$bin[i], protein = m, genus = taxa[m, "genus"],
        annotation = substr(taxa[m, "annotation"], 1, 60),
        KO = if (is.null(ko)) NA else unname(ko[m]),
        category = if (is.null(category)) NA else unname(category[m]),
        mean_Casein = round(colMeans(counts_raw[md$diet == "Casein", m,
            drop = FALSE]), 1),
        mean_top = round(colMeans(counts_raw[md$diet == picks$diet[i], m,
            drop = FALSE]), 1), row.names = NULL)
}))
if (nrow(picks)) {
    fig("bin", plot_bin_similarity(b, picks$bin[1], S, x), width = 6, height = 5)
    hl <- merge(tsne, res$members[c("protein", "bin")])
    fig("tsne", ggplot(tsne, aes(tsne1, tsne2)) +
        geom_point(colour = "grey70", alpha = 0.3, size = 0.3) +
        geom_point(data = hl, aes(colour = bin), size = 2) +
        labs(x = "t-SNE 1", y = "t-SNE 2", colour = "responsive bin") +
        theme_minimal() + theme(axis.text = element_blank()),
        width = 6.5, height = 5.5)
}

# Last, the same floor without the coverage gate: its largest component is
# the biggest dense MCL problem here.
S_u <- protein_similarity(x, min_coverage = 0)
rm(S)
invisible(gc())
b_u <- attempt(timed("bin_proteins, similarity only", bin_proteins(S_u)))
rm(S_u)
res$bins <- rbind(if (is.list(b_u)) summarise_bins(b_u, "similarity only",
    perf[["bin_proteins, similarity only"]]$seconds), res$bins)
res$similarity_only_mcl <- if (is.character(b_u)) b_u else "ok"

res$perf <- do.call(rbind, perf)
res$finished <- format(Sys.time())
res$session <- utils::capture.output(utils::sessionInfo())
saveRDS(res, file.path(out, "results.rds"))
message("done")
