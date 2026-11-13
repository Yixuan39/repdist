# Build curated/simulation_kofam.rds: a small universe of Bacteria/Archaea
# Swiss-Prot proteins whose KEGG Orthology (KO) label comes from official
# KEGG REST lookups at build time. The treatment group is every KO that
# KEGG's ko-rn link table assigns to one shared reaction (interchangeable
# alternative enzymes for that reaction, not required to be structurally
# similar); background KOs carry effect size 0 by design. See
# KOFAM_SIM_PLAN.md for the full algorithm (sections 2-6) this implements.
#
# Run from the root of this worktree:
#   Rscript data-raw/simulation_kofam.R          # sizes the full download, then stops
#   KOFAM_SMOKE=1 Rscript data-raw/simulation_kofam.R
#
# Needs network access to rest.kegg.jp (rate-limited to <= 3 calls/s, per
# KEGG's academic API policy) and rest.uniprot.org, and the `mmseqs` binary
# on PATH. KEGG responses are cached, untracked, in data/kegg/ (gitignored);
# KEGG is not an open-licensed database (kofam-official-sources.md), so raw
# KEGG files and the complete ko-rn/BRITE tables are never committed. In
# KOFAM_SMOKE=1 mode, embeddings are random unit vectors, not real TM-Vec
# (the plan explicitly allows this for smoke runs; it is not a result).
#
# Mapping route (Section 2.1): Swiss-Prot first, not KO first. An earlier
# version of this script fetched each KO's genes from KEGG and mapped them
# to UniProt; most of a KO's genes across all organisms have no reviewed
# Swiss-Prot entry (they resolve to unreviewed TrEMBL accessions), so random
# per-KO sampling almost always came up empty. This version instead pulls
# reviewed Bacteria/Archaea Swiss-Prot entries with any KEGG cross-reference
# and maps their genes to KO with the official `link/ko/<genes>` endpoint, so
# every candidate/background check is a lookup in one local, already-reviewed
# table instead of a new pair of network round trips.

smoke <- identical(Sys.getenv("KOFAM_SMOKE", "0"), "1")

out <- file.path("curated", "simulation_kofam.rds")
kegg_dir <- file.path("data", "kegg")
dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)

# Parameters fixed before looking at any KEGG or UniProt result ---------------
length_range <- c(50L, 1000L)          # the TM-Vec head's training range
min_proteins_per_ko <- 3L              # Step E: a KO needs this many qualifying proteins
n_per_ko <- if (smoke) 10L else 20L    # Step 5 cap per KO, proposed default
global_cap <- if (smoke) 60L else 3000L # proposed hard cap on the whole universe
max_probe <- if (smoke) 60L else 25L   # candidate reactions actually inspected (Steps D-F);
                                        # cheap in smoke, see Step 7 (a local table lookup
                                        # unless a candidate clears the protein-count bar)
k_target <- 1L                         # primary treatment groups to accept (budget: see notes)
smoke_sample_size <- 15000L            # Swiss-Prot entries pulled in smoke mode (see Step 5)
seed <- 20261112L                      # date of KOFAM_SIM_PLAN.md

set.seed(seed)

# KEGG REST, rate-limited to <= 3 calls/s --------------------------------------
.kegg_state <- new.env()
.kegg_state$last <- 0
kegg_wait <- function() {
    gap <- (.kegg_state$last + 1 / 3) - as.numeric(Sys.time())
    if (gap > 0) Sys.sleep(gap)
    .kegg_state$last <- as.numeric(Sys.time())
}
# `path` is the part after https://rest.kegg.jp/. Returns NULL (not an error)
# on a non-200 or empty response so one bad id does not abort the whole build;
# the required bulk endpoints below check the result and stop if it is NULL.
kegg_fetch <- function(path) {
    kegg_wait()
    res <- tryCatch(curl::curl_fetch_memory(paste0("https://rest.kegg.jp/", path)),
        error = function(e) NULL)
    if (is.null(res) || res$status_code != 200L || length(res$content) == 0L) {
        return(NULL)
    }
    rawToChar(res$content)
}
kegg_lines <- function(path) {
    txt <- kegg_fetch(path)
    if (is.null(txt)) return(character(0))
    l <- strsplit(txt, "\n")[[1]]
    l[nzchar(l)]
}
# Cache a one-time bulk download under data/kegg/; never committed (.gitignore).
kegg_cached <- function(path, file) {
    file <- file.path(kegg_dir, file)
    if (file.exists(file)) return(readLines(file, warn = FALSE))
    txt <- kegg_fetch(path)
    if (is.null(txt)) stop("KEGG fetch failed (not HTTP 200 or empty): ", path, call. = FALSE)
    writeLines(txt, file)
    strsplit(txt, "\n")[[1]]
}
# Batched `get`, up to 10 ids per call (KEGG API limit); returns one string.
kegg_get_batch <- function(prefix, ids) {
    ids <- unique(ids)
    if (!length(ids)) return("")
    chunks <- split(ids, ceiling(seq_along(ids) / 10L))
    txt <- vapply(chunks, function(g) {
        t <- kegg_fetch(paste0("get/", paste0(prefix, ":", g, collapse = "+")))
        if (is.null(t)) "" else t
    }, "")
    paste(txt, collapse = "\n")
}
# Batched `link/ko/<genes>` (Section 2.1 item 2): KEGG gene id -> KO. 150 ids
# per call is a measured safe batch size (300 hit an HTTP 403, likely a URL
# length limit rather than a KEGG API limit); well under the per-call cost of
# going one gene at a time at the same <= 3 calls/s rate.
kegg_genes_to_ko <- function(genes) {
    genes <- unique(genes)
    out <- lapply(split(genes, ceiling(seq_along(genes) / 150L)), function(g) {
        l <- kegg_lines(paste0("link/ko/", paste(g, collapse = "+")))
        if (!length(l)) return(NULL)
        data.frame(gene = sub("\t.*", "", l), ko = sub(".*\tko:", "", l),
            stringsAsFactors = FALSE)
    })
    do.call(rbind, out)
}

# MMseqs2 all-vs-all + clustering, on a small named-character-vector universe ----
# Same settings as data-raw/simulation_mmseqs2.R (Step 4.1): -s 7.5, 80%
# coverage of both sequences; 50%/30% identity clusters and E <= 1e-3 hits.
mmseqs_run <- function(seqs) {
    tmp <- tempfile("mmseqs")
    dir.create(tmp)
    fasta <- file.path(tmp, "universe.fasta")
    Biostrings::writeXStringSet(Biostrings::AAStringSet(seqs), fasta)
    run <- function(...) {
        args <- c(...)
        stopifnot(system2("mmseqs", args, stdout = FALSE, stderr = FALSE) == 0L)
        paste("mmseqs", paste(args, collapse = " "))
    }
    cluster <- function(id, tag) {
        prefix <- file.path(tmp, tag)
        cmd <- run("easy-cluster", fasta, prefix, file.path(tmp, paste0("w", tag)),
            "--min-seq-id", id, "-c", 0.8, "--cov-mode", 0, "--cluster-mode", 0, "-s", 7.5)
        m <- utils::read.delim(paste0(prefix, "_cluster.tsv"), header = FALSE,
            col.names = c("representative", "member"))
        list(membership = setNames(m$representative, m$member),
            command = gsub(tmp, "<tmp>", cmd, fixed = TRUE))
    }
    c50 <- cluster(0.5, "c50")
    c30 <- cluster(0.3, "c30")
    m8 <- file.path(tmp, "hits.m8")
    search <- run("easy-search", fasta, fasta, m8, file.path(tmp, "wsearch"),
        "-s", 7.5, "-e", "1e-3", "--max-seqs", 10000,
        "--format-output", "query,target,fident,alnlen,qcov,tcov,evalue")
    h <- utils::read.delim(m8, header = FALSE, col.names = c("query", "target",
        "fident", "alnlen", "qcov", "tcov", "evalue"))
    h <- h[h$query != h$target, ]
    if (nrow(h)) {
        h$p1 <- pmin(h$query, h$target)
        h$p2 <- pmax(h$query, h$target)
        h <- h[order(-h$fident), ]
        hits <- h[!duplicated(h[c("p1", "p2")]), c("p1", "p2", "fident", "alnlen",
            "qcov", "tcov", "evalue")]
    } else {
        hits <- h[, c("query", "target", "fident", "alnlen", "qcov", "tcov", "evalue")]
        names(hits)[1:2] <- c("p1", "p2")
    }
    rownames(hits) <- NULL
    version <- system2("mmseqs", "version", stdout = TRUE)
    unlink(tmp, recursive = TRUE)
    list(c50 = c50$membership, c30 = c30$membership, hits = hits,
        provenance = list(cluster_50 = c50$command, cluster_30 = c30$command,
            search = gsub(tmp, "<tmp>", search, fixed = TRUE), version = version))
}

message("1. KEGG info endpoints")
# Required per plan Section 2.2: validate the service and record it.
info <- list(kegg = kegg_fetch("info/kegg"), ko = kegg_fetch("info/ko"),
    reaction = kegg_fetch("info/reaction"))
if (any(vapply(info, is.null, TRUE))) {
    stop("KEGG info/* endpoints did not return HTTP 200; KEGG REST may be ",
        "unavailable.", call. = FALSE)
}

message("2. KO-reaction links (Step A)")
# The literal "rn" alias in the plan's documented URL (link/rn/ko) returned
# HTTP 400 when checked live; the "reaction" name works and returns the same
# two-column ko/reaction mapping. Try the documented alias first, in case a
# later KEGG release restores it, and record which one actually worked.
ko_rn_source <- NULL
ko_rn_lines <- NULL
for (name in c("rn", "reaction")) {
    f <- paste0("ko_rn_", name, ".tsv")
    cached <- file.path(kegg_dir, f)
    if (file.exists(cached)) {
        ko_rn_lines <- readLines(cached, warn = FALSE)
        ko_rn_source <- paste0("link/", name, "/ko")
        break
    }
    txt <- kegg_fetch(paste0("link/", name, "/ko"))
    if (!is.null(txt) && nzchar(trimws(txt))) {
        writeLines(txt, cached)
        ko_rn_lines <- strsplit(txt, "\n")[[1]]
        ko_rn_source <- paste0("link/", name, "/ko")
        break
    }
}
stopifnot(!is.null(ko_rn_lines))
ko_rn_lines <- ko_rn_lines[nzchar(ko_rn_lines)]
ko_rn_all <- read.delim(text = paste(ko_rn_lines, collapse = "\n"), header = FALSE,
    col.names = c("ko", "reaction"), stringsAsFactors = FALSE)
well_formed <- grepl("^ko:K[0-9]+$", ko_rn_all$ko) & grepl("^rn:R[0-9]+$", ko_rn_all$reaction)
steps <- c(ko_rn_rows = nrow(ko_rn_all), ko_rn_well_formed = sum(well_formed))
ko_rn_all <- ko_rn_all[well_formed, ]
ko_rn_all$ko <- sub("^ko:", "", ko_rn_all$ko)
ko_rn_all$reaction <- sub("^rn:", "", ko_rn_all$reaction)
ko_rn_all <- unique(ko_rn_all)

message("3. BRITE ko00001 (auxiliary; ko01000 is never fetched or used)")
brite_json <- kegg_cached("get/br:ko00001/json", "ko00001.json")
brite_root <- jsonlite::fromJSON(paste(brite_json, collapse = "\n"), simplifyVector = FALSE)
ko_top <- local({
    rows <- list()
    recurse <- function(node, top) {
        if (is.null(node$children)) {
            m <- regmatches(node$name, regexpr("^K[0-9]+", node$name))
            if (length(m) && nzchar(m)) {
                rows[[length(rows) + 1]] <<- c(ko = m, top = top)
            }
            return(invisible(NULL))
        }
        for (ch in node$children) recurse(ch, top)
    }
    for (ch in brite_root$children) recurse(ch, ch$name)
    d <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    tapply(d$top, d$ko, function(v) paste(sort(unique(v)), collapse = "; "))
})

message("4. Candidate groups (Steps B, C, G): reactions with >= 2 distinct KOs")
rxn_kos <- lapply(split(ko_rn_all$ko, ko_rn_all$reaction), unique)
cand_kos <- rxn_kos[lengths(rxn_kos) >= 2L]
steps[["candidate_reactions"]] <- length(cand_kos)
cand_rank <- data.frame(reaction = names(cand_kos), n_ko = lengths(cand_kos),
    stringsAsFactors = FALSE)
cand_rank$metabolism <- vapply(cand_kos, function(ko) {
    any(grepl("Metabolism", ko_top[ko], fixed = TRUE), na.rm = TRUE)
}, logical(1))
# Prefer Metabolism-annotated candidates (Step C); within that, and within the
# rest, rank by number of KOs desc, reaction id asc (Step G items i and iii;
# item ii, total qualifying proteins, is only known after Step E, so it only
# orders candidates that tie on i after probing -- not reordered here).
cand_rank <- cand_rank[order(-cand_rank$metabolism, -cand_rank$n_ko, cand_rank$reaction), ]

message("5. UniProtKB/Swiss-Prot, Bacteria/Archaea, any KEGG cross-reference (Section 2.1)")
uniprot_query <- paste(
    "(reviewed:true) AND (fragment:false)",
    "AND (taxonomy_id:2 OR taxonomy_id:2157)",
    "AND (xref:kegg-*)",
    sprintf("AND (length:[%d TO %d])", length_range[1], length_range[2])
)
uniprot_fields <- paste0("accession,protein_name,organism_name,organism_id,",
    "lineage,length,xref_kegg,sequence")
uniprot_page <- function(cursor = NULL) {
    url <- paste0("https://rest.uniprot.org/uniprotkb/search?format=json",
        "&query=", curl::curl_escape(uniprot_query), "&fields=", uniprot_fields,
        "&size=500", if (!is.null(cursor)) paste0("&cursor=", cursor) else "")
    res <- curl::curl_fetch_memory(url)
    stopifnot(res$status_code == 200L)
    list(results = jsonlite::fromJSON(rawToChar(res$content), simplifyVector = FALSE)$results,
        link = curl::parse_headers_list(res$headers)[["link"]])
}
if (!smoke) {
    # One cheap probe to size it, then stop: this query matches essentially
    # every reviewed Bacteria/Archaea Swiss-Prot entry (measured live on
    # 2026-11-12/13: 334,791 of 355,745 total), which is a large bulk download
    # with sequences, not a small probe. Per KOFAM_SIM_PLAN.md Section 10,
    # stop rather than start it; see IMPLEMENTATION_NOTES.md.
    probe <- curl::curl_fetch_memory(paste0("https://rest.uniprot.org/uniprotkb/search?",
        "format=json&query=", curl::curl_escape(uniprot_query),
        "&fields=accession&size=1"))
    total <- curl::parse_headers_list(probe$headers)[["x-total-results"]]
    stop("Full build: the Swiss-Prot query matches ", total, " entries (reviewed, ",
        "non-fragment, Bacteria/Archaea, length ", length_range[1], "-",
        length_range[2], ", any KEGG cross-reference), a large download with ",
        "sequences. Stopping per KOFAM_SIM_PLAN.md Section 10; see ",
        "IMPLEMENTATION_NOTES.md. Run with KOFAM_SMOKE=1 for the pipeline check ",
        "that was actually run, or edit this stop deliberately to proceed.",
        call. = FALSE)
}
# Steps 5-6 are the expensive, purely network-bound part of a run (UniProt
# pagination, then rate-limited KEGG gene->KO batches); cache the result
# (untracked, data/kegg/) so iterating on what follows does not repeat them.
pool_cache <- file.path(kegg_dir, sprintf("pool_%s_%d.rds", if (smoke) "smoke" else "full",
    smoke_sample_size))
if (file.exists(pool_cache)) {
    message("  using cached pool: ", pool_cache)
    cached <- readRDS(pool_cache)
    pool <- cached$pool
    steps <- c(steps, cached$steps)
} else {
    entries <- list()
    cursor <- NULL
    repeat {
        p <- uniprot_page(cursor)
        entries <- c(entries, p$results)
        message("  fetched ", length(entries), " Swiss-Prot entries")
        if (length(entries) >= smoke_sample_size || !length(p$results)) break
        m <- regmatches(p$link, regexpr("cursor=[^&>]+", p$link))
        if (!length(m) || !nzchar(m)) break
        cursor <- sub("^cursor=", "", m)
    }
    pool <- do.call(rbind, lapply(entries, function(e) {
        kegg_genes <- vapply(Filter(function(x) identical(x$database, "KEGG"),
            e$uniProtKBCrossReferences), `[[`, "", "id")
        data.frame(
            protein = e$primaryAccession,
            name = if (is.null(e$proteinDescription$recommendedName$fullName$value)) {
                NA_character_
            } else e$proteinDescription$recommendedName$fullName$value,
            organism = e$organism$scientificName, taxon_id = e$organism$taxonId,
            superkingdom = e$organism$lineage[[1]], length = e$sequence$length,
            sequence = e$sequence$value,
            kegg_genes = paste(kegg_genes, collapse = ";"),
            stringsAsFactors = FALSE
        )
    }))
    steps[["swissprot_entries_fetched"]] <- nrow(pool)
    pool <- pool[!grepl("[^ACDEFGHIKLMNPQRSTVWY]", pool$sequence), ]
    steps[["swissprot_standard_residues"]] <- nrow(pool)

    message("6. Map KEGG genes to KO (Section 2.1 item 2: link/ko/<genes>)")
    all_genes <- unique(unlist(strsplit(pool$kegg_genes, ";")))
    gene_ko <- kegg_genes_to_ko(all_genes)
    steps[["kegg_genes_total"]] <- length(all_genes)
    steps[["kegg_genes_with_ko"]] <- if (is.null(gene_ko)) 0L else nrow(gene_ko)
    ko_of_gene <- if (is.null(gene_ko)) character(0) else setNames(gene_ko$ko, gene_ko$gene)
    pool$ko <- vapply(strsplit(pool$kegg_genes, ";"), function(g) {
        k <- unique(ko_of_gene[g])
        k <- k[!is.na(k)]
        if (length(k) == 1L) k else NA_character_ # 0 or >1 distinct KO: unannotated/ambiguous
    }, NA_character_)
    steps[["ambiguous_or_unannotated_KO"]] <- sum(is.na(pool$ko))
    pool <- pool[!is.na(pool$ko), c("protein", "name", "organism", "taxon_id",
        "superkingdom", "length", "sequence", "ko")]
    pool <- pool[order(pool$protein), ]
    pool <- pool[!duplicated(pool$sequence), ] # one accession per unique sequence
    steps[["swissprot_with_one_KO_unique_sequence"]] <- nrow(pool)
    saveRDS(list(pool = pool, steps = steps), pool_cache, compress = "xz")
}

message("7. Probe loop (Steps D-F): inspect ranked candidates until ", k_target,
    " accepted or ", max_probe, " probed, against the local Swiss-Prot/KO table")
probe_log <- list()
accepted <- NULL
has_subunit_language <- function(reaction, kos) {
    txt <- paste(kegg_fetch(paste0("get/rn:", reaction)),
        kegg_get_batch("ko", kos))
    grepl("subunit", txt, ignore.case = TRUE)
}
for (i in seq_len(min(nrow(cand_rank), max_probe))) {
    if (length(accepted) >= k_target) break
    reaction <- cand_rank$reaction[i]
    kos <- cand_kos[[reaction]]
    members <- pool[pool$ko %in% kos, ]
    message("  probing ", reaction, " (", length(kos), " KOs, ", nrow(members),
        " proteins in the local table)")
    if (!nrow(members)) {
        probe_log[[reaction]] <- list(reaction = reaction, kos = kos, n_ko = length(kos),
            status = "rejected", reason = "no_proteins_in_fetched_sample")
        next
    }
    qualifying_kos <- names(which(table(members$ko) >= min_proteins_per_ko))
    if (length(qualifying_kos) < 2L) {
        message("    rejected: only ", length(qualifying_kos), " KO(s) with >= ",
            min_proteins_per_ko, " proteins")
        probe_log[[reaction]] <- list(reaction = reaction, kos = kos, n_ko = length(kos),
            status = "rejected", reason = "fewer_than_2_KOs_with_enough_proteins",
            n_proteins = nrow(members))
        next
    }
    if (has_subunit_language(reaction, qualifying_kos)) {
        message("    rejected: subunit language")
        probe_log[[reaction]] <- list(reaction = reaction, kos = kos, n_ko = length(kos),
            status = "rejected", reason = "subunit_language_in_ko_or_reaction_entry")
        next
    }
    members <- members[members$ko %in% qualifying_kos, ]
    members <- do.call(rbind, lapply(split(members, members$ko), function(d) {
        if (nrow(d) > n_per_ko) d <- d[sample.int(nrow(d), n_per_ko), ]
        d
    }))
    # Step F: no cross-KO hit or shared 30%-identity cluster within the group.
    mm <- mmseqs_run(setNames(members$sequence, members$protein))
    cross <- mm$hits[mm$hits$p1 %in% members$protein & mm$hits$p2 %in% members$protein, ]
    ko_of <- setNames(members$ko, members$protein)
    cross_ko_hit <- nrow(cross) && any(ko_of[cross$p1] != ko_of[cross$p2])
    c30_tab <- table(mm$c30[members$protein], ko_of[members$protein])
    cross_ko_cluster <- any(rowSums(c30_tab > 0) > 1)
    if (cross_ko_hit || cross_ko_cluster) {
        message("    rejected: cross-KO sequence homology among ", nrow(members), " proteins")
        probe_log[[reaction]] <- list(reaction = reaction, kos = kos, n_ko = length(kos),
            status = "rejected", reason = "cross_KO_sequence_homology",
            n_proteins = nrow(members))
        next
    }
    message("    accepted: ", length(qualifying_kos), " KOs, ", nrow(members), " proteins")
    probe_log[[reaction]] <- list(reaction = reaction, kos = qualifying_kos,
        n_ko = length(qualifying_kos), status = "accepted", reason = NA_character_,
        n_proteins = nrow(members))
    members$role <- "treatment"
    members$treatment_group <- reaction
    accepted <- members
}
candidates_table <- do.call(rbind, lapply(probe_log, function(p) data.frame(
    reaction = p$reaction, n_ko = p$n_ko, kos = paste(p$kos, collapse = ";"),
    status = p$status, reason = p$reason,
    n_proteins = if (is.null(p$n_proteins)) NA_integer_ else p$n_proteins,
    stringsAsFactors = FALSE
)))
rownames(candidates_table) <- NULL
if (is.null(accepted)) {
    print(candidates_table)
    stop("No candidate reaction survived Steps D-F within max_probe = ", max_probe,
        " candidates against the ", nrow(pool), "-protein fetched sample. Raise ",
        "max_probe or smoke_sample_size.", call. = FALSE)
}
treatment <- accepted
message("Accepted treatment reaction ", treatment$treatment_group[1], " (",
    length(unique(treatment$ko)), " KOs, ", nrow(treatment), " proteins).")

message("8. Background: every other KO in the fetched sample, effect size 0 by design")
background <- pool[!pool$ko %in% treatment$ko, ]
background <- do.call(rbind, lapply(split(background, background$ko), function(d) {
    if (nrow(d) > n_per_ko) d <- d[sample.int(nrow(d), n_per_ko), ]
    d
}))
background$role <- "background"
background$treatment_group <- NA_character_
steps[["background_KOs_available"]] <- length(unique(background$ko))
steps[["background_proteins_before_cap"]] <- nrow(background)
# Cap background to the global budget *before* the leakage check below, so
# that step (and the final mmseqs object) scales with global_cap, not with
# every KO touched by the fetched Swiss-Prot sample.
max_background <- global_cap - nrow(treatment)
if (nrow(background) > max_background) {
    background <- background[sample.int(nrow(background), max_background), ]
}
steps[["background_KOs_kept"]] <- length(unique(background$ko))
steps[["background_proteins_after_cap"]] <- nrow(background)

message("9. Combine, remove leakage, apply the global cap")
universe <- rbind(treatment, background)
mm <- mmseqs_run(setNames(universe$sequence, universe$protein))
treat_ids <- universe$protein[universe$role == "treatment"]
bg_ids <- universe$protein[universe$role == "background"]
h <- mm$hits
leak_hit <- union(h$p1[h$p1 %in% bg_ids & h$p2 %in% treat_ids],
    h$p2[h$p2 %in% bg_ids & h$p1 %in% treat_ids])
leak_cluster <- bg_ids[mm$c30[bg_ids] %in% unique(mm$c30[treat_ids])]
leak <- union(leak_hit, leak_cluster)
steps[["background_leakage_removed"]] <- length(leak)
universe <- universe[!universe$protein %in% leak, ]
mm$c50 <- mm$c50[names(mm$c50) %in% universe$protein]
mm$c30 <- mm$c30[names(mm$c30) %in% universe$protein]
mm$hits <- mm$hits[mm$hits$p1 %in% universe$protein & mm$hits$p2 %in% universe$protein, ]

if (sum(universe$role == "treatment") > global_cap) {
    stop("Treatment alone (", sum(universe$role == "treatment"), " proteins) exceeds ",
        "the global cap (", global_cap, "); stopping rather than thinning the ",
        "treatment set. Lower n_per_ko or the number of treatment KOs.", call. = FALSE)
}
if (nrow(universe) > global_cap) {
    bg <- universe[universe$role == "background", ]
    keep_bg <- bg$protein[sample.int(nrow(bg), global_cap - sum(universe$role == "treatment"))]
    universe <- universe[universe$role == "treatment" | universe$protein %in% keep_bg, ]
}
steps[["final_universe"]] <- nrow(universe)
rownames(universe) <- NULL

message("10. Reaction/BRITE annotation for every retained KO")
ko_reactions <- tapply(ko_rn_all$reaction[ko_rn_all$ko %in% universe$ko],
    ko_rn_all$ko[ko_rn_all$ko %in% universe$ko],
    function(v) paste(sort(unique(v)), collapse = ";"))
universe$reactions <- unname(ko_reactions[universe$ko])
universe$brite_top <- unname(ko_top[universe$ko])

message("11. Embeddings")
if (smoke) {
    # Explicitly sanctioned for smoke only (KOFAM_SIM_PLAN.md Section 10): not
    # a result, just enough structure for the pipeline to run end to end.
    Z <- matrix(rnorm(nrow(universe) * 512), nrow(universe), 512)
    Z <- Z / sqrt(rowSums(Z^2))
    rownames(Z) <- universe$protein
    embed_note <- "KOFAM_SMOKE=1: random unit vectors, not real TM-Vec embeddings"
} else {
    library(repdist)
    basilisk::setBasiliskShared(FALSE)
    basilisk::setBasiliskFork(FALSE)
    x0 <- Biostrings::AAStringSet(setNames(universe$sequence, universe$protein))
    Z <- S4Vectors::mcols(embed_proteins(x0, model = "tmvec"))$embedding
    Z <- round(Z, 6)
    embed_note <- "embed_proteins(model = 'tmvec')"
}
x <- Biostrings::AAStringSet(setNames(universe$sequence, universe$protein))
S4Vectors::mcols(x)$embedding <- Z[names(x), ]
universe$sequence <- NULL

saveRDS(list(
    proteins = x,
    annotation = universe,
    ko_rn = ko_rn_all[ko_rn_all$ko %in% universe$ko, ],
    candidates = candidates_table,
    filters = list(
        min_proteins_per_ko = min_proteins_per_ko, n_per_ko = n_per_ko,
        global_cap = global_cap, max_probe = max_probe,
        mmseqs = mm$provenance, steps = steps
    ),
    provenance = list(
        sources = c(
            kegg = "https://rest.kegg.jp (Kanehisa Laboratories; academic use, <=3 req/s; see kofam-official-sources.md)",
            uniprot = "https://rest.uniprot.org (CC BY 4.0)"
        ),
        ko_rn_endpoint = ko_rn_source,
        kegg_info = info,
        uniprot_query = uniprot_query,
        accessed = format(Sys.Date()),
        length_range = length_range,
        seed = seed,
        smoke = smoke,
        embedding = embed_note,
        mmseqs_version = mm$provenance$version,
        repdist = tryCatch(as.character(utils::packageVersion("repdist")),
            error = function(e) NA_character_),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")

print(steps)
print(candidates_table)
cat("Universe:", nrow(universe), "proteins,", length(unique(universe$ko)), "KOs (",
    sum(universe$role == "treatment"), "treatment,", sum(universe$role == "background"),
    "background )\n")
