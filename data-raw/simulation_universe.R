# Build curated/simulation_universe.rds: microbial Swiss-Prot enzymes with one
# experimentally supported, complete EC number, embedded from their sequences
# with the default repdist model.
#
# Run from the root of the gh-pages branch with repdist installed:
#   Rscript data-raw/simulation_universe.R
# Needs network access to rest.uniprot.org, the curl and jsonlite packages, and
# the default `tmvec` model (downloaded on first use). A later UniProt release
# returns a different pool; the release used is recorded in `provenance`.
library(repdist)

out <- file.path("curated", "simulation_universe.rds")
length_range <- c(50L, 1000L) # the TM-Vec head's training range

# 1. Candidates -----------------------------------------------------------------
# Reviewed, non-fragment Bacteria/Archaea entries with at least one catalytic
# activity carrying experimental evidence.
query <- paste(
    "(reviewed:true) AND (fragment:false)",
    "AND (taxonomy_id:2 OR taxonomy_id:2157)",
    "AND (cc_catalytic_activity_exp:*)"
)
fields <- paste0(
    "accession,protein_name,organism_name,organism_id,lineage,length,",
    "ec,cc_catalytic_activity,sequence"
)
res <- curl::curl_fetch_memory(paste0(
    "https://rest.uniprot.org/uniprotkb/stream?format=json&compressed=true",
    "&query=", curl::curl_escape(query), "&fields=", curl::curl_escape(fields)
))
stopifnot(res$status_code == 200L)
headers <- curl::parse_headers_list(res$headers)
entries <- jsonlite::fromJSON(rawToChar(memDecompress(res$content, "gzip")),
    simplifyVector = FALSE)$results

# 2. One complete, experimental EC number -----------------------------------------
# Every EC number on the entry -- on any catalytic-activity reaction, whatever
# its evidence, and in the protein name including multi-domain parts -- must
# be the same single number, specified to all four levels (no "-" and no
# preliminary "n" numbers), and at least one reaction carrying it must have
# experimental evidence (ECO:0000269).
reaction_ec <- function(e) {
    ca <- Filter(function(c) identical(c$commentType, "CATALYTIC ACTIVITY"),
        e$comments)
    ec <- vapply(ca, function(c) {
        if (is.null(c$reaction$ecNumber)) NA_character_ else c$reaction$ecNumber
    }, character(1))
    exp <- vapply(ca, function(c) any(vapply(c$reaction$evidences,
        function(v) identical(v$evidenceCode, "ECO:0000269"), logical(1))),
        logical(1))
    list(any = unique(ec[!is.na(ec)]), exp = unique(ec[exp & !is.na(ec)]))
}
name_ec <- function(pd) {
    found <- character()
    walk <- function(o) {
        if (!is.list(o)) return()
        if (!is.null(o$ecNumbers)) {
            found <<- c(found, vapply(o$ecNumbers, `[[`, character(1), "value"))
        }
        for (k in o) walk(k)
    }
    walk(pd)
    unique(found)
}

pool <- do.call(rbind, lapply(entries, function(e) {
    rx <- reaction_ec(e)
    ec <- unique(c(rx$any, name_ec(e$proteinDescription)))
    ok <- length(ec) == 1L && grepl("^[0-9]+\\.[0-9]+\\.[0-9]+\\.[0-9]+$", ec) &&
        ec %in% rx$exp
    rec <- e$proteinDescription$recommendedName$fullName$value
    data.frame(
        protein = e$primaryAccession,
        ec4 = if (ok) ec else NA_character_,
        name = if (is.null(rec)) NA_character_ else rec,
        organism = e$organism$scientificName,
        taxon_id = e$organism$taxonId,
        superkingdom = e$organism$lineage[[1]],
        length = e$sequence$length,
        sequence = e$sequence$value
    )
}))
steps <- c(candidates = nrow(pool))
pool <- pool[!is.na(pool$ec4), ]
steps[["single_complete_ec4"]] <- nrow(pool)
pool <- pool[pool$length >= length_range[1] & pool$length <= length_range[2], ]
steps[["length_50_1000"]] <- nrow(pool)
# The coverage alignment in protein_similarity() scores with BLOSUM62, which
# has no selenocysteine (U) or pyrrolysine (O).
pool <- pool[!grepl("[UO]", pool$sequence), ]
steps[["no_U_or_O"]] <- nrow(pool)

# 3. One protein per UniRef50 cluster ----------------------------------------------
# No two retained proteins share more than ~50% sequence identity, so bins
# cannot be explained by near-identical sequences alone. The representative is
# the lowest accession in each cluster.
h <- curl::new_handle()
curl::handle_setform(h, ids = paste(pool$protein, collapse = ","),
    from = "UniProtKB_AC-ID", to = "UniRef50")
job <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory(
    "https://rest.uniprot.org/idmapping/run", handle = h)$content))$jobId
repeat {
    st <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory(paste0(
        "https://rest.uniprot.org/idmapping/status/", job))$content))
    if (!is.null(st$results) || identical(st$jobStatus, "FINISHED")) break
    stopifnot(!identical(st$jobStatus, "ERROR"))
    Sys.sleep(5)
}
map <- utils::read.delim(text = rawToChar(curl::curl_fetch_memory(paste0(
    "https://rest.uniprot.org/idmapping/stream/", job, "?format=tsv"))$content))
pool$uniref50 <- map$To[match(pool$protein, map$From)]
pool <- pool[!is.na(pool$uniref50), ]
pool <- pool[order(pool$protein), ]
pool <- pool[!duplicated(pool$uniref50), ]
steps[["uniref50_representatives"]] <- nrow(pool)
pool$ec3 <- sub("\\.[0-9]+$", "", pool$ec4)
rownames(pool) <- NULL

# 4. Embed every sequence with the default model -------------------------------------
# In chunks of 500, each in a fresh worker process: on Apple GPUs (MPS) one
# long run can exhaust device memory and fall back to the CPU for most
# sequences. Chunking changes nothing but run time.
embed_chunked <- function(seqs, chunk = 500L) {
    basilisk::setBasiliskShared(FALSE)
    basilisk::setBasiliskFork(FALSE)
    parts <- lapply(split(seq_along(seqs), ceiling(seq_along(seqs) / chunk)),
        function(i) {
            message("embedding ", max(i), " / ", length(seqs))
            S4Vectors::mcols(embed_proteins(seqs[i], model = "tmvec"))$embedding
        })
    x <- Biostrings::AAStringSet(seqs)
    S4Vectors::mcols(x)$embedding <- do.call(rbind, unname(parts))
    x
}
x <- embed_chunked(setNames(pool$sequence, pool$protein))
stopifnot(identical(names(x), pool$protein),
    identical(rownames(S4Vectors::mcols(x)$embedding), pool$protein))
# Six decimals keep the file small; the cosine error is below 1e-5.
S4Vectors::mcols(x)$embedding <- round(S4Vectors::mcols(x)$embedding, 6)
pool$sequence <- NULL

saveRDS(list(
    proteins = x,
    annotation = pool,
    provenance = list(
        source = "UniProtKB/Swiss-Prot via https://rest.uniprot.org (CC BY 4.0)",
        uniprot_release = headers[["x-uniprot-release"]],
        uniprot_release_date = headers[["x-uniprot-release-date"]],
        accessed = format(Sys.Date()),
        query = query,
        ec_rule = paste("all EC numbers on the entry are one complete EC4,",
            "carried by a reaction with ECO:0000269 evidence"),
        length_range = length_range,
        redundancy = "one protein per UniRef50 cluster, lowest accession",
        steps = steps,
        model = c(key = "tmvec", known_models()[["tmvec"]]),
        repdist = as.character(utils::packageVersion("repdist")),
        session = utils::capture.output(utils::sessionInfo())
    )
), out, compress = "xz")
print(steps)
print(table(table(pool$ec3) >= 50))
