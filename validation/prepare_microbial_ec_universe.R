# Build validation/microbial_ec_universe.rds: a fixed universe of microbial
# enzymes with experimentally supported EC numbers, embedded with TM-Vec.
# The eligibility rules are the whole design, so they are all stated here and
# recorded in the RDS `provenance`; see DATA_PROVENANCE.md.
#
# Run from the source root after installing repdist. Needs network access to
# rest.uniprot.org, the curl and jsonlite packages, and the default `tmvec`
# model (several GB on first use). The UniProt release is recorded: a later
# release returns a different pool, so the shipped RDS, not this script, is
# the reference universe.
library(repdist)

# Two groups of `n_per_group` samples each upregulate `n_up` private
# proteins, so one target EC3 must hold n_per_group * n_up proteins.
n_per_group <- 10L
n_up <- 5L
per_ec3 <- n_per_group * n_up
length_range <- c(50L, 1000L) # TM-Vec's head was trained through 1000 aa
seed <- 20260924L

# 1. Candidate entries -------------------------------------------------------
# cc_catalytic_activity_exp selects entries holding at least one catalytic-
# activity comment with experimental evidence (ECO:0000269). It is an
# evidence filter, not a review-status filter; reviewed status is recorded,
# never required. Bacteria (taxon 2) and Archaea (2157) define "microbial".
query <- paste(
    "(cc_catalytic_activity_exp:*) AND (fragment:false)",
    "AND (taxonomy_id:2 OR taxonomy_id:2157)"
)
fields <- paste(
    "accession,reviewed,protein_name,organism_name,organism_id,lineage,",
    "length,ec,cc_catalytic_activity,sequence,fragment",
    sep = ""
)
h <- curl::new_handle()
res <- curl::curl_fetch_memory(paste0(
    "https://rest.uniprot.org/uniprotkb/stream?format=json&compressed=true",
    "&query=", curl::curl_escape(query),
    "&fields=", curl::curl_escape(fields)
), handle = h)
stopifnot(res$status_code == 200L)
headers <- curl::parse_headers_list(res$headers)
release <- headers[["x-uniprot-release"]]
entries <- jsonlite::fromJSON(
    rawToChar(memDecompress(res$content, "gzip")),
    simplifyVector = FALSE
)$results

# 2. Experimental EC and a single EC3 ----------------------------------------
# Experimental EC: the EC number on a CATALYTIC ACTIVITY reaction whose own
# evidence includes ECO:0000269. EC numbers in the protein name, and reactions
# with only curator-inferred or automatic evidence, never make an EC
# experimental.
has_exp <- function(evidences) {
    any(vapply(evidences, function(v) identical(v$evidenceCode, "ECO:0000269"),
        logical(1)))
}
reaction_ec <- function(e) {
    ca <- Filter(function(c) identical(c$commentType, "CATALYTIC ACTIVITY"),
        e$comments)
    ec <- vapply(ca, function(c) {
        if (is.null(c$reaction$ecNumber)) NA_character_ else c$reaction$ecNumber
    }, character(1))
    exp <- vapply(ca, function(c) has_exp(c$reaction$evidences), logical(1))
    list(exp = unique(ec[exp & !is.na(ec)]), any = unique(ec[!is.na(ec)]))
}
# Every EC anywhere in the protein name, including multi-domain
# `includes`/`contains` parts.
name_ec <- function(pd) {
    out <- character()
    walk <- function(o) {
        if (!is.list(o)) return()
        if (!is.null(o$ecNumbers)) {
            out <<- c(out, vapply(o$ecNumbers, `[[`, character(1), "value"))
        }
        for (k in o) walk(k)
    }
    walk(pd)
    unique(out)
}
ec_levels <- function(ec) strsplit(ec, ".", fixed = TRUE)[[1]]
# EC3 of a number specified to at least three levels, else NA ("3.1.-.-").
ec3_of <- function(ec) {
    p <- ec_levels(ec)
    if (length(p) < 3L || !all(grepl("^[0-9]+$", p[1:3]))) return(NA_character_)
    paste(p[1:3], collapse = ".")
}
# An EC agrees with an EC3 on every level it specifies ("3.1.-.-" agrees
# with 3.1.3; "2.7.-.-" does not).
agrees <- function(ec, ec3) {
    p <- ec_levels(ec)
    q <- ec_levels(ec3)
    k <- which(!grepl("^[0-9]+$", p))[1]
    lev <- min(3L, if (is.na(k)) length(p) else k - 1L)
    identical(p[seq_len(lev)], q[seq_len(lev)])
}

pool <- do.call(rbind, lapply(entries, function(e) {
    rx <- reaction_ec(e)
    all_ec <- unique(c(rx$any, name_ec(e$proteinDescription)))
    ec3 <- unique(stats::na.omit(vapply(rx$exp, ec3_of, character(1))))
    # One EC3 from the experimental reactions, and every other EC on the
    # entry, experimental or not, consistent with it.
    ok <- length(ec3) == 1L && all(vapply(all_ec, agrees, logical(1), ec3))
    data.frame(
        protein = e$primaryAccession,
        reviewed = grepl("reviewed", e$entryType, fixed = TRUE) &&
            !grepl("unreviewed", e$entryType, fixed = TRUE),
        exp_ec = paste(rx$exp, collapse = ";"),
        ec3 = if (ok) ec3 else NA_character_,
        superkingdom = e$organism$lineage[[1]],
        organism = e$organism$scientificName,
        taxon_id = e$organism$taxonId,
        name = if (is.null(e$proteinDescription$recommendedName)) {
            NA_character_
        } else {
            e$proteinDescription$recommendedName$fullName$value
        },
        length = e$sequence$length,
        sequence = e$sequence$value,
        stringsAsFactors = FALSE
    )
}))
counts <- c(candidates = nrow(pool), exp_reaction_ec = sum(nzchar(pool$exp_ec)))
pool <- pool[!is.na(pool$ec3), ]
counts[["single_consistent_ec3"]] <- nrow(pool)
pool <- pool[pool$length >= length_range[1] & pool$length <= length_range[2], ]
counts[["length_filter"]] <- nrow(pool)

# 3. Redundancy: one protein per UniRef50 cluster -----------------------------
# UniRef50 groups sequences at 50% identity, so the universe holds no
# near-duplicates that sequence identity alone would already group.
run <- curl::new_handle()
curl::handle_setform(run,
    ids = paste(pool$protein, collapse = ","),
    from = "UniProtKB_AC-ID", to = "UniRef50"
)
job <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory(
    "https://rest.uniprot.org/idmapping/run", handle = run
)$content))$jobId
repeat {
    st <- jsonlite::fromJSON(rawToChar(curl::curl_fetch_memory(paste0(
        "https://rest.uniprot.org/idmapping/status/", job
    ))$content))
    if (!is.null(st$results) || identical(st$jobStatus, "FINISHED")) break
    Sys.sleep(3)
}
map <- utils::read.delim(text = rawToChar(curl::curl_fetch_memory(paste0(
    "https://rest.uniprot.org/idmapping/stream/", job, "?format=tsv"
))$content))
pool$uniref50 <- map$To[match(pool$protein, map$From)]
pool <- pool[!is.na(pool$uniref50), ]
counts[["uniref50_mapped"]] <- nrow(pool)
pool <- pool[order(pool$protein), ]
pool <- pool[!duplicated(pool$uniref50), ]
counts[["uniref50_dedup"]] <- nrow(pool)

# 4. Universe: every EC3 large enough to be a target --------------------------
# Eligibility is a count only, fixed before any embedding or clustering, so no
# EC3 is in or out because of how TM-Vec treats it. Each eligible EC3 is
# subsampled to exactly `per_ec3` proteins, so every EC3 in the universe can be
# a target and none dominates.
ec3_size <- table(pool$ec3)
eligible <- sort(names(ec3_size)[ec3_size >= per_ec3])
set.seed(seed)
members <- do.call(rbind, lapply(eligible, function(e3) {
    m <- pool[pool$ec3 == e3, ]
    m[sort(sample(nrow(m), per_ec3)), ]
}))
rownames(members) <- NULL
counts[["universe"]] <- nrow(members)

# 5. Embed with the default model ---------------------------------------------
seqs <- stats::setNames(members$sequence, members$protein)
Z <- S4Vectors::mcols(embed_proteins(seqs, model = "tmvec"))$embedding
stopifnot(identical(rownames(Z), members$protein))

members$sequence <- NULL
u <- list(
    members = members,
    embeddings = round(Z, 4),
    sequences = seqs,
    provenance = list(
        source = "UniProtKB REST API, https://rest.uniprot.org",
        uniprot_release = release,
        accessed = as.character(Sys.Date()),
        query = query,
        experimental_ec = paste(
            "EC number on a CATALYTIC ACTIVITY reaction whose evidence",
            "includes ECO:0000269; protein-name EC numbers never qualify"
        ),
        ec3_rule = paste(
            "exactly one EC3 among the experimental ECs, and every EC on the",
            "entry agrees with it on every level it specifies"
        ),
        length_range = length_range,
        redundancy = "one protein per UniRef50 cluster, lowest accession",
        eligibility = sprintf(
            "EC3 with at least %d proteins after deduplication", per_ec3
        ),
        n_per_group = n_per_group, n_up = n_up, per_ec3 = per_ec3,
        seed = seed,
        counts = counts,
        ec3_sizes_after_dedup = c(ec3_size),
        embedding_model = known_models()[["tmvec"]],
        repdist = as.character(utils::packageVersion("repdist")),
        session = utils::capture.output(utils::sessionInfo())
    )
)
saveRDS(u, file.path("validation", "microbial_ec_universe.rds"),
    compress = "xz"
)
print(counts)
