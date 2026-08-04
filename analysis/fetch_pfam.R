# Downloads the protein universe for analysis/simulation_study.Rmd.
#
# Two functional groups (role A = transferase, role B = hydrolase), each
# pooling candidates from two Pfam families so family identity and group are
# not the same variable -- see simulation_study.Rmd for why that matters.
# 600 proteins are drawn per role (well above the 240 the simulation actually
# uses: 30 samples/group x 8 proteins/sample), pooled across a role's two
# families and randomly sampled -- no organism-matching between families.
# Length isn't matched by construction; simulation_study.Rmd checks the
# resulting distribution instead of engineering it.
#
#   PF00162  PGK               role A, transferase, EC 2.7.2.3
#   PF01546  Peptidase_M20     role B, hydrolase,   EC 3.4.-.-
#   PF00406  Adenylate_kinase  role A, transferase, EC 2.7.4.3
#   PF00293  Nudix_hydrolase   role B, hydrolase,   EC 3.6.1.-
#
#   Rscript analysis/fetch_pfam.R

N_PER_ROLE <- 600L
OUT <- file.path("data", "simulated", "pfam_universe.rds")

FAMILIES <- rbind(
  data.frame(pfam = "PF00162", family = "PGK", role = "A", ec = NA,
             len_min = 150, len_max = 400, note = "phosphoglycerate kinase (transferase)"),
  data.frame(pfam = "PF01546", family = "Peptidase_M20", role = "B", ec = NA,
             len_min = 150, len_max = 400, note = "M20 zinc peptidase (hydrolase)"),
  data.frame(pfam = "PF00406", family = "Adenylate_kinase", role = "A", ec = NA,
             len_min = 150, len_max = 260, note = "adenylate kinase (transferase)"),
  data.frame(pfam = "PF00293", family = "Nudix_hydrolase", role = "B", ec = "3.6.1.*",
             len_min = 150, len_max = 260, note = "Nudix pyrophosphohydrolase (hydrolase)")
)

# TSV rather than FASTA so the source organism comes back in the same request.
# `ec` narrows Pfam families (like Nudix) that cover several EC numbers to the
# hydrolase-annotated members only -- PF00293 also contains a reannotated
# isomerase subfamily (EC 5.3.3.2) that would blur the transferase/hydrolase
# contrast if left in.
uniprot_tsv <- function(pfam, len_min, len_max, ec = NA, size = 500L) {   # 500 is UniProt's per-page maximum
  q <- sprintf("(xref:pfam-%s) AND (reviewed:true) AND (taxonomy_id:2) AND (length:[%d TO %d])",
              pfam, len_min, len_max)
  if (!is.na(ec)) q <- sprintf("%s AND (ec:%s)", q, ec)
  url <- sprintf(paste0("https://rest.uniprot.org/uniprotkb/search?query=%s&format=tsv",
                        "&fields=accession,organism_name,organism_id,sequence&size=%d"),
                 utils::URLencode(q, reserved = TRUE), size)
  tmp <- tempfile(fileext = ".tsv")
  utils::download.file(url, tmp, quiet = TRUE)
  d <- utils::read.delim(tmp, stringsAsFactors = FALSE, quote = "")
  file.remove(tmp)
  names(d) <- c("accession", "organism", "taxid", "seq")
  d <- d[!duplicated(d$accession) & nzchar(d$seq), ]
  d[!duplicated(d$taxid), ]          # one protein per organism per family
}

dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
set.seed(1)

per_family <- lapply(seq_len(nrow(FAMILIES)), function(i) {
  f <- FAMILIES[i, ]
  d <- uniprot_tsv(f$pfam, f$len_min, f$len_max, f$ec)
  message(sprintf("%-8s %-18s role %s: %4d candidates, length median %d",
                  f$pfam, f$family, f$role, nrow(d), as.integer(median(nchar(d$seq)))))
  data.frame(protein = paste0(f$pfam, "_", d$accession), accession = d$accession,
             pfam = f$pfam, family = f$family, role = f$role, note = f$note,
             organism = d$organism, taxid = d$taxid, seq = d$seq, stringsAsFactors = FALSE)
})

universe <- do.call(rbind, lapply(split(FAMILIES, FAMILIES$role), function(fs) {
  pool <- do.call(rbind, per_family[match(fs$family, FAMILIES$family)])
  if (nrow(pool) < N_PER_ROLE)
    stop("role ", fs$role[1], ": only ", nrow(pool), " candidates (need ", N_PER_ROLE, ")")
  pool[sample(nrow(pool), N_PER_ROLE), ]
}))

stopifnot(!anyDuplicated(universe$protein))
saveRDS(universe, OUT)

message("\nlength by role after random selection:")
print(tapply(nchar(universe$seq), universe$role, summary))
message(sprintf("\n%d proteins, %d families, written to %s",
                nrow(universe), nrow(FAMILIES), OUT))
