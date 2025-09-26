#!/usr/bin/env Rscript
suppressPackageStartupMessages({ library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
usage <- function() {
  cat("
Usage:
  Rscript gene_id_sanity.R --bulk <bulk.tsv> --sig <signature.tsv[.gz]>
                           [--out <dir>] [--write_noversion TRUE|FALSE]

Notes:
  * bulk.tsv and signature.tsv must have 'gene_id' as the first column name
  * Writes a text report and (optionally) 'noversion' copies if Ensembl versions are present
\n"); quit(save="no", status=1)
}
opt <- list(bulk=NULL, sig=NULL, out=".", write_noversion="TRUE")
kv <- strsplit(args, "=")
for (x in kv) {
  if (length(x) != 2) next
  k <- sub("^--", "", x[[1]]); v <- x[[2]]
  if (k %in% names(opt)) opt[[k]] <- v
}
if (is.null(opt$bulk) || is.null(opt$sig)) usage()
write_noversion <- toupper(opt$write_noversion) %in% c("TRUE","T","1","YES","Y")

dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
report_file <- file.path(opt$out, "gene_id_overlap_report.txt")

is_ens  <- function(x) grepl("^ENS[A-Z]*G\\d+(\\.\\d+)?$", x)
has_ver <- function(x) grepl("^ENS[A-Z]*G\\d+\\.\\d+$", x)
strip_ver <- function(x) sub("\\.\\d+$", "", x)

cat("[INFO] Loading bulk:", opt$bulk, "\n")
bulk_dt <- fread(opt$bulk); setnames(bulk_dt, 1, "gene_id")
cat("[INFO] Loading signature:", opt$sig, "\n")
sig_dt  <- fread(opt$sig);  setnames(sig_dt, 1, "gene_id")

bulk_ids <- bulk_dt$gene_id; sig_ids <- sig_dt$gene_id

bulk_ens_frac <- mean(is_ens(bulk_ids)); sig_ens_frac <- mean(is_ens(sig_ids))
bulk_ver_frac <- mean(has_ver(bulk_ids)); sig_ver_frac <- mean(has_ver(sig_ids))

overlap0 <- length(intersect(bulk_ids, sig_ids))

bulk_ids2 <- if (bulk_ver_frac > 0.05) strip_ver(bulk_ids) else bulk_ids
sig_ids2  <- if (sig_ver_frac  > 0.05) strip_ver(sig_ids)  else sig_ids
overlap1  <- length(intersect(bulk_ids2, sig_ids2))

only_bulk <- head(setdiff(unique(bulk_ids2), unique(sig_ids2)), 10)
only_sig  <- head(setdiff(unique(sig_ids2),  unique(bulk_ids2)), 10)

cat("==== Gene ID style guess ====\n")
cat(sprintf("Bulk:  %.1f%% Ensembl-like (%.1f%% with versions)\n", 100*bulk_ens_frac, 100*bulk_ver_frac))
cat(sprintf("Sig :  %.1f%% Ensembl-like (%.1f%% with versions)\n", 100*sig_ens_frac,  100*sig_ver_frac))
if (bulk_ens_frac > 0.7 && sig_ens_frac < 0.3) cat("⚠️  Likely: bulk=Ensembl, signature=SYMBOLs\n")
if (sig_ens_frac > 0.7 && bulk_ens_frac < 0.3) cat("⚠️  Likely: signature=Ensembl, bulk=SYMBOLs\n")
if (bulk_ver_frac > 0.3 || sig_ver_frac > 0.3) cat("ℹ️  Ensembl version suffixes (e.g., .12) detected\n")
cat(sprintf("Overlap (as-is): %d genes\n", overlap0))
cat(sprintf("Overlap (after stripping versions): %d genes\n", overlap1))
cat("Bulk-only examples: ", paste(only_bulk, collapse=", "), "\n")
cat("Sig-only  examples: ", paste(only_sig,  collapse=", "), "\n")

# Write report
cat(
  paste0(
    "Bulk Ensembl-like: ", sprintf("%.1f%%", 100*bulk_ens_frac),
    " (with version: ", sprintf("%.1f%%", 100*bulk_ver_frac), ")\n",
    "Sig  Ensembl-like: ", sprintf("%.1f%%", 100*sig_ens_frac),
    " (with version: ", sprintf("%.1f%%", 100*sig_ver_frac), ")\n",
    "Overlap (as-is): ", overlap0, "\n",
    "Overlap (strip versions): ", overlap1, "\n",
    "Bulk-only top10: ", paste(only_bulk, collapse=", "), "\n",
    "Sig-only  top10: ", paste(only_sig,  collapse=", "), "\n"
  ),
  file = report_file
)
cat("[WRITE] Report ->", report_file, "\n")

# Optionally write 'noversion' copies if they help
if (write_noversion && overlap1 > overlap0) {
  out_bulk <- file.path(opt$out, sub("\\.tsv(\\.gz)?$", ".noversion.tsv", basename(opt$bulk)))
  out_sig  <- file.path(opt$out,  sub("\\.tsv(\\.gz)?$", ".noversion.tsv", basename(opt$sig)))
  bulk_dt$gene_id <- bulk_ids2; sig_dt$gene_id <- sig_ids2
  fwrite(bulk_dt, out_bulk, sep="\t"); fwrite(sig_dt, out_sig, sep="\t")
  cat("✅ Wrote normalized (no-version) copies:\n  ", out_bulk, "\n  ", out_sig, "\n")
} else {
  cat("[INFO] No benefit from stripping versions or disabled; originals unchanged.\n")
}
