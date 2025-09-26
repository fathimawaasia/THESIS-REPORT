#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# ---------- tiny arg parser ----------
argv <- commandArgs(trailingOnly = TRUE)
args <- list()
if (length(argv)) {
  for (i in seq(1, length(argv), by = 2)) {
    key <- sub("^--", "", argv[i]); val <- argv[i + 1]; args[[key]] <- val
  }
}

# ---------- auto-detect defaults if missing ----------
pick_first_existing <- function(cands) {
  for (p in cands) if (file.exists(p)) return(p)
  return(NA_character_)
}

if (is.null(args$dc2_auc))
  args$dc2_auc <- pick_first_existing(c("DC2_auc_cells_as_cols.csv","DC2_auc.csv"))
if (is.null(args$dp_auc))
  args$dp_auc  <- pick_first_existing(c("DP_auc_cells_as_cols.csv","DP_auc.csv"))
if (is.null(args$dc2_meta))
  args$dc2_meta <- pick_first_existing(c("DC2_meta.tsv","DC2_meta.csv"))
if (is.null(args$dp_meta))
  args$dp_meta  <- pick_first_existing(c("DP_meta.tsv","DP_meta.csv"))
if (is.null(args$outdir))
  args$outdir <- "SCENIC_SUMMARY_R"

needed <- c("dc2_auc","dp_auc","dc2_meta","dp_meta","outdir")
missing <- needed[is.na(unlist(args[needed]))]
if (length(missing)) {
  cat("\n[USAGE] Rscript 42_attach_auc_and_plots_final.R",
      "--dc2_auc DC2_auc_cells_as_cols.csv",
      "--dp_auc DP_auc_cells_as_cols.csv",
      "--dc2_meta DC2_meta.tsv",
      "--dp_meta DP_meta.tsv",
      "--outdir SCENIC_SUMMARY_R\n\n")
  stop(sprintf("Missing inputs (not found in CWD): %s", paste(missing, collapse=", ")))
}

outdir <- args$outdir
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] Inputs chosen:\n",
    "  dc2_auc : ", args$dc2_auc,  "\n",
    "  dp_auc  : ", args$dp_auc,   "\n",
    "  dc2_meta: ", args$dc2_meta, "\n",
    "  dp_meta : ", args$dp_meta,  "\n",
    "  outdir  : ", outdir,        "\n", sep="")

# ---------- helpers ----------
read_meta <- function(path) {
  # works for tsv or csv
  ext <- tools::file_ext(path)
  if (tolower(ext) == "csv") {
    read_csv(path, show_col_types = FALSE) |> rename(cell_id = 1) |> mutate(cell_id = as.character(cell_id))
  } else {
    read_tsv(path, show_col_types = FALSE) |> rename(cell_id = 1) |> mutate(cell_id = as.character(cell_id))
  }
}

read_auc_any <- function(path) {
  x <- data.table::fread(path)
  # Case A: already cells-as-columns (first col is regulon name, subsequent headers look like cell IDs)
  first2 <- names(x)[1:2]
  looks_reg_firstcol <- grepl("\\(\\+\\)$", first2[1])  # "FOXP3(+)" etc
  looks_reg_second   <- grepl("\\(\\+\\)$", first2[2])

  if (looks_reg_firstcol && !looks_reg_second) {
    # headers are wrong; user gave AUC with regulon in col1 and regulon also in col2 -> not expected
    # fall back to transpose path below
    looks_reg_firstcol <- FALSE
  }

  if (!looks_reg_second) {
    # assume cells are columns already
    reg <- x[[1]]; x[[1]] <- NULL
    m <- as.matrix(x); rownames(m) <- reg
    return(m) # rows=regulons, cols=cells
  }

  # Case B: rows=cells -> transpose
  cells <- x[[1]]; x[[1]] <- NULL
  m <- as.matrix(x); rownames(m) <- cells
  m <- t(m)
  return(m)           # rows=regulons, cols=cells
}

melt_auc <- function(M) {
  df <- as.data.frame(M)
  df$regulon <- rownames(M)
  df |>
    relocate(regulon) |>
    pivot_longer(-regulon, names_to = "cell_id", values_to = "AUC") |>
    mutate(cell_id = as.character(cell_id))
}

qc_write <- function(path, ...) {
  capture.output(list(...), file = path, append = TRUE)
}

# ---------- load ----------
cat(sprintf("[INFO] Output dir: %s\n", outdir))
meta_dc2 <- read_meta(args$dc2_meta)
meta_dp  <- read_meta(args$dp_meta)
M_dc2    <- read_auc_any(args$dc2_auc)
M_dp     <- read_auc_any(args$dp_auc)

# ---------- long tables + joins ----------
L_dc2 <- melt_auc(M_dc2) |> inner_join(meta_dc2, by = "cell_id")
L_dp  <- melt_auc(M_dp)  |> inner_join(meta_dp,  by = "cell_id")

diag_file <- file.path(outdir, "_DIAG.txt")
unlink(diag_file)
qc_write(diag_file,
         sprintf("[DIAG] DC2: %d regulons, %d cells (pre-join)", nrow(M_dc2), ncol(M_dc2)),
         sprintf("[DIAG]  DP: %d regulons, %d cells (pre-join)", nrow(M_dp),  ncol(M_dp)),
         "[DIAG] DC2 head:", utils::head(L_dc2),
         "[DIAG]  DP head:", utils::head(L_dp))

# ---------- summaries ----------
p_dc2_cond <- L_dc2 |>
  ggplot(aes(x = disease_group, y = AUC)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~state_label, scales = "free_y") +
  theme_bw() + labs(title = "DC2 AUCell by condition", x = "Condition")

p_dp_cond <- L_dp |>
  ggplot(aes(x = disease_group, y = AUC)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~state_label, scales = "free_y") +
  theme_bw() + labs(title = "DP AUCell by condition", x = "Condition")

p_dc2_tis <- L_dc2 |>
  ggplot(aes(x = tissue, y = AUC)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~state_label, scales = "free_y") +
  theme_bw() + labs(title = "DC2 AUCell by tissue", x = "Tissue")

p_dp_tis <- L_dp |>
  ggplot(aes(x = tissue, y = AUC)) +
  geom_boxplot(outlier.size = 0.2) +
  facet_wrap(~state_label, scales = "free_y") +
  theme_bw() + labs(title = "DP AUCell by tissue", x = "Tissue")

med_dc2_tis <- L_dc2 |>
  group_by(regulon, tissue) |>
  summarise(median_AUC = median(AUC, na.rm = TRUE), .groups = "drop")

med_dp_tis <- L_dp |>
  group_by(regulon, tissue) |>
  summarise(median_AUC = median(AUC, na.rm = TRUE), .groups = "drop")

top_by_var <- function(df, n=40) {
  df |>
    group_by(regulon) |>
    summarise(varA = var(median_AUC, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(varA)) |>
    slice_head(n = n) |>
    pull(regulon)
}
top_dc2 <- top_by_var(med_dc2_tis)
top_dp  <- top_by_var(med_dp_tis)

hm_dc2 <- med_dc2_tis |> filter(regulon %in% top_dc2) |>
  mutate(regulon = factor(regulon, levels = unique(regulon))) |>
  ggplot(aes(tissue, regulon, fill = median_AUC)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_minimal() + labs(title = "DC2 | median AUCell by tissue")

hm_dp <- med_dp_tis |> filter(regulon %in% top_dp) |>
  mutate(regulon = factor(regulon, levels = unique(regulon))) |>
  ggplot(aes(tissue, regulon, fill = median_AUC)) +
  geom_tile() + scale_fill_viridis_c() +
  theme_minimal() + labs(title = "DP | median AUCell by tissue")

# ---------- outputs ----------
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
write_csv(med_dc2_tis, file.path(outdir, "S2_DC2_medianAUC_by_tissue.csv"))
write_csv(med_dp_tis,  file.path(outdir, "S4_DP_medianAUC_by_tissue.csv"))

ggsave(file.path(outdir, "01_DC2_AUC_box_by_condition.png"), p_dc2_cond, width = 12, height = 8, dpi = 200)
ggsave(file.path(outdir, "02_DC2_AUC_box_by_tissue.png"),     p_dc2_tis,  width = 12, height = 8, dpi = 200)
ggsave(file.path(outdir, "05_DP_AUC_box_by_condition.png"),   p_dp_cond,  width = 12, height = 8, dpi = 200)
ggsave(file.path(outdir, "06_DP_AUC_box_by_tissue.png"),      p_dp_tis,   width = 12, height = 8, dpi = 200)
ggsave(file.path(outdir, "03_DC2_heatmap_state_by_tissue.png"), hm_dc2,   width = 10, height = 12, dpi = 200)
ggsave(file.path(outdir, "07_DP_heatmap_state_by_tissue.png"),  hm_dp,    width = 10, height = 12, dpi = 200)

cat("[OK] All plots in:", outdir, "\n")



