#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr)
  library(stringr); library(ggplot2); library(forcats)
})

## ---------- CLI (fixed) ----------
args <- commandArgs(trailingOnly = TRUE)

auc_csv  <- if (length(args) >= 1) args[1] else "DP_auc_cells_as_cols.csv"
meta_tsv <- if (length(args) >= 2) args[2] else "DP_meta.tsv"
med_csv  <- if (length(args) >= 3) args[3] else "SCENIC_SUMMARY_R/S4_DP_medianAUC_by_tissue.csv"
out_dir  <- if (length(args) >= 4) args[4] else "SCENIC_SUMMARY_R/HEADLINE_FROM_MEDIANS_1"
TOP_N    <- if (length(args) >= 5) as.integer(args[5]) else 5L


message("[INFO] Using medians: ", med_csv)
message("[INFO] AUC:          ", auc_csv)
message("[INFO] META:         ", meta_tsv)
message("[INFO] Output dir:   ", out_dir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---------- headless device ----------
use_ragg <- requireNamespace("ragg", quietly = TRUE)
save_png <- function(file, width=10, height=7, dpi=300) {
  if (use_ragg) ragg::agg_png(file, width=width, height=height, units="in", res=dpi)
  else if (capabilities("cairo")) png(file, width=width, height=height, units="in", res=dpi, type="cairo")
  else { options(bitmapType="cairo"); png(file, width=width, height=height, units="in", res=dpi) }
}

## ---------- read data ----------
auc  <- read_csv(auc_csv, show_col_types = FALSE)
stopifnot("regulon" %in% names(auc))
meta <- read.delim(meta_tsv, sep="\t", check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("cell_id","disease_group","tissue") %in% names(meta)))
meta$disease_group <- factor(meta$disease_group, levels=c("Control","UC","CD"))
meta$tissue        <- factor(meta$tissue,        levels=c("Colon","Ileum","Rectum"))

## overlap cells only
cells_auc <- setdiff(names(auc), "regulon")
cells_ok  <- intersect(cells_auc, meta$cell_id)
stopifnot(length(cells_ok) > 0)

long_all <- auc |>
  pivot_longer(cols = all_of(cells_ok), names_to="cell_id", values_to="AUCell") |>
  left_join(meta[,c("cell_id","disease_group","tissue")], by="cell_id") |>
  filter(!is.na(disease_group), !is.na(tissue))

## ---------- choose TOP-N TFs from medians ----------
med <- read_csv(med_csv, show_col_types = FALSE,
                col_types = cols(regulon=col_character(),
                                 tissue=col_character(),
                                 median_AUC=col_double()))
top_tf <- med |>
  group_by(regulon) |>
  summarise(max_median = max(median_AUC, na.rm=TRUE), .groups="drop") |>
  arrange(desc(max_median)) |>
  slice_head(n = TOP_N) |>
  pull(regulon)

message("[INFO] Top-", TOP_N, " TFs by max tissue median: ", paste(top_tf, collapse=", "))

## ---------- theme & palettes ----------
BASE_TEXT <- 18
THEME_PUB <- theme_bw(base_size = BASE_TEXT) +
  theme(
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    strip.background  = element_rect(fill="grey92", colour="grey60"),
    strip.text        = element_text(size = BASE_TEXT, face="bold"),
    axis.title        = element_text(size = BASE_TEXT, face="bold"),
    axis.text.x       = element_text(size = BASE_TEXT-4, angle = 25, vjust = 1, hjust = 1),
    axis.text.y       = element_text(size = BASE_TEXT-3),
    legend.position   = "right",
    plot.title        = element_text(size = BASE_TEXT+1, face="bold", hjust = 0)
  )

COLS_COND  <- c(Control="#4D9221", UC="#35978F", CD="#C51B7D")
COLS_TISS  <- c(Colon="#1B9E77", Ileum="#D95F02", Rectum="#7570B3")

## Utility: thin jitter points for huge groups (max 5k)
thin_jitter <- function(df, max_n = 5000) {
  if (nrow(df) <= max_n) return(df)
  set.seed(42); df[sample.int(nrow(df), max_n), , drop=FALSE]
}

## ---------- plotting helpers (per TF) ----------
p_by_condition <- function(df_tf, tf_label) {
  ggplot(thin_jitter(df_tf), aes(disease_group, AUCell, fill=disease_group)) +
    geom_violin(trim = TRUE, colour="grey30", linewidth=0.3) +
    geom_boxplot(width=.15, outlier.shape = NA, fill="white", colour="grey20") +
    geom_point(position = position_jitter(width=.12, height=0), size=.5, alpha=.25, colour="grey20") +
    facet_wrap(~ tissue, nrow = 1) +
    scale_fill_manual(values = COLS_COND, guide = "none") +
    labs(x="Condition", y="AUCell",
         title = paste0(tf_label, " | headline TFs by condition")) +
    THEME_PUB
}

p_by_tissue <- function(df_tf, tf_label) {
  ggplot(thin_jitter(df_tf), aes(tissue, AUCell, fill=tissue)) +
    geom_violin(trim = TRUE, colour="grey30", linewidth=0.3) +
    geom_boxplot(width=.15, outlier.shape = NA, fill="white", colour="grey20") +
    geom_point(position = position_jitter(width=.12, height=0), size=.5, alpha=.25, colour="grey20") +
    facet_wrap(~ disease_group, nrow = 1) +
    scale_fill_manual(values = COLS_TISS, guide = "none") +
    labs(x="Tissue", y="AUCell",
         title = paste0(tf_label, " | headline TFs by tissue")) +
    THEME_PUB
}

p_grid_cond_by_tiss <- function(df_tf, tf_label) {
  ggplot(thin_jitter(df_tf), aes(disease_group, AUCell, fill=disease_group)) +
    geom_violin(trim = TRUE, colour="grey30", linewidth=0.25, scale="width") +
    geom_boxplot(width=.12, outlier.shape = NA, fill="white", colour="grey20", linewidth=0.25) +
    facet_grid(tissue ~ ., scales="free_y") +
    scale_fill_manual(values = COLS_COND, guide = "none") +
    labs(x="Condition", y="AUCell",
         title = paste0(tf_label, " | by condition within each tissue")) +
    THEME_PUB
}

## ---------- per-TF plots ----------
panel_dir <- file.path(out_dir, "PANELS"); dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

unique_tfs <- sort(unique(long_all$regulon))
for (tf in unique_tfs) {
  d <- filter(long_all, regulon == tf)

  # by condition (facet tissue)
  f1 <- file.path(panel_dir, paste0("DP_", gsub("[^A-Za-z0-9_+]+","_", tf), "_by_condition.png"))
  save_png(f1, width=12, height=6, dpi=300); print(p_by_condition(d, tf)); dev.off()

  # by tissue (facet condition)
  f2 <- file.path(panel_dir, paste0("DP_", gsub("[^A-Za-z0-9_+]+","_", tf), "_by_tissue.png"))
  save_png(f2, width=12, height=6, dpi=300); print(p_by_tissue(d, tf)); dev.off()

  # grid condition within tissue
  f3 <- file.path(panel_dir, paste0("DP_", gsub("[^A-Za-z0-9_+]+","_", tf), "_grid_condition_by_tissue.png"))
  save_png(f3, width=8.5, height=10, dpi=300); print(p_grid_cond_by_tiss(d, tf)); dev.off()
}

## ---------- TOP-N summary panels ----------
# build a tidy table limited to top TFs
long_top <- long_all |>
  filter(regulon %in% top_tf) |>
  mutate(regulon = fct_relevel(regulon, top_tf))

# per condition (facet over TFs, 2 rows)
p_top_by_cond <- ggplot(thin_jitter(long_top), aes(disease_group, AUCell, fill=disease_group)) +
  geom_violin(trim=TRUE, colour="grey30", linewidth=0.3) +
  geom_boxplot(width=.15, outlier.shape=NA, fill="white", colour="grey20") +
  facet_wrap(~ regulon, ncol = min(5, length(top_tf))) +
  scale_fill_manual(values = COLS_COND) +
  labs(x="Condition", y="AUCell", title=paste0("Top-", TOP_N, " headline TFs | by condition")) +
  THEME_PUB

f_top_cond <- file.path(out_dir, paste0("DP_top", TOP_N, "_headline_TFs_by_condition.png"))
save_png(f_top_cond, width=16, height=9, dpi=300); print(p_top_by_cond); dev.off()

# per tissue (facet over TFs)
p_top_by_tiss <- ggplot(thin_jitter(long_top), aes(tissue, AUCell, fill=tissue)) +
  geom_violin(trim=TRUE, colour="grey30", linewidth=0.3) +
  geom_boxplot(width=.15, outlier.shape=NA, fill="white", colour="grey20") +
  facet_wrap(~ regulon, ncol = min(5, length(top_tf))) +
  scale_fill_manual(values = COLS_TISS) +
  labs(x="Tissue", y="AUCell", title=paste0("Top-", TOP_N, " headline TFs | by tissue")) +
  THEME_PUB

f_top_tiss <- file.path(out_dir, paste0("DP_top", TOP_N, "_headline_TFs_by_tissue.png"))
save_png(f_top_tiss, width=16, height=9, dpi=300); print(p_top_by_tiss); dev.off()

message("[OK] Panels in: ", panel_dir)
message("[OK] Top-N summaries: ", f_top_cond, " | ", f_top_tiss)


