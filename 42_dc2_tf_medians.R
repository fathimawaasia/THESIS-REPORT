#!/usr/bin/env Rscript
# Summarise AUCell medians for headline TFs by tissue & disease group (DC2)
# and generate publication-ready plots (headless PNG).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(scales)
})

# ---------- Headless PNG helpers ----------
use_ragg <- requireNamespace("ragg", quietly = TRUE)
open_png <- function(file, width = 1800, height = 1200, res = 150) {
  if (use_ragg) {
    ragg::agg_png(filename = file, width = width/72, height = height/72,
                  units = "in", res = res)
  } else if (capabilities("cairo")) {
    grDevices::png(filename = file, width = width, height = height,
                   res = res, type = "cairo")
  } else {
    options(bitmapType = "cairo")
    grDevices::png(filename = file, width = width, height = height, res = res)
  }
}
close_png <- function() grDevices::dev.off()
# ------------------------------------------

# Args: AUC CSV, META TSV, OUTDIR
args <- commandArgs(trailingOnly = TRUE)
auc_csv <- ifelse(length(args) >= 1, args[1], "DC2_auc_cells_as_cols.csv")
meta_tsv <- ifelse(length(args) >= 2, args[2], "DC2_meta.tsv")
out_dir  <- ifelse(length(args) >= 3, args[3], "SCENIC_SUMMARY_R")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

message("[INFO] AUC file:  ", auc_csv)
message("[INFO] Meta file: ", meta_tsv)
message("[INFO] Out dir:   ", out_dir)
stopifnot(file.exists(auc_csv), file.exists(meta_tsv))

# Headline TFs to plot
tfs <- c("ARID3A(+)", "ATF3(+)", "BATF(+)", "IRF4(+)", "KLF4(+)")

# Read data
auc <- read_csv(auc_csv, show_col_types = FALSE, progress = FALSE)
if (!"regulon" %in% names(auc)) {
  stop("[ERROR] 'regulon' column not found in AUC CSV.")
}
meta <- read.delim(meta_tsv, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

req_meta <- c("cell_id", "disease_group", "tissue")
if (!all(req_meta %in% names(meta))) {
  stop(sprintf("[ERROR] Meta file must contain: %s", paste(req_meta, collapse = ", ")))
}

# Overlapping cells
auc_cells <- setdiff(names(auc), "regulon")
ok_cells  <- intersect(auc_cells, meta$cell_id)
if (length(ok_cells) == 0) stop("[ERROR] No overlapping cell IDs between AUC and meta.")

# Diagnostics
diag_path <- file.path(out_dir, "_DC2_DIAG.txt")
cat(
  sprintf("[DIAG] Total AUC cells: %d\n", length(auc_cells)),
  sprintf("[DIAG] Total META rows: %d\n", nrow(meta)),
  sprintf("[DIAG] Overlap cells used: %d\n", length(ok_cells)),
  file = diag_path
)

# Tidy factor levels
desired_disease_levels <- c("Control", "UC", "CD")
present_levels <- intersect(desired_disease_levels, unique(meta$disease_group))
if (length(present_levels) == 0) present_levels <- unique(meta$disease_group)
meta <- meta %>% mutate(
  disease_group = factor(disease_group, levels = present_levels),
  tissue = as.factor(tissue)
)

# Long table (only headline TFs + overlapping cells)
long <- auc %>%
  filter(regulon %in% tfs) %>%
  pivot_longer(cols = all_of(ok_cells), names_to = "cell_id", values_to = "AUCell") %>%
  left_join(meta[, req_meta], by = "cell_id") %>%
  filter(!is.na(disease_group), !is.na(tissue))

# Counts for sanity
counts <- long %>%
  distinct(cell_id, tissue, disease_group) %>%
  count(tissue, disease_group, name = "n_cells") %>%
  arrange(tissue, disease_group)
write.table(counts, file = file.path(out_dir, "DC2_counts_by_tissue_disease.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Medians table
summary_tbl <- long %>%
  group_by(tissue, disease_group, regulon) %>%
  summarise(n = dplyr::n(),
            median = median(AUCell, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(regulon, tissue, disease_group)
write.table(summary_tbl,
            file = file.path(out_dir, "DC2_headline_TF_medians_by_tissue_condition.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---------- PLOTS ----------
# 1) Per-TF boxplots: AUCell ~ disease_group, faceted by tissue
plot_dir <- file.path(out_dir, "plots_DC2_headline_TFs")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

bp_theme <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))

for (tf in tfs) {
  df <- long %>% filter(regulon == tf)
  if (nrow(df) == 0) next
  p <- ggplot(df, aes(x = disease_group, y = AUCell)) +
    geom_boxplot(outlier.size = 0.7) +
    facet_wrap(~ tissue, scales = "free_y") +
    labs(title = paste0("DC2: ", tf, " AUCell by disease (faceted by tissue)"),
         x = NULL, y = "AUCell (0–1)") +
    bp_theme +
    coord_cartesian(ylim = c(0, 1))
  f <- file.path(plot_dir, paste0("BOX_DC2_", gsub("[^A-Za-z0-9]+","_", tf), ".png"))
  open_png(f, width = 1600, height = 900, res = 160)
  print(p)
  close_png()
}

# 2) Combined multi-TF boxplot: facet by regulon and tissue
p_all <- ggplot(long, aes(x = disease_group, y = AUCell)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(regulon ~ tissue, scales = "free_y") +
  labs(title = "DC2: headline regulons — AUCell by disease × tissue",
       x = NULL, y = "AUCell (0–1)") +
  bp_theme +
  coord_cartesian(ylim = c(0, 1))
f_all <- file.path(plot_dir, "BOX_DC2_all_TFs.png")
open_png(f_all, width = 2200, height = 1400, res = 150)
print(p_all)
close_png()

# 3) Heatmap of **median** AUCell per tissue × disease, faceted by TF
hm <- summary_tbl %>%
  mutate(disease_group = factor(disease_group, levels = present_levels))
p_hm <- ggplot(hm, aes(x = disease_group, y = tissue, fill = median)) +
  geom_tile(color = "grey85", linewidth = 0.2) +
  scale_fill_gradientn(colours = c("#2166AC", "white", "#B2182B"),
                       values  = rescale(c(0, 0.5, 1)),
                       limits  = c(0, 1),
                       oob = scales::squish,
                       name = "Median AUCell") +
  facet_wrap(~ regulon, ncol = 3) +
  labs(title = "DC2: Median AUCell by tissue × disease (headline regulons)",
       x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5))
f_hm <- file.path(plot_dir, "HEATMAP_DC2_median_AUCell.png")
open_png(f_hm, width = 1600, height = 1100, res = 160)
print(p_hm)
close_png()

message("[OK] Outputs written to: ", out_dir)
