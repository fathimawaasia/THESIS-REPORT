# =======================
# 60_ileum_TFprograms.R  (final patched)
# White background, exact palette, VIOLIN plots
# Headless PNG device + robust file reads + safe correlation handling
# =======================

options(bitmapType = "cairo")  # headless PNG (no X11)

suppressPackageStartupMessages({
  library(tidyverse)
  library(janitor)
  library(pheatmap)
  library(scales)
  library(cowplot)
  library(igraph)
  library(ggraph)
  library(ggrepel)
})

# ---- paths ----
scores_fp <- "meanz_scores.tsv"              # TSV: sample_id + TF program z-scores
meta_fp   <- "sample_metadata_core.tsv"      # TSV: sample_id, tissue, condition (may have trailing tab)
frax_fp   <- "cibersortx_fractions_core.csv" # CSV: sample_id + cell-state fractions
outdir    <- "ILEUM_TFPROGRAMS_OUTPUT"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ---- palette ----
pal_heat <- colorRampPalette(c("#2166AC","#FFFFFF","#B2182B"))(101)
col_tissue <- c("Colon"="#00A88F","Ileum"="#E07B00","Rectum"="#6A51A3")
col_condition <- c("Control"="#9E9E9E","UC"="#2E7D32","CD"="#C2185B")
col_pos <- "#B2182B"; col_neg <- "#2166AC"

# ---- read data (robust, quiet) ----
scores <- readr::read_tsv(scores_fp, progress = FALSE, show_col_types = FALSE) %>% clean_names()

# keep only first 3 columns even if file has a trailing tab/blank column
meta <- readr::read_tsv(meta_fp, progress = FALSE, show_col_types = FALSE, col_select = 1:3) %>%
        clean_names()
names(meta)[1:3] <- c("sample_id","tissue","condition")

frax <- readr::read_csv(frax_fp, progress = FALSE, show_col_types = FALSE) %>% clean_names()

stopifnot(all(c("sample_id") %in% names(scores)))
stopifnot(all(c("sample_id","tissue","condition") %in% names(meta)))
stopifnot(all(c("sample_id") %in% names(frax)))

# Ileum only
meta_ileum   <- meta %>% filter(str_to_lower(tissue) == "ileum")
scores_ileum <- meta_ileum %>% select(sample_id, tissue, condition) %>% inner_join(scores, by="sample_id")
frax_ileum   <- meta_ileum %>% select(sample_id, tissue, condition) %>% inner_join(frax,   by="sample_id")

tf_cols    <- setdiff(names(scores_ileum), c("sample_id","tissue","condition"))
state_cols <- setdiff(names(frax_ileum),   c("sample_id","tissue","condition"))
conds      <- intersect(c("UC","CD","HC","Control"), unique(scores_ileum$condition))

# =======================
# PART 1 — Top-20 TF programs (violins)
# =======================
top20_activity <- scores_ileum |>
  group_by(condition) |>
  summarise(across(all_of(tf_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop") |>
  pivot_longer(cols = -condition, names_to = "tf", values_to = "mean_z") |>
  mutate(abs_mean = abs(mean_z)) |>
  group_by(condition) |>
  slice_max(order_by = abs_mean, n = 20, with_ties = FALSE) |>
  ungroup()

write_csv(top20_activity, file.path(outdir, "Ileum_top20_TFprograms_activity_per_condition.csv"))

make_violin <- function(cond) {
  keep <- top20_activity %>% filter(condition == cond) %>% pull(tf)
  df   <- scores_ileum %>% filter(condition == cond) %>% select(sample_id, all_of(keep))
  long <- df %>% pivot_longer(-sample_id, names_to="tf", values_to="z") %>%
    mutate(tf = fct_reorder(tf, z, .fun = median, .desc = TRUE))

  ggplot(long, aes(x = tf, y = z)) +
    geom_violin(fill = "white", color = "grey40", linewidth = 0.4) +
    stat_summary(fun = "median", geom = "point", size = 1.8, color = "black") +
    geom_jitter(width = 0.12, alpha = 0.35, size = 0.7) +
    coord_flip() +
    labs(title = paste0("Ileum ", cond, ": Top 20 TF programs"),
         x = NULL, y = "TF program activity (z-score)") +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill="white", colour=NA))
}

plots_v <- lapply(conds, make_violin); names(plots_v) <- conds
for (cn in names(plots_v)) {
  ggsave(file.path(outdir, paste0("Ileum_", cn, "_top20_TFprograms_violin.png")),
         plots_v[[cn]], width = 8, height = 9, dpi = 300, bg = "white")
}

# Faceted heatmap: union of top-20 TFs across conditions (means by condition)
tf_union <- unique(top20_activity$tf)
mat_fac  <- scores_ileum |>
  group_by(condition) |>
  summarise(across(all_of(tf_union), ~ mean(.x, na.rm = TRUE)), .groups="drop") |>
  pivot_longer(-condition, names_to = "tf", values_to = "mean_z") |>
  pivot_wider(names_from = condition, values_from = mean_z) |>
  column_to_rownames("tf") %>% as.matrix()
mat_fac <- mat_fac[order(apply(abs(mat_fac),1,max), decreasing = TRUE), , drop = FALSE]

pheatmap(mat_fac,
         color = pal_heat, cluster_rows = FALSE, cluster_cols = FALSE,
         border_color = NA, fontsize_row = 7,
         main = "Ileum: Top TF programs across conditions",
         filename = file.path(outdir, "Ileum_TFprograms_faceted_heatmap.png"),
         width = 6, height = 12)

# =======================
# PART 2 — TF-program ↔ Cell-state correlation networks
# =======================
thr <- 0.5  # |rho| threshold

compute_edges <- function(cond) {
  s <- scores_ileum %>% filter(condition == cond)
  f <- frax_ileum   %>% filter(condition == cond)

  ids <- intersect(s$sample_id, f$sample_id)
  s <- s %>% filter(sample_id %in% ids) %>% arrange(sample_id)
  f <- f %>% filter(sample_id %in% ids) %>% arrange(sample_id)

  if (nrow(s) < 6) return(tibble())   # need enough samples

  s_mat <- as.matrix(s[, tf_cols, drop = FALSE])
  f_mat <- as.matrix(f[, state_cols, drop = FALSE])

  suppressWarnings({
    c_all <- cor(s_mat, f_mat, method = "spearman", use = "pairwise.complete.obs")
  })
  if (is.null(dim(c_all)) || any(dim(c_all) == 0)) return(tibble())

  rank_tbl <- tibble(tf = rownames(c_all),
                     mean_abs_r = rowMeans(abs(c_all), na.rm = TRUE)) %>%
              mutate(mean_abs_r = if_else(is.na(mean_abs_r), 0, mean_abs_r)) %>%
              arrange(desc(mean_abs_r)) %>%
              slice_head(n = 20)                 # <-- constant n (fixed)

  tf_rank <- rank_tbl$tf
  if (!length(tf_rank)) return(tibble())

  suppressWarnings({
    c_sub <- cor(s_mat[, tf_rank, drop = FALSE], f_mat,
                 method = "spearman", use = "pairwise.complete.obs")
  })
  if (is.null(dim(c_sub)) || any(dim(c_sub) == 0)) return(tibble())

  as.data.frame(as.table(c_sub)) %>%
    as_tibble() %>%
    rename(tf = Var1, state = Var2, rho = Freq) %>%  # Var1/Var2/Freq are from as.table()
    filter(!is.na(rho), abs(rho) >= thr) %>%
    mutate(condition = cond)
}

edges_list <- lapply(conds, compute_edges); names(edges_list) <- conds
for (cn in names(edges_list)) {
  ed <- edges_list[[cn]]
  if (nrow(ed)) write_csv(ed, file.path(outdir, paste0("Ileum_", cn, "_TFprogram_state_edges_absrho_ge_0p5.csv")))
}

plot_network <- function(ed, title_txt) {
  if (nrow(ed) == 0) return(ggplot() + theme_void() + ggtitle(paste0(title_txt, " (no edges ≥ 0.5)")))
  nodes_tf <- unique(ed$tf); nodes_st <- unique(ed$state)
  nodes <- tibble(name = c(nodes_tf, nodes_st),
                  type = c(rep("TF_program", length(nodes_tf)), rep("Cell_state", length(nodes_st))))
  g <- graph_from_data_frame(d = transmute(ed, from = tf, to = state, rho),
                             vertices = nodes, directed = FALSE)
  E(g)$color <- ifelse(E(g)$rho >= 0, col_pos, col_neg)
  E(g)$width <- rescale(abs(E(g)$rho), to = c(0.6, 3.2))
  V(g)$shape <- ifelse(V(g)$type == "TF_program", "square", "circle")

  set.seed(42)
  ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_colour = I(E(g)$color), edge_width = I(E(g)$width)), alpha = 0.95) +
    geom_node_point(aes(shape = type), size = 3, colour = "grey20", fill = "white", stroke = 0.7) +
    geom_node_text(aes(label = name), size = 2.8, repel = TRUE) +
    scale_shape_manual(values = c(TF_program = 22, Cell_state = 21)) +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          legend.position = "none") +
    ggtitle(title_txt)
}

nets <- lapply(names(edges_list), function(cn) {
  plot_network(edges_list[[cn]],
               paste0("Ileum ", cn, ": TF programs \u2194 cell states (|", "\u03C1", "| \u2265 0.5)"))
})
names(nets) <- names(edges_list)

for (cn in names(nets)) {
  ggsave(file.path(outdir, paste0("Ileum_", cn, "_TFprogram_state_network.png")),
         nets[[cn]], width = 8, height = 6, dpi = 300, bg = "white")
}
panel <- plot_grid(plotlist = nets, nrow = 1, labels = "AUTO")
ggsave(file.path(outdir, "Ileum_TFprogram_state_networks_UC_CD_HC_panel.png"),
       panel, width = 18, height = 6, dpi = 300, bg = "white")

message("Done. Outputs written to: ", normalizePath(outdir))



