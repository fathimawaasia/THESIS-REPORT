#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(vroom)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(purrr)
  library(scales)
  library(cowplot)
})

cat("[SIMPLE] bulk TF validation + figure builder\n")

# ------------------------ CLI ------------------------
option_list <- list(
  make_option("--expr", type="character", help="Bulk expression matrix (TSV) with gene ID column + samples"),
  make_option("--meta", type="character", help="Sample metadata (TSV). Should contain sample/tissue/group info"),
  make_option("--gmt",  type="character", help="Comma-separated GMT paths (e.g., DC2_unsigned.gmt,DP_unsigned.gmt)"),
  make_option("--fractions-core", type="character", default=NA,
              help="CIBERSORTx core fractions CSV (DC2/DP cell-state fractions) [optional]"),
  make_option("--outdir", type="character", default="SS_SIMPLE_DC2_DP",
              help="Output directory [default %default]"),
  make_option("--expr-id-col", type="character", default=NA,
              help="Name of gene column in --expr (auto if NA) [default %default]"),
  make_option("--sample-col", type="character", default="sample_id",
              help="Sample column name in metadata [default %default]"),
  make_option("--tissue-col", type="character", default="tissue",
              help="Tissue column name in metadata [default %default]"),
  make_option("--disease-col", type="character", default=NA,
              help="Disease/condition column name (auto if NA) [default %default]"),
  make_option("--tissue-name", type="character", default="Ileum",
              help="Focus tissue for headline panels [default %default]"),
  make_option("--min-overlap", type="integer", default=3,
              help="Minimum overlap genes between regulon & matrix [default %default]"),
  make_option("--log2p1", action="store_true", default=TRUE,
              help="Apply log2(x+1) [default %default]"),
  make_option("--symbol-mode", type="character", default="auto",
              help="Symbol handling: auto|upper|none [default %default]"),
  make_option("--topk", type="integer", default=10,
              help="Top K regulons for headline panels [default %default]"),
  make_option("--alpha", type="double", default=0.05,
              help="FDR threshold for stars [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$expr) || is.null(opt$meta) || is.null(opt$gmt)) {
  stop("Missing --expr / --meta / --gmt")
}

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cat("[INFO] Using:\n",
    "  expr:   ", opt$expr, "\n",
    "  meta:   ", opt$meta, "\n",
    "  gmt(s): ", opt$gmt, "\n",
    "  outdir: ", opt$outdir, "\n", sep = "")

# ------------------------ helpers ------------------------

harmonize_symbols <- function(x, mode=c("auto","upper","none")) {
  mode <- match.arg(mode)
  if (mode == "none") return(x)
  x2 <- x
  if (mode %in% c("auto","upper")) x2 <- toupper(x2)
  if (mode == "auto") {
    # Use upper-case; preserve as-is if needed. No external deps.
    # If you want exact HGNC normalization, switch to HGNChelper here.
  }
  x2
}

# tuple-aware GMT reader (handles tuple weights or plain GMT)
read_gmt <- function(path) {
  lines <- readr::read_lines(path)
  out <- lapply(lines, function(L){
    xs <- strsplit(L, "\t", fixed = TRUE)[[1]]
    if (length(xs) < 2) return(NULL)
    set  <- xs[1]
    desc <- xs[2]
    rest <- xs[-c(1,2)]

    genes <- character(0)
    if (length(rest) >= 1) {
      # tuple style: tokens like "('GENE' 1.23)"
      if (any(grepl("\\(['\"][A-Za-z0-9._-]+['\"]\\s+", rest)) ||
          any(grepl("^\\[?\\(['\"]", rest))) {
        m <- stringr::str_match(rest, "['\"]([A-Za-z0-9._-]+)['\"]")
        genes <- m[,2]
        genes <- genes[!is.na(genes)]
      } else {
        genes <- rest
        if (length(genes) == 1 && grepl("[,; ]", genes)) {
          genes <- unlist(strsplit(genes, "[,; ]+"))
        }
      }
    }
    genes <- stringr::str_trim(genes)
    genes <- genes[nzchar(genes)]
    if (!length(genes)) return(NULL)
    tibble::tibble(set = set, gene = genes)
  })
  dplyr::bind_rows(out)
}

find_first_present <- function(df, candidates, fallback = NULL) {
  for (c in candidates) if (c %in% names(df)) return(c)
  fallback
}

score_sets_meanZ <- function(expr_mat, sets, min_overlap = 3, outdir=".") {
  stopifnot(is.matrix(expr_mat))
  z <- t(scale(t(expr_mat)))
  z[!is.finite(z)] <- 0

  ov <- tibble(
    set = names(sets),
    overlap = vapply(sets, function(g) sum(rownames(expr_mat) %in% g), numeric(1))
  )
  readr::write_tsv(ov, file.path(outdir, "overlap_diagnostics.tsv"))
  keep <- ov$set[ov$overlap >= min_overlap]
  if (!length(keep)) {
    stop("No gene sets have overlap >= ", min_overlap, ". See overlap_diagnostics.tsv")
  }
  sets <- sets[keep]

  score_list <- lapply(names(sets), function(sname){
    g <- intersect(sets[[sname]], rownames(expr_mat))
    colMeans(z[g, , drop=FALSE])
  })
  S <- do.call(rbind, score_list)
  rownames(S) <- names(sets)
  S
}

fdr_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*", "")))
}

nice_num <- function(x) sprintf("%.2f", x)

# ------------------------ load expression ------------------------
cat("[INFO] Loading expression…\n")
expr_raw <- suppressWarnings(readr::read_tsv(opt$expr, show_col_types = FALSE, progress = FALSE))
if (!is.na(opt$`expr-id-col`)) {
  gene_col <- opt$`expr-id-col`
} else {
  gene_col <- find_first_present(expr_raw,
              c("SYMBOL","symbol","gene","Gene","gene_symbol","gene_id","id", names(expr_raw)[1]),
              names(expr_raw)[1])
}
if (!gene_col %in% names(expr_raw)) stop("Gene column not found: ", gene_col)

# move gene col to SYMBOL
expr_tbl <- expr_raw %>%
  mutate(SYMBOL = .data[[gene_col]]) %>%
  select(SYMBOL, where(is.numeric))

if (isTRUE(opt$log2p1) || max(expr_tbl[,-1], na.rm = TRUE) > 50) {
  cat("[INFO] Applying log2(x+1)…\n")
  expr_tbl[,-1] <- log2(expr_tbl[,-1, drop=FALSE] + 1)
}

expr_tbl <- expr_tbl %>%
  mutate(SYMBOL = harmonize_symbols(SYMBOL, opt$`symbol-mode`)) %>%
  filter(!is.na(SYMBOL) & nzchar(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarise(across(where(is.numeric), \(x) max(x, na.rm = TRUE)), .groups="drop")

expr_mat <- expr_tbl %>% column_to_rownames("SYMBOL") %>% as.matrix()
cat("[INFO] Matrix: ", nrow(expr_mat), " genes × ", ncol(expr_mat), " samples\n", sep="")

# ------------------------ metadata ------------------------
cat("[INFO] Loading metadata…\n")
meta <- readr::read_tsv(opt$meta, show_col_types = FALSE, progress = FALSE)

sample_col  <- if (opt$sample_col %in% names(meta)) opt$sample_col else
  find_first_present(meta, c("sample_id","sample","Sample","GSM","id"), "sample_id")
tissue_col  <- if (opt$tissue_col %in% names(meta)) opt$tissue_col else
  find_first_present(meta, c("tissue","Tissue"), NA)
disease_col <- if (!is.na(opt$`disease-col`) && opt$`disease-col` %in% names(meta)) opt$`disease-col` else
  find_first_present(meta, c("disease_group","disease","condition","Condition","group"), NA)

meta <- meta %>%
  rename(sample_id = all_of(sample_col))
if (!is.na(tissue_col))  meta <- meta %>% rename(tissue = all_of(tissue_col))
if (!is.na(disease_col)) meta <- meta %>% rename(disease_group = all_of(disease_col))

# ------------------------ read GMTs ------------------------
cat("[INFO] Reading GMTs…\n")
gmt_paths <- strsplit(opt$gmt, ",", fixed = TRUE)[[1]] %>% trimws()
gmt_df <- map_dfr(gmt_paths, \(p){
  df <- read_gmt(p)
  if (!nrow(df)) return(df)
  pref <- tools::file_path_sans_ext(basename(p))
  df$set <- paste0(pref, "::", df$set)
  df
})
gmt_df <- gmt_df %>%
  mutate(gene = harmonize_symbols(gene, opt$`symbol-mode`)) %>%
  filter(!is.na(gene) & nzchar(gene))
gs <- split(gmt_df$gene, gmt_df$set)

# ------------------------ scoring ------------------------
cat("[INFO] Scoring regulons (mean Z)…\n")
S <- score_sets_meanZ(expr_mat, gs, min_overlap = opt$`min-overlap`, outdir = opt$outdir)
scores <- as.data.frame(t(S), stringsAsFactors = FALSE) %>% rownames_to_column("sample_id")
readr::write_tsv(scores, file.path(opt$outdir, "meanz_scores.tsv"))

# ------------------------ effect sizes (Ileum) ------------------------
cat("[INFO] Computing Ileum effect sizes…\n")
df_long <- scores %>%
  pivot_longer(-sample_id, names_to="set", values_to="score") %>%
  left_join(meta, by = "sample_id")

df_ileum <- df_long %>% filter(!is.na(tissue), tissue == opt$`tissue-name`)
if (!"disease_group" %in% names(df_ileum)) {
  stop("No disease_group column present in metadata (after rename).")
}

summ_one <- function(sdf, grp, ctrl_label="Control") {
  # return delta, ci, p, n per set
  out <- sdf %>%
    group_by(set) %>%
    summarise(
      m_grp = mean(score[disease_group==grp], na.rm=TRUE),
      m_ctl = mean(score[disease_group==ctrl_label], na.rm=TRUE),
      sd_grp = sd(score[disease_group==grp], na.rm=TRUE),
      sd_ctl = sd(score[disease_group==ctrl_label], na.rm=TRUE),
      n_grp  = sum(disease_group==grp & !is.na(score)),
      n_ctl  = sum(disease_group==ctrl_label & !is.na(score)),
      p      = tryCatch(wilcox.test(score[disease_group==grp],
                                    score[disease_group==ctrl_label])$p.value,
                        error=function(e) NA_real_),
      .groups="drop"
    ) %>%
    mutate(
      delta = m_grp - m_ctl,
      se = sqrt((sd_grp^2/n_grp) + (sd_ctl^2/n_ctl)),
      ci = 1.96*se
    )
  out
}

eff_uc <- summ_one(df_ileum, "UC")
eff_cd <- summ_one(df_ileum, "CD")

effects <- eff_uc %>%
  select(set, delta_uc = delta, ci_uc = ci, p_uc = p, n_uc = n_grp, n_ctl_uc = n_ctl) %>%
  full_join(
    eff_cd %>% select(set, delta_cd = delta, ci_cd = ci, p_cd = p, n_cd = n_grp, n_ctl_cd = n_ctl),
    by="set"
  ) %>%
  mutate(
    fdr_uc = p.adjust(p_uc, method="BH"),
    fdr_cd = p.adjust(p_cd, method="BH"),
    star_uc = fdr_stars(fdr_uc),
    star_cd = fdr_stars(fdr_cd)
  )

readr::write_tsv(effects, file.path(opt$outdir, "effect_sizes_ileum.tsv"))

# ------------------------ correlation with fractions (core) ------------------------
corr_tbl <- NULL
if (!is.na(opt$`fractions-core`) && file.exists(opt$`fractions-core`)) {
  cat("[INFO] Correlating regulon scores with core fractions…\n")
  frac <- suppressWarnings(readr::read_csv(opt$`fractions-core`, show_col_types=FALSE))
  # pick sample id column
  cand <- c("Mixture","sample_id","sample","Sample","ID","SubjectID","rowname")
  sid <- find_first_present(frac, cand, cand[1])
  if (!(sid %in% names(frac))) {
    # maybe sample IDs are rownames
    if (!is.null(rownames(frac))) {
      frac <- frac %>% rownames_to_column("sample_id")
      sid <- "sample_id"
    } else stop("Cannot find sample id column in fractions-core.")
  }
  frac <- frac %>% rename(sample_id = all_of(sid))
  # keep numeric fraction columns
  num_cols <- names(frac)[sapply(frac, is.numeric)]
  # sometimes file includes confidence cols; keep 0..1 columns
  frac_mat <- frac %>% select(sample_id, all_of(num_cols))

  # join to scores (samples)
  sc <- scores %>% inner_join(frac_mat, by="sample_id")
  state_cols <- setdiff(names(frac_mat), "sample_id")
  set_cols   <- setdiff(names(scores), "sample_id")

  corr_tbl <- tidyr::crossing(set = set_cols, state = state_cols) %>%
    mutate(r = purrr::map2_dbl(set, state, ~{
      x <- sc[[.x]]; y <- sc[[.y]]
      if (sd(x, na.rm=TRUE)==0 || sd(y, na.rm=TRUE)==0) return(NA_real_)
      suppressWarnings(cor(x, y, use="complete.obs"))
    })) %>%
    filter(!is.na(r))

  readr::write_tsv(corr_tbl, file.path(opt$outdir, "regulon_fraction_correlations.tsv"))
}

# ------------------------ select headline regulons (topK) ------------------------
topk <- opt$topk
ord_uc <- effects %>% arrange(desc(abs(delta_uc))) %>% slice_head(n=topk)
ord_cd <- effects %>% arrange(desc(abs(delta_cd))) %>% slice_head(n=topk)
headline <- union(ord_uc$set, ord_cd$set) %>% unique()

# ------------------------ Panel A: Ileum dumbbell (annotated) ------------------------
plotA_df <- effects %>%
  filter(set %in% headline) %>%
  transmute(
    set,
    group = "UC-Control",  delta = delta_uc,  ci = ci_uc,  star = star_uc
  ) %>%
  bind_rows(
    effects %>% filter(set %in% headline) %>%
      transmute(set, group = "CD-Control", delta = delta_cd, ci = ci_cd, star = star_cd)
  ) %>%
  mutate(set = factor(set, levels = (effects %>% arrange(desc(pmax(abs(delta_uc),abs(delta_cd)))) %>%
                                      pull(set) %>% intersect(headline)))) %>%
  drop_na(delta)

# label map from top correlated cell state per regulon (if we have corr_tbl)
lab_map <- NULL
if (!is.null(corr_tbl)) {
  top_state <- corr_tbl %>%
    filter(set %in% levels(plotA_df$set)) %>%
    group_by(set) %>%
    slice_max(order_by = abs(r), n=1, with_ties=FALSE) %>%
    ungroup() %>%
    mutate(nice = paste0(set, " — ", state, " (r=", nice_num(r), ")"))
  lab_map <- setNames(top_state$nice, top_state$set)
}
lab_fun <- function(vals) {
  if (is.null(lab_map)) return(vals)
  out <- lab_map[vals]; out[is.na(out)] <- vals[is.na(out)]; unname(out)
}

pA <- ggplot(plotA_df, aes(x=delta, y=set)) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey60") +
  geom_errorbar(aes(xmin=delta-ci, xmax=delta+ci, colour=group),
                width=0, alpha=0.6) +
  geom_point(aes(colour=group), size=3) +
  geom_line(aes(group=set), colour="grey60") +
  geom_text(aes(label=star),
            nudge_y=0.26, size=5, show.legend=FALSE) +
  scale_color_manual(values=c("UC-Control"="#D62728","CD-Control"="#1F77B4")) +
  scale_y_discrete(labels = lab_fun) +
  labs(title="A. Ileum TF-program rewiring (Top regulons)",
       subtitle="Mean(Z) difference vs Control ±95% CI; stars = FDR (Wilcoxon)\nY-labels show top correlated DC2/DP state (r)",
       x="Mean(Z) difference", y=NULL, colour="Contrast") +
  theme_bw(base_size=11)
ggsave(file.path(opt$outdir, "PanelA_Ileum_dumbbell_topK.png"),
       pA, width=9, height=6.5, dpi=300)

# ------------------------ Panel B: UC vs CD rank scatter ------------------------
pick_rank <- effects %>% select(set, delta_uc, delta_cd) %>% filter(set %in% headline)
pB <- ggplot(pick_rank, aes(x=delta_cd, y=delta_uc, label=set)) +
  geom_hline(yintercept=0, linetype="dashed", colour="grey70") +
  geom_vline(xintercept=0, linetype="dashed", colour="grey70") +
  geom_point(size=2.8, colour="#444444") +
  geom_abline(slope=1, intercept=0, colour="steelblue", alpha=0.5) +
  labs(title="B. UC vs CD effect size (Ileum)",
       x="CD − Control", y="UC − Control") +
  theme_bw(base_size=11)
ggsave(file.path(opt$outdir, "PanelB_UC_vs_CD_scatter.png"),
       pB, width=5.2, height=4.8, dpi=300)

# ------------------------ Panel C: Tissue × group heatmap (headline regulons) ------------------------
mean_by_tg <- df_long %>%
  filter(set %in% headline) %>%
  group_by(set, tissue, disease_group) %>%
  summarise(mean_score = mean(score, na.rm=TRUE), .groups="drop")

# pick only 4–6 tissues most populated to keep it readable
top_tissues <- mean_by_tg %>%
  count(tissue, sort=TRUE) %>% slice_head(n=6) %>% pull(tissue)

hmC <- mean_by_tg %>% filter(tissue %in% top_tissues)
pC <- ggplot(hmC, aes(x=disease_group, y=set, fill=mean_score)) +
  geom_tile() +
  facet_wrap(~ tissue, nrow=1, scales="free_x") +
  scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0) +
  labs(title="C. Mean regulon scores by tissue and group (headline regulons)",
       x="Group", y=NULL, fill="Mean(Z)") +
  theme_bw(base_size=10) +
  theme(panel.spacing.x = unit(6, "pt"))
ggsave(file.path(opt$outdir, "PanelC_tissue_group_heatmap.png"),
       pC, width=max(8, 3*length(top_tissues)), height=6, dpi=300)

# ------------------------ Panel D: Correlation heatmap (headline × top states) ------------------------
pD <- NULL
if (!is.null(corr_tbl)) {
  # pick top M states overall by |r|
  M <- 10
  top_states <- corr_tbl %>% filter(set %in% headline) %>%
    mutate(a = abs(r)) %>%
    group_by(state) %>% summarise(a = max(a, na.rm=TRUE), .groups="drop") %>%
    arrange(desc(a)) %>% slice_head(n=M) %>% pull(state)

  corr_h <- corr_tbl %>% filter(set %in% headline, state %in% top_states)
  pD <- ggplot(corr_h, aes(x=state, y=set, fill=r)) +
    geom_tile() +
    scale_fill_gradient2(low="#313695", mid="white", high="#A50026", midpoint=0, limits=c(-1,1)) +
    labs(title="D. Correlation of regulon scores with DC2/DP cell-state fractions",
         x="Cell state (CIBERSORTx core)", y=NULL, fill="r") +
    theme_bw(base_size=10) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
  ggsave(file.path(opt$outdir, "PanelD_correlations_heatmap.png"),
         pD, width=max(8, 0.7*length(top_states)+4), height=6, dpi=300)
}

# ------------------------ Panel E: Ileum small-multiples (top headline) ------------------------
head6 <- {
  top_uc <- effects %>% arrange(desc(abs(delta_uc))) %>% slice_head(n=3) %>% pull(set)
  top_cd <- effects %>% arrange(desc(abs(delta_cd))) %>% slice_head(n=3) %>% pull(set)
  unique(c(top_uc, top_cd))
}
df_e <- df_ileum %>% filter(set %in% head6)
pE <- ggplot(df_e, aes(x=disease_group, y=score, fill=disease_group)) +
  geom_violin(trim=FALSE, alpha=0.5, colour=NA) +
  geom_boxplot(width=0.18, outlier.shape=NA, alpha=0.8) +
  facet_wrap(~ set, scales="free_y", ncol=3) +
  scale_fill_manual(values=c("Control"="#4daf4a","UC"="#e41a1c","CD"="#377eb8","NA"="#999999")) +
  labs(title="E. Ileum distributions for headline regulons",
       x="Group", y="Mean(Z)") +
  theme_bw(base_size=11) + theme(legend.position="none")
ggsave(file.path(opt$outdir, "PanelE_Ileum_violin_headline6.png"),
       pE, width=9, height=6.5, dpi=300)

# ------------------------ Composite figure ------------------------
# Assemble available panels
panels_row1 <- plot_grid(pA, pB, ncol=2, rel_widths=c(1.3, 0.9), labels=NULL, align="h")
if (!is.null(pD)) {
  panels_row2 <- plot_grid(pC, pD, ncol=2, rel_widths=c(1,1), labels=NULL, align="h")
} else {
  panels_row2 <- plot_grid(pC, NULL, ncol=2, rel_widths=c(1,0.001))
}
final_fig <- plot_grid(panels_row1, pE, panels_row2, ncol=1, rel_heights=c(1.1, 0.9, 1.1))

ggsave(file.path(opt$outdir, "FIG_main_bulk_validation.png"),
       final_fig, width=16, height=12, dpi=300)

cat("\n[HINT] If you ever encounter 'No gene sets have overlap ≥ K', try:\n",
    "  --symbol-mode none   (uses symbols exactly as-is)\n",
    "or reduce --min-overlap slightly for tuple-sourced GMTs.\n", sep="")
cat("[DONE] Outputs in: ", normalizePath(opt$outdir), "\n", sep="")
