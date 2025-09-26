#!/usr/bin/env python3
# Detect & rank marker genes per Leiden cluster + enrichment-ready lists & plots.
import os, argparse, shutil
import numpy as np
import pandas as pd
import scanpy as sc

def parse_args():
    p = argparse.ArgumentParser(description="Marker detection by Leiden clusters.")
    p.add_argument("--input", default="qc_outputs/adata_scanorama_clustered_umap_1.h5ad",
                   help="Integrated h5ad with UMAP/Leiden.")
    p.add_argument("--groupby", default="leiden",
                   help="Clustering key in .obs (e.g., 'leiden' or 'leiden_post').")
    p.add_argument("--method", default="wilcoxon",
                   choices=["wilcoxon","t-test","t-test_overestim_var","logreg"],
                   help="DE method.")
    p.add_argument("--n_top", type=int, default=100,
                   help="Top N genes to rank per group (for plotting).")
    p.add_argument("--outdir", default="results", help="Where to save CSVs.")
    p.add_argument("--figdir", default="figures_markers", help="Where to save figures.")
    p.add_argument("--enrichdir", default="enrichment_lists", help="Where to save enrichment TXT lists.")
    p.add_argument("--top_summary", type=int, default=10, help="Top N for summary/TXT.")
    p.add_argument("--filter_adj_p", type=float, default=0.05, help="Adj p-value cutoff for TXT lists.")
    p.add_argument("--min_logfc", type=float, default=0.0, help="Min logFC for 'UP' TXT lists (if available).")
    return p.parse_args()

def ensure_dirs(*paths):
    for p in paths:
        os.makedirs(p, exist_ok=True)

def copy_plot(src_name, dst_name, figdir, outdir):
    src = os.path.join(figdir, src_name)
    dst = os.path.join(outdir, dst_name)
    if os.path.exists(src):
        shutil.copy2(src, dst)

def main():
    args = parse_args()
    ensure_dirs(args.outdir, args.figdir, args.enrichdir)

    print(f"🔹 Loading: {args.input}")
    adata = sc.read(args.input)
    if args.groupby not in adata.obs.columns:
        raise ValueError(f"'{args.groupby}' not in .obs. Available: {list(adata.obs.columns)[:30]} ...")

    use_raw = adata.raw is not None
    print(f"🔹 rank_genes_groups by '{args.groupby}' (method={args.method}, use_raw={use_raw})")
    sc.tl.rank_genes_groups(
        adata,
        groupby=args.groupby,
        method=args.method,
        use_raw=use_raw,
        n_genes=args.n_top,
        tie_correct=True,
        corr_method="benjamini-hochberg",
        key_added="rank_genes_leiden"
    )

    # --- CSV outputs ---
    markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_leiden").copy()
    # Normalize column names across Scanpy versions
    if "names" in markers_df.columns and "gene" not in markers_df.columns:
        markers_df.rename(columns={"names": "gene"}, inplace=True)

    print(f"🔎 Columns in markers_df: {markers_df.columns.tolist()[:10]}")
    required_cols = {"group", "gene", "pvals_adj"}
    missing = required_cols - set(markers_df.columns)
    if missing:
        raise RuntimeError(f"Missing expected columns in DE table: {missing}")

    combined_csv = os.path.join(args.outdir, f"scanorama_{args.groupby}_{args.method}_markers_combined.csv")
    markers_df.to_csv(combined_csv, index=False)
    print(f"✅ DEGs (combined): {combined_csv}")

    top_rows, gene_list_blocks = [], []
    # silence FutureWarning by passing observed=False
    for g, sub in markers_df.groupby("group", sort=False, observed=False):
        per_csv = os.path.join(args.outdir, f"scanorama_{args.groupby}_{args.method}_markers_{g}.csv")
        sub.to_csv(per_csv, index=False)

        top_g = sub.head(args.top_summary).copy()
        top_g["rank"] = range(1, len(top_g) + 1)
        top_rows.append(top_g)

        # Build enrichment list
        if "logfoldchanges" in sub.columns and not sub["logfoldchanges"].isna().all():
            up = sub[(sub["pvals_adj"] < args.filter_adj_p) & (sub["logfoldchanges"] >= args.min_logfc)]
        else:
            up = sub[sub["pvals_adj"] < args.filter_adj_p]
        up_list = up["gene"].dropna().astype(str).tolist()
        with open(os.path.join(args.enrichdir, f"{args.groupby}_{g}_UP.txt"), "w") as f:
            f.write("\n".join(up_list))
        gene_list_blocks.append(f">{g}\n" + "\n".join(up_list[:args.top_summary]) + "\n")

    top_summary = pd.concat(top_rows, ignore_index=True) if top_rows else pd.DataFrame()
    top_csv = os.path.join(args.outdir, f"scanorama_{args.groupby}_{args.method}_top{args.top_summary}_summary.csv")
    top_summary.to_csv(top_csv, index=False)
    print(f"✅ Top-{args.top_summary} summary: {top_csv}")

    small_txt = os.path.join(args.outdir, f"{args.groupby}_{args.method}_top{args.top_summary}_gene_lists.txt")
    with open(small_txt, "w") as f:
        f.writelines(gene_list_blocks)
    print(f"✅ Compact gene lists: {small_txt}")

    # --- Plots (saved in figdir and copied to results) ---
    old_figdir = sc.settings.figdir
    sc.settings.figdir = args.figdir

    top20_name   = f"rank_genes_groups_scanorama_{args.groupby}_{args.method}_top20.png"
    heatmap_name = f"rank_genes_groups_scanorama_{args.groupby}_{args.method}_heatmap.png"
    dotplot_name = f"rank_genes_groups_scanorama_{args.groupby}_{args.method}_dotplot.png"
    violin_name  = f"rank_genes_groups_scanorama_{args.groupby}_{args.method}_stackedviolin.png"
    matrix_name  = f"rank_genes_groups_scanorama_{args.groupby}_{args.method}_matrixplot.png"

    sc.pl.rank_genes_groups(adata, n_genes=min(args.n_top, 20), sharey=False,
                            key="rank_genes_leiden",
                            save=f"_scanorama_{args.groupby}_{args.method}_top20.png",
                            show=False)
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=min(args.n_top, 10), groupby=args.groupby,
                                    key="rank_genes_leiden", show_gene_labels=True, dendrogram=False,
                                    swap_axes=True,
                                    save=f"_scanorama_{args.groupby}_{args.method}_heatmap.png", show=False)
    sc.pl.rank_genes_groups_dotplot(adata, n_genes=min(args.n_top, 10), groupby=args.groupby,
                                    key="rank_genes_leiden", standard_scale="var",
                                    save=f"_scanorama_{args.groupby}_{args.method}_dotplot.png", show=False)
    sc.pl.rank_genes_groups_stacked_violin(adata, n_genes=min(args.n_top, 10), groupby=args.groupby,
                                           key="rank_genes_leiden",
                                           save=f"_scanorama_{args.groupby}_{args.method}_stackedviolin.png",
                                           show=False)
    sc.pl.rank_genes_groups_matrixplot(adata, n_genes=min(args.n_top, 10), groupby=args.groupby,
                                       key="rank_genes_leiden", standard_scale="var",
                                       save=f"_scanorama_{args.groupby}_{args.method}_matrixplot.png",
                                       show=False)
    sc.settings.figdir = old_figdir

    copy_plot(top20_name,   "top20_leiden.png",         args.figdir, args.outdir)
    copy_plot(heatmap_name, "heatmap_leiden.png",       args.figdir, args.outdir)
    copy_plot(dotplot_name, "dotplot_leiden.png",       args.figdir, args.outdir)
    copy_plot(violin_name,  "stackedviolin_leiden.png", args.figdir, args.outdir)
    copy_plot(matrix_name,  "matrixplot_leiden.png",    args.figdir, args.outdir)

    print(f"🎉 Done. Figures copied to {args.outdir}/ and also in {args.figdir}/")
    print(f"🧬 Enrichment lists: {args.enrichdir}/")

if __name__ == "__main__":
    main()

