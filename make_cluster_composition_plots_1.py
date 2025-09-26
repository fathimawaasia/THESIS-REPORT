#!/usr/bin/env python3
# make_cluster_composition_plots_1.py
# Creates stacked bar charts (counts & %) for:
#   1) cluster × CellTypist label
#   2) cluster × Sample ID
#   3) cluster × Group
# and writes the contingency tables to CSVs.

import os
import argparse
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def pick_first_present(obs, candidates, required=False, what=""):
    for c in candidates:
        if c in obs.columns:
            return c
    if required:
        raise KeyError(
            f"Missing required column for {what}. "
            f"Looked for any of: {', '.join(candidates)}. "
            f"Use --{what.replace(' ', '-')}-key to set explicitly "
            f"or add the column to your .h5ad."
        )
    return None

def stacked_bar(df, title, ylabel, out_png, legend_cols=None):
    # choose columns for legend automatically
    n_items = len(df.columns)
    if legend_cols is None:
        if n_items <= 18:
            legend_cols = 1
        elif n_items <= 36:
            legend_cols = 2
        else:
            legend_cols = 3

    fig, ax = plt.subplots(figsize=(16, 7))  # wider canvas
    df.plot(kind='bar', stacked=True, ax=ax)
    ax.set_title(title)
    ax.set_xlabel("Cluster")
    ax.set_ylabel(ylabel)

    ax.legend(
        bbox_to_anchor=(1.02, 1),
        loc='upper left',
        frameon=False,
        ncol=legend_cols,
        fontsize=8,
        borderaxespad=0.0,
    )
    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_tables_and_plots(adata, cluster_key, col_key, tag, outdir):
    # contingency tables
    table = pd.crosstab(adata.obs[cluster_key], adata.obs[col_key])
    table_pct = table.div(table.sum(axis=1), axis=0) * 100

    # save tables
    table.to_csv(os.path.join(outdir, f"composition_{tag}_counts_1.csv"))
    table_pct.to_csv(os.path.join(outdir, f"composition_{tag}_percent_1.csv"))

    # plots
    stacked_bar(
        table,
        title=f"Cluster composition by {tag} (Counts)",
        ylabel="Number of cells",
        out_png=os.path.join(outdir, f"composition_{tag}_counts_1.png"),
    )
    stacked_bar(
        table_pct,
        title=f"Cluster composition by {tag} (%)",
        ylabel="Percentage of cells",
        out_png=os.path.join(outdir, f"composition_{tag}_percent_1.png"),
    )

def main():
    p = argparse.ArgumentParser(description="Make composition-per-cluster plots from a CellTypist-annotated AnnData.")
    p.add_argument("--input", required=True, help="Path to .h5ad (e.g. qc_outputs/adata_celltypist_annotated_1.h5ad)")
    p.add_argument("--outdir", default="composition_outputs", help="Directory to save CSVs & PNGs")
    p.add_argument("--cluster-key", default=None, help="Cluster column in .obs (default: auto-detect)")
    p.add_argument("--celltype-key", default=None, help="Cell type column in .obs (default: auto-detect)")
    p.add_argument("--sample-key", default=None, help="Sample ID column in .obs (default: auto-detect)")
    p.add_argument("--group-key", default=None, help="Group/condition column in .obs (default: auto-detect)")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Load data
    adata = sc.read_h5ad(args.input)

    # Auto-detect columns if not provided
    cluster_key = args.cluster_key or pick_first_present(
        adata.obs,
        candidates=["leiden", "leiden_post", "louvain", "cluster", "clusters"],
        required=True,
        what="cluster key"
    )

    celltype_key = args.celltype_key or pick_first_present(
        adata.obs,
        candidates=["cell_typist_labels", "predicted_labels", "cell_type", "celltype", "CellTypist"],
        required=True,
        what="celltype key"
    )

    sample_key = args.sample_key or pick_first_present(
        adata.obs,
        candidates=["sample_id", "Sample", "sample", "GSM_ID", "GSM", "library_id", "sample_name"],
        required=False,
        what="sample key"
    )

    group_key = args.group_key or pick_first_present(
        adata.obs,
        candidates=["Group", "group", "condition", "disease_group", "diagnosis", "status"],
        required=False,
        what="group key"
    )

    print(f"[INFO] Using cluster_key = {cluster_key}")
    print(f"[INFO] Using celltype_key = {celltype_key}")
    if sample_key:
        print(f"[INFO] Using sample_key   = {sample_key}")
    else:
        print("[WARN] No sample_key found; sample composition plots will be skipped.")
    if group_key:
        print(f"[INFO] Using group_key    = {group_key}")
    else:
        print("[WARN] No group_key found; group composition plots will be skipped.")

    # 1) cluster × CellTypist label
    make_tables_and_plots(adata, cluster_key, celltype_key, tag="celltype", outdir=args.outdir)

    # 2) cluster × Sample ID (optional)
    if sample_key:
        make_tables_and_plots(adata, cluster_key, sample_key, tag="sample", outdir=args.outdir)

    # 3) cluster × Group (optional)
    if group_key:
        make_tables_and_plots(adata, cluster_key, group_key, tag="group", outdir=args.outdir)

    print(f"[DONE] CSVs and PNGs written to: {args.outdir}")

if __name__ == "__main__":
    main()