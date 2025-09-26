#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Batch enrichment runner for UC vs Healthy (and other group) UP gene lists.

Inputs  : TXT files with one gene symbol per line (e.g., Group_UC_UP.txt)
Outputs : CSVs + barplots per gene set library

Requires: gseapy, pandas, numpy, matplotlib
Install : conda install -c conda-forge gseapy matplotlib pandas
"""

import argparse
import os
import glob
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless
import matplotlib.pyplot as plt

import gseapy as gp


GENE_SETS = [
    "GO_Biological_Process_2021",
    "KEGG_2021_Human",
    "Reactome_2022",
    "WikiPathways_2023_Human",
    "MSigDB_Hallmark_2020",
]


def read_gene_list(path):
    """Read a 1-col text file of gene symbols."""
    genes = []
    with open(path, "r") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    # de-duplicate, keep case but strip spaces
    genes = pd.Series(genes).str.strip().dropna().unique().tolist()
    return genes


def sanitize(name):
    """Safe filename chunk."""
    return re.sub(r"[^A-Za-z0-9_.\-]+", "_", str(name))


def run_enrichr(gene_list, gene_sets, out_csv, out_png, title, cutoff=0.05):
    """
    Run Enrichr (gseapy) and write CSV + a horizontal barplot of top 20 terms.
    """
    if len(gene_list) == 0:
        print(f"[WARN] Empty gene list for {title}; skipping.")
        return None

    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism="Human",
        cutoff=cutoff,          # adjusted P-value cutoff
        no_plot=True,           # we'll plot ourselves
        outdir=None,            # we handle outputs
        verbose=False,
    )

    # gseapy returns one dataframe per gene set in enr.results (concatenated)
    res = enr.results.copy()
    if res is None or res.empty:
        print(f"[INFO] No enriched terms passed cutoff for {title}.")
        return None

    # Sort nicely: by Adjusted P-value then Combined Score
    sort_cols = [c for c in ["Adjusted P-value", "P-value", "Combined Score"] if c in res.columns]
    res = res.sort_values(sort_cols, ascending=[True] + [False]*(len(sort_cols)-1))
    res.to_csv(out_csv, index=False)

    # Build a compact barplot (top 20)
    topN = res.head(20).copy()
    # Choose a score column for the bars (prefer Combined Score if present)
    score_col = "Combined Score" if "Combined Score" in topN.columns else (
        "-log10(adjP)" if "-log10(adjP)" in topN.columns else None
    )
    if score_col is None:
        # Fallback: use -log10(Adjusted P-value)
        topN["-log10(adjP)"] = -np.log10(topN["Adjusted P-value"].astype(float))
        score_col = "-log10(adjP)"

    # Shorten long term names for readability
    topN["Term_short"] = topN["Term"].str.replace("_", " ").str.slice(0, 70)

    plt.figure(figsize=(8, max(3, 0.35*len(topN))))
    plt.barh(topN["Term_short"][::-1], topN[score_col][::-1])
    plt.xlabel(score_col)
    plt.title(title, pad=10)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    return res


def main():
    ap = argparse.ArgumentParser(description="Run GO/KEGG pathway enrichment on group UP gene lists.")
    ap.add_argument("--input-dir", type=str, default="enrichment_lists_by_group",
                    help="Folder containing Group_*_UP.txt files.")
    ap.add_argument("--pattern", type=str, default="Group_*_UP.txt",
                    help="Glob pattern to match gene list files.")
    ap.add_argument("--outdir", type=str, default="enrichment_results",
                    help="Where to write CSVs and plots.")
    ap.add_argument("--cutoff", type=float, default=0.05,
                    help="Adjusted p-value cutoff.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    files = sorted(glob.glob(os.path.join(args.input_dir, args.pattern)))
    if not files:
        print(f"[ERROR] No files found matching {os.path.join(args.input_dir, args.pattern)}")
        return

    print(f"[INFO] Found {len(files)} gene lists.")
    print(f"[INFO] Gene sets: {', '.join(GENE_SETS)}")

    for f in files:
        base = os.path.splitext(os.path.basename(f))[0]      # e.g., Group_UC_UP
        print(f"\n[RUN] {base}")

        genes = read_gene_list(f)
        print(f"[INFO] {len(genes)} genes in {base}")

        # Run once across all requested libraries (combined result comes back)
        out_csv  = os.path.join(args.outdir, f"{sanitize(base)}_enrichment.csv")
        out_plot = os.path.join(args.outdir, f"{sanitize(base)}_enrichment_top20.png")
        title    = f"{base} (GO/KEGG/Reactome/WikiPathways/Hallmark)"

        run_enrichr(
            gene_list=genes,
            gene_sets=GENE_SETS,
            out_csv=out_csv,
            out_png=out_plot,
            title=title,
            cutoff=args.cutoff,
        )

    print("\n[OK] Enrichment complete. Results in:", args.outdir)


if __name__ == "__main__":
    main()
