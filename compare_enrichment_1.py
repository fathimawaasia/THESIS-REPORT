#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Columns gseapy.enrichr typically returns
CORE_COLS = [
    "Gene_set", "Term", "Overlap", "P-value",
    "Adjusted P-value", "Old P-value", "Old Adjusted P-value",
    "Z-score", "Combined Score", "Genes"
]

def read_enrich(path, label):
    df = pd.read_csv(path)
    # Keep only the columns we need and add label prefix
    for c in ["Term", "Gene_set", "Adjusted P-value", "Combined Score"]:
        if c not in df.columns:
            raise ValueError(f"Missing column '{c}' in {path}")
    # Clean term text (nicer merging)
    df["Term"] = (df["Term"].astype(str)
                  .str.replace("_", " ")
                  .str.replace(r"\s+", " ", regex=True)
                  .str.strip())
    # Build compact numeric score for plotting: -log10 padj
    df[f"-log10padj_{label}"] = -np.log10(df["Adjusted P-value"].astype(float).clip(lower=1e-300))
    df = df.rename(columns={
        "Adjusted P-value": f"padj_{label}",
        "Combined Score":   f"combined_{label}",
        "Gene_set":         "Library"
    })
    keep = ["Term", "Library", f"padj_{label}", f"combined_{label}", f"-log10padj_{label}"]
    return df[keep].copy()

def main():
    ap = argparse.ArgumentParser(description="Compare UC vs HC enrichment CSVs.")
    ap.add_argument("--uc", required=True, help="Path to Group_UC_UP_enrichment.csv")
    ap.add_argument("--hc", required=True, help="Path to Group_HC_UP_enrichment.csv")
    ap.add_argument("--outdir", default="results/compare_UC_vs_HC",
                    help="Output directory for merged CSV and plots.")
    ap.add_argument("--padj", type=float, default=0.05, help="Adjusted p-value cutoff.")
    ap.add_argument("--topn", type=int, default=20, help="Top N to plot for unique terms.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    df_uc = read_enrich(args.uc, "UC")
    df_hc = read_enrich(args.hc, "HC")

    # Merge on Library+Term (most robust)
    merged = pd.merge(df_uc, df_hc, how="outer", on=["Library", "Term"])

    # Significance flags
    merged["sig_UC"] = merged["padj_UC"].lt(args.padj)
    merged["sig_HC"] = merged["padj_HC"].lt(args.padj)

    # Category: UC_only / HC_only / Common / None
    conds = [
        merged["sig_UC"] & ~merged["sig_HC"],
        ~merged["sig_UC"] & merged["sig_HC"],
        merged["sig_UC"] & merged["sig_HC"]
    ]
    cats = ["UC_only", "HC_only", "Common"]
    merged["enrichment_category"] = np.select(conds, cats, default="None")

    # Delta scores for ranking (where present)
    merged["delta_-log10padj_UC_minus_HC"] = merged["-log10padj_UC"].fillna(0) - merged["-log10padj_HC"].fillna(0)
    merged["delta_combined_UC_minus_HC"]   = merged["combined_UC"].fillna(0)   - merged["combined_HC"].fillna(0)

    # Save merged table
    out_csv = os.path.join(args.outdir, "UC_vs_HC_enrichment_merged.csv")
    merged.sort_values(["Library", "enrichment_category", "delta_-log10padj_UC_minus_HC"], ascending=[True, True, False]).to_csv(out_csv, index=False)
    print(f"[OK] Wrote {out_csv}")

    # --- Plots ---

    # 1) Scatter for Common terms
    common = merged.query("enrichment_category == 'Common'").copy()
    if not common.empty:
        plt.figure(figsize=(6.5, 6))
        plt.scatter(common["-log10padj_HC"], common["-log10padj_UC"], s=10, alpha=0.6)
        lim = max(common["-log10padj_HC"].max(), common["-log10padj_UC"].max()) * 1.05
        plt.plot([0, lim], [0, lim], ls="--", lw=1, color="gray")
        plt.xlabel("-log10(adjP) HC")
        plt.ylabel("-log10(adjP) UC")
        plt.title("Common enriched terms (significant in both)")
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, "scatter_common_terms_UC_vs_HC.png"), dpi=300)
        plt.close()

    # 2) Barplots: UC-only and HC-only (top N by -log10padj)
    for label, colflag, scorecol, fname, ttl in [
        ("UC", "sig_UC", "-log10padj_UC", "bar_UC_only_top.png",  "UC-only enriched terms"),
        ("HC", "sig_HC", "-log10padj_HC", "bar_HC_only_top.png",  "HC-only enriched terms"),
    ]:
        sub = merged.query(f"enrichment_category == '{label}_only'").copy()
        if sub.empty:
            continue
        sub["term_label"] = (sub["Library"] + ": " + sub["Term"]).str.slice(0, 70)
        sub = sub.sort_values(scorecol, ascending=False).head(args.topn)
        plt.figure(figsize=(8, max(3, 0.35*len(sub))))
        plt.barh(sub["term_label"][::-1], sub[scorecol][::-1])
        plt.xlabel("-log10(adjP)")
        plt.title(ttl)
        plt.tight_layout()
        plt.savefig(os.path.join(args.outdir, fname), dpi=300)
        plt.close()

    # 3) Summary counts
    counts = (merged["enrichment_category"]
              .value_counts()
              .rename_axis("category")
              .reset_index(name="n"))
    counts.to_csv(os.path.join(args.outdir, "summary_counts.csv"), index=False)
    print(counts)

if __name__ == "__main__":
    main()