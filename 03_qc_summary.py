#!/usr/bin/env python3
"""
06_qc_summary.py  (robust to column names)
- For each subset H5AD, writes:
  *_counts_by_group.csv, *_counts_by_tissue.csv, *_counts_by_leiden.csv
  *_top10_markers_per_leiden.csv (scanpy rank_genes_groups on log1p if needed)
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np

SUBSETS = {
    "DC2": "DC2_subset_reclustered_fullX.h5ad",
    "DP":  "DP_subset_reclustered_fullX.h5ad",
}
OUTDIR = "QC_SUMMARY"
os.makedirs(OUTDIR, exist_ok=True)

def pick_col(ad, candidates, required=False):
    for c in candidates:
        if c in ad.obs.columns:
            return c
    if required:
        raise KeyError(f"None of the columns found in .obs: {candidates}")
    return None

def write_counts(ad, tag):
    ad.obs["n"] = 1

    # group (prefer 'group', fallback 'disease_group')
    gcol = pick_col(ad, ["group", "disease_group"])
    if gcol is not None:
        ad.obs.groupby(gcol, observed=False)["n"].sum().to_csv(f"{OUTDIR}/{tag}_counts_by_group.csv")
    else:
        # write an empty file with a note
        pd.Series(dtype=int, name="n").to_csv(f"{OUTDIR}/{tag}_counts_by_group.csv")

    # tissue (if present)
    tcol = pick_col(ad, ["tissue"])
    if tcol is not None:
        ad.obs.groupby(tcol, observed=False)["n"].sum().to_csv(f"{OUTDIR}/{tag}_counts_by_tissue.csv")
    else:
        pd.Series(dtype=int, name="n").to_csv(f"{OUTDIR}/{tag}_counts_by_tissue.csv")

    # leiden_sub (required for clustering summaries)
    lcol = pick_col(ad, ["leiden_sub"], required=True)
    ad.obs.groupby(lcol, observed=False)["n"].sum().to_csv(f"{OUTDIR}/{tag}_counts_by_leiden.csv")

def safe_markers(ad, tag, groupby="leiden_sub", topn=20):
    # Try to detect if X is counts-ish; if so, normalize/log before ranking on a temp copy.
    tmp = ad.copy()
    try:
        xmax = tmp.X.max() if hasattr(tmp.X, "max") else np.max(tmp.X)
    except Exception:
        xmax = 1.0
    if xmax is not None and xmax > 50:
        sc.pp.normalize_total(tmp, target_sum=1e4)
        sc.pp.log1p(tmp)
    sc.tl.rank_genes_groups(tmp, groupby=groupby, method="wilcoxon", use_raw=False, pts=True)
    res = sc.get.rank_genes_groups_df(tmp, group=None)
    res = res.sort_values(["group","scores"], ascending=[True, False])
    out = res.groupby("group").head(topn)
    out.to_csv(f"{OUTDIR}/{tag}_top20_markers_per_leiden.csv", index=False)

def main():
    for tag, fn in SUBSETS.items():
        ad = sc.read(fn)
        write_counts(ad, tag)
        safe_markers(ad, tag)
        print(f"[OK] wrote QC summaries for {tag}")

if __name__ == "__main__":
    main()



