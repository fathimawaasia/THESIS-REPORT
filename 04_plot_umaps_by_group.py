#!/usr/bin/env python3
import os
import pandas as pd
import scanpy as sc

# --- inputs ---
SUBSETS = {
    "DC2": "DC2_subset_reclustered_fullX.h5ad",
    "DP":  "DP_subset_reclustered_fullX.h5ad",
}

# --- outputs ---
OUTDIR = "FIGS_UMAP"
os.makedirs(OUTDIR, exist_ok=True)
sc.settings.figdir = OUTDIR          # force all scanpy saves here
sc.settings.autoshow = False

# columns we want if present (first available per group is chosen)
CANDIDATE_COLS = [
    ["leiden_sub"],              # clusters
    ["group", "disease_group"],  # cohort label
    ["tissue"],
    ["dataset_id"],
]

def pick_existing_cols(ad, candidate_lists):
    cols = []
    for cand in candidate_lists:
        for c in cand:
            if c in ad.obs.columns:
                cols.append(c)
                break
    # dedupe, keep order
    seen, ordered = set(), []
    for c in cols:
        if c not in seen:
            seen.add(c); ordered.append(c)
    return ordered

def _clean_palette_if_needed(ad, col):
    """
    Scanpy caches palettes in ad.uns[f'{col}_colors'].
    If cached palette length != number of categories, drop it to let scanpy rebuild.
    """
    key = f"{col}_colors"
    if col not in ad.obs:
        return
    ser = ad.obs[col]
    if pd.api.types.is_categorical_dtype(ser.dtype):
        cats = list(ser.cat.categories)
        if key in ad.uns:
            colors = ad.uns[key]
            if not isinstance(colors, (list, tuple)) or len(colors) != len(cats):
                del ad.uns[key]
    else:
        if key in ad.uns:
            del ad.uns[key]

def save_single(ad, color, tag):
    _clean_palette_if_needed(ad, color)
    legend_loc = "on data" if color == "leiden_sub" else "right"
    sc.pl.umap(
        ad, color=color, frameon=False, legend_loc=legend_loc,
        show=False, save=f"_{tag}_{color}.png"
    )

def save_multi(ad, colors, tag):
    for c in colors:
        _clean_palette_if_needed(ad, c)
    sc.pl.umap(
        ad, color=colors, ncols=2, frameon=False, show=False,
        save=f"_{tag}_multi.png"
    )

def main():
    for tag, path in SUBSETS.items():
        if not os.path.exists(path):
            print(f"[SKIP] {tag}: {path} not found")
            continue

        ad = sc.read(path)
        cols = pick_existing_cols(ad, CANDIDATE_COLS)
        if not cols:
            print(f"[WARN] {tag}: no known obs columns to plot")
            continue

        # individual panels
        for c in cols:
            save_single(ad, c, tag)

        # combined 2x2 (or 2xn) panel
        save_multi(ad, cols, tag)

        print(f"[OK] Saved UMAPs for {tag} -> {OUTDIR}")

if __name__ == "__main__":
    main()



