#!/usr/bin/env python3
# Export BayesPrism CPM reference from a finished scANVI-refined AnnData

import os, sys, numpy as np, pandas as pd
import scanpy as sc
from pathlib import Path
from scipy import sparse
import matplotlib.pyplot as plt

# ----------- CONFIG (adjust paths if needed) -----------
INP_AD  = "adata_all_scanvi_refined_qc.h5ad"     # your saved scANVI-refined file
OUTD    = "scanvi_for_bayesprism"                # output folder
FINAL_LABEL = "cell_type_final"                  # scANVI-refined labels
SEED_LABEL  = "ct_seed"                          # CellTypist seed labels (optional for plots)
BATCH_KEY   = "dataset_id"                       # for composition plot
COUNTS_CANDIDATES = ["counts","raw_counts","X_counts"]

# QC thresholds (start lenient, tighten later if desired)
CONF_MIN_SEED  = 0.50
CONF_MIN_FINAL = 0.60
MIN_DATASETS   = 2
MIN_PER_CT_TOTAL   = 100
MIN_PER_CT_DATASET = 30
CORR_MIN       = 0.85
CPM_GENE_MIN   = 1.0
# -------------------------------------------------------

Path(OUTD).mkdir(parents=True, exist_ok=True)
SIG_TSV_FILTERED   = f"{OUTD}/reference_scanvi_meanCPM.filtered.tsv.gz"
SIG_TSV_UNFILTERED = f"{OUTD}/reference_scanvi_meanCPM.unfiltered.tsv.gz"
QC_REPORT          = f"{OUTD}/qc_report_scanvi.csv"

print(f"[LOAD] {INP_AD}")
ad = sc.read(INP_AD)
print(ad)

for key in [FINAL_LABEL, BATCH_KEY]:
    if key not in ad.obs:
        sys.exit(f"ERROR: obs['{key}'] not found in {INP_AD}")

if "scanvi_confidence" not in ad.obs:
    ad.obs["scanvi_confidence"] = 1.0  # fallback

# -------- find counts matrix ----------
Xc = None
for cand in COUNTS_CANDIDATES:
    if cand in ad.layers:
        Xc = ad.layers[cand]
        print(f"[COUNTS] Using layers['{cand}']")
        break
if Xc is None:
    if ad.raw is not None:
        print("[COUNTS] Using ad.raw.X (ensure these are counts, not log)")
        Xc = ad.raw.X
    else:
        sys.exit("ERROR: No counts layer found. Attach counts to ad.layers['counts'] and re-run.")

if not sparse.issparse(Xc):
    Xc = sparse.csr_matrix(Xc)

# ------------- CPM transform -------------
lib = np.asarray(Xc.sum(axis=1)).ravel()
lib[lib == 0] = 1.0
X_cpm = Xc.multiply((1.0 / (lib / 1e6))[:, None]).tocsr()

labels  = ad.obs[FINAL_LABEL].astype(str).values
batches = ad.obs[BATCH_KEY].astype(str).values
conf    = ad.obs["scanvi_confidence"].astype(float).values
genes   = ad.var_names.to_list()

def mean_cpm(mask):
    if mask.sum() == 0: return None
    v = X_cpm[mask,:].mean(axis=0)
    return v.A1 if sparse.issparse(v) else np.asarray(v).ravel()

def safe_corr(a,b):
    if a is None or b is None: return np.nan
    if a.size != b.size: return np.nan
    if np.allclose(a, a[0]) or np.allclose(b, b[0]): return np.nan
    return np.corrcoef(a,b)[0,1]

cts = sorted(set(labels) - {"Unknown"})
qc_rows, kept, ref_unf_cols, ref_unf_mat = [], [], [], []

for ct in cts:
    idx_all = (labels == ct)
    n_total = int(idx_all.sum())

    # unfiltered mean (useful fallback)
    mu_all = mean_cpm(idx_all)
    if mu_all is not None:
        ref_unf_cols.append(ct); ref_unf_mat.append(mu_all)

    if n_total < MIN_PER_CT_TOTAL:
        qc_rows.append([ct, n_total, 0, np.nan, False, "total_cells_lt_min"]); continue

    batches_ct, counts_ct = np.unique(batches[idx_all], return_counts=True)
    valid_batches = list(batches_ct[counts_ct >= MIN_PER_CT_DATASET])
    if len(valid_batches) < MIN_DATASETS:
        qc_rows.append([ct, n_total, len(valid_batches), np.nan, False, "datasets_lt_min"]); continue

    means = {}
    for b in valid_batches:
        idx = (labels == ct) & (batches == b) & (conf >= CONF_MIN_FINAL)
        means[b] = mean_cpm(idx)

    corrs = []
    for i in range(len(valid_batches)):
        for j in range(i+1, len(valid_batches)):
            r = safe_corr(means[valid_batches[i]], means[valid_batches[j]])
            if not np.isnan(r): corrs.append(r)

    min_corr = float(np.nanmin(corrs)) if corrs else np.nan
    passed = (len(corrs) > 0) and (min_corr >= CORR_MIN)
    qc_rows.append([ct, n_total, len(valid_batches), min_corr, passed,
                    "ok" if passed else "low_pairwise_corr"])
    if passed: kept.append(ct)

pd.DataFrame(qc_rows, columns=["cell_type","n_cells_total","n_datasets_valid","min_pairwise_corr","pass","reason"])\
  .to_csv(QC_REPORT, index=False)
print(f"[QC] Wrote: {QC_REPORT}")
print(f"[QC] Passed CTs: {len(kept)} / {len(cts)}")

# ---- UNFILTERED reference (always) ----
if ref_unf_cols:
    ref_unf = pd.DataFrame(np.vstack(ref_unf_mat).T, index=genes, columns=ref_unf_cols)
    ref_unf = ref_unf.loc[ref_unf.max(axis=1) >= CPM_GENE_MIN]
    ref_unf.to_csv(SIG_TSV_UNFILTERED, sep="\t", compression="gzip")
    print(f"[REF] UNFILTERED: {SIG_TSV_UNFILTERED}")

# ---- FILTERED reference (if any kept) ----
if kept:
    cols, mats = [], []
    for ct in kept:
        idx = (labels == ct) & (conf >= CONF_MIN_FINAL)
        v = mean_cpm(idx)
        if v is None: continue
        cols.append(ct); mats.append(v)
    ref = pd.DataFrame(np.vstack(mats).T, index=genes, columns=cols)
    ref = ref.loc[ref.max(axis=1) >= CPM_GENE_MIN]
    ref.to_csv(SIG_TSV_FILTERED, sep="\t", compression="gzip")
    print(f"[REF] FILTERED: {SIG_TSV_FILTERED}")
else:
    print("[REF] No CTs passed QC for FILTERED file (use UNFILTERED or relax thresholds).")

# ================== CLEAN SINGLE-PANEL PLOTS + COLOR HARMONIZATION ==================
import re, numpy as np, pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from pathlib import Path

# expected to exist already from your script:
# - AnnData object: ad
# - output directory: OUTD
# - keys in .obs: ct_seed (seed), celltypist_majority, celltypist_labels, cell_type_final, dataset_id
SEED_LABEL  = locals().get("SEED_LABEL", "ct_seed")
FINAL_LABEL = locals().get("FINAL_LABEL", "cell_type_final")
BATCH_KEY   = locals().get("BATCH_KEY", "dataset_id")
Path(OUTD).mkdir(parents=True, exist_ok=True)

# --- ensure a UMAP is available (prefer scANVI latent) ---
if "X_umap" not in ad.obsm:
    rep = "X_scANVI" if "X_scANVI" in ad.obsm else None
    if rep is not None:
        sc.pp.neighbors(ad, use_rep=rep, n_neighbors=15)
    else:
        sc.pp.neighbors(ad, n_neighbors=15)
    sc.tl.umap(ad, min_dist=0.3)

# --- helpers to get palettes and draw single UMAP with bottom legend ---
def _extract_palette(field, enforce=False):
    """Return (categories, colors) for a given obs field; if missing, create via a silent draw."""
    if field not in ad.obs:
        return [], []
    cats = list(ad.obs[field].astype("category").cat.categories)
    key  = f"{field}_colors"
    cols = ad.uns.get(key, None)
    if (cols is None or len(cols) < len(cats)) or enforce:
        sc.pl.umap(ad, color=field, legend_loc="none", show=False, save=None)  # silent; just to populate colors
        plt.close()
        cols = ad.uns.get(key, None)
    if cols is None:
        cols = [None] * len(cats)  # fallback (matplotlib default)
    return cats, list(cols)[:len(cats)]

def single_umap(field, out_png, palette_dict=None, order=None, title=None):
    if field not in ad.obs:
        print(f"[SKIP] {field} not in .obs"); return
    if order is not None:
        # keep only categories that exist in this field, ordered by 'order'
        existing = list(ad.obs[field].astype("category").cat.categories)
        ordered  = [c for c in order if c in existing] + [c for c in existing if c not in order]
        ad.obs[field] = ad.obs[field].astype("category").cat.set_categories(ordered)
    # colors (optionally override with provided palette)
    cats, cols = _extract_palette(field)
    if palette_dict:
        cols = [palette_dict.get(c, cols[i] if i < len(cols) else None) for i, c in enumerate(cats)]
        ad.uns[f"{field}_colors"] = cols

    sc.pl.umap(ad, color=field, legend_loc="none", show=False, save=None, title=title or field)
    handles = [mpatches.Patch(color=cols[i], label=str(cats[i])) for i in range(len(cats))]
    ncol = min(8, max(2, int(np.ceil(len(cats) / 10))))  # spread legend across rows
    plt.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.14),
               ncol=ncol, frameon=False, fontsize=6, title=field, title_fontsize=8)
    plt.subplots_adjust(bottom=0.25)
    plt.savefig(f"{OUTD}/{out_png}", dpi=300, bbox_inches="tight")
    plt.close()

# --- get final palette (authoritative colors) ---
final_cats, final_cols = _extract_palette(FINAL_LABEL, enforce=True)
final_palette = dict(zip(final_cats, final_cols)) if final_cols else {}

# --- harmonize CellTypist/seed colors to match cell_type_final ---
DEFAULT_GREY = "#D3D3D3"
JACCARD_MIN = 0.60  # token-overlap threshold for fuzzy matches

# add aliases for any common name mismatches you notice
ALIASES = {
    # "TFH": "Follicular helper T cells",
    # "pDC": "pDC",
    # "NK cells": "NK cells",
}

def _norm_tokens(s: str):
    s = s.lower().replace("+", " plus ")
    s = re.sub(r"[^a-z0-9]+", " ", s)
    toks = [t for t in s.split() if t not in {"cell", "cells"}]
    return set(toks)

def _best_match(label, ref_labels):
    t = _norm_tokens(label)
    best, best_score = None, 0.0
    for r in ref_labels:
        r_t = _norm_tokens(r)
        if not t or not r_t:
            continue
        inter = len(t & r_t); union = len(t | r_t)
        score = inter / union if union else 0.0
        if score > best_score:
            best, best_score = r, score
    return best, best_score

def make_palette_like_final(field, out_map_name):
    if field not in ad.obs:
        return None
    cats = list(ad.obs[field].astype("category").cat.categories)
    pal, rows = {}, []
    for c in cats:
        if c in ALIASES and ALIASES[c] in final_palette:
            pal[c] = final_palette[ALIASES[c]]
            rows.append([field, c, ALIASES[c], "alias", 1.0, pal[c]]); continue
        if c in final_palette:
            pal[c] = final_palette[c]
            rows.append([field, c, c, "exact", 1.0, pal[c]]); continue
        match, score = _best_match(c, final_cats)
        if match is not None and score >= JACCARD_MIN:
            pal[c] = final_palette[match]
            rows.append([field, c, match, f"jaccard≥{JACCARD_MIN}", score, pal[c]])
        else:
            pal[c] = DEFAULT_GREY
            rows.append([field, c, "(unmatched)", "grey", 0.0, pal[c]])
    ad.uns[f"{field}_colors"] = [pal[c] for c in cats]
    pd.DataFrame(rows, columns=["field","label","mapped_to_final","method","score","color"]).to_csv(
        f"{OUTD}/{out_map_name}", sep="\t", index=False
    )
    return pal

seed_pal = make_palette_like_final(SEED_LABEL, "palette_mapping_ct_seed.tsv") if SEED_LABEL in ad.obs else None
maj_pal  = make_palette_like_final("celltypist_majority", "palette_mapping_celltypist_majority.tsv") \
           if "celltypist_majority" in ad.obs else None
lab_pal  = make_palette_like_final("celltypist_labels", "palette_mapping_celltypist_labels.tsv") \
           if "celltypist_labels" in ad.obs else None

# --- draw single-panel UMAPs (harmonized colors) ---
if SEED_LABEL in ad.obs and seed_pal:
    single_umap(SEED_LABEL, "umap_ct_seed.png", palette_dict=seed_pal)
if "celltypist_majority" in ad.obs and maj_pal:
    single_umap("celltypist_majority", "umap_celltypist_majority.png", palette_dict=maj_pal)
if "celltypist_labels" in ad.obs and lab_pal:
    single_umap("celltypist_labels", "umap_celltypist_labels.png", palette_dict=lab_pal)
single_umap(FINAL_LABEL, "umap_cell_type_final.png", palette_dict=final_palette, order=final_cats, title="cell_type_final")

# --- composition by dataset (bottom, multi-column legend) + color key ---
comp = (
    ad.obs[[FINAL_LABEL, BATCH_KEY]]
    .assign(n=1)
    .groupby([BATCH_KEY, FINAL_LABEL])["n"].sum()
    .groupby(level=0).apply(lambda s: s / s.sum())
    .unstack(fill_value=0)
)
# same order/colors as final UMAP
comp = comp.reindex(columns=[c for c in final_cats if c in comp.columns])
colors = [final_palette[c] for c in comp.columns] if final_palette else None

# cell type -> color key
pd.Series(final_palette).rename("color").to_csv(f"{OUTD}/celltype_color_map.tsv", sep="\t")

n_datasets, n_ct = comp.shape
fig_w = max(10, 1.2 * n_datasets + 5); fig_h = 6
ax = comp.plot(kind="bar", stacked=True, figsize=(fig_w, fig_h), legend=False, color=colors, width=0.85)
ax.set_ylabel("Fraction"); ax.set_xlabel("dataset_id"); ax.set_title("scANVI composition by dataset")
ax.tick_params(axis="x", rotation=45, ha="right")

handles = [mpatches.Patch(color=(colors[i] if colors else None), label=ct) for i, ct in enumerate(comp.columns)]
ncol_bottom = min(8, max(2, int(np.ceil(n_ct / 10))))
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18),
          frameon=False, ncol=ncol_bottom, fontsize=8, title="Cell type", title_fontsize=9)
plt.subplots_adjust(bottom=0.28)
plt.savefig(f"{OUTD}/composition_by_dataset_scanvi_with_legend.png", dpi=300, bbox_inches="tight")
plt.close()
# =================================================================================================
