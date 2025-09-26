#!/usr/bin/env python
# coding: utf-8

# ==============================
# GSE290695_FINAL
# ==============================

# ==============================
# STEP 01: import libraries
# ==============================

import os
import glob
import pickle
import zipfile

import numpy as np
import pandas as pd

import scanpy as sc
import anndata as ad
from scipy import sparse
from scipy.sparse import csr_matrix
from scipy.stats import median_abs_deviation

import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

import scrublet as scr
# import scanorama as scrma   # (loaded if you need later)

# ==============================
# STEP 02: Setting variables
# ==============================

wd = "/rds/projects/e/elhamsak-ibd-single-cell/"
proj = "/THESIS/GSE290695/"
dir_raws = "FW_GSE290695_raws/"
dir_res_plots = "plots/"
dir_res_h5ad = "h5ad/"
raw_file_suffix = ".zip"

# ==============================
# STEP 03: Settings (logging/figures)
# ==============================

sc.settings.verbosity = 3  # errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()

mpl.rcParams["figure.dpi"] = 80
mpl.rcParams["figure.facecolor"] = "white"
mpl.rcParams["savefig.facecolor"] = "white"
mpl.rcParams["axes.facecolor"] = "white"
mpl.rcParams["axes.edgecolor"] = "black"
mpl.rcParams["axes.grid"] = False

sc.settings.figdir = wd + proj + dir_res_plots
os.makedirs(sc.settings.figdir, exist_ok=True)

# ==============================
# STEP 04: Paths sanity check
# ==============================

print("Current working directory:", os.getcwd())
print("Files in directory:", os.listdir())

proj_dir_raws = os.path.join(wd, proj.strip("/"), dir_raws.strip("/"))
zip_files_raw = [f for f in os.listdir(proj_dir_raws) if f.endswith(".zip") and not f.startswith(".")]
zip_paths_raw = [os.path.join(proj_dir_raws, f) for f in zip_files_raw]
print("Found zip files:", zip_paths_raw)

# ==============================
# STEP 05: Unzip raw archive
# ==============================

zip_path = "GSE290695_RAW.zip"
output_dir = "FW_GSE290695_raws_1"

with zipfile.ZipFile(zip_path, "r") as zip_ref:
    zip_ref.extractall(output_dir)

print("Folders in output_dir:", os.listdir(output_dir))
print("ZIP extracted successfully.")

# ==============================
# STEP 06: Sample metadata
# ==============================

metadata_path = "sample_to_donor.csv"
meta_df = pd.read_csv(metadata_path)

meta_df.columns = meta_df.columns.str.strip()
for c in ["GSM_ID", "tissue", "donor_id"]:
    assert c in meta_df.columns, f"Missing column {c} in {metadata_path}"

meta_df["GSM_ID"] = meta_df["GSM_ID"].astype(str).str.strip()
meta_df = meta_df.set_index("GSM_ID")

print(meta_df.head())
print(meta_df["donor_id"].value_counts())

# ==============================
# STEP 07: List sample folders
# ==============================

print("output_dir:", output_dir)
print("Contents of output_dir:", os.listdir(output_dir))
for item in os.listdir(output_dir):
    print("*", item, "->", os.path.join(output_dir, item))

adatas = []
sample_root = os.path.join(output_dir, "GSE290695_RAW")
sample_folders = sorted(glob.glob(os.path.join(sample_root, "GSM*")))
print("Sample folders found:", len(sample_folders))

for folder in sample_folders:
    sample_id = os.path.basename(folder)
    print(f" Reading: {sample_id}")

    X = sc.read_mtx(os.path.join(folder, "matrix.mtx")).T  # AnnData with .X (sparse)
    genes = pd.read_csv(os.path.join(folder, "features.tsv"), header=None, sep="\t")
    barcodes = pd.read_csv(os.path.join(folder, "barcodes.tsv"), header=None)

    ad_obj = sc.AnnData(X.X)
    ad_obj.var_names = genes[1].astype(str)
    ad_obj.obs_names = barcodes[0].astype(str)
    ad_obj.var_names_make_unique()
    ad_obj.obs["GSM_ID"] = sample_id

    if sample_id in meta_df.index:
        ad_obj.obs["donor_id"] = meta_df.loc[sample_id, "donor_id"]
        ad_obj.obs["tissue"] = meta_df.loc[sample_id, "tissue"]
    else:
        ad_obj.obs["donor_id"] = "Unknown"
        ad_obj.obs["tissue"] = "Unknown"

    ad_obj.obs["donor_id"] = ad_obj.obs["donor_id"].astype("category")
    ad_obj.obs["tissue"] = ad_obj.obs["tissue"].astype("category")

    adatas.append(ad_obj)
    print(f" Loaded {len(adatas)} samples")

loaded = {a.obs["GSM_ID"].iloc[0] for a in adatas}
missing = set(meta_df.index) - loaded
print("Missing GSMs from data:", missing)

print(pd.concat([a.obs[["GSM_ID", "donor_id", "tissue"]].head(1) for a in adatas]))

# ==============================
# STEP 08: QC plotting filenames
# ==============================

files_pl_vp_pct_cmt = []
files_pl_vp_tc = []
files_pl_vp_ngc = []
files_pl_vp_pct_cercc = []
files_pl_sc_col_pct_mt = []
files_pl_sc_col_pct_ercc = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    files_pl_vp_pct_cmt.append(f"{sample_id}_QC_plot_pct_counts_mt.png")
    files_pl_vp_tc.append(f"{sample_id}_QC_plot_total_counts.png")
    files_pl_vp_ngc.append(f"{sample_id}_QC_plot_n_genes_by_counts.png")
    files_pl_vp_pct_cercc.append(f"{sample_id}_QC_plot_pct_counts_ercc.png")
    files_pl_sc_col_pct_mt.append(f"{sample_id}_QC_scatter_col_pct_counts_mt.png")
    files_pl_sc_col_pct_ercc.append(f"{sample_id}_QC_scatter_col_pct_counts_ercc.png")

# ==============================
# STEP 09: QC + filtering function
# ==============================

from pathlib import Path
import scipy.io

def pp_qc_fun(raw_paths, pid, donor_id, tissue, vp1, vp2, vp3, vp4, sc1, sc2):
    raw_path = Path(raw_paths)

    matrix = scipy.io.mmread(raw_path / "matrix.mtx").tocsr().T
    genes = pd.read_csv(raw_path / "features.tsv", header=None, sep="\t")
    if genes.shape[1] == 3:
        genes.columns = ["gene_id", "gene_symbols", "feature_type"]
    elif genes.shape[1] == 2:
        genes.columns = ["gene_id", "gene_symbols"]
    else:
        raise ValueError(f"Unexpected number of columns in features.tsv: {genes.shape[1]}")
    barcodes = pd.read_csv(raw_path / "barcodes.tsv", header=None)[0].tolist()

    adata = sc.AnnData(X=matrix)
    adata.var_names = genes["gene_symbols"]
    adata.obs_names = barcodes
    adata.var_names_make_unique()

    adata.obs["sample_id"] = pid
    adata.obs["donor_id"] = donor_id
    adata.obs["tissue"] = tissue
    adata.obs["donor_id"] = adata.obs["donor_id"].astype("category")
    adata.obs["tissue"] = adata.obs["tissue"].astype("category")

    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ercc"] = adata.var_names.str.upper().str.startswith("ERCC")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ercc", "ribo", "hb"], inplace=True, log1p=True)

    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)

    unwanted = adata.var["mt"] | adata.var["ribo"] | adata.var["hb"] | adata.var["ercc"]
    adata = adata[:, ~unwanted].copy()

    sc.pl.violin(adata, ["total_counts"], jitter=0.4, save=f"_{vp1}")
    sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.4, save=f"_{vp2}")
    sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.4, save=f"_{vp3}")
    sc.pl.violin(adata, ["pct_counts_ercc"], jitter=0.4, save=f"_{vp4}")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", save=f"_{sc1}")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_ercc", save=f"_{sc2}")

    return adata

# ==============================
# STEP 10: Run QC per sample
# ==============================

adata_lst = []

for i in range(len(sample_folders)):
    raw_path = sample_folders[i]
    sample_id = os.path.basename(raw_path)

    print(f"Processing sample: {sample_id}")
    print(f" Raw path: {raw_path}")
    print(f" VP1: {files_pl_vp_tc[i]}")
    print(f" VP2: {files_pl_vp_ngc[i]}")
    print(f" VP3: {files_pl_vp_pct_cmt[i]}")
    print(f" VP4: {files_pl_vp_pct_cercc[i]}")
    print(f" SC1: {files_pl_sc_col_pct_mt[i]}")
    print(f" SC2: {files_pl_sc_col_pct_ercc[i]}")

    donor = meta_df.loc[sample_id, "donor_id"] if sample_id in meta_df.index else "Unknown"
    tiss = meta_df.loc[sample_id, "tissue"] if sample_id in meta_df.index else "Unknown"

    adata_ps = pp_qc_fun(
        raw_path,
        sample_id,
        donor,
        tiss,
        files_pl_vp_tc[i],
        files_pl_vp_ngc[i],
        files_pl_vp_pct_cmt[i],
        files_pl_vp_pct_cercc[i],
        files_pl_sc_col_pct_mt[i],
        files_pl_sc_col_pct_ercc[i],
    )

    adata_lst.append(adata_ps)

# quick previews
adata = adata_lst[0]
print(adata.obs["sample_id"].unique(), adata.n_obs, adata.n_vars)
for ad in adata_lst:
    print(ad.shape)
for i, ad in enumerate(adata_lst):
    sid = ad.obs["sample_id"].iloc[0]
    print(f" Sample: {sid} | Cells: {ad.n_obs} | Genes: {ad.n_vars}")

# save individual preprocessed
output_dir = "GSE290695_preprocessed_1/h5ad_preprocessed"
os.makedirs(output_dir, exist_ok=True)
for ad in adata_lst:
    sid = ad.obs["sample_id"].iloc[0]
    ad.write(os.path.join(output_dir, f"{sid}_preprocessed.h5ad"))

with open("adata_4th_preprocessed.pkl", "wb") as f:
    pickle.dump(adata_lst, f)

# ==============================
# STEP 11: MAD outlier filtering
# ==============================

def is_outlier(adata, metric: str, nmads: int = 5):
    M = adata.obs[metric]
    med = np.median(M)
    mad = median_abs_deviation(M)
    return (M < med - nmads * mad) | (M > med + nmads * mad)

with open("adata_4th_preprocessed.pkl", "rb") as f:
    adata_lst = pickle.load(f)

meta_df_csv = pd.read_csv("sample_to_donor.csv")
gsm_ids = meta_df_csv["GSM_ID"].tolist()

os.makedirs("filtered_2_h5ad", exist_ok=True)
os.makedirs("outlier_2_summary", exist_ok=True)

summary = []
for i, ad in enumerate(adata_lst):
    gsm_id = gsm_ids[i]
    ad.obs["outlier"] = (
        is_outlier(ad, "log1p_total_counts", 5)
        | is_outlier(ad, "log1p_n_genes_by_counts", 5)
        | is_outlier(ad, "pct_counts_mt", 3)
        | is_outlier(ad, "pct_counts_in_top_50_genes", 5)
    )
    total_cells = ad.n_obs
    n_outliers = ad.obs["outlier"].sum()
    n_retained = total_cells - n_outliers
    percent_outliers = round(n_outliers / total_cells * 100, 2)

    filtered = ad[~ad.obs["outlier"]].copy()
    filtered.write(f"filtered_2_h5ad/{gsm_id}_filtered_2_.h5ad")

    summary.append([gsm_id, total_cells, n_outliers, n_retained, percent_outliers])

summary_df = pd.DataFrame(
    summary,
    columns=["GSM_ID", "total_cells", "outliers_removed", "cells_retained", "percent_outliers"],
)
summary_df.to_csv("outlier_2_summary/mad_outlier_2_summary.csv", index=False)
print("✅ Outlier filtering complete with GSM-based filenames.")

# sanity check
folder = "filtered_2_h5ad"
files = sorted([f for f in os.listdir(folder) if f.endswith(".h5ad")])
for f in files:
    ad_check = sc.read(os.path.join(folder, f))
    sid = ad_check.obs["sample_id"].unique()[0]
    print(f"Sample: {sid} | Cells: {ad_check.n_obs} | Genes: {ad_check.n_vars}")

# ==============================
# STEP 12: Scrublet doublet removal
# ==============================

input_dir = "filtered_2_h5ad"
out_dir = "filtered_2_h5ad_scrublet"
metadata_csv = "sample_to_donor.csv"
os.makedirs(out_dir, exist_ok=True)

meta = pd.read_csv(metadata_csv)

filtered_adatas = []
scrublet_summary = []
results = {}

file_list = [f for f in os.listdir(input_dir) if f.endswith(".h5ad")]

for fname in file_list:
    fpath = os.path.join(input_dir, fname)
    ad = sc.read(fpath)
    sample_id = ad.obs["sample_id"][0] if "sample_id" in ad.obs.columns else fname.replace(".h5ad", "")
    print(f"\n✅ Running Scrublet for {sample_id}...")

    try:
        if ad.raw is not None:
            counts_matrix = ad.raw.X.toarray() if not isinstance(ad.raw.X, np.ndarray) else ad.raw.X
        else:
            counts_matrix = ad.X.toarray() if not isinstance(ad.X, np.ndarray) else ad.X

        scrub = scr.Scrublet(counts_matrix)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()

        ad.obs["doublet_score"] = doublet_scores
        ad.obs["predicted_doublet"] = predicted_doublets
        ad.uns["scrublet_threshold"] = scrub.threshold_

        ad_filtered = ad[~ad.obs["predicted_doublet"], :].copy()
        filtered_adatas.append(ad_filtered)

        out_file = os.path.join(out_dir, fname.replace(".h5ad", "_scrublet_filtered_2.h5ad"))
        ad_filtered.write(out_file)

        total_cells = ad.n_obs
        doublets_detected = int(np.sum(predicted_doublets))
        scrublet_summary.append(
            {
                "sample_id": sample_id,
                "total_cells": total_cells,
                "doublets_detected": doublets_detected,
                "cells_retained": int(total_cells - doublets_detected),
                "doublet_rate (%)": round(100 * doublets_detected / total_cells, 2),
            }
        )

        results[sample_id] = {"scores": doublet_scores, "threshold": scrub.threshold_}

    except Exception as e:
        print(f" Scrublet error in {sample_id}: {e}")
        scrublet_summary.append(
            {
                "sample_id": sample_id,
                "total_cells": ad.n_obs,
                "doublets_detected": "ERROR",
                "cells_retained": "ERROR",
                "doublet_rate (%)": "ERROR",
            }
        )

summary_df = pd.DataFrame(scrublet_summary)
summary_df.to_csv("scrublet_summary_per_sample_2.csv", index=False)
print("\n Summary saved to 'scrublet_summary_per_sample_2.csv'")

# quick shapes
scrublet_dir = "filtered_2_h5ad_scrublet"
files = sorted([f for f in os.listdir(scrublet_dir) if f.endswith(".h5ad")])
for f in files:
    adx = sc.read(os.path.join(scrublet_dir, f))
    sid = adx.obs["sample_id"].unique()[0] if "sample_id" in adx.obs else f
    print(f"Sample: {sid} | Cells: {adx.n_obs} | Genes: {adx.n_vars}")

# ==============================
# STEP 13: Scrublet histograms
# ==============================

cwd = os.getcwd()
input_dir_hist = "filtered_2_h5ad_scrublet"
summary_csv = "scrublet_summary_per_sample_2.csv"
output_dir_hist = "scrublet_histograms_2"
os.makedirs(output_dir_hist, exist_ok=True)

summary_df = pd.read_csv(summary_csv)

all_scores = []
for fname in os.listdir(input_dir_hist):
    if not fname.endswith(".h5ad"):
        continue
    sample_id = fname.replace(".h5ad", "")
    fpath = os.path.join(input_dir_hist, fname)
    try:
        ad = sc.read(fpath)
        scores = ad.obs["doublet_score"]
        threshold = ad.uns.get("scrublet_threshold", None)
        all_scores.extend(scores.tolist())

        plt.figure(figsize=(6, 4))
        sns.histplot(scores, bins=50, kde=False, color="skyblue", edgecolor="black")
        if threshold is not None:
            plt.axvline(threshold, color="red", linestyle="--", label=f"Threshold: {threshold:.2f}")
            plt.legend()
        plt.title(f"Scrublet Doublet Scores: {sample_id}")
        plt.xlabel("Doublet Score")
        plt.ylabel("Cell Count")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_hist, f"{sample_id}_doublet_score_hist.png"), dpi=150)
        plt.close()
    except Exception as e:
        print(f" Error in {sample_id}: {e}")

plt.figure(figsize=(7, 5))
sns.histplot(all_scores, bins=100, kde=False, color="darkgreen", edgecolor="black")
plt.title("Global Scrublet Doublet Score Distribution")
plt.xlabel("Doublet Score")
plt.ylabel("Cell Count")
plt.tight_layout()
plt.savefig(os.path.join(output_dir_hist, "global_doublet_score_hist.png"), dpi=150)
plt.close()
print(" All histograms saved to:", output_dir_hist)

# ==============================
# STEP 14: Concatenate, save pre-HVG
# ==============================

filtered_files = sorted([f for f in os.listdir(scrublet_dir) if f.endswith(".h5ad")])
filtered_adatas = [sc.read(os.path.join(scrublet_dir, f)) for f in filtered_files]
adata_all = sc.concat(filtered_adatas, label="sample_id", index_unique="-")

print("Final concatenated shape:", adata_all.shape)

os.makedirs("qc_outputs", exist_ok=True)
adata_all.write("qc_outputs/adata_all_postqc_preHVG.h5ad")

# ==============================
# STEP 15: HVG filtering
# ==============================

print("▶ Loading data...")
adata = sc.read("qc_outputs/adata_all_postqc_preHVG.h5ad")

print("▶ Converting matrix to float32...")
if sparse.issparse(adata.X):
    adata.X = adata.X.astype(np.float32)
    print(" Sparse matrix dtype:", adata.X.dtype)
else:
    adata.X = adata.X.astype(np.float32)
    print(" Dense matrix dtype:", adata.X.dtype)

print(" Applying log1p transformation manually to avoid overflow...")
sc.pp.log1p(adata)

if sparse.issparse(adata.X):
    X_data = adata.X.data
    print("Inf (sparse):", np.isinf(X_data).sum())
    print("NaN (sparse):", np.isnan(X_data).sum())
else:
    print("Inf (dense):", np.isinf(adata.X).sum())
    print("NaN (dense):", np.isnan(adata.X).sum())

print(" Running HVG filtering...")
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000, batch_key="sample_id")
adata = adata[:, adata.var["highly_variable"]]
print(" Number of HVGs selected:", adata.shape[1])

print("▶ Saving filtered .h5ad...")
adata.write("qc_outputs/adata_all_HVGfiltered.h5ad")
print(" HVG filtering complete and saved.")

adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ercc"],
    jitter=0.4,
    multi_panel=True,
    save="_hvg_filtered_violin.png",
)

# quick info
adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")
print(" HVG-filtered AnnData shape:", adata.shape)
print("\nTop 10 HVG gene names:")
print(adata.var_names[:10])

# scale
adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")
sc.pp.scale(adata, zero_center=True, max_value=10)
adata.write("qc_outputs/adata_all_HVGscaled.h5ad")
print(" Scaling complete. Scaled data saved to 'adata_all_HVGscaled.h5ad'")
print("Scaled AnnData shape:", adata.shape)

# ensure metadata types + save meta h5ad
adata_all = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")
for c in ["sample_id", "donor_id", "tissue"]:
    assert c in adata_all.obs, f"Missing {c} in obs"
    adata_all.obs[c] = adata_all.obs[c].astype("category")

adata_all.write("qc_outputs/adata_all_HVGfiltered_meta.h5ad")
print("Saved qc_outputs/adata_all_HVGfiltered_meta.h5ad")

# ==============================
# STEP 16: PCA/UMAP/Leiden + plots
# ==============================

adata_all = sc.read("qc_outputs/adata_all_HVGfiltered_meta.h5ad")
print(adata_all.obs["donor_id"].value_counts(dropna=False).head())
print(adata_all.obs["tissue"].value_counts(dropna=False).head())

for c in ["donor_id", "tissue"]:
    if c not in adata_all.obs:
        raise KeyError(f"Missing '{c}' in adata_all.obs")
    adata_all.obs[c] = adata_all.obs[c].astype("category")

sc.settings.figdir = "figures"
os.makedirs(sc.settings.figdir, exist_ok=True)

print("Running PCA...")
sc.tl.pca(adata_all, svd_solver="arpack")
sc.pl.pca(adata_all, color="donor_id", save="_pca_by_donor_id.png")

print("Computing neighbors...")
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=40)

print("Computing UMAP...")
sc.tl.umap(adata_all)

print("Running Leiden clustering...")
sc.tl.leiden(adata_all, resolution=0.5, flavor="igraph", n_iterations=2, directed=False)

print("Saving UMAP plots...")
sc.pl.umap(
    adata_all,
    color="leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    frameon=False,
    save="_leiden.png",
    title="Leiden Clusters",
)
sc.pl.umap(
    adata_all,
    color="donor_id",
    legend_loc="right margin",
    frameon=False,
    save="_by_donor_id.png",
    title="UMAP by donor_id",
)
sc.pl.umap(
    adata_all,
    color="tissue",
    legend_loc="right margin",
    frameon=False,
    save="_by_tissue.png",
    title="UMAP by tissue",
)
sc.pl.umap(
    adata_all,
    color=["leiden", "donor_id", "tissue"],
    ncols=1,
    wspace=0.3,
    legend_loc="right margin",
    frameon=False,
    save="_leiden_donor_tissue.png",
    title=["Leiden", "donor_id", "tissue"],
)

adata_all.write("qc_outputs/adata_all_clustered_umap.h5ad")
print(" Analysis complete. Output saved.")

# ==============================
# STEP 17: Expression matrix CSV
# ==============================

adata_all = sc.read("qc_outputs/adata_all_clustered_umap.h5ad")
expression_df = pd.DataFrame(
    adata_all.X.toarray() if sparse.issparse(adata_all.X) else adata_all.X,
    index=adata_all.obs_names,
    columns=adata_all.var_names,
)
expression_df.to_csv("qc_outputs/adata_all_hvg_filtered_expression_matrix_1.csv")
print(" Expression matrix saved as 'adata_all_hvg_filtered_expression_matrix_1.csv'.")

