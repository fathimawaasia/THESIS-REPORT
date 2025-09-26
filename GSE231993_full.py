#!/usr/bin/env python
# coding: utf-8
#==============================
# GSE231993_FINAL
#==============================

#==============================
# STEP 01: import libraries////
#==============================

import os
import numpy as np
import pandas as pd
import random as rd
import scanpy as sc
import anndata as ad
import seaborn as sns
import glob
import matplotlib.pyplot as pl
from scipy.sparse import csr_matrix
import scanorama as scrma
import pickle
from scipy.stats import median_abs_deviation
import scrublet as scr
from scipy import sparse
import matplotlib.pyplot as plt
import matplotlib as mpl
import zipfile
import scipy.io
from pathlib import Path

#==============================
# STEP 02: Setting variables
#==============================

wd = '/rds/projects/e/elhamsak-ibd-single-cell/'
proj = '/THESIS/GSE231993/'
dir_raws = 'FW_GSE231993_raws/'
dir_res_plots = 'plots/'
dir_res_h5ad = 'h5ad/'
raw_file_suffix = '.zip'

#==============================
# STEP 03: Settings
#==============================
#Set verbosity and logging
sc.settings.verbosity = 3 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
# Set global matplotlib figure parameters
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['figure.facecolor'] = 'white'
mpl.rcParams['savefig.facecolor'] = 'white'
mpl.rcParams['axes.facecolor'] = 'white'
mpl.rcParams['axes.edgecolor'] = 'black'
mpl.rcParams['axes.grid'] = False

#Set Scanpy figure output directory
sc.settings.figdir = wd + proj + dir_res_plots

#==========================================
# STEP 04: Check the current directory path
#==========================================
import os
print("Current working directory:", os.getcwd())
print("Files in directory:", os.listdir())

#==========================================
# STEP 05: Check the zip file path
#==========================================

# Raw directory
proj_dir_raws = os.path.join(wd, proj.strip('/'), dir_raws.strip('/'))

#List of zip files in the raw directory
zip_files_raw =[f for f in os.listdir(proj_dir_raws) if f.endswith('.zip') and not f.startswith('.')]

# Indices of all samples
sample_idcs = range(len(zip_files_raw))

# List of full paths to zip files
zip_paths_raw = [os.path.join(proj_dir_raws, f) for f in zip_files_raw]

print("Found zip files:", zip_paths_raw)

#==========================================
# STEP 06: Unzipping
#==========================================

zip_path = "GSE231993_RAW.zip"
output_dir = "FW_GSE231993_raws_1"

#Extract zip contents
with zipfile.ZipFile(zip_path, 'r') as zip_ref:
    zip_ref.extractall(output_dir)

# Confirm extraction
print("Folders found inside output_dir:")
print(os.listdir(output_dir))
print("ZIP extracted successfully.")

#==========================================
# STEP 07: Defining Sample Metadata
#==========================================

metadata_path = "sample_metadata.csv"

meta_df = pd.read_csv(metadata_path)
meta_df = meta_df.set_index("GSM_ID")

print("Metadata preview:")
print(meta_df.head())

print("output_dir:", output_dir)
print("Contents of output_dir:",os.listdir(output_dir))

for item in os.listdir(output_dir):
    print("*", item, "->", os.path.join(output_dir, item))

#======================================================================
# STEP 08: Loading 10X matrices and building Anndata objects per sample
#======================================================================

adatas = []
sample_root = os.path.join(output_dir, "GSE231993_RAW")
sample_folders = sorted(glob.glob(os.path.join(sample_root, "GSM*")))
print("Sample folders found:", len(sample_folders))
for folder in sample_folders:
    sample_id = os.path.basename(folder)
    print(f" Reading: {sample_id}")

    # Manually load files
    X = sc.read_mtx(os.path.join(folder, "matrix.mtx")).T
    genes = pd.read_csv(os.path.join(folder, "features.tsv"), header=None, sep="\t")
    barcodes = pd.read_csv(os.path.join(folder, "barcodes.tsv"), header=None)

    ad = sc.AnnData(X.X)   
    ad.var_names = genes[1].astype(str)
    ad.obs_names = barcodes[0].astype(str)

    ad.var_names_make_unique()
    ad.obs['GSM_ID'] = sample_id
    group_values = meta_df.loc[sample_id, 'Group'] if sample_id in meta_df.index else 'Unknown'
    ad.obs['group'] = [group_values] * ad.n_obs

    adatas.append(ad)

    print(f"\n Loaded {len(adatas)} samples")

#======================================================================
# STEP 09: Check sample and group metadata in one sample
#======================================================================

# preview metadata of first sample
adatas[0].obs.head()

#===========================================
# STEP 10: Summarize sample counts per group
#===========================================

all_obs = pd.concat([ad.obs for ad in adatas])
print(all_obs['group'].value_counts())

meta_df['Group'].value_counts()

# check group assignments for all loaded samples are correct
# Note: there should be no output, if all is correct.
for ad in adatas:
    sample = ad.obs['GSM_ID'].iloc[0]
    groups = ad.obs['group'].unique()
    if len(groups) >1:
        print(f"Sample {sample} has multiple group labels: {groups}")

# Recheck cell count distribution per group label
all_obs = pd.concat([ad.obs for ad in adatas])
print(all_obs['group'].value_counts())

# check each sample's assigned group and cell count
for ad in adatas:
    sample = ad.obs['GSM_ID'].iloc[0]
    group = ad.obs['group'].iloc[0]
    n_cells = ad.n_obs
    print(f"Sample: {sample} | Group: {group} | Cells: {n_cells}")

all_obs = pd.concat([ad.obs for ad in adatas])
print(all_obs['group'].value_counts())

#===========================================
# STEP 11 : QC-Preprocessing
#===========================================
# Step:A- Defining files for Voilin and scatter plots:

# Violin plot QC_plot_pct_counts_mt
files_pl_vp_pct_cmt = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    vp_pct_cmt = f"{sample_id}_QC_plot_pct_counts_mt.png"
    files_pl_vp_pct_cmt.append(vp_pct_cmt)

# Violin plot QC_plot_total_counts
files_pl_vp_tc = []
for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    vp_tc = f"{sample_id}_QC_plot_total_counts.png"
    files_pl_vp_tc.append(vp_tc)

# Violin plot QC_plot_n_genes_by_counts
files_pl_vp_ngc = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    vp_ngc = f"{sample_id}_QC_plot_n_genes_by_counts.png"
    files_pl_vp_ngc.append(vp_ngc)

# violin plot QC_plot_pct_counts_ercc
files_pl_vp_pct_cercc = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    vp_pct_cercc = f"{sample_id}_QC_plot_pct_counts_ercc.png"
    files_pl_vp_pct_cercc.append(vp_pct_cercc)

# Scatter plot QC_scatter_col_pct_counts_mt
files_pl_sc_col_pct_mt = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    sc_col_pct_mt = f"{sample_id}_QC_scatter_col_pct_counts_mt.png"
    files_pl_sc_col_pct_mt.append(sc_col_pct_mt)

# Scatter plot QC_scatter_col_pct_ercc
files_pl_sc_col_pct_ercc = []

for raw_path in sample_folders:
    sample_id = os.path.basename(raw_path)
    sc_col_pct_ercc = f"{sample_id}_QC_scatter_col_pct_counts_ercc.png"
    files_pl_sc_col_pct_ercc.append(sc_col_pct_ercc)

# Step:B- Filtering of mt,ribosomal,hb genes and low quality cells:
# First Normalization happens in this step.
import pandas as pd
import scipy.io
import scanpy as sc
from pathlib import Path

def pp_qc_fun(raw_paths, pid, vp1, vp2, vp3, vp4, sc1, sc2):
    raw_path = Path(raw_paths)

    # ======================
    # Step 1: Load 10x data
    # ======================
    matrix = scipy.io.mmread(raw_path / "matrix.mtx").tocsr().T
    genes = pd.read_csv(raw_path / "features.tsv", header=None, sep="\t")
    genes.columns = ['gene_id', 'gene_symbols', 'feature_type']
    barcodes = pd.read_csv(raw_path / "barcodes.tsv", header=None)[0].tolist()

    # ======================
    # Step 2: Create AnnData
    # ======================
    adata = sc.AnnData(X=matrix)
    adata.var_names = genes["gene_symbols"]
    adata.obs_names = barcodes
    adata.var_names_make_unique()
    adata.obs['sample_id'] = pid

    # ======================
    # Step 3: QC Annotations
    # ======================
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    adata.var["ercc"] = adata.var_names.str.upper().str.startswith("ERCC")
    adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]", regex=True)

    # =============================
    # Step 4: Calculate QC metrics
    # =============================
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ercc", "ribo", "hb"], inplace=True, log1p=True)

    # ======================
    # Step 5: Filter cells
    # ======================
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_genes=200)

    #==============================
    # Step 6: Remove unwanted genes
    #==============================
    unwanted = adata.var["mt"] | adata.var["ribo"] | adata.var["hb"] | adata.var["ercc"]
    adata = adata[:, ~unwanted].copy()

    # ======================
    # Step 7: QC plots
    # ======================
    sc.pl.violin(adata, ["total_counts"], jitter=0.4, save=f"_{vp1}")
    sc.pl.violin(adata, ["n_genes_by_counts"], jitter=0.4, save=f"_{vp2}")
    sc.pl.violin(adata, ["pct_counts_mt"], jitter=0.4, save=f"_{vp3}")
    sc.pl.violin(adata, ["pct_counts_ercc"], jitter=0.4, save=f"_{vp4}")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", save=f"_{sc1}")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_ercc", save=f"_{sc2}")
    return adata

# Create an empty list to store processed AnnData objects
adata_lst = []

for i in range(len(sample_folders)):
    raw_path = sample_folders[i]                      
    sample_id = os.path.basename(raw_path)            

    print(f"Processing sample: {sample_id}")
    print(f" Raw path: {raw_path}")
    print(f" VP1: {files_pl_vp_tc}")
    print(f" VP2: {files_pl_vp_ngc[i]}")
    print(f" VP3: {files_pl_vp_pct_cmt[i]}")
    print(f" VP4: {files_pl_vp_pct_cercc[i]}")
    print(f" SC1: {files_pl_sc_col_pct_mt[i]}")
    print(f" SC2: {files_pl_sc_col_pct_ercc[i]}")

    adata_ps = pp_qc_fun(
        raw_path, 
        sample_id, 
        files_pl_vp_tc[i], 
        files_pl_vp_ngc[i], 
        files_pl_vp_pct_cmt[i], 
        files_pl_vp_pct_cercc[i], 
        files_pl_sc_col_pct_mt[i], 
        files_pl_sc_col_pct_ercc[i]
    )     

    adata_lst.append(adata_ps)

adata = adata_lst[0]
adata.obs['sample_id'].unique(), adata.n_obs, adata.n_vars

for adata in adata_lst:
    print(adata.shape)  # (n_cells, n_genes)

for i, adata in enumerate(adata_lst):
    sample_id = adata.obs['sample_id'].iloc[0]  # Proper way to extract scalar
    print(f" Sample: {sample_id} | Cells: {adata.n_obs} | Genes: {adata.n_vars}")

output_dir = "GSE231993_preprocessed_1/h5ad_preprocessed"
os.makedirs(output_dir, exist_ok=True)

for adata in adata_lst:
    sample_id = adata.obs['sample_id'].iloc[0]
    adata.write(os.path.join(output_dir, f"{sample_id}_preprocessed.h5ad"))
    
with open("adata_2nd_preprocessed.pkl", "wb") as f:
    pickle.dump(adata_lst, f)
# Save the list of AnnData objects to a pickle file

### Step:C-  Filter outlier cells using the MAD(Median Absolute Deviation)
#Note: 2nd Normalization of preprocessing is in this step.

# ========== Function to identify outliers ==========
def is_outlier(adata, metric: str, nmads: int = 5):
    M = adata.obs[metric]
    med = np.median(M)
    mad = median_abs_deviation(M)
    outlier = (M < med - nmads * mad) | (M > med + nmads * mad)
    return outlier

# ========== Load adata_lst ==========
with open("adata_2nd_preprocessed.pkl", "rb") as f:
    adata_lst = pickle.load(f)

# ========== Load GSM ID mapping ==========
# Ensure this file has exactly 12 rows in order: GSM_ID, group
meta_df = pd.read_csv("sample_metadata.csv")  # adjust path if needed
gsm_ids = meta_df["GSM_ID"].tolist()  # Assumes it's in correct order

# ========== Output folders ==========
os.makedirs("filtered_2_h5ad", exist_ok=True)
os.makedirs("outlier_2_summary", exist_ok=True)
# ========== Initialize summary ==========
summary = []
# ========== Loop over each AnnData ==========
for i, adata in enumerate(adata_lst):
    gsm_id = gsm_ids[i]
    adata.obs["outlier"] = (
        is_outlier(adata, "log1p_total_counts", 5)
        | is_outlier(adata, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata, "pct_counts_mt", 3)
        | is_outlier(adata, "pct_counts_in_top_50_genes", 5)
    )

    total_cells = adata.n_obs
    n_outliers = adata.obs["outlier"].sum()
    n_retained = total_cells - n_outliers
    percent_outliers = round(n_outliers / total_cells * 100, 2)
 # Save filtered AnnData
    filtered = adata[~adata.obs["outlier"]].copy()
    filtered.write(f"filtered_2_h5ad/{gsm_id}_filtered_2_.h5ad")
# Append to summary
    summary.append([gsm_id, total_cells, n_outliers, n_retained, percent_outliers])

# ========== Save CSV ==========
summary_df = pd.DataFrame(
    summary,
    columns=["GSM_ID", "total_cells", "outliers_removed", "cells_retained", "percent_outliers"]
)
summary_df.to_csv("outlier_2_summary/mad_outlier_2_summary.csv", index=False)

print(" Outlier filtering complete with GSM-based filenames.")

# check gene count after MAD:
import os
import scanpy as sc
folder = "filtered_2_h5ad"

# Optional: Sort files by GSM ID order
files = sorted([f for f in os.listdir(folder) if f.endswith(".h5ad")])

for f in files:
    adata = sc.read(os.path.join(folder, f))
    sample_id = adata.obs["sample_id"].unique()[0]
    print(f"Sample: {sample_id} | Cells: {adata.n_obs} | Genes: {adata.n_vars}")

# Step:D- Doublet removal by Scrublet package using adaptive threshold as per each sample.
#Note: Adaptive thresholds are used to remove the bais of having one threshold of 0.2 or 0.25 which will not give correct outputs.

# Path to filtered .h5ad files
input_dir = "filtered_2_h5ad"
out_dir = "filtered_2_h5ad_scrublet"
metadata_csv = "sample_metadata.csv"  # Must have columns: GSM_ID, Sample_ID
os.makedirs(out_dir, exist_ok=True)

# ========= Load Metadata =========
meta = pd.read_csv(metadata_csv)
sample_to_gsm = dict(zip(meta["GSM_ID"], meta["Group"]))

# ========= Init containers =========
filtered_adatas = []
scrublet_summary = []
results = {}

# Get list of input .h5ad files
file_list = [f for f in os.listdir(input_dir) if f.endswith(".h5ad")]

for fname in file_list:
    fpath = os.path.join(input_dir, fname)
    ad = sc.read(fpath)
    sample_id = ad.obs['PID'][0] if 'PID' in ad.obs.columns else fname.replace(".h5ad", "")
    print(f"\n Running Scrublet for {sample_id}...")

    try:
        # Extract raw counts
        if ad.raw is not None:
            counts_matrix = ad.raw.X.toarray() if not isinstance(ad.raw.X, np.ndarray) else ad.raw.X
        else:
            counts_matrix = ad.X.toarray() if not isinstance(ad.X, np.ndarray) else ad.X

        # Run Scrublet
        scrub = scr.Scrublet(counts_matrix)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()

        # Annotate and filter
        ad.obs['doublet_score'] = doublet_scores
        ad.obs['predicted_doublet'] = predicted_doublets
        ad.uns['scrublet_threshold'] = scrub.threshold_ 

        ad_filtered = ad[~ad.obs['predicted_doublet'], :].copy()
        filtered_adatas.append(ad_filtered)

        # Save filtered AnnData
        out_file = os.path.join(out_dir, fname.replace(".h5ad", "_scrublet_filtered_2.h5ad"))
        ad_filtered.write(out_file)

        # Summary
        total_cells = ad.n_obs
        doublets_detected = np.sum(predicted_doublets)
        scrublet_summary.append({
            "sample_id": sample_id,
            "total_cells": total_cells,
            "doublets_detected": int(doublets_detected),
            "cells_retained": int(total_cells - doublets_detected),
            "doublet_rate (%)": round(100 * doublets_detected / total_cells, 2)
        })

        results[sample_id] = {
            "scores": doublet_scores,
            "threshold": scrub.threshold_
        }

    except Exception as e:
        print(f" Scrublet error in {sample_id}: {e}")
        scrublet_summary.append({
            "sample_id": sample_id,
            "total_cells": ad.n_obs,
            "doublets_detected": "ERROR",
            "cells_retained": "ERROR",
            "doublet_rate (%)": "ERROR"
        })

# Save summary table
summary_df = pd.DataFrame(scrublet_summary)
summary_df.to_csv("scrublet_summary_per_sample_2.csv", index=False)
print("\n Summary saved to 'scrublet_summary_per_sample_2.csv'")

# Path to Scrublet-filtered h5ad files
scrublet_dir = "filtered_2_h5ad_scrublet"  

# List and sort all h5ad files
files = sorted([f for f in os.listdir(scrublet_dir) if f.endswith(".h5ad")])

# Print shape (cells, genes) per sample
for f in files:
    adata = sc.read(os.path.join(scrublet_dir, f))
    sample_id = adata.obs['sample_id'].unique()[0] if 'sample_id' in adata.obs else f
    print(f"Sample: {sample_id} | Cells: {adata.n_obs} | Genes: {adata.n_vars}")

# Step:E- Histograms for the final doublet removed -preprocessed adatas of all samples
# ============================
#  Set Robust Paths
# ============================

cwd = os.getcwd()  # Already in: FW_GSE231993_raws

input_dir = "filtered_2_h5ad_scrublet"
summary_csv = "scrublet_summary_per_sample_2.csv"  # Now confirmed in current folder
output_dir = "scrublet_histograms_2"

os.makedirs(output_dir, exist_ok=True)

# ============================
#  Load Scrublet Summary
# ============================
summary_df = pd.read_csv(summary_csv)

# ============================
#  Plot per-sample histograms
# ============================
all_scores = []

for fname in os.listdir(input_dir):
    if not fname.endswith(".h5ad"):
        continue

    sample_id = fname.split("_filtered_2_scrublet_filtered_2.h5ad")[0]
    fpath = os.path.join(input_dir, fname)

    try:
        ad = sc.read(fpath)
        scores = ad.obs["doublet_score"]
        threshold = ad.uns.get("scrublet_threshold", None)
        all_scores.extend(scores.tolist())
        # Plot
        plt.figure(figsize=(6, 4))
        sns.histplot(scores, bins=50, kde=False, color="skyblue", edgecolor="black")
        if threshold is not None:
            plt.axvline(threshold, color="red", linestyle="--", label=f"Threshold: {threshold:.2f}")
            plt.legend()
        plt.title(f"Scrublet Doublet Scores: {sample_id}")
        plt.xlabel("Doublet Score")
        plt.ylabel("Cell Count")
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f"{sample_id}_doublet_score_hist.png"), dpi=150)
        plt.close()
    except Exception as e:
        print(f" Error in {sample_id}: {e}")

# ============================
#  Global Histogram
# ============================
plt.figure(figsize=(7, 5))
sns.histplot(all_scores, bins=100, kde=False, color="darkgreen", edgecolor="black")
plt.title("Global Scrublet Doublet Score Distribution")
plt.xlabel("Doublet Score")
plt.ylabel("Cell Count")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "global_doublet_score_hist.png"), dpi=150)
plt.close()

print(" All histograms saved to:", output_dir)

# Step 1: Define scrublet-filtered directory
scrublet_dir = "filtered_2_h5ad_scrublet"

# Step 2: List files
filtered_files = sorted([f for f in os.listdir(scrublet_dir) if f.endswith(".h5ad")])
print(f"Total samples found: {len(filtered_files)}")

# Step 3: Load files into a list
filtered_adatas = []

for f in filtered_files:
    ad = sc.read(os.path.join(scrublet_dir, f))
    filtered_adatas.append(ad)
    print(f" Loaded: {f} | Cells: {ad.n_obs} | Genes: {ad.n_vars}")

# ====================================
#  Step12- Concatenate filtered adatas
# ====================================
# concatenate all filtered adatas into one AnnData object

adata_all = sc.concat(filtered_adatas, label="sample_id", index_unique="-")

print("Final concatenated shape:", adata_all.shape)

os.makedirs("qc_outputs", exist_ok=True)
adata_all.write("qc_outputs/adata_all_postqc_preHVG.h5ad")

#===============================================
# STEP 12: HVG Filtering
#===============================================

print(" Loading data...")
adata = sc.read("qc_outputs/adata_all_postqc_preHVG.h5ad")

# === Convert matrix to float32 to avoid overflow during expm1 ===
print(" Converting matrix to float32...")
if sparse.issparse(adata.X):
    adata.X = adata.X.astype(np.float32)
    print("Sparse matrix dtype:", adata.X.dtype)
else:
    adata.X = adata.X.astype(np.float32)
    print("Dense matrix dtype:", adata.X.dtype)

# === Apply log1p transformation explicitly to stabilize values ===
print("Applying log1p transformation manually to avoid overflow...")
sc.pp.log1p(adata)

# === Check for Inf/NaN post-log1p ===
if sparse.issparse(adata.X):
    X_data = adata.X.data
    print("Inf (sparse):", np.isinf(X_data).sum())
    print("NaN (sparse):", np.isnan(X_data).sum())
else:
    print("Inf (dense):", np.isinf(adata.X).sum())
    print("NaN (dense):", np.isnan(adata.X).sum())

# === HVG Filtering ===
print("▶ Running HVG filtering...")
sc.pp.highly_variable_genes(
    adata, flavor="seurat", n_top_genes=2000, batch_key="sample_id"
)

# === Subset to HVGs ===
adata = adata[:, adata.var["highly_variable"]]
print(" Number of HVGs selected:", adata.shape[1])

# === Save output ===
print(" Saving filtered .h5ad...")
adata.write("qc_outputs/adata_all_HVGfiltered.h5ad")
print(" HVG filtering complete and saved.")

#  generate volin plots of ncounts, ngenes, pct_counts_mt, pct_counts_ercc
adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")
sc.pl.violin(
    adata, 
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ercc"], 
    jitter=0.4, 
    multi_panel=True,
    save="_hvg_filtered_violin.png"
)

# Load HVG-filtered AnnData
adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")

# Print shape: (n_cells, n_genes)
print(" HVG-filtered AnnData shape:", adata.shape)

# Show top 10 gene names
print("\nTop 10 HVG gene names:")
print(adata.var_names[:10])

#===============================================
# STEP 13: Scaling HVG Filtered Data
#===============================================
# Scale and regress out confounders

# === Load HVG-filtered AnnData ===
adata = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")

# === Scale the data (zero-mean, unit-variance) ===
# Max value is clipped to 10 to reduce impact of extreme outliers
sc.pp.scale(adata, zero_center=True, max_value=10)

# === Save the scaled object ===
adata.write("qc_outputs/adata_all_HVGscaled.h5ad")

# === Confirm ===
print(" Scaling complete. Scaled data saved to 'adata_all_HVGscaled.h5ad'")
print("Scaled AnnData shape:", adata.shape)

#========================================================================
# STEP 14: Merging Metadata -HVG Filtered Data for visualization by plots
#========================================================================
# Load HVG-filtered data
adata_all = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")

# Load metadata (must match order of loading samples!)
meta_df = pd.read_csv("sample_metadata.csv")  # Ensure this has GSM_ID and Group

# Extract sample index from obs_names like 'AAAC...-1-0'
adata_all.obs["sample_index"] = adata_all.obs_names.str.split("-").str[-1].astype(int)

# Map sample index to GSM_ID
gsm_ids = meta_df["GSM_ID"].tolist()
adata_all.obs["GSM_ID"] = adata_all.obs["sample_index"].map(lambda i: gsm_ids[i])

# Merge metadata
adata_all.obs = adata_all.obs.merge(meta_df, on="GSM_ID", how="left")

# Make group categorical
adata_all.obs["Group"] = adata_all.obs["Group"].astype("category")

# Save with metadata attached
adata_all.write("qc_outputs/adata_all_HVGfiltered_meta.h5ad")

#===========================================
# STEP 15: PCA, UMAP, and Leiden clustering
#===========================================

# ========= 1. Load HVG + metadata filtered AnnData =========
adata_all = sc.read("qc_outputs/adata_all_HVGfiltered_meta.h5ad")

# ========= 2. PCA =========
print("Running PCA...")
sc.tl.pca(adata_all, svd_solver="arpack")
sc.pl.pca(adata_all, color="Group", save="_pca_by_Group.png")

# ========= 3. Compute neighbors =========
print("Computing neighbors...")
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=40)

# ========= 4. Run UMAP =========
print("Computing UMAP...")
sc.tl.umap(adata_all)

# ========= 5. Leiden clustering =========
print("Running Leiden clustering...")
sc.tl.leiden(adata_all, resolution=0.5)

# ========= 6. Save UMAP plots =========
print("Saving UMAP plots...")
sc.pl.umap(adata_all, color="leiden", save="_leiden.png", title="Leiden Clusters")
sc.pl.umap(adata_all, color="Group", save="_by_Group.png", title="UMAP by Group")
sc.pl.umap(adata_all, color="GSM_ID", save="_by_GSMID.png", title="UMAP by GSM ID")

# ========= 7. Save updated AnnData =========
adata_all.write("qc_outputs/adata_all_clustered_umap.h5ad")
print(" Analysis complete. Output saved.")

#generate expression matrix as csv file
# Load the clustered AnnData
adata_all = sc.read("qc_outputs/adata_all_clustered_umap.h5ad")
# Convert to DataFrame
expression_df = pd.DataFrame(
    adata_all.X.toarray(),  # Convert sparse matrix to dense
    index=adata_all.obs_names,
    columns=adata_all.var_names
)
# Save to CSV
expression_df.to_csv("qc_outputs/adata_all_hvg_filtered_expression_matrix.csv")
print(" Expression matrix saved as 'adata_all_hvg_filtered_expression_matrix.csv'.")

# Load the HVG-filtered AnnData object
adata_all = sc.read("qc_outputs/adata_all_HVGfiltered.h5ad")

# Now check the observation names
print("Obs names example:", adata_all.obs_names[:5])

#==================================================
# STEP 16: Save pre-integration log-normalized Data
#==================================================

# Load clustered + normalized AnnData
adata_all = sc.read("qc_outputs/adata_all_clustered_umap.h5ad")

# Save log-normalized values to .raw before integration
adata_all.raw = adata_all

# Save updated object with .raw
adata_all.write("qc_outputs/adata_all_raw_set.h5ad")

print(" .raw data stored and saved in: qc_outputs/adata_all_raw_set.h5ad")

#==================================================
# STEP 17: Integration by Scanorama
#==================================================

# ========== Step 1: Load preprocessed AnnData ==========
adata_all = sc.read("qc_outputs/adata_all_raw_set.h5ad")
# ========== Step 2: Split AnnData by sample_id ==========
adatas_split = [
    adata_all[adata_all.obs["sample_id"] == sid].copy()
    for sid in adata_all.obs["sample_id"].unique()
]
# ========== Step 3: Run Scanorama integration ==========
corrected_adatas = scanorama.correct_scanpy(adatas_split, return_dimred=True)
# ========== Step 4: Concatenate back into a single AnnData ==========
# Keep track of sample names
sample_ids = [a.obs["sample_id"][0] for a in corrected_adatas]
adata_scanorama = sc.concat(corrected_adatas, label="sample_id", keys=sample_ids)
adata_scanorama.obs_names_make_unique(join="_")
# ========== Step 5: Save integrated expression matrix ==========
# Convert expression matrix to dense format if needed
df_expr = pd.DataFrame(
    adata_scanorama.X.toarray() if not isinstance(adata_scanorama.X, np.ndarray) else adata_scanorama.X,
    index=adata_scanorama.obs_names,
    columns=adata_scanorama.var_names
)
df_expr.to_csv("scanorama_integrated_expression_matrix.csv")

# Save AnnData object itself
adata_scanorama.write("qc_outputs/adata_scanorama_integrated.h5ad")

print("\n Scanorama integration complete!")
print("Saved integrated expression matrix: scanorama_integrated_expression_matrix.csv")
print("Saved integrated AnnData object: qc_outputs/adata_scanorama_integrated.h5ad")

# saving expression matrix with header:
import scanpy as sc
import pandas as pd
import numpy as np

# Load the integrated Scanorama AnnData object
adata = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")

# Convert the expression matrix to a dense DataFrame with gene and cell names
df_expr = pd.DataFrame(
    adata.X.toarray() if not isinstance(adata.X, np.ndarray) else adata.X,
    index=adata.obs_names,
    columns=adata.var_names
)

# Save to CSV with row index (cells) and column headers (genes)
df_expr.to_csv("scanorama_integrated_expression_matrix.csv", index=True, header=True)

print("Saved: scanorama_integrated_expression_matrix.csv with row and column labels")

# Rows are cells and columns are genes.
# Scanorama integrated expression values are -transformed gene expression per cell.
print(df_expr.head())

adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")
print(adata_scanorama.obs.columns.tolist())

#======================================================================
# STEP 17: PCA, UMAP and Leiden clustering on Scanorama integrated data
#======================================================================

# Load integrated data
adata_scanorama = sc.read("qc_outputs/adata_scanorama_integrated.h5ad")

# ========== Step 1: Assign Scanorama embeddings as PCA ==========
adata_scanorama.obsm["X_pca"] = adata_scanorama.obsm["X_scanorama"]

# PCA plot
sc.pl.pca(
    adata_scanorama,
    color=["Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    save="_scanorama_Group_PCA.png"
)

print("PCA plot saved to figures/_scanorama_Group_PCA.png")

# ========== Step 2: Neighbors, UMAP ==========
sc.pp.neighbors(adata_scanorama, use_rep="X_pca", n_neighbors=15, n_pcs=40)
sc.tl.umap(adata_scanorama)

# ========== Step 3: Leiden Clustering ==========
sc.tl.leiden(adata_scanorama, resolution=0.5, flavor="igraph", directed=False)

# ========== Step 4: UMAP Plotting ==========
sc.set_figure_params(dpi=150)

# Grouping overview
sc.pl.umap(
    adata_scanorama,
    color=["Group", "sample_id", "leiden", "GSM_ID"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    save="_scanorama_integrated_UMAP.png"
)

# Cluster overview
sc.pl.umap(
    adata_scanorama,
    color=["leiden", "Group"],
    legend_loc="right margin",
    legend_fontsize=10,
    ncols=2,
    wspace=0.4,
    save="_scanorama_leiden_Group.png"
)

# ========== Step 5: Save clustered object ==========
adata_scanorama.write("qc_outputs/adata_scanorama_clustered_umap.h5ad")

print("UMAP & clustering complete. Plots saved to figures/. AnnData saved to qc_outputs/adata_scanorama_clustered_umap.h5ad") 

#===============================================
# STEP 18: Detection and Ranking of Marker genes
#===============================================

# A: Rank genes by Leiden: Ran as sbatch using run_marker_leiden_detection.py: This is for genes distinguishing clusters.

# B:Rank genes by Group(Condition-UC,UCA,HC): Ran as sbatch using marker_genes_by_group.py. :This gives genes that define condition-specific responses across all cells or per cell type/cluster. Also DEGs between groups.
#==========================================================================================================================================================================================================================

# A: Rank genes by Leiden: Ran as sbatch using run_marker_leiden_detection.py: This is for genes distinguishing clusters.
# Load clustered data
adata = sc.read("qc_outputs/adata_scanorama_clustered_umap.h5ad")

# Identify marker genes for each Leiden cluster
sc.tl.rank_genes_groups(
    adata,
    groupby="leiden",
    method="wilcoxon",
    key_added="rank_genes_leiden"
)

# Create results directory
os.makedirs("results", exist_ok=True)

# Save top marker table
markers_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_leiden")
markers_df.to_csv("results/scanorama_leiden_markers.csv", index=False)

# Plot top 20 genes per cluster
sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False,
    key="rank_genes_leiden",
    save="_scanorama_leiden_top20.png"
)

print(" Marker gene detection complete: results/scanorama_leiden_markers.csv + top20 plot saved")

# B:Rank genes by Group(Condition-UC,UCA,HC): Ran as sbatch using marker_genes_by_group.py. :This gives genes that define condition-specific responses across all cells or per cell type/cluster. Also DEGs between groups.
# Load clustered data
adata = sc.read("qc_outputs/adata_scanorama_clustered_umap.h5ad")
# Identify DEGs by condition group (e.g., UC vs CD vs HC)
sc.tl.rank_genes_groups(
    adata,
    groupby="Group",  # must match your adata.obs['Group']
    method="wilcoxon",
    key_added="rank_genes_Group"
)
# Create results folder
os.makedirs("result", exist_ok=True)

# Save all group comparisons
markers_group_df = sc.get.rank_genes_groups_df(adata, group=None, key="rank_genes_Group")
markers_group_df.to_csv("results/condition_DEGs_by_Group.csv", index=False)

# Plot top 20 per group
sc.pl.rank_genes_groups(
    adata,
    n_genes=20,
    sharey=False,
    key="rank_genes_Group",
    save="_top20_Group.png"
)

print(" Marker gene detection by Group complete.")

#===========================================
# STEP 18: Feature Plots for specific genes
#===========================================
adata_path = "qc_outputs/adata_scanorama_clustered_umap.h5ad"
output_dir = "umap_feature_plots"
os.makedirs(output_dir, exist_ok=True)
# Load AnnData
adata = sc.read(adata_path)
# Set output directory and figure parameters
sc.settings.figdir = output_dir
sc.set_figure_params(dpi=150, format="png")

# Feature genes
genes_to_plot = ["LYZ", "MKI67"]

# Plot each gene
for gene in genes_to_plot:
    sc.pl.umap(adata, color=gene, show=False, save=f"_feature_{gene}.png")

print("Feature plots saved in:", output_dir)

# # end===============================================================================================================================
