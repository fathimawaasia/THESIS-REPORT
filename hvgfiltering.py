#!/usr/bin/env python
# coding: utf-8

"""
GSE125527 Post-QC Integration and Clustering
- Input: scrublet-filtered h5ad files
- Output: HVG-filtered, scaled, clustered AnnData with metadata attached
"""

import os
import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse

# ========== Paths ==========
qc_dir = "qc_outputs"
filtered_dir = f"{qc_dir}/scrublet_filtered_h5ad"
meta_file = "sample_metadata_7.csv"
os.makedirs(qc_dir, exist_ok=True)

# ========== Step 1: Load all scrublet-filtered AnnData files ==========
print("▶ Loading scrublet-filtered h5ad files...")
filtered_files = sorted([f for f in os.listdir(filtered_dir) if f.endswith(".h5ad")])
adata_list = [sc.read(os.path.join(filtered_dir, f)) for f in filtered_files]

# ========== Step 2: Concatenate ==========
adata_all = sc.concat(adata_list, label="sample_id", index_unique="-")
print(f"✅ Concatenated shape: {adata_all.shape}")

# Save pre-HVG filtered file
adata_all.write(f"{qc_dir}/adata_all_postqc_preHVG.h5ad")

# ========== Step 3: HVG Filtering ==========
sc.pp.highly_variable_genes(
    adata_all,
    flavor="seurat",
    n_top_genes=2000,
    batch_key="sample_id"
)
adata_all = adata_all[:, adata_all.var.highly_variable]
print(f"✅ HVG filtering done. Shape: {adata_all.shape}")

# Save HVG-filtered data
adata_all.write(f"{qc_dir}/adata_all_HVGfiltered.h5ad")

# ========== Step 4: Scaling ==========
adata_all = adata_all.copy()  # avoid view warning
sc.pp.scale(adata_all, zero_center=True, max_value=10)
adata_all.write(f"{qc_dir}/adata_all_HVGscaled.h5ad")
print("✅ Scaling complete.")

# ========== Step 5: Attach sample metadata ==========
meta_df = pd.read_csv(meta_file)
meta_df.columns = [col.strip() for col in meta_df.columns]  # remove any trailing spaces

# Try to identify 'Group'-like column
possible_group_cols = [col for col in meta_df.columns if col.lower().strip() in ["group", "condition"]]
if possible_group_cols:
    print(f"🧠 Found group column: {possible_group_cols[0]}")
    meta_df.rename(columns={possible_group_cols[0]: "Group"}, inplace=True)
else:
    raise ValueError("❌ No suitable 'Group' column found in metadata. Columns available: " + ", ".join(meta_df.columns))

# Ensure GSM_ID is in obs
if "GSM_ID" not in adata_all.obs.columns:
    print("🛠 Adding GSM_ID from sample_index...")
    adata_all.obs["sample_index"] = adata_all.obs_names.str.split("-").str[-1].astype(int)
    gsm_ids = meta_df["GSM_ID"].tolist()
    adata_all.obs["GSM_ID"] = adata_all.obs["sample_index"].map(lambda i: gsm_ids[i])

# Show sample before merge
print("🧪 Before merge:")
print(adata_all.obs[["GSM_ID"]].head())

# Merge metadata
adata_all.obs = adata_all.obs.merge(meta_df, on="GSM_ID", how="left")

# Check and set Group as category
has_group = "Group" in adata_all.obs.columns and adata_all.obs["Group"].notna().any()
if has_group:
    adata_all.obs["Group"] = adata_all.obs["Group"].astype("category")
    print("✅ Metadata attached successfully.")
    # Log failed merges
    unmerged = adata_all.obs[adata_all.obs["Group"].isna()].copy()
    if not unmerged.empty:
        unmerged.to_csv("qc_outputs/unmatched_GSM_rows.csv")
        print(f"⚠️ Warning: {len(unmerged)} cells could not match Group info. See unmatched_GSM_rows.csv")
else:
    print("⚠️ 'Group' column not found or all values are NaN after merge. Skipping category conversion.")

# Save HVG filtered + metadata attached
adata_all.write(f"{qc_dir}/adata_all_HVGfiltered_meta.h5ad")

# ========== Step 6: PCA, neighbors, UMAP, Leiden ==========
print("▶ Running PCA...")
sc.tl.pca(adata_all, svd_solver="arpack")
if has_group:
    sc.pl.pca(adata_all, color="Group", save="_GSE125527_PCA_Group.png")
else:
    sc.pl.pca(adata_all, save="_GSE125527_PCA.png")

print("▶ Computing neighbors...")
sc.pp.neighbors(adata_all, n_neighbors=15, n_pcs=40)

print("▶ Running UMAP...")
sc.tl.umap(adata_all)
sc.pl.umap(adata_all, color=["sample_id"] + (["Group"] if has_group else []), save="_GSE125527_UMAP.png")

print("▶ Running Leiden clustering...")
sc.tl.leiden(adata_all, resolution=0.5)
sc.pl.umap(adata_all, color=["leiden"], save="_GSE125527_UMAP_leiden.png")

# ========== Step 7: Violin plots ==========
sc.pl.violin(
    adata_all,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
    save="_GSE125527_violin_postqc.png"
)

# ========== Step 8: Save clustered file ==========
adata_all.write(f"{qc_dir}/adata_all_clustered_umap.h5ad")
print("✅ Clustering complete. All files saved to qc_outputs/")
