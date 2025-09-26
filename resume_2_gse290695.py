#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#==============================
# GSE290695_resume_1
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

# Detection and Ranking of Marker genes:
# A: Rank genes by Leiden: Ran as sbatch using run_marker_leiden_detection.py: This is for genes distinguishing clusters.

# B:Rank genes by Group(Condition-UC,UCA,HC): Ran as sbatch using marker_genes_by_group.py. :This gives genes that define condition-specific responses across all cells or per cell type/cluster. Also DEGs between groups.


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


# In[ ]:


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


# In[ ]:


# feature plots
adata_path = "qc_outputs/adata_scanorama_clustered_umap.h5ad"
output_dir = "umap_feature_plots"
os.makedirs(output_dir, exist_ok=True)
# Load AnnData
adata = sc.read(adata_path)
# Set output directory and figure parameters
sc.settings.figdir = output_dir
sc.set_figure_params(dpi=150, format="png")

# Feature genes
genes_to_plot = ["CLCN2", "IL5RA"]

# Plot each gene
for gene in genes_to_plot:
    sc.pl.umap(adata, color=gene, show=False, save=f"_feature_{gene}.png")

print("Feature plots saved in:", output_dir)