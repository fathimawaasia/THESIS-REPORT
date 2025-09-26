#!/usr/bin/env python
# coding: utf-8

"""
Marker Gene Detection for GSE125527
- Input: Scanorama clustered AnnData (h5ad)
- Output: Marker genes ranked by Leiden clusters and Group, plus feature plots
"""

import os
import scanpy as sc

# ========== Load clustered Scanorama data ==========
adata = sc.read("qc_outputs/adata_scanorama_clustered_umap_1.h5ad")

# ========== Step C: Feature plots for selected genes ==========
print("▶ Generating UMAP feature plots for selected genes...")
output_dir = "umap_feature_plots"
os.makedirs(output_dir, exist_ok=True)
sc.settings.figdir = output_dir
sc.set_figure_params(dpi=150, format="png")

genes_to_plot = ["CD69", "IGHG3"]
for gene in genes_to_plot:
    sc.pl.umap(adata, color=gene, show=False, save=f"_feature_{gene}.png")

print("✅ Feature plots saved to umap_feature_plots/")