#!/usr/bin/env python3
"""
Prepare pooled AnnData with joint HVG selection for scVI (counts preserved).
Joins ALL metadata by GSM_ID.

Example:
  python prepare_joint_hvg_for_scvi.py \
    --inputs  adata_all_postqc_preHVG_GSE125527_dataset_1.h5ad \
              adata_all_postqc_preHVG_GSE148837_dataset_2.h5ad \
              adata_all_postqc_preHVG_GSE182270_dataset_3.h5ad \
              adata_all_postqc_preHVG_GSE231993_dataset_4.h5ad \
              adata_all_postqc_preHVG_GSE290695_dataset_5.h5ad \
    --metas   sample_metadata_GSE125527_dataset_1.csv \
              sample_metadata_GSE148837_dataset_2.csv \
              sample_metadata_GSE182270_dataset_3.csv \
              sample_metadata_GSE231993_dataset_4.csv \
              sample_metadata_GSE290695_dataset_5.csv \
    --names   GSE125527 GSE148837 GSE182270 GSE231993 GSE290695 \
    --n_hvg   2000 \
    --out     ../ALL5/qc_outputs/adata_all_HVGfiltered_for_scvi.h5ad
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
from scipy.sparse import issparse

# ---------- helpers ----------
def is_int_like(mat) -> bool:
    if mat is None:
        return False
    if issparse(mat):
        d = mat.data
        return d.size == 0 or (d.min() >= 0 and np.allclose(d, np.round(d), atol=1e-8))
    arr = np.asarray(mat)
    return arr.size == 0 or (arr.min() >= 0 and np.allclose(arr, np.round(arr), atol=1e-8))

def ensure_counts_layer(ad, ds_name: str):
    """Ensure ad.layers['counts'] contains raw/integer-like counts."""
    if "counts" in ad.layers and is_int_like(ad.layers["counts"]):
        return
    if ad.raw is not None and ad.raw.X is not None and is_int_like(ad.raw.X):
        ad.layers["counts"] = ad.raw.X
        return
    if is_int_like(ad.X):
        ad.layers["counts"] = ad.X
        return
    raise ValueError(
        f"{ds_name}: no integer-like counts in layers['counts'], .raw.X or .X."
    )
def attach_metadata_gsm(ad, meta_path: Path, ds_name: str):
    """Join metadata strictly by GSM_ID; add only columns that don't already exist."""
    if not meta_path.exists():
        print(f"[WARN] Missing metadata CSV: {meta_path}")
        return
    mdf = pd.read_csv(meta_path)
    mdf.columns = [c.strip() for c in mdf.columns]

    # Ensure AnnData.obs has GSM_ID (rename sample_id -> GSM_ID if needed)
    if "GSM_ID" not in ad.obs.columns and "sample_id" in ad.obs.columns:
        ad.obs = ad.obs.rename(columns={"sample_id": "GSM_ID"})

    if "GSM_ID" not in ad.obs.columns or "GSM_ID" not in mdf.columns:
        print(f"[WARN] {ds_name}: GSM_ID not present in both obs and CSV; skipping metadata.")
        return

    # Deduplicate metadata by GSM_ID (keep first)
    mdf = mdf.drop_duplicates(subset=["GSM_ID"], keep="first").copy()

    # String keys for stable matching
    ad.obs["GSM_ID"] = ad.obs["GSM_ID"].astype(str)
    mdf["GSM_ID"] = mdf["GSM_ID"].astype(str)
    mdf = mdf.set_index("GSM_ID")

    # Add only NEW columns (avoid overlaps like 'Group')
    cols_to_add = [c for c in mdf.columns if c not in ad.obs.columns]
    added = 0
    for c in cols_to_add:
        ad.obs[c] = ad.obs["GSM_ID"].map(mdf[c])
        added += 1

    matched = ad.obs["GSM_ID"].isin(mdf.index).mean()
    print(f"[INFO] {ds_name}: GSM_ID mapping — matched {matched*100:.1f}% cells, added {added} new col(s): {cols_to_add}")

# ---------- main ----------
def parse_args():
    ap = argparse.ArgumentParser("Prepare pooled HVGs for scVI (counts preserved, GSM_ID join)")
    ap.add_argument("--inputs", nargs=5, required=True, help="5 postQC, pre-HVG .h5ad files")
    ap.add_argument("--metas",  nargs=5, required=True, help="5 metadata CSVs (each with GSM_ID)")
    ap.add_argument("--names",  nargs=5, required=True, help="Short dataset names (e.g., GSE...)")
    ap.add_argument("--n_hvg",  type=int, default=2000, help="Number of HVGs")
    ap.add_argument("--out",    required=True, help="Output .h5ad path")
    return ap.parse_args()

def main():
    args = parse_args()

    adatas = []
    for p_h5, p_meta, ds in zip(args.inputs, args.metas, args.names):
        p_h5 = Path(p_h5); p_meta = Path(p_meta)
        if not p_h5.exists():
            raise FileNotFoundError(f"Missing input: {p_h5}")
        print(f"[LOAD] {ds}: {p_h5}")
        ad = sc.read(p_h5)

        # 1) counts layer
        ensure_counts_layer(ad, ds)

        # 2) dataset id
        ad.obs["dataset_id"] = ds

        # 3) metadata by GSM_ID (handles rename sample_id->GSM_ID)
        attach_metadata_gsm(ad, p_meta, ds)

        # Keep .X = counts for downstream tools expecting counts in .X
        ad.X = ad.layers["counts"]
        adatas.append(ad)

    # 4) concatenate and make obs names unique
    ad_all = an.concat(adatas, join="inner", label="internal_batch", keys=args.names, index_unique=None)
    ad_all.obs_names_make_unique()

    if "counts" not in ad_all.layers:
        ad_all.layers["counts"] = ad_all.X

    # 5) HVG on counts (Seurat v3), batched by dataset
    print(f"[HVG] Selecting {args.n_hvg} genes (flavor='seurat_v3') on layer='counts', batch_key='dataset_id'.")
    sc.pp.highly_variable_genes(
        ad_all,
        n_top_genes=args.n_hvg,
        flavor="seurat_v3",
        layer="counts",
        batch_key="dataset_id",
        inplace=True,
    )

    mask = ad_all.var["highly_variable"].values
    kept = int(mask.sum())
    print(f"[HVG] kept {kept} genes.")

    ad_hvg = ad_all[:, mask].copy()
    outp = Path(args.out)
    outp.parent.mkdir(parents=True, exist_ok=True)
    ad_hvg.write(outp, compression="gzip")
    print(f"[DONE] Wrote pooled HVG file for scVI: {outp}")

if __name__ == "__main__":
    main()


