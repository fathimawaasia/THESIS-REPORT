#!/usr/bin/env python3
"""
Preprocess GSE125527 from dense cell×gene TSVs while PRESERVING raw counts.

Input  :
  - A folder of *.tsv.gz (rows=cells, cols=genes, values=UMIs)
  - A sample metadata CSV with at least a GSM_ID column (and optional Group)
Output :
  - Per-sample postQC files   : <outdir>/<sample>_postqc_counts.h5ad
  - Pooled postQC (pre-HVG)   : <outdir>/adata_all_postqc_preHVG_GSE125527_counts.h5ad
Notes  :
  - Raw counts are stashed in layers["counts"] and kept in .X as well
  - QC metrics are computed from counts (no log/normalize in place)
  - Scrublet runs on counts
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
from scipy.sparse import issparse, csr_matrix
from scipy.stats import median_abs_deviation
import scrublet as scr
import matplotlib.pyplot as plt


def parse_args():
    ap = argparse.ArgumentParser(description="Preprocess GSE125527 with counts preserved.")
    ap.add_argument("--input_dir", required=True, help="Folder of *.tsv or *.tsv.gz (cell×gene UMI tables).")
    ap.add_argument("--metadata_csv", required=True, help="CSV with at least GSM_ID column; optional Group column.")
    ap.add_argument("--outdir", required=True, help="Output directory.")
    ap.add_argument("--min_cells", type=int, default=3, help="Filter genes expressed in fewer than this many cells.")
    ap.add_argument("--min_genes", type=int, default=200, help="Filter cells with fewer detected genes.")
    ap.add_argument("--nmads_counts", type=float, default=5.0, help="MAD threshold for total_counts and n_genes.")
    ap.add_argument("--nmads_mt", type=float, default=3.0, help="MAD threshold for pct_counts_mt.")
    ap.add_argument("--nmads_top50", type=float, default=5.0, help="MAD threshold for pct_counts_in_top_50_genes.")
    ap.add_argument("--write_per_sample", action="store_true", help="Write per-sample postQC .h5ad files.")
    return ap.parse_args()


def is_outlier(x: np.ndarray, nmads: float) -> np.ndarray:
    med = np.median(x)
    mad = median_abs_deviation(x, scale='normal')
    if mad == 0:
        return np.zeros_like(x, dtype=bool)
    return (x < med - nmads * mad) | (x > med + nmads * mad)


def ensure_counts_in_layers(ad):
    """Put raw counts into layers['counts'] and optionally mirror to .X."""
    if "counts" not in ad.layers:
        if issparse(ad.X):
            ad.layers["counts"] = ad.X.copy()
        else:
            ad.layers["counts"] = csr_matrix(ad.X)
    # For clarity, keep .X == counts
    ad.X = ad.layers["counts"]


def load_sample_tables(input_dir: Path, meta_csv: Path):
    # metadata
    meta_df = pd.read_csv(meta_csv)
    meta_df.columns = [c.strip() for c in meta_df.columns]
    if "GSM_ID" not in meta_df.columns:
        raise ValueError("metadata_csv must contain a 'GSM_ID' column")
    group_map = dict(zip(meta_df["GSM_ID"].astype(str), meta_df.get("Group", pd.Series(["Unknown"]*len(meta_df))).astype(str)))

    # data tables
    tsvs = sorted(list(input_dir.glob("*.tsv")) + list(input_dir.glob("*.tsv.gz")))
    if not tsvs:
        raise FileNotFoundError(f"No .tsv or .tsv.gz files found in {input_dir}")

    adatas = []
    for f in tsvs:
        df = pd.read_csv(f, sep="\t", index_col=0)
        # df: rows=cells, cols=genes, values=counts (integers)
        ad = sc.AnnData(df.values)
        ad.obs_names = df.index.astype(str)
        ad.var_names = df.columns.astype(str)
        ad.var_names_make_unique()

        # sample id from filename prefix before first underscore
        sid = f.name.split("_")[0]
        ad.obs["sample_id"] = sid
        ad.obs["GSM_ID"] = sid
        ad.obs["Group"] = group_map.get(sid, "Unknown")

        ensure_counts_in_layers(ad)
        adatas.append(ad)
    return adatas


def qc_and_filter(ad: sc.AnnData,
                  min_genes: int,
                  min_cells: int,
                  nmads_counts: float,
                  nmads_mt: float,
                  nmads_top50: float) -> sc.AnnData:
    # gene flags
    ad.var["mt"] = ad.var_names.str.upper().str.startswith("MT-")
    ad.var["ercc"] = ad.var_names.str.upper().str.startswith("ERCC")
    ad.var["ribo"] = ad.var_names.str.upper().str.startswith(("RPS", "RPL"))
    ad.var["hb"] = ad.var_names.str.contains("^HB(?!P)", regex=True)

    # QC metrics from counts (no log1p!)
    sc.pp.calculate_qc_metrics(
        ad,
        qc_vars=["mt", "ercc", "ribo", "hb"],
        layer="counts",
        log1p=False,
        inplace=True,
    )

    # remove unwanted genes (applies to X and layers)
    unwanted = ad.var["mt"] | ad.var["ribo"] | ad.var["hb"] | ad.var["ercc"]
    ad = ad[:, ~unwanted].copy()

    # basic filters
    sc.pp.filter_cells(ad, min_genes=min_genes)
    sc.pp.filter_genes(ad, min_cells=min_cells)

    # MAD-based outlier filtering
    out = (
        is_outlier(ad.obs["total_counts"].values, nmads_counts) |
        is_outlier(ad.obs["n_genes_by_counts"].values, nmads_counts) |
        is_outlier(ad.obs["pct_counts_mt"].values, nmads_mt) |
        is_outlier(ad.obs["pct_counts_in_top_50_genes"].values, nmads_top50)
    )
    ad = ad[~out].copy()

    # keep counts intact
    ensure_counts_in_layers(ad)
    return ad


def run_scrublet_on_counts(ad: sc.AnnData, figpath: Path = None) -> sc.AnnData:
    Xc = ad.layers["counts"]
    Xc = Xc.toarray() if issparse(Xc) else np.asarray(Xc)
    scru = scr.Scrublet(Xc)
    scores, preds = scru.scrub_doublets()
    ad.obs["doublet_score"] = scores
    ad.obs["predicted_doublet"] = preds
    ad = ad[~preds].copy()

    if figpath is not None:
        figpath.parent.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(6, 4))
        plt.hist(scores, bins=50)
        plt.axvline(scru.threshold_, linewidth=2)
        plt.title(f"Scrublet scores")
        plt.xlabel("score"); plt.ylabel("cells")
        plt.tight_layout(); plt.savefig(figpath, dpi=120); plt.close()
    ensure_counts_in_layers(ad)
    return ad


def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "scrublet_histograms").mkdir(exist_ok=True)

    # 1) Load all samples as AnnData with counts in layers
    adatas = load_sample_tables(input_dir, Path(args.metadata_csv))

    # 2) QC + Scrublet per-sample
    postqc = []
    for ad in adatas:
        sid = ad.obs["sample_id"].unique()[0]
        ad = qc_and_filter(
            ad,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            nmads_counts=args.nmads_counts,
            nmads_mt=args.nmads_mt,
            nmads_top50=args.nmads_top50,
        )
        ad = run_scrublet_on_counts(ad, figpath=outdir / "scrublet_histograms" / f"{sid}_scrublet_hist.png")
        if args.write_per_sample:
            ad.write(outdir / f"{sid}_postqc_counts.h5ad", compression="gzip")
        postqc.append(ad)

    # 3) Concatenate all samples to a single postQC, pre-HVG object
    ad_all = an.concat(postqc, join="inner", label="sample_id_from",
                       keys=[a.obs["sample_id"].unique()[0] for a in postqc],
                       index_unique=None)
    ad_all.obs_names_make_unique()

    # ensure counts present and mirrored in .X
    if "counts" not in ad_all.layers:
        ad_all.layers["counts"] = ad_all.X
    ad_all.X = ad_all.layers["counts"]

    outpath = outdir / "adata_all_postqc_preHVG_GSE125527_counts.h5ad"
    ad_all.write(outpath, compression="gzip")
    print(f"[DONE] Wrote: {outpath} | cells={ad_all.n_obs} | genes={ad_all.n_vars}")
    print("[NEXT] Use this file as the dataset_1 input for pooled HVG + scVI (layer='counts').")


if __name__ == "__main__":
    main()
