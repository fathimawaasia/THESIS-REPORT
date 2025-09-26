import scrublet as scr
import scanpy as sc
import numpy as np
import pandas as pd
import os

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
    print(f"\n✅ Running Scrublet for {sample_id}...")

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
        print(f"❌ Scrublet error in {sample_id}: {e}")
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
print("\n✅ Summary saved to 'scrublet_summary_per_sample_2.csv'")
