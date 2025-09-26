#!/usr/bin/env bash
#SBATCH --job-name=scenic_dp_boost
#SBATCH --partition=emerald-shared
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --chdir=/rds/projects/e/elhamsak-ibd-single-cell/THESIS/pySCENIC_TFs
#SBATCH --output=logs/dp_grn_%j.out
#SBATCH --error=logs/dp_grn_%j.err

set -euo pipefail
mkdir -p logs

# ---- env ----
eval "$(conda shell.bash hook)"
conda activate scenic_env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# ---- inputs ----
TF_LIST=resources/hs_hgnc_tfs.txt
DB1=resources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
DB2=resources/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
MOT=resources/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl
MTX=SCENIC_INPUT/DP_expr.tsv.gz

for f in "$TF_LIST" "$DB1" "$DB2" "$MOT" "$MTX"; do
  [[ -s "$f" ]] || { echo "[ERR] Missing: $f"; exit 1; }
done

echo "===== [DP] GRN (GRNBoost2) ====="
[[ -s DP_adjacencies.tsv ]] || pyscenic grn "$MTX" "$TF_LIST" \
  -o DP_adjacencies.tsv --method grnboost2 --num_workers 8

echo "===== [DP] ctx ====="
[[ -s DP_regulons.csv ]] || pyscenic ctx DP_adjacencies.tsv "$DB1" "$DB2" \
  --annotations_fname "$MOT" \
  --expression_mtx_fname "$MTX" \
  --mode custom_multiprocessing \
  --chunk_size 20000 \
  --output DP_regulons.csv \
  --num_workers 8

echo "===== [DP] AUCell ====="
[[ -s DP_auc.csv ]] || pyscenic aucell "$MTX" DP_regulons.csv \
  --num_workers 8 --output DP_auc.csv

echo "[OK] DP done."
