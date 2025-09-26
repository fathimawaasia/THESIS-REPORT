#!/usr/bin/env bash
#SBATCH --job-name=scenic_dc2_boost
#SBATCH --partition=emerald-shared
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --chdir=/rds/projects/e/elhamsak-ibd-single-cell/THESIS/pySCENIC_TFs
#SBATCH --output=logs/dc2_grn_%j.out
#SBATCH --error=logs/dc2_grn_%j.err

set -euo pipefail
mkdir -p logs

# ---- env ----
eval "$(conda shell.bash hook)"
conda activate scenic_env

# avoid thread oversubscription
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1

# ---- inputs ----
TF_LIST=resources/hs_hgnc_tfs.txt
DB1=resources/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather
DB2=resources/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather
MOT=resources/motifs-v10-nr.hgnc-m0.00001-o0.0.tbl
MTX=SCENIC_INPUT/DC2_expr.tsv.gz

# sanity
for f in "$TF_LIST" "$DB1" "$DB2" "$MOT" "$MTX"; do
  [[ -s "$f" ]] || { echo "[ERR] Missing: $f"; exit 1; }
done

echo "===== [DC2] GRN (GRNBoost2) ====="
[[ -s DC2_adjacencies.tsv ]] || pyscenic grn "$MTX" "$TF_LIST" \
  -o DC2_adjacencies.tsv --method grnboost2 --num_workers 8

echo "===== [DC2] ctx ====="
[[ -s DC2_regulons.csv ]] || pyscenic ctx DC2_adjacencies.tsv "$DB1" "$DB2" \
  --annotations_fname "$MOT" \
  --expression_mtx_fname "$MTX" \
  --mode custom_multiprocessing \
  --chunk_size 20000 \
  --output DC2_regulons.csv \
  --num_workers 8

echo "===== [DC2] AUCell ====="
[[ -s DC2_auc.csv ]] || pyscenic aucell "$MTX" DC2_regulons.csv \
  --num_workers 8 --output DC2_auc.csv

echo "[OK] DC2 done."
