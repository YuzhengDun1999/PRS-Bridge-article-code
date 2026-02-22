#!/usr/bin/env bash
#SBATCH --job-name=PRSBridge
#SBATCH --output=PRSBridge_%A_%a.out
#SBATCH --error=PRSBridge_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --array=0-1407%200

set -euo pipefail
trait="${1:?Usage: sbatch PRSBridge.sh TRAIT}"

# -----------------------------
# Parameter grids
# -----------------------------
alpha_list=("auto" "0.5" "0.25" "0.125")
percent_list=("0.2" "0.4" "0.6" "0.8")
ref_paths=("data/ref_1kg_large/chr" "data/ref_ukbb/chr" "data/ref_ukbb_large/chr" "data/ref_1kg/chr")
out_subdirs=("1kg/Bridge_large" "ukbb/Bridge_small" "ukbb/Bridge_large" "1kg/Bridge_small")

task_id="${SLURM_ARRAY_TASK_ID}"

# 22 chr × 4 alpha × 4 percent = 352 tasks per ref setting
tasks_per_ref=352

ref_idx=$(( task_id / tasks_per_ref ))   # 0..3
within=$(( task_id % tasks_per_ref ))    # 0..351

# Within each ref setting:
# 16 combos per chr (4 alpha × 4 percent)
chr=$(( within / 16 + 1 ))               # 1..22
combo=$(( within % 16 ))                 # 0..15
alpha_idx=$(( combo / 4 ))               # 0..3
percent_idx=$(( combo % 4 ))             # 0..3

alpha="${alpha_list[$alpha_idx]}"
percent="${percent_list[$percent_idx]}"
REF_PATH="${ref_paths[$ref_idx]}"
OUT_BASE="${out_subdirs[$ref_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/hm3_sumdat.txt"

outdir="${trait}/${OUT_BASE}/result/chr${chr}/percent${percent}/alpha${alpha}"
mkdir -p "${outdir}"

echo "Trait      : ${trait}"
echo "Ref idx    : ${ref_idx}"
echo "REF_PATH   : ${REF_PATH}"
echo "OUT_BASE   : ${OUT_BASE}"
echo "Chr        : ${chr}"
echo "Alpha      : ${alpha}"
echo "Percent    : ${percent}"
echo "Sumdat     : ${sumdat}"
echo "Outdir     : ${outdir}"

python PRSBridge/PRSBridge.py \
  --percent "${percent}" \
  --chr "${chr}" \
  --alpha "${alpha}" \
  --ref "${REF_PATH}" \
  --sumdat "${sumdat}" \
  --h2 0 \
  --h2_se 0 \
  --method cg \
  --output "${outdir}"