#!/usr/bin/env bash
#SBATCH --job-name=PRScs_reg
#SBATCH --output=PRScs_reg_%A_%a.out
#SBATCH --error=PRScs_reg_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-219%200
# Total tasks = 2 refs × 22 chr × 5 phi = 220
# %200 caps concurrency (adjust to your cluster)

set -euo pipefail

trait="${1:?Usage: sbatch PRScs_regularized.sh TRAIT}"

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p "${trait}/1kg/PRScs_regularized/result" "${trait}/ukbb/PRScs_regularized/result"

# -----------------------------
# Parameter grids
# -----------------------------
a="1"
regularized="1"

phi_list=("0" "1" "0.01" "0.0001" "0.000001")
phi_name_list=("_auto" "_1" "_e2" "_e4" "_e6")  # kept for logging

ref_names=("1kg" "ukbb")
ref_dirs=("data/ldblk_1kg_eur" "data/ldblk_ukbb_eur")

# -----------------------------
# Task mapping
# -----------------------------
task_id="${SLURM_ARRAY_TASK_ID}"

# Per ref: 22 * 5 = 110 tasks
tasks_per_ref=110
ref_idx=$(( task_id / tasks_per_ref ))   # 0..1
within=$(( task_id % tasks_per_ref ))    # 0..109

# Per chr: 5 tasks (phi)
tasks_per_chr=5
chr=$(( within / tasks_per_chr + 1 ))    # 1..22
phi_idx=$(( within % tasks_per_chr ))    # 0..4

phi="${phi_list[$phi_idx]}"
phi_tag="${phi_name_list[$phi_idx]}"

ref="${ref_names[$ref_idx]}"
ref_dir="${ref_dirs[$ref_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/PRScs_sumdat.txt"
ldsc="${trait}/data/chr${chr}/ldsc-hm3.txt"
outdir="${trait}/${ref}/PRScs_regularized/result/regularized${regularized}/chr${chr}"

mkdir -p "${outdir}"

# Compute N = floor(mean(N)) from ldsc-hm3.txt
N=$(awk 'NR>1 {sum+=$4; n++} END {print sum/n}' "${ldsc}")


echo "Trait       : ${trait}"
echo "Ref         : ${ref} (${ref_dir})"
echo "Chr         : ${chr}"
echo "a           : ${a}"
echo "phi         : ${phi} (${phi_tag})"
echo "regularized : ${regularized}"
echo "N_gwas      : ${N}"
echo "sumdat      : ${sumdat}"
echo "ldsc        : ${ldsc}"
echo "outdir      : ${outdir}"

# PRS-CS-Regularized command
# NOTE: --bim_prefix should point to validation BIM prefix; your R script used chr${chr}.
# Adjust to a full path prefix if needed.
if [[ "${phi}" == "0" ]]; then
  python PRS-CS-proj/PRScs_Regularized.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --regularized="${regularized}"
else
  python PRS-CS-proj/PRScs_Regularized.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --phi="${phi}" \
    --regularized="${regularized}"
fi