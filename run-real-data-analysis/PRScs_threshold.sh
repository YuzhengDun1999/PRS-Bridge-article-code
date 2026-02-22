#!/usr/bin/env bash
#SBATCH --job-name=PRScs_thr
#SBATCH --output=PRScs_thr_%A_%a.out
#SBATCH --error=PRScs_thr_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-1099%200
# Total tasks = 2 refs × 22 chr × 5 phi × 5 threshold = 1100
# %200 caps concurrency (adjust to your cluster)

set -euo pipefail

trait="${1:?Usage: sbatch PRScs_threshold.sh TRAIT}"

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p "${trait}/1kg/PRScs_threshold/result" "${trait}/ukbb/PRScs_threshold/result"
# -----------------------------
# Parameter grids
# -----------------------------
a="1"
phi_list=("0" "1" "0.01" "0.0001" "0.000001")
phi_name_list=("_auto" "_1" "_e2" "_e4" "_e6")  # kept for logging
threshold_list=("0.01" "0.1" "1" "10" "100")
ref_names=("ukbb" "1kg")
ref_dirs=("data/ldblk_ukbb_eur" "data/ldblk_1kg_eur")

# -----------------------------
# Task mapping
# -----------------------------
task_id="${SLURM_ARRAY_TASK_ID}"

# Per ref: 22 * 5 * 5 = 550
tasks_per_ref=550
ref_idx=$(( task_id / tasks_per_ref ))   # 0..1
within=$(( task_id % tasks_per_ref ))    # 0..549

# Per chr: 5 phi * 5 threshold = 25
tasks_per_chr=25
chr=$(( within / tasks_per_chr + 1 ))    # 1..22
combo=$(( within % tasks_per_chr ))      # 0..24

# Let threshold vary fastest, phi vary slowest (matches your loops: phi -> threshold -> a)
phi_idx=$(( combo / 5 ))                 # 0..4
thr_idx=$(( combo % 5 ))                 # 0..4

phi="${phi_list[$phi_idx]}"
phi_tag="${phi_name_list[$phi_idx]}"
threshold="${threshold_list[$thr_idx]}"

ref="${ref_names[$ref_idx]}"
ref_dir="${ref_dirs[$ref_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/PRScs_sumdat.txt"
ldsc="${trait}/data/chr${chr}/ldsc-hm3.txt"
outdir="${trait}/${ref}/PRScs_threshold/result/threshold${threshold}/chr${chr}"

mkdir -p "${outdir}"

N=$(awk 'NR>1 {sum+=$4; n++} END {print sum/n}' "${ldsc}")

echo "Trait     : ${trait}"
echo "Ref       : ${ref} (${ref_dir})"
echo "Chr       : ${chr}"
echo "a         : ${a}"
echo "phi       : ${phi} (${phi_tag})"
echo "threshold : ${threshold}"
echo "N_gwas    : ${N}"
echo "sumdat    : ${sumdat}"
echo "ldsc      : ${ldsc}"
echo "outdir    : ${outdir}"

# PRS-CS-threshold command
if [[ "${phi}" == "0" ]]; then
  python PRS-CS-proj/PRScs_threshold.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --threshold="${threshold}"
else
  python PRS-CS-proj/PRScs_threshold.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --phi="${phi}" \
    --threshold="${threshold}"
fi