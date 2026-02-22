#!/usr/bin/env bash
#SBATCH --job-name=PRScs_proj
#SBATCH --output=PRScs_proj_%A_%a.out
#SBATCH --error=PRScs_proj_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-879%200
# Total tasks = 2 refs × 22 chr × 5 phi × 4 percent = 880
# %200 caps concurrency (adjust to your cluster)

set -euo pipefail

trait="${1:?Usage: sbatch PRScs_proj.sh TRAIT}"

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p "${trait}/1kg/PRScs_proj/result" "${trait}/ukbb/PRScs_proj/result"

# -----------------------------
# Parameter grids
# -----------------------------
phi_list=("0" "1" "0.01" "0.0001" "0.000001")
phi_name_list=("_auto" "_1" "_e2" "_e4" "_e6")  # kept for consistency/logging
percent_list=("0.2" "0.4" "0.6" "0.8")

ref_names=("1kg" "ukbb")
ref_dirs=("data/ldblk_1kg_eur" "data/ldblk_ukbb_eur")

# -----------------------------
# Task mapping
# -----------------------------
task_id="${SLURM_ARRAY_TASK_ID}"

# Per ref: 22 * 5 * 4 = 440 tasks
tasks_per_ref=440
ref_idx=$(( task_id / tasks_per_ref ))     # 0..1
within=$(( task_id % tasks_per_ref ))      # 0..439

# Per chromosome: 5 * 4 = 20 tasks
tasks_per_chr=20
chr=$(( within / tasks_per_chr + 1 ))      # 1..22
combo=$(( within % tasks_per_chr ))        # 0..19

# Let percent vary fastest, phi vary slowest (matches your loop nesting: phi -> percent -> a)
phi_idx=$(( combo / 4 ))                   # 0..4
percent_idx=$(( combo % 4 ))               # 0..3

phi="${phi_list[$phi_idx]}"
phi_tag="${phi_name_list[$phi_idx]}"
percent="${percent_list[$percent_idx]}"
a="1"

ref="${ref_names[$ref_idx]}"
ref_dir="${ref_dirs[$ref_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/PRScs_sumdat.txt"
ldsc="${trait}/data/chr${chr}/ldsc-hm3.txt"
outdir="${trait}/${ref}/PRScs_proj/result/percent${percent}/chr${chr}"

mkdir -p "${outdir}"

N=$(awk 'NR>1 {sum+=$4; n++} END {print sum/n}' "${ldsc}")

echo "Trait       : ${trait}"
echo "Ref         : ${ref} (${ref_dir})"
echo "Chr         : ${chr}"
echo "a           : ${a}"
echo "phi         : ${phi} (${phi_tag})"
echo "eigenval_rm : ${percent}"
echo "N_gwas      : ${N}"
echo "sumdat      : ${sumdat}"
echo "ldsc        : ${ldsc}"
echo "outdir      : ${outdir}"

# PRS-CS-Projection command
if [[ "${phi}" == "0" ]]; then
  python PRS-CS-proj/PRScs_proj.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --eigenval_rm="${percent}"
else
  python PRS-CS-proj/PRScs_proj.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --phi="${phi}" \
    --eigenval_rm="${percent}"
fi