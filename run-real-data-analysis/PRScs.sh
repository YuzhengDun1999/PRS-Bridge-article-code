#!/usr/bin/env bash
#SBATCH --job-name=PRScs
#SBATCH --output=PRScs_%A_%a.out
#SBATCH --error=PRScs_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-659%200
# Total tasks = 2 refs × 22 chr × 3 a × 5 phi = 660
# %200 caps concurrency (adjust for your cluster)

set -euo pipefail

trait="${1:?Usage: sbatch PRScs.sh TRAIT}"

# -----------------------------
# Create output directories (like the R script)
# -----------------------------
mkdir -p "${trait}/1kg/PRScs/result" "${trait}/ukbb/PRScs/result"

# -----------------------------
# Parameter grids
# -----------------------------
a_list=("0.5" "1" "1.5")
phi_list=("0" "1" "0.01" "0.0001" "0.000001")
phi_name_list=("_auto" "_1" "_e2" "_e4" "_e6")

# Reference settings mapping
ref_names=("ukbb" "1kg")
ref_dirs=("data/ldblk_ukbb_eur" "data/ldblk_1kg_eur")

# -----------------------------
# Task mapping
# -----------------------------
task_id="${SLURM_ARRAY_TASK_ID}"

# Per reference: 22 * 3 * 5 = 330 tasks
tasks_per_ref=330
ref_idx=$(( task_id / tasks_per_ref ))     # 0..1
within=$(( task_id % tasks_per_ref ))      # 0..329

# Per chromosome: 3 * 5 = 15 tasks
tasks_per_chr=15
chr=$(( within / tasks_per_chr + 1 ))      # 1..22
combo=$(( within % tasks_per_chr ))        # 0..14

# combo encodes (phi_idx, a_idx):
# Let phi vary slowest (as in your loops): phi_i outside, a inside
# -> a_idx cycles fastest
phi_idx=$(( combo / 3 ))                   # 0..4
a_idx=$(( combo % 3 ))                     # 0..2

phi="${phi_list[$phi_idx]}"
phi_tag="${phi_name_list[$phi_idx]}"
a="${a_list[$a_idx]}"

ref_name="${ref_names[$ref_idx]}"
ref_dir="${ref_dirs[$ref_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/PRScs_sumdat.txt"
ldsc="${trait}/data/chr${chr}/ldsc-hm3.txt"
outdir="${trait}/${ref_name}/PRScs/result/chr${chr}"

mkdir -p "${outdir}"
N=$(awk 'NR>1 {sum+=$4; n++} END {print sum/n}' "${ldsc}")

echo "Trait   : ${trait}"
echo "Ref     : ${ref_name} (${ref_dir})"
echo "Chr     : ${chr}"
echo "a       : ${a}"
echo "phi     : ${phi}"
echo "N_gwas  : ${N}"
echo "sumdat  : ${sumdat}"
echo "ldsc    : ${ldsc}"
echo "outdir  : ${outdir}"

# PRS-CS command
if [[ "${phi}" == "0" ]]; then
  python PRScs/PRScs.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}"
else
  python PRScs/PRScs.py \
    --ref_dir="${ref_dir}" \
    --bim_prefix="chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas="${N}" \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --phi="${phi}"
fi