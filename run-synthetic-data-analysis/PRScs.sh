#!/usr/bin/env bash
#SBATCH --job-name=PRScs_sim
#SBATCH --output=PRScs_sim_%A_%a.out
#SBATCH --error=PRScs_sim_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --array=0-109%100
# Total tasks = 22 chr × 5 phi = 110
# %100 caps concurrency (adjust to your cluster)

set -euo pipefail

RHO="${1:?Usage: sbatch PRScs.sh RHO GA}"
GA="${2:?Usage: sbatch PRScs.sh RHO GA}"

trait="rho${RHO}GA${GA}"

# -----------------------------
# Create output directories
# -----------------------------
mkdir -p "${trait}/1kg/PRScs/result"

# -----------------------------
# Parameter grids
# -----------------------------
a="1"
phi_list=("0" "1" "0.01" "0.0001" "0.000001")
phi_name_list=("_auto" "_1" "_e2" "_e4" "_e6")  # for logging only
phi_name_prs=("auto" "1e+00" "1e-02" "1e-04" "1e-06")  # for optional downstream naming

# -----------------------------
# Task mapping
# -----------------------------
task_id="${SLURM_ARRAY_TASK_ID}"

# Per chromosome: 5 phi settings
tasks_per_chr=5
chr=$(( task_id / tasks_per_chr + 1 ))    # 1..22
phi_idx=$(( task_id % tasks_per_chr ))    # 0..4

phi="${phi_list[$phi_idx]}"
phi_tag="${phi_name_list[$phi_idx]}"
phi_out="${phi_name_prs[$phi_idx]}"

# -----------------------------
# Inputs / outputs
# -----------------------------
sumdat="${trait}/data/chr${chr}/PRScs_sumdat.txt"
outdir="${trait}/1kg/PRScs/result/chr${chr}"

mkdir -p "${outdir}"

echo "Trait   : ${trait}"
echo "Chr     : ${chr}"
echo "a       : ${a}"
echo "phi     : ${phi} (${phi_tag})"
echo "N_gwas  : ${N}"
echo "sumdat  : ${sumdat}"
echo "outdir  : ${outdir}"

# PRS-CS command
# NOTE: --bim_prefix should point to validation BIM prefix; your R script used chr${chr}.
# Adjust to a full path prefix if needed.
REF_DIR="data/ldblk_1kg_eur"

if [[ "${phi}" == "0" ]]; then
  python PRScs/PRScs.py \
    --ref_dir="${REF_DIR}" \
    --bim_prefix="data/EUR_chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas=100000 \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}"
else
  python PRScs/PRScs.py \
    --ref_dir="${REF_DIR}" \
    --bim_prefix="data/chr${chr}" \
    --sst_file="${sumdat}" \
    --n_gwas=100000 \
    --out_dir="${outdir}" \
    --a="${a}" \
    --chrom="${chr}" \
    --phi="${phi}"
fi
