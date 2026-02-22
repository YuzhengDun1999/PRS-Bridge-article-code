#!/usr/bin/env bash
#SBATCH --job-name=PRSBridge_sim
#SBATCH --output=PRSBridge_sim_%A_%a.out
#SBATCH --error=PRSBridge_sim_%A_%a.err
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=0-263%100
# 264 tasks total:
# 22 chromosomes × 3 alpha × 4 percent = 264
# %100 caps concurrency at 100 running tasks (adjust as needed)

set -euo pipefail

RHO="${1:?Usage: sbatch PRSBridge.sh RHO GA}"
GA="${2:?Usage: sbatch PRSBridge.sh RHO GA}"

trait="rho${RHO}GA${GA}"

# -----------------------------
# Parameter grids
# -----------------------------
alpha_list=("0.5" "0.25" "0.125")
percent_list=("0.2" "0.4" "0.6" "0.8")

# Per chromosome: 3 × 4 = 12 combos
combos_per_chr=12

task_id="${SLURM_ARRAY_TASK_ID}"

chr=$(( task_id / combos_per_chr + 1 ))   # 1..22
combo=$(( task_id % combos_per_chr ))     # 0..11

alpha_idx=$(( combo / 4 ))                # 0..2
percent_idx=$(( combo % 4 ))              # 0..3

alpha="${alpha_list[$alpha_idx]}"
percent="${percent_list[$percent_idx]}"

REF_PATH="data/ref_1kg/chr"
sumdat="${trait}/sumdat/chr${chr}/hm3_sumdat.txt"
outdir="${trait}/1kg/Bridge_small/chr${chr}/percent${percent}/alpha${alpha}"

mkdir -p "${outdir}"

echo "Trait   : ${trait}"
echo "Chr     : ${chr}"
echo "Alpha   : ${alpha}"
echo "Percent : ${percent}"
echo "Sumdat  : ${sumdat}"
echo "Outdir  : ${outdir}"

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