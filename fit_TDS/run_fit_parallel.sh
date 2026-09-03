#!/bin/bash

#SBATCH --job-name=boron-tds-fit
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=300:00:00
#SBATCH --output=fit_%j.out
#SBATCH --error=fit_%j.err

set -euo pipefail

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONUNBUFFERED=1

cd "${SLURM_SUBMIT_DIR}"

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate boron-films-env

N_TRAPS="$1"
CACHE_ROOT="${SLURM_TMPDIR:-/tmp/${USER}}/festim_cache_${SLURM_JOB_ID}"

mkdir -p "${CACHE_ROOT}/xdg"
export XDG_CACHE_HOME="${CACHE_ROOT}/xdg"
trap 'rm -rf "${CACHE_ROOT}"' EXIT

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Model order: ${N_TRAPS}"

srun --mpi=none --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores python -u fit_TDS/fit.py "${N_TRAPS}"