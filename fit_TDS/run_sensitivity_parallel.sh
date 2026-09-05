#!/bin/bash
#SBATCH --job-name=boron-tds-sensitivity
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=20G
#SBATCH --time=300:00:00
#SBATCH --array=0-8
#SBATCH --output=sensitivity_%A_%a.out
#SBATCH --error=sensitivity_%A_%a.err

set -euo pipefail

k0_values=(
    1.53187207e-17
    1.53187207e-17
    1.53187207e-17
    1.53187207e-16
    1.53187207e-16
    1.53187207e-16
    1.53187207e-15
    1.53187207e-15
    1.53187207e-15
)

p0_values=(
    1.0e12
    1.0e13
    1.0e14
    1.0e12
    1.0e13
    1.0e14
    1.0e12
    1.0e13
    1.0e14
)

k0="${k0_values[$SLURM_ARRAY_TASK_ID]}"
p0="${p0_values[$SLURM_ARRAY_TASK_ID]}"

cd "${SLURM_SUBMIT_DIR}"

source "$HOME/anaconda3/etc/profile.d/conda.sh"
conda activate boron-films-env

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export PYTHONUNBUFFERED=1

cache_root="${SLURM_TMPDIR:-/tmp/${USER}}/festim_cache_${SLURM_JOB_ID}"
mkdir -p "${cache_root}/xdg"
export XDG_CACHE_HOME="${cache_root}/xdg"
trap 'rm -rf "${cache_root}"' EXIT

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "k0: ${k0} m3/s"
echo "p0: ${p0} s^-1"

srun --mpi=none --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" --cpu-bind=cores python -u fit_TDS/sensitivity.py "${k0}" "${p0}"
