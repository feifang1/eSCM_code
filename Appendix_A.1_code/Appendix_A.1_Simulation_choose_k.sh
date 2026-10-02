#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -J eSCM_sim
#SBATCH -o logs/log.%A_%a.out
#SBATCH -e logs/log.%A_%a.err
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --time=14:00:00
#SBATCH --array=1-36

# Before submission: place this script, the R script, and params_v7.txt in
# the same project directory, then run mkdir -p logs before calling sbatch.
# Rename the attached params_v7(1).txt to params_v7.txt.

echo '-------------------------------'
cd "${SLURM_SUBMIT_DIR}"
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo Script Simulation_choose_k_lambda_v8.R
echo SLURM_ARRAY_JOB_ID is ${SLURM_ARRAY_JOB_ID}
echo SLURM_ARRAY_TASK_ID is ${SLURM_ARRAY_TASK_ID}
echo '-------------------------------'
echo -e '\n\n'

mkdir -p logs
# Replace with the output directory specified in the R script.
mkdir -p /path/to/your/output

export TOTALPROCS=${SLURM_CPUS_ON_NODE}

# Adjust the R module name to match your cluster.
module load R/4.4.2-gfbf-2024a

# Read one row from params_v7.txt according to the array task ID
# Expected column order: n p shape lambda quantile_val
PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" params_v7.txt)

if [[ -z "${PARAMS}" ]]; then
  echo "ERROR: No parameter row found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 1
fi

read -r n p shape lambda quantile_val <<< "${PARAMS}"

if [[ -z "${quantile_val}" ]]; then
  echo "ERROR: params_v7.txt must contain five columns: n p shape lambda quantile_val" >&2
  exit 1
fi

echo '-------------------------------'
echo "Parameters for this job:"
echo "n      = ${n}"
echo "p      = ${p}"
echo "shape  = ${shape}"
echo "lambda       = ${lambda}"
echo "quantile_val = ${quantile_val}"
echo '-------------------------------'

Rscript ./Simulation_choose_k_lambda_v8.R "${n}" "${p}" "${shape}" "${lambda}" "${quantile_val}"






