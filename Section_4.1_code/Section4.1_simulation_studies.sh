#!/usr/bin/env bash

#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH -J eSCM_sim
#SBATCH -o logs/log.%A_%a.out
#SBATCH -e logs/log.%A_%a.err
#SBATCH --mail-type=begin
#SBATCH --mail-type=end
#SBATCH --time=23:00:00
#SBATCH --array=1-12

# Place this script, Section4.1_simulation_studies.R, and param_simulation.txt
# in the same project directory. Before job submission, run mkdir -p logs.
# Submit from that directory with sbatch Section4.1_simulation_studies.sh.
# The 12 array tasks correspond to the 12 parameter rows.

echo '-------------------------------'
cd "${SLURM_SUBMIT_DIR}"
echo ${SLURM_SUBMIT_DIR}
echo Running on host $(hostname)
echo Time is $(date)
echo SLURM_NODES are $(echo ${SLURM_NODELIST})
echo Script Section4.1_simulation_studies.R
echo SLURM_ARRAY_JOB_ID is ${SLURM_ARRAY_JOB_ID}
echo SLURM_ARRAY_TASK_ID is ${SLURM_ARRAY_TASK_ID}
echo '-------------------------------'
echo -e '\n\n'

export TOTALPROCS=${SLURM_CPUS_ON_NODE}

# Replace with the output directory specified in the R script.
mkdir -p /path/to/your/output

# Adjust the R module name to match your cluster.
module load R/4.4.2-gfbf-2024a

################################################################################
# Read one parameter row per array task: p shape (no header).
# p = DAG size; shape = Pareto shape parameter alpha_0.
################################################################################

# This file is read from the directory where the job was submitted.
PARAM_FILE="param_simulation.txt"

PARAMS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${PARAM_FILE}")

if [[ -z "${PARAMS}" ]]; then
    echo "ERROR: No parameter row found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

read -r p shape <<< "${PARAMS}"

if [[ -z "${p}" || -z "${shape}" ]]; then
    echo "ERROR: param_simulation.txt must contain two columns: p shape" >&2
    exit 1
fi

echo '-------------------------------'
echo "Parameters for this job:"
echo "p     = ${p}"
echo "shape = ${shape}"
echo '-------------------------------'

Rscript ./Section4.1_simulation_studies.R "${p}" "${shape}"






