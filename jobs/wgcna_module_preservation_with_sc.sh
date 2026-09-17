#!/bin/bash

#SBATCH --job-name=module_preservation
#SBATCH --array=1-6
#SBATCH --output=/projects/b1169/thomas/als_multitissue/Visium/hdWGCNA/CNS_Consensus/module_preservation/logs/%x_%A_%a.log
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=b1169
#SBATCH --partition=b1169

cd /projects/b1169/thomas/als_multitissue/Visium/hdWGCNA/CNS_Consensus/module_preservation

module load R/4.4.0

PARAMS_FILE="module_preservation_params.txt"

# this is for the one that failed for no reason
PARAM1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)
PARAM2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f2 -d,)
NAME1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f3 -d,)
NAME2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f4 -d,)
TYPE1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f5 -d,)
TYPE2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f6 -d,)

echo "Assessing preservation of modules from ${PARAM1} in ${PARAM2}"

Rscript module_preservation_analysis.R "${PARAM1}" "${PARAM2}" "${NAME1}" "${NAME2}" "${TYPE1}" "${TYPE2}"