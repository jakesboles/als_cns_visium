#!/bin/bash

#SBATCH --job-name=module_preservation
#SBATCH --array=1-2
#SBATCH --output=/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --time=48:00:00
#SBATCH --mem=300G
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --account=b1169
#SBATCH --partition=b1169

# --array matches jobs/wgcna_module_preservation_with_sc_params.txt's
# current 2 rows (Microglia, both directions) -- widen back to 1-6 (and
# add the Oligodendrocyte/Astrocyte rows back to that file) if/when this
# expands beyond Microglia.

cd /projects/b1169/boles/als_cns_visium

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

# Filename typo fixed (was "..._params.txt.txt", which doesn't exist) --
# the real file has a single .txt extension.
PARAMS_FILE="/projects/b1169/boles/als_cns_visium/jobs/wgcna_module_preservation_with_sc_params.txt"

PARAM1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f1 -d,)
PARAM2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f2 -d,)
NAME1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f3 -d,)
NAME2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f4 -d,)
TYPE1=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f5 -d,)
TYPE2=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $PARAMS_FILE | cut -f6 -d,)

echo "Assessing preservation of modules from ${NAME1} in ${NAME2}"

# Was "module_preservation_analysis.R", a filename that doesn't exist in
# this repo -- the actual script this job runs.
Rscript /projects/b1169/boles/als_cns_visium/scripts/wgcna_module_preservation_with_scrna.R "${PARAM1}" "${PARAM2}" "${NAME1}" "${NAME2}" "${TYPE1}" "${TYPE2}"