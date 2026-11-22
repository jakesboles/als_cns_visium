#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --job-name 07a_c2l
#SBATCH --array=1-4
#SBATCH --gres gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches COMBOS = [(tissue, label_set) for tissue in
# ["mcx","sc"] for label_set in ["discrete","microglia_subtypes"]] in
# scripts/07a_c2l.py: task 1 = mcx/discrete, 2 = mcx/microglia_subtypes,
# 3 = sc/discrete, 4 = sc/microglia_subtypes. All 4 tasks are independent
# and run in parallel.
#
# Must run (and finish) before jobs/07b_c2l.sh's matching combo task,
# since that step loads this one's saved reference model directly.
# Submit with a dependency so each combo's mapping step only waits on
# that same combo's own regression step:
#   id1=$(sbatch --parsable jobs/07a_c2l.sh)
#   sbatch --dependency=aftercorr:$id1 jobs/07b_c2l.sh

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/07a_c2l.py
