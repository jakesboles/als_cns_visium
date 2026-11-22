#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 06b_c2l_object_prep
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 100GB
#SBATCH --time 4:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

# Must run after jobs/06a_c2l_object_prep.sh (reads its .mtx/.csv
# output) -- submit as a dependent job, e.g.:
#   id1=$(sbatch --parsable jobs/06a_c2l_object_prep.sh)
#   sbatch --dependency=afterok:$id1 jobs/06b_c2l_object_prep.sh

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/06b_c2l_object_prep.py
