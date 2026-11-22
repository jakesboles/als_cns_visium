#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --job-name 07b_c2l
#SBATCH --array=1-2
#SBATCH --gres gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 300GB
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches TISSUES = ["mcx", "sc"] in scripts/07b_c2l.py, same
# indexing as jobs/07a_c2l.sh -- see that script's header for the
# --dependency=aftercorr submission pattern this depends on (this step
# reads 07a's saved reference model for the same tissue directly).

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/07b_c2l.py
