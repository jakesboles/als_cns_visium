#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 03c_make_halo_feature_gdfs
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 100GB
#SBATCH --time 4:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

# Load python
module purge
module load python-miniconda3/4.12.0
module load mamba/23.1.0

# Load environment
source activate "/projects/b1169/thomas/StardistEnv/env/Stardist"

python3 "/gpfs/projects/b1169/boles/als_cns_visium/scripts/03c_make_halo_feature_gdfs.py"
