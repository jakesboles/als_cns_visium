#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name NMF_MCXw
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 350GB
#SBATCH --time 48:00:00
#SBATCH --output /gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/logs/%x_oe%j.log
#SBATCH --verbose

# Load python
module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

# Load environment
source activate "/projects/b1169/thomas/Cell2LocationEnv/env/Cell2Location100626"

# Run script
python3 "/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/NMF_WM.py"