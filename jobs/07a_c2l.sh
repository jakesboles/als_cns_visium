#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --job-name MCX_S1
#SBATCH --gres=gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/logs/%x_oe%j.log
#SBATCH --verbose
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.watson1@northwestern.edu

# State date and name
date
echo 'Thomas Watson'

# Load python
module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

# Load environment
source activate "/projects/b1169/thomas/Cell2LocationEnv/env/Cell2Location1024"

# Run script
python3 "/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx/mcx_step1.py"