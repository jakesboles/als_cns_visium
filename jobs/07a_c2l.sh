#!/bin/bash
#SBATCH --account b1042
#SBATCH --partition genomics-gpu
#SBATCH --job-name 07a_c2l
#SBATCH --array=1-2
#SBATCH --gres gpu:a100:1
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 180GB
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches TISSUES = ["mcx", "sc"] in scripts/07a_c2l.py (task 1 =
# mcx, task 2 = sc) -- both tasks are independent and run in parallel.
#
# Must run (and finish) before jobs/07b_c2l.sh's matching tissue task,
# since that step loads this one's saved reference model directly.
# Submit with a dependency so each tissue's mapping step only waits on
# that same tissue's own regression step, not both:
#   id1=$(sbatch --parsable jobs/07a_c2l.sh)
#   sbatch --dependency=aftercorr:$id1 jobs/07b_c2l.sh

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/07a_c2l.py
