#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 08_nmf
#SBATCH --array=1-6
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 500GB
#SBATCH --time 48:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches COMBOS = [("mcx","GM"), ("mcx","WM"), ("sc","GM"),
# ("sc","WM"), ("mcx","whole"), ("sc","whole")] in scripts/08_nmf.py --
# the last two are whole-tissue runs (GM+WM pooled, not split), added
# per explicit user request. All 6 tasks are independent and run in
# parallel. No GPU requested: ArchetypalAnalysis (py_pcha) is a CPU-only
# (numpy) algorithm, same as 07a/07b (also CPU now, not GPU-trained --
# see that pair's own job scripts/CLAUDE.md section).
#
# Must run after jobs/07b_c2l.sh's matching tissue task has finished
# (reads that step's saved adata_vis_after.h5ad) -- e.g.:
#   id2=$(sbatch --parsable jobs/07b_c2l.sh)
#   sbatch --dependency=afterok:$id2 jobs/08_nmf.sh
# (afterok here, not aftercorr: mcx tasks (array indices 1,2,5 of this
# job) depend on 07b's mcx task, and sc tasks (3,4,6) on 07b's sc task
# -- not a clean 1:1 index match even before the whole-tissue combos
# made it less contiguous, so wait on the whole 07b array instead of
# trying to line up individual indices.)

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/08_nmf.py
