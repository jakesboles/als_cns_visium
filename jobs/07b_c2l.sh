#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07b_c2l
#SBATCH --array=1-4
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 300GB
#SBATCH --time 48:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches COMBOS in scripts/07b_c2l.py, built the same way (and
# in the same order) as 07a_c2l.py's own COMBOS -- see that script's
# header for the --dependency=aftercorr submission pattern this depends
# on (this step reads 07a's saved reference model for the same combo
# directly).
#
# No GPU request here anymore (dropped --gres gpu:a100:1, switched off
# b1042/genomics-gpu) -- scripts/07b_c2l.py was switched to
# accelerator='cpu', same change and same reasoning as 07a_c2l.py's own
# CPU switch (see that script's/job's header), per explicit user
# request after confirming 07a's switch. Account/partition moved to
# b1169/b1169, matching 07a_c2l.sh and 08_nmf.sh's own CPU-only
# convention. --cpus-per-task 16/--time 48:00:00 are flagged, unmeasured
# starting points, not confirmed-correct numbers -- this script's Visium
# query is much smaller than 07a's scRNAseq reference, so CPU wall-clock
# here is plausibly shorter, but that's not benchmarked in this session
# (no cluster access); check `seff <jobid>` after the first real run and
# adjust.

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

python3 /projects/b1169/boles/als_cns_visium/scripts/07b_c2l.py
