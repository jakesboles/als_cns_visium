#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 06a_c2l_object_prep
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 4
#SBATCH --mem 150GB
#SBATCH --time 12:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

# Single task -- writes both tissues' .mtx/.csv triplets in one run (the
# script's own loop over tissues), not an array job. Sizing is an
# estimate (holding one tissue's filtered scRNAseq counts matrix in
# memory at a time is the most memory-hungry step) -- check `seff
# <jobid>` after the first run and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_visium/scripts/06a_c2l_object_prep.R
