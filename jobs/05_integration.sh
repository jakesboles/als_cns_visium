#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 05_integration
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 300GB
#SBATCH --time 48:00:00
#SBATCH --output /projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript /projects/b1169/boles/als_cns_visium/scripts/05_integration.R
