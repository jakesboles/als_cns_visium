#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name deseq2_by_compartment
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem 32GB
#SBATCH --time 2:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

# --mem/--time are an unmeasured estimate, not a measurement -- this loads
# the full CNS BPCells matrix once, then loops over 4 compartments
# pseudobulking and running DESeq2 per compartment. Sized well below
# als_cns_scrnaseq's per-tissue deseq2.R job (64G/12:00:00 per tissue,
# looping over a variable number of cell types) since there are only 4
# compartments here and far fewer genes captured per Visium spot than per
# cell. Check `seff <jobid>` after this runs and adjust.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript "/gpfs/projects/b1169/boles/als_cns_visium/scripts/deseq2_by_compartment.R"
