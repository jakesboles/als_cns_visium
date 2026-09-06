#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name wgcna_consensus
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 16
#SBATCH --mem 200G
#SBATCH --time 24:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%j.log
#SBATCH --verbose

# Single task, no --array -- wgcna_consensus.R processes all 4 anatomical
# compartments (mcx GM, mcx WM, sc GM, sc WM) in one run via hdWGCNA's own
# consensus mechanism (SetMultiExpr()/TestSoftPowersConsensus()), unlike
# als_cns_scrnaseq's wgcna_single.R/wgcna_consensus.R, which array over a
# separate target (cell type) per task.
#
# --ntasks-per-node matches enableWGCNAThreads(nThreads = 16) in the script.
#
# --mem/--time carried over from als_cns_scrnaseq's WGCNA job scripts --
# wgcna_single.sh, wgcna_consensus.sh, and wgcna_consensus_muscle.sh all
# use this same 200G/24:00:00 sizing for the same kind of metacell +
# (consensus) network construction work. An unmeasured estimate there
# too, not a measurement specific to this dataset -- TOM/network
# construction cost scales mainly with the number of selected genes
# rather than cell/spot count, so it's a reasonable starting point despite
# Visium having far fewer spots than that project has cells. Check
# `seff <jobid>` once this runs and adjust.
#
# NOTE: all three of those sibling-repo WGCNA jobs submit under
# --account b1042 --partition genomics, not b1169/b1169 like every other
# job in this repo -- possibly because b1169's partition isn't sized for
# a 200G/16-core request. Worth confirming b1169 can actually serve that
# before relying on this; switch to b1042/genomics (or whatever
# allocation covers this scale) if b1169 rejects the job or queues it
# indefinitely.

module load R/4.4.0
module load hdf5/1.14.1-2-gcc-12.3.0

Rscript "/gpfs/projects/b1169/boles/als_cns_visium/scripts/wgcna_consensus.R"
