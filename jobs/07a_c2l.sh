#!/bin/bash
#SBATCH --account b1169
#SBATCH --partition b1169
#SBATCH --job-name 07a_c2l
#SBATCH --array=1-4
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 180GB
#SBATCH --time 48:00:00
#SBATCH --output /gpfs/projects/b1169/boles/als_cns_visium/logs/%x_%A_%a.log
#SBATCH --verbose

# --array matches COMBOS = [(tissue, label_set) for tissue in
# ["mcx","sc"] for label_set in ["discrete","microglia_subtypes"]] in
# scripts/07a_c2l.py: task 1 = mcx/discrete, 2 = mcx/microglia_subtypes,
# 3 = sc/discrete, 4 = sc/microglia_subtypes. All 4 tasks are independent
# and run in parallel.
#
# No GPU request here anymore (dropped --gres gpu:a100:1, switched off
# b1042/genomics-gpu) -- scripts/07a_c2l.py was switched to
# accelerator='cpu' per explicit user request, since this cluster's GPU
# queue times are very long and this step (unlike 07b's deconvolution
# model) has no actual hard GPU requirement (confirmed directly against
# scvi-tools/cell2location source -- see that script's header). Account/
# partition moved to b1169/b1169, matching every other CPU-only job
# script in this repo (e.g. 08_nmf.sh, also scvi/cell2location-adjacent
# CPU work). --cpus-per-task 16 and --time 48:00:00 are both flagged,
# unmeasured starting points, not confirmed-correct numbers -- CPU
# training wall-clock for this step (hundreds of thousands of
# scRNAseq cells) hasn't been benchmarked anywhere in this session (no
# cluster access); check `seff <jobid>` after the first real run and
# adjust both, same as this repo's other resource-sizing caveats (e.g.
# wgcna_consensus.sh).
#
# Must run (and finish) before jobs/07b_c2l.sh's matching combo task,
# since that step loads this one's saved reference model directly.
# Submit with a dependency so each combo's mapping step only waits on
# that same combo's own regression step:
#   id1=$(sbatch --parsable jobs/07a_c2l.sh)
#   sbatch --dependency=aftercorr:$id1 jobs/07b_c2l.sh

module purge
module load python-miniconda3/4.10.3
module load mamba/23.1.0

source activate "/projects/b1169/boles/als_cns_visium/envs/cell2location"

# Make sure PyTorch's BLAS backend (and any other OMP/MKL-threaded
# library, e.g. numpy) actually uses the --cpus-per-task above, rather
# than whatever it would otherwise detect on its own -- a well-known
# PyTorch-on-SLURM gotcha. scripts/07a_c2l.py also pins
# torch.set_num_threads() to this same value directly, as a
# belt-and-suspenders measure (see that script's header).
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK

python3 /projects/b1169/boles/als_cns_visium/scripts/07a_c2l.py
