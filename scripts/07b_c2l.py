# Step 2 of Cell2Location: maps 07a_c2l.py's trained reference cell-type
# signatures onto the Visium data itself, estimating absolute cell-type
# abundance at each spot.
#
# Runs as a 4-task SLURM array (jobs/07b_c2l.sh): 2 tissues x 2 label
# sets, matching 07a_c2l.py's own COMBOS exactly (see that script's
# header for what the two label sets -- "discrete"/cell_type3 and
# "microglia_subtypes"/cell_type4 -- actually are). Must run after
# 07a_c2l.py's array job has finished for the SAME combo (see
# jobs/07b_c2l.sh for the --dependency note), since this reads that
# step's saved reference model directly. The Visium query itself
# (06b_c2l_object_prep.py's {tissue}_query.h5ad) doesn't depend on label
# set at all -- only the scRNAseq reference side does -- so the same
# query file is read for both of a tissue's label-set combos.
#
# Cleaned up from a colleague's stale draft the same way as 07a_c2l.py
# (see that script's header for the general approach). Specific fixes
# here, both confirmed rather than guessed:
# - batch_key="sample_id" doesn't exist anywhere in this project's Visium
#   metadata -- the real per-Visium-section id column is "code" (see
#   CLAUDE.md's cohort/sample-naming notes; confirmed there's no
#   "sample_id" column anywhere in 04_spot_annotation.R's saved
#   metadata). This would have raised a KeyError the moment
#   setup_anndata() tried to read it. Replaced with "code" throughout,
#   including the per-section adata export loop at the bottom, which
#   read the same nonexistent column.
# - Added the output directory creation this script's save calls need --
#   same gap as 07a_c2l.py.
#
# Reproducibility of the diagnostic plots (mod.plot_history()'s training
# curve, mod.plot_QC()'s reconstruction-accuracy scatter), per the user:
# same fix and same reasoning as 07a_c2l.py's header -- see that script
# for the full explanation (checked directly against scvi-tools 1.1.6
# source for what mod.save()/.load() actually preserves). Short version:
# mod.history_["elbo_train"] is saved to CSV right after training in
# case a reload doesn't restore it, and mod.plot_QC() is reimplemented
# inline (not called directly) so the (observed, predicted) values
# behind c2l_training_curve.pdf/reconstruction_plot.pdf are saved as
# plain data too -- self.samples (needed by plot_QC()) is confirmed NOT
# restored by mod.load() at all (no trailing underscore, so it's not in
# save()'s attr_dict), and plot_QC()'s own random use_n_obs subsample
# isn't seeded, so there'd be no way to recreate the exact same plot
# later even after rerunning export_posterior().
#
# Directory rework, per explicit user request ("I hate the way Thomas
# saved data... I don't want many nested folders"): data objects (both
# scvi model save dirs, all the adata_vis_*.h5ad files) now live under
# data/07_c2l/ as flat, combo-prefixed files/dirs, matching
# 06a/06b_c2l_object_prep's own {tissue}_ref.h5ad-style naming -- see
# 07a_c2l.py's header for the full reasoning (same rework, same script
# pair). Diagnostic plots/CSVs stay one level deep under
# results/07_c2l/<combo>/, flat inside (no more models/adata_objects
# subfolders). This script's own inf_aver.csv write is now the ONLY one
# -- 07a_c2l.py used to also write a file by that same name into the
# same shared folder, recomputing the identical table from the same
# saved data; since this script always runs after 07a and the content
# was identical either way, 07a's copy was dropped rather than kept as
# a guaranteed-stale duplicate (see 07a_c2l.py's header).
#
# Real cluster failure found in 07a_c2l.py's first production run (see
# that script's header for the full writeup): mod.view_anndata_setup()
# (called below too) crashes with "AttributeError: module 'rich' has no
# attribute 'pretty'" under envs/cell2location.yml's pinned rich==15.0.0,
# because scvi-tools 1.1.6's own scvi/data/_manager.py calls
# rich.pretty.pprint() without ever importing rich.pretty itself, relying
# on some other import having already registered it as a side effect --
# which this rich version doesn't do. Fixed the same way as 07a_c2l.py:
# an explicit `import rich.pretty` below, which permanently registers
# that attribute on the `rich` package for the rest of the process.
#
# Switched to CPU training (accelerator='cpu' in both mod.train() and
# export_posterior()'s sample_kwargs), per explicit user request to
# avoid this cluster's long GPU queue times -- same change, same
# reasoning as 07a_c2l.py (confirmed via source that neither
# scvi-tools' Pyro mixins nor cell2location's own source hardcode a
# CUDA-only device anywhere; accelerator is a plain passthrough to
# PyTorch Lightning's Trainer). jobs/07b_c2l.sh's resources were
# reworked to match -- see that file. This script's Visium query is
# much smaller than 07a's full scRNAseq reference (spots per section
# vs. hundreds of thousands of cells), so CPU wall-clock here is
# plausibly less painful than for 07a, but that's not measured in this
# session either (no cluster access) -- watch the first real run.
#
# CPU thread count, per explicit follow-up request ("can more
# cores/workers be assigned to speed things up") -- same change, same
# reasoning as 07a_c2l.py: torch.set_num_threads() pinned to
# jobs/07b_c2l.sh's --cpus-per-task, plus OMP_NUM_THREADS/MKL_NUM_THREADS
# exported in that job script too, so PyTorch's (and any other BLAS-
# backed library's) intra-op parallelism actually uses the cores the job
# requests. Multi-process CPU parallelism (scvi-tools' `device`/
# `devices` argument spawning several processes under a DDP strategy,
# rather than more threads within one process) was deliberately NOT
# turned on here either -- see 07a_c2l.py's header for why: nothing
# found confirming this is a validated path for Pyro-based models like
# Cell2location specifically, and a wrong DDP setup could silently
# distort the ELBO's minibatch-rescaling math rather than crash
# outright, which isn't a risk worth taking without being able to test
# it first.

import os
import sys
import rich.pretty
import torch
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import cell2location
from scvi import REGISTRY_KEYS

sys.stdout.reconfigure(line_buffering=True)

os.chdir("/projects/b1169/boles/als_cns_visium")

# Pin PyTorch's intra-op thread pool to the job's actual CPU allocation
# -- see header note above. Falls back to os.cpu_count() so this still
# behaves sensibly if run outside SLURM (SLURM_CPUS_PER_TASK unset).
n_threads = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))
torch.set_num_threads(n_threads)
print(f"INFO: Using {n_threads} CPU threads for PyTorch")

TISSUES = ["mcx", "sc"]

# Must match 07a_c2l.py's LABEL_SETS exactly -- only the human-readable
# name is used here (for locating 07a's output directory); the real
# column name doesn't matter on this side since the Visium query object
# never carries a cell-type label column at all.
LABEL_SETS = ["discrete", "microglia_subtypes"]

COMBOS = [(tissue, label_set) for tissue in TISSUES for label_set in LABEL_SETS]

task_id_str = os.environ.get("SLURM_ARRAY_TASK_ID")
if task_id_str is None:
    sys.exit("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run "
              "as a SLURM job array (see jobs/07b_c2l.sh), one task per "
              "(tissue, label set) combo in COMBOS, not as a standalone "
              "invocation.")
task_id = int(task_id_str)
if task_id < 1 or task_id > len(COMBOS):
    sys.exit(f"SLURM_ARRAY_TASK_ID ({task_id}) is out of range for "
             f"{len(COMBOS)} combos -- check the --array range in "
             f"jobs/07b_c2l.sh.")
tissue, label_set = COMBOS[task_id - 1]
combo = f"{tissue}_{label_set}"

print(f"INFO: Started! Combo: {combo} (task {task_id}/{len(COMBOS)})")

input_dir = "data/06_c2l_object_prep/"
data_dir = "data/07_c2l/"
results_folder = f"results/07_c2l/{combo}"

os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_folder, exist_ok=True)

# Load post-regression reference data and Visium data -- the query file
# is tissue-only (not label-set-specific), see header note above.
adata_ref = sc.read(f"{data_dir}{combo}_ref_after_training.h5ad")
adata_vis = sc.read(f"{input_dir}{tissue}_query.h5ad")

# Create "SYMBOL" variable for gene names
adata_vis.var['SYMBOL'] = adata_vis.var.index

# Remove MT and HB genes, store original expression in obsm
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var_names]
adata_vis.var['HB_gene'] = [gene.startswith('HB') for gene in adata_vis.var_names]
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]
adata_vis.obsm['HB'] = adata_vis[:, adata_vis.var['HB_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['HB_gene'].values]

# Load regression model
mod = cell2location.models.RegressionModel.load(f"{data_dir}{combo}_regression_model", adata_ref)

# Export estimated expression in each cell type
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.to_csv(f"{results_folder}/inf_aver.csv")

# Subset to shared genes
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
inf_aver.to_csv(f"{results_folder}/inf_aver_genes_used.csv")

# Prepare Visium data for cell2location model -- batch_key is "code"
# (per-Visium-section id), not "sample_id" (doesn't exist here). See
# header note above.
cell2location.models.Cell2location.setup_anndata(
    adata=adata_vis, batch_key="code",
    continuous_covariate_keys=["age"], categorical_covariate_keys=["sex"]
)

# Initialize model: https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # Expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology
    N_cells_per_location=7,
    # Hyperparameter controlling normalisation of within-experiment variation in RNA detection
    detection_alpha=20  # Recommended to test 20 and 200 (see tutorial for details) - lower for more within-sample technical effects
)

# View anndata_setup as a sanity check (only prints in log, doesn't save to pdf)
mod.view_anndata_setup()

# Train model
print("INFO: Training!")
mod.train(max_epochs=250,
          batch_size=10000,  # Set to None to train using all data
          train_size=1,  # Use all data points in training because we need to estimate cell abundance at all locations
          accelerator='cpu')

# Plot ELBO loss history for all training epochs -- also saved as plain
# CSV (see header note above).
mod.history_["elbo_train"].to_csv(f"{results_folder}/c2l_training_curve.csv")

mod.plot_history(0)
plt.legend(labels=['full data training'])
plt.savefig(f"{results_folder}/c2l_training_curve.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

# Export the estimated cell abundance (summary of the posterior distribution)
# https://github.com/BayraktarLab/cell2location/issues/316
print("INFO: Exporting posterior!")
adata_vis = mod.export_posterior(
    adata_vis, use_quantiles=True,
    add_to_obsm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'accelerator': 'cpu'}
)

# Save model
mod.save(f"{data_dir}{combo}_c2l_model", overwrite=True)

expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results -- one h5ad per Visium section ("code"),
# plus the combined object.
for section in adata_vis.obs['code'].unique():
    adat = adata_vis[adata_vis.obs['code'] == section].copy()
    adata_file = f"{data_dir}{combo}_adata_vis_{section}.h5ad"
    adat.write(adata_file)

adata_vis.write_h5ad(f"{data_dir}{combo}_adata_vis_after.h5ad")

# plot_QC works for c2l model - only outputs one plot -- reimplemented
# inline from PltExportMixin.plot_QC() (not called directly anymore) so
# the underlying (observed, predicted) values are saved as a CSV
# alongside the PDF, not just rendered and discarded. use_n_obs=1000
# matches plot_QC()'s own default (not passed explicitly below since
# this call never overrode it either). Produces the identical plot --
# see 07a_c2l.py's header for the full reasoning.
use_n_obs = 1000
summary_name = "q50"
ind_x = np.random.choice(mod.adata_manager.adata.n_obs, np.min((use_n_obs, mod.adata.n_obs)), replace=False)
expected_nb_param = mod.module.model.compute_expected(
    mod.samples[f"post_sample_{summary_name}"], mod.adata_manager, ind_x=ind_x
)
x_data = mod.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[ind_x, :]
if hasattr(x_data, "toarray"):
    x_data = np.asarray(x_data.toarray())

pd.DataFrame({
    "observed": x_data.flatten(),
    "predicted_mu": expected_nb_param["mu"].flatten(),
}).to_csv(f"{results_folder}/reconstruction_plot_data.csv", index=False)

plt.figure(figsize=(5.5, 5.5))
plt.hist2d(
    np.log10(x_data.flatten() + 1),
    np.log10(expected_nb_param["mu"].flatten() + 1),
    bins=50,
    norm=mpl.colors.LogNorm(),
)
plt.gca().set_aspect("equal", adjustable="box")
plt.xlabel("Data, log10")
plt.ylabel("Posterior expected value, log10")
plt.title("Reconstruction accuracy")
plt.tight_layout()
plt.savefig(f"{results_folder}/reconstruction_plot.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

print("INFO: Done")
