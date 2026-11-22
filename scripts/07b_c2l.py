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

import os
import sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import cell2location

sys.stdout.reconfigure(line_buffering=True)

os.chdir("/projects/b1169/boles/als_cns_visium")

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

data_dir = "data/06_c2l_object_prep/"
ref_folder = f"results/07_c2l/{combo}"
output_folder = f"results/07_c2l/{combo}"

os.makedirs(output_folder, exist_ok=True)

# Load post-regression reference data and Visium data -- the query file
# is tissue-only (not label-set-specific), see header note above.
adata_ref = sc.read(f"{ref_folder}/adata_objects/ref_after_training.h5ad")
adata_vis = sc.read(f"{data_dir}{tissue}_query.h5ad")

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
mod = cell2location.models.RegressionModel.load(f"{ref_folder}/models/model/", adata_ref)

# Export estimated expression in each cell type
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.to_csv(f"{ref_folder}/inf_aver.csv")

# Subset to shared genes
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()
inf_aver.to_csv(f"{ref_folder}/inf_aver_genes_used.csv")

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
          accelerator='gpu')

# Plot ELBO loss history for all training epochs
mod.plot_history(0)
plt.legend(labels=['full data training'])
plt.savefig(f"{output_folder}/c2l_training_curve.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

# Export the estimated cell abundance (summary of the posterior distribution)
# https://github.com/BayraktarLab/cell2location/issues/316
print("INFO: Exporting posterior!")
adata_vis = mod.export_posterior(
    adata_vis, use_quantiles=True,
    add_to_obsm=["q05", "q50", "q95", "q0001"],
    sample_kwargs={'batch_size': 2500, 'accelerator': 'gpu'}
)

# Save model
mod.save(f"{output_folder}", overwrite=True)

expected_dict = mod.module.model.compute_expected_per_cell_type(
    mod.samples["post_sample_q05"], mod.adata_manager
)

for i, n in enumerate(mod.factor_names_):
    adata_vis.layers[n] = expected_dict['mu'][i]

# Save anndata object with results -- one h5ad per Visium section ("code"),
# plus the combined object.
for section in adata_vis.obs['code'].unique():
    adat = adata_vis[adata_vis.obs['code'] == section].copy()
    adata_file = f"{output_folder}/adata_vis_{section}.h5ad"
    adat.write(adata_file)

adata_vis.write_h5ad(f"{output_folder}/adata_vis_after.h5ad")

# plot_QC works for c2l model - only outputs one plot
mod.plot_QC(summary_name="q50")
plt.savefig(f"{output_folder}/reconstruction_plot.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

print("INFO: Done")
