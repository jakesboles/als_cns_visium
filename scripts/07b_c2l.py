import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
os.environ["THEANO_FLAGS"] = 'device=cuda,floatX=float32,force_device=True'
import cell2location
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
import torch
from scvi import REGISTRY_KEYS

# Path to reference output folder
ref_folder = '/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx'

# Path to output folder 
output_folder = '/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx'

# Load post-regression reference data and Visium data
adata_ref = sc.read(f"{ref_folder}/adata_objects/ref_after_training.h5ad")
adata_vis = sc.read("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects/mcx_query.h5ad")

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

# Prepare Visium data for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key = "sample_id", continuous_covariate_keys = ["age"], categorical_covariate_keys = ["sex"])

# Initialize model: https://github.com/BayraktarLab/cell2location/blob/master/docs/images/Note_on_selecting_hyperparameters.pdf
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # Expected average cell abundance: tissue-dependent hyper-prior which can be estimated from paired histology
    N_cells_per_location=7, 
    # Hyperparameter controlling normalisation of within-experiment variation in RNA detection
    detection_alpha=20 # Recommended to test 20 and 200 (see tutorial for details) - lower for more within-sample technical effects
)

# View anndata_setup as a sanity check (only prints in log, doesn't save to pdf)
mod.view_anndata_setup()

# Train model 
mod.train(max_epochs=250,
          batch_size=10000, # Set to None to train using all data 
          train_size=1, # Use all data points in training because we need to estimate cell abundance at all locations
          accelerator='gpu'
         )

# Plot ELBO loss history for all training epochs
mod.plot_history(0)
plt.legend(labels=['full data training'])
plt.savefig(f"{output_folder}/c2l_training_curve.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Export the estimated cell abundance (summary of the posterior distribution)
# https://github.com/BayraktarLab/cell2location/issues/316
adata_vis = mod.export_posterior(
    # adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
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

# Save anndata object with results
for sample in adata_vis.obs['sample_id'].unique():
    adat = adata_vis[adata_vis.obs['sample_id'] == sample].copy()
    adata_file = f"{output_folder}/adata_vis_{sample}.h5ad"
    adat.write(adata_file)
    
adata_vis.write_h5ad(f"{output_folder}/adata_vis_after.h5ad")

# plot_QC works for c2l model - only outputs one plot
mod.plot_QC(summary_name = "q50")
plt.savefig(f"{output_folder}/reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()
