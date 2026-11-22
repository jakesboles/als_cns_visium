import os
import scanpy as sc
import numpy as np
import pandas as pd
import scipy as sp
import torch
import cell2location as c2l
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from matplotlib import rcParams
import matplotlib.pyplot as plt
import matplotlib as mpl
from scvi import REGISTRY_KEYS
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

print("INFO: Started!")

sc_ref = sc.read_h5ad('/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects/mcx_ref.h5ad')
results_folder = "/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/mcx"

# permissive gene filtering per C2L's recommendation
selected = filter_genes(sc_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

# save plot
plt.savefig(f"{results_folder}/plots/filter_genes.pdf", bbox_inches = "tight")
plt.close()

# filter the object
sc_ref = sc_ref[:, selected].copy()

# set up NB regression
c2l.models.RegressionModel.setup_anndata(adata=sc_ref,
                                         labels_key='cell_type4',
                                         batch_key = "batch"
                                         )
                                         
# NB reg
mod = RegressionModel(sc_ref)
# view anndata_setup as a sanity check
mod.view_anndata_setup()
print("INFO: Training!")
# train model
mod.train(max_epochs=250, batch_size=5000, train_size=1, lr=0.002, use_gpu=True)

mod.plot_history(0)
plt.savefig(f"{results_folder}/plots/training_loss.pdf")
plt.close()
plt.cla()
plt.clf()

# Export posterior
sc_ref = mod.export_posterior(
    sc_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 5000, 'use_gpu': True}
)

# export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in sc_ref.varm.keys():
    inf_aver = sc_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}' 
                                    for i in sc_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = sc_ref.var[[f'means_per_cluster_mu_fg_{i}' 
                                    for i in sc_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = sc_ref.uns['mod']['factor_names']

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(sc_ref.var_names, inf_aver.index)
inf_aver = inf_aver.loc[intersect, :].copy()

# save predictions to CSV
inf_aver.to_csv(f"{results_folder}/inf_aver.csv")

summary_name = "means"
use_n_obs = 1000
scale_average_detection = True

# First plot
super(type(mod), mod).plot_QC(summary_name = summary_name, use_n_obs = use_n_obs)
plt.savefig(f"{results_folder}/plots/first_reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Second plot
inf_aver = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
  inf_aver = inf_aver * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]

plt.hist2d(
  np.log10(aver.values.flatten() + 1),
  np.log10(inf_aver.flatten() + 1),
  bins=50,
  norm=mpl.colors.LogNorm(),
)

plt.xlabel("Mean expression for every gene in every cluster")
plt.ylabel("Estimated expression for every gene in every cluster")

plt.savefig(f"{results_folder}/plots/second_reconstruction_plot.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# Training curve
mod.plot_history(0)
plt.savefig(f"{results_folder}/plots/training_loss.pdf", format = "pdf", bbox_inches = "tight")
plt.close()
plt.cla()
plt.clf()

# print to log file
print("INFO: Saving!")

# Save model
mod.save(f"{results_folder}/models/model", overwrite=True)

# Save anndata objects with results
adata_file = f"{results_folder}/adata_objects/ref_after_training.h5ad"
sc_ref.write(adata_file)
