# Step 1 of Cell2Location: trains cell2location's negative-binomial
# regression model on the scRNAseq reference to get per-cell-type
# average expression signatures, for later use (07b_c2l.py) as the
# reference cell_state_df in the actual spatial deconvolution model.
#
# Runs as a 2-task SLURM array (jobs/07a_c2l.sh), one task per tissue
# (motor cortex, cervical spinal cord) -- each tissue gets its own
# reference, trained only against that tissue's own scRNAseq cells
# (06a_c2l_object_prep.R already split the reference by tissue for this
# reason). Both tasks are independent and run in parallel.
#
# Cleaned up from a colleague's stale draft (absolute paths under a
# different project's /gpfs/projects/b1169/thomas/... tree; hardcoded to
# motor cortex only, with no spinal cord equivalent at all):
# - Parametrized by SLURM_ARRAY_TASK_ID over TISSUES, matching this
#   project's own established array-job convention (e.g.
#   als_cns_scrnaseq's wgcna_consensus_cns.R) rather than introducing a
#   new commandArgs/params-file convention for just a 2-element list.
# - mod.train()'s `use_gpu=True` is gone: scvi-tools 1.1.6 (pinned in
#   envs/cell2location.yml, checked directly against that exact tagged
#   source) removed that kwarg from PyroSviTrainMixin.train() entirely in
#   favor of PyTorch Lightning's accelerator/device API -- as written,
#   this would have raised a TypeError before training ever started.
#   Replaced with accelerator='gpu' (also in export_posterior()'s own
#   sample_kwargs, which hit the same deprecated kwarg), matching
#   07b_c2l.py's own already-correct call.
# - Added the plots/models/adata_objects output directories this script
#   writes into -- the draft never created them, so the very first
#   plt.savefig() call would have failed outright unless someone had
#   already mkdir'd them by hand.
# - Dropped a duplicate mod.plot_history(0) call (the draft called it
#   once right after training and again, identically, right before
#   saving -- both wrote the same training_loss.pdf, so the first call
#   was pure wasted recomputation).
#
# NOTE: the second-reconstruction-plot diagnostic block (mod.
# _compute_cluster_averages(), a private method) is carried over as-is
# from the draft -- unlike the three fixes above, this one wasn't
# independently re-verified against scvi-tools 1.1.6's actual source, so
# keep an eye on it specifically if this step errors deep inside that
# plotting block.

import os
import sys
import scanpy as sc
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cell2location as c2l
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
from matplotlib import rcParams
from scvi import REGISTRY_KEYS

sys.stdout.reconfigure(line_buffering=True)
rcParams['pdf.fonttype'] = 42

os.chdir("/projects/b1169/boles/als_cns_visium")

TISSUES = ["mcx", "sc"]

task_id_str = os.environ.get("SLURM_ARRAY_TASK_ID")
if task_id_str is None:
    sys.exit("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run "
              "as a SLURM job array (see jobs/07a_c2l.sh), one task per "
              "tissue in TISSUES, not as a standalone invocation.")
task_id = int(task_id_str)
if task_id < 1 or task_id > len(TISSUES):
    sys.exit(f"SLURM_ARRAY_TASK_ID ({task_id}) is out of range for "
             f"{len(TISSUES)} tissues -- check the --array range in "
             f"jobs/07a_c2l.sh.")
tissue = TISSUES[task_id - 1]

print(f"INFO: Started! Tissue: {tissue} (task {task_id}/{len(TISSUES)})")

data_dir = "data/06_c2l_object_prep/"
results_folder = f"results/07_c2l/{tissue}"

for sub in ["plots", "models", "adata_objects"]:
    os.makedirs(f"{results_folder}/{sub}", exist_ok=True)

sc_ref = sc.read_h5ad(f"{data_dir}{tissue}_ref.h5ad")

# Permissive gene filtering per Cell2Location's own recommendation.
selected = filter_genes(sc_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

plt.savefig(f"{results_folder}/plots/filter_genes.pdf", bbox_inches="tight")
plt.close()

sc_ref = sc_ref[:, selected].copy()

c2l.models.RegressionModel.setup_anndata(adata=sc_ref,
                                         labels_key='cell_type4',
                                         batch_key="batch")

mod = RegressionModel(sc_ref)
mod.view_anndata_setup()

print("INFO: Training!")
mod.train(max_epochs=250, batch_size=5000, train_size=1, lr=0.002, accelerator='gpu')

sc_ref = mod.export_posterior(
    sc_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 5000, 'accelerator': 'gpu'}
)

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

inf_aver.to_csv(f"{results_folder}/inf_aver.csv")

summary_name = "means"
use_n_obs = 1000
scale_average_detection = True

# First plot
super(type(mod), mod).plot_QC(summary_name=summary_name, use_n_obs=use_n_obs)
plt.savefig(f"{results_folder}/plots/first_reconstruction_plot.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

# Second plot
inf_aver_plot = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
  inf_aver_plot = inf_aver_plot * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()
aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]

plt.hist2d(
  np.log10(aver.values.flatten() + 1),
  np.log10(inf_aver_plot.flatten() + 1),
  bins=50,
  norm=mpl.colors.LogNorm(),
)

plt.xlabel("Mean expression for every gene in every cluster")
plt.ylabel("Estimated expression for every gene in every cluster")

plt.savefig(f"{results_folder}/plots/second_reconstruction_plot.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

# Training curve
mod.plot_history(0)
plt.savefig(f"{results_folder}/plots/training_loss.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

print("INFO: Saving!")

mod.save(f"{results_folder}/models/model", overwrite=True)

adata_file = f"{results_folder}/adata_objects/ref_after_training.h5ad"
sc_ref.write(adata_file)

print("INFO: Done")
