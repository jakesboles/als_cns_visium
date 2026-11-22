# Step 1 of Cell2Location: trains cell2location's negative-binomial
# regression model on the scRNAseq reference to get per-cell-type
# average expression signatures, for later use (07b_c2l.py) as the
# reference cell_state_df in the actual spatial deconvolution model.
#
# Runs as a 4-task SLURM array (jobs/07a_c2l.sh): 2 tissues (motor
# cortex, cervical spinal cord) x 2 label sets. Per the user: this
# project has actually run cell2location against two different scRNAseq
# annotation columns -- "discrete" (cell_type3, one label per cell, no
# further splitting) and "microglia_subtypes" (cell_type4, identical to
# cell_type3 except microglia cells carry finer Microglia0-N subcluster
# labels instead of one combined "Microglia" label) -- and both need to
# keep running, not just the microglia_subtypes one. Checked the
# colleague's original pre-cleanup draft directly for this: it only ever
# referenced labels_key='cell_type4', with no cell_type3 branch anywhere
# -- that second label set was real work that had been done outside
# what got committed to this repo, not something dropped by a previous
# edit here.
#
# No change needed in 06a_c2l_object_prep.R/06b_c2l_object_prep.py for
# this: als_cns_scrnaseq's 20_subcluster_phenotyping.R (already read in
# detail for the previous cleanup -- see that script's lines ~258-271)
# keeps BOTH cell_type3 and cell_type4 as separate columns on the same
# full_metadata.rds object, so 06's existing per-tissue *_ref.h5ad files
# already carry both label columns in .obs; only which column gets used
# as labels_key needed to become a second parametrized axis here.
#
# Each (tissue, label set) combo gets its own reference, trained only
# against that tissue's own scRNAseq cells with that label set as the
# cell-type key. All 4 tasks are independent and run in parallel.
#
# Cleaned up from a colleague's stale draft (absolute paths under a
# different project's /gpfs/projects/b1169/thomas/... tree; hardcoded to
# motor cortex only, with no spinal cord equivalent at all, and only the
# microglia_subtypes label set, as described above):
# - Parametrized by SLURM_ARRAY_TASK_ID over TISSUES x LABEL_SETS,
#   matching this project's own established array-job convention (e.g.
#   als_cns_scrnaseq's wgcna_consensus_cns.R) rather than introducing a
#   new commandArgs/params-file convention for just a 4-element combo
#   list.
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
#
# Reproducibility of the diagnostic plots, per the user: these are
# currently throwaway diagnostics, but some may end up as supplemental
# figures, and re-running the full GPU training just to redraw a plot
# would be wasteful. Checked directly against scvi-tools 1.1.6's source
# what mod.save()/.load() actually preserves (not guessed):
# - self.history_ (what mod.plot_history() reads, i.e. the ELBO training
#   loss curve) IS an attribute whose name ends in "_", which is exactly
#   what save()'s own _get_user_attributes() filter looks for and
#   load()'s _initialize_model() restores onto the reloaded model via
#   setattr() -- so in principle it round-trips. But save()'s own
#   docstring explicitly says "Neither the trainer optimizer state nor
#   the trainer history are saved", which directly contradicts that. This
#   session has no way to actually run the save/train/reload cycle to
#   settle which claim is right for these specific Pyro-based model
#   classes, so rather than gamble on it, mod.history_["elbo_train"] is
#   now also written to a plain CSV right after training, while it's
#   definitely still live in memory -- this guarantees the training curve
#   is reproducible later regardless of how that ambiguity resolves.
# - self.samples (export_posterior()'s result, needed by plot_QC() --
#   it explicitly raises RuntimeError without it) has NO trailing
#   underscore, so it's definitely NOT swept into save()'s attr_dict --
#   confirmed, not ambiguous. A freshly mod.load()-ed model's plot_QC()
#   would fail immediately unless export_posterior() is rerun first (a
#   GPU posterior-sampling pass -- much cheaper than retraining, but not
#   free, and the first_reconstruction_plot's use_n_obs=1000 random
#   subsample isn't seeded, so even a rerun wouldn't reproduce the exact
#   same plot). So plot_QC()'s few lines are reimplemented inline below
#   (not called via super() anymore) purely so the actual (observed,
#   predicted) values it plots can be saved as a CSV alongside the PDF --
#   the rendered plot is identical either way.
# - The second reconstruction plot and the training-curve plot already
#   have their underlying arrays as plain local variables (aver,
#   inf_aver_plot, mod.history_) -- those just needed a .to_csv() added,
#   no reimplementation.
# - filter_genes.pdf needs no new saving at all: it's a pure function of
#   data/06_c2l_object_prep/{tissue}_ref.h5ad (already saved, untouched
#   by this script), so rerunning filter_genes() on that file later
#   reproduces it exactly, with no GPU and no model involved.
#
# Real cluster failure found and fixed after the first actual production
# run of this script (everything above this had only ever been verified
# by reading source, never executed): mod.view_anndata_setup() crashed
# with "AttributeError: module 'rich' has no attribute 'pretty'". Traced
# directly into scvi-tools 1.1.6's own source
# (scvi/data/_manager.py:509, inside view_setup_method_args() ->
# view_anndata_setup()): that file does `import rich` and
# `from rich import box`/`from rich.console import Console`, but never
# `import rich.pretty` anywhere, even though it calls
# `rich.pretty.pprint(setup_args)` -- it's relying on `rich.pretty`
# already being registered as an attribute of the top-level `rich`
# package by some unrelated import earlier in the process, which is true
# for some rich versions and not others. With envs/cell2location.yml's
# pinned rich==15.0.0, it isn't true, so the attribute lookup fails.
# Fixed with the one-line, no-risk workaround below rather than touching
# the rich pin (which could have other transitive-dependency effects):
# `import rich.pretty` on its own, run anywhere before
# view_anndata_setup() is called, is enough -- Python's import system
# binds a submodule as an attribute of its parent package as a side
# effect of importing it (see the Python docs on the import system), and
# that binding persists for the rest of the process regardless of what
# scvi-tools/rich do internally afterward. 07b_c2l.py calls the identical
# scvi-tools code path (its own mod.view_anndata_setup() call) and got
# the same fix for the same reason.

import os
import sys
import rich.pretty
import scanpy as sc
import numpy as np
import pandas as pd
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

# Maps a human-readable name (used for output directory naming) to the
# real scRNAseq metadata column used as labels_key -- see header note
# above for what each one is.
LABEL_SETS = {
    "discrete": "cell_type3",
    "microglia_subtypes": "cell_type4",
}

COMBOS = [(tissue, label_set) for tissue in TISSUES for label_set in LABEL_SETS]

task_id_str = os.environ.get("SLURM_ARRAY_TASK_ID")
if task_id_str is None:
    sys.exit("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run "
              "as a SLURM job array (see jobs/07a_c2l.sh), one task per "
              "(tissue, label set) combo in COMBOS, not as a standalone "
              "invocation.")
task_id = int(task_id_str)
if task_id < 1 or task_id > len(COMBOS):
    sys.exit(f"SLURM_ARRAY_TASK_ID ({task_id}) is out of range for "
             f"{len(COMBOS)} combos -- check the --array range in "
             f"jobs/07a_c2l.sh.")
tissue, label_set = COMBOS[task_id - 1]
labels_key = LABEL_SETS[label_set]
combo = f"{tissue}_{label_set}"

print(f"INFO: Started! Combo: {combo} (task {task_id}/{len(COMBOS)})")

data_dir = "data/06_c2l_object_prep/"
results_folder = f"results/07_c2l/{combo}"

for sub in ["plots", "models", "adata_objects"]:
    os.makedirs(f"{results_folder}/{sub}", exist_ok=True)

sc_ref = sc.read_h5ad(f"{data_dir}{tissue}_ref.h5ad")

# Permissive gene filtering per Cell2Location's own recommendation.
selected = filter_genes(sc_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)

plt.savefig(f"{results_folder}/plots/filter_genes.pdf", bbox_inches="tight")
plt.close()

sc_ref = sc_ref[:, selected].copy()

c2l.models.RegressionModel.setup_anndata(adata=sc_ref,
                                         labels_key=labels_key,
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

# First plot -- reimplemented inline from PltExportMixin.plot_QC()
# (bayraktarlab/cell2location/cell2location/models/base/_pyro_mixin.py)
# rather than called via super(), purely so the actual (observed,
# predicted) values get saved, not just rendered -- see header note
# above. Produces the identical plot; the only addition is the .csv.
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
}).to_csv(f"{results_folder}/plots/first_reconstruction_plot_data.csv", index=False)

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
plt.savefig(f"{results_folder}/plots/first_reconstruction_plot.pdf", format="pdf", bbox_inches="tight")
plt.close()
plt.cla()
plt.clf()

# Second plot -- aver/inf_aver_plot were already plain local variables,
# just needed a .to_csv() alongside the existing plot (see header note
# above).
aver = mod._compute_cluster_averages(key=REGISTRY_KEYS.LABELS_KEY)
aver = aver[mod.factor_names_]
inf_aver_plot = mod.samples[f"post_sample_{summary_name}"]["per_cluster_mu_fg"].T
if scale_average_detection and ("detection_y_c" in list(mod.samples[f"post_sample_{summary_name}"].keys())):
  inf_aver_plot = inf_aver_plot * mod.samples[f"post_sample_{summary_name}"]["detection_y_c"].mean()

pd.DataFrame({
    "cluster_mean_expression": aver.values.flatten(),
    "estimated_expression": inf_aver_plot.flatten(),
}).to_csv(f"{results_folder}/plots/second_reconstruction_plot_data.csv", index=False)

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

# Training curve -- mod.history_["elbo_train"] saved to CSV alongside
# the PDF, since this is the one plot whose underlying data might not
# survive a later mod.load() at all (see header note above).
mod.history_["elbo_train"].to_csv(f"{results_folder}/plots/training_loss.csv")

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
