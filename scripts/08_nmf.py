# NMF-style co-location analysis (Cell2Location's downstream
# ArchetypalAnalysis model, run via run_colocation()) on 07b_c2l.py's
# per-spot cell-type abundance estimates: finds "archetypal" tissue zones
# -- spatially coherent combinations of co-located cell types -- run
# separately per tissue AND per anatomical compartment (gray vs. white
# matter), since these are expected to have genuinely different cell-type
# compositions and shouldn't be pooled into one factorization.
#
# Runs as a 6-task SLURM array (jobs/08_nmf.sh): mcx GM, mcx WM, sc GM,
# sc WM, plus two whole-tissue combos (mcx GM+WM together, sc GM+WM
# together) added per explicit user request, so the archetypal
# factorization can also be compared against one run per tissue that
# doesn't split gray/white matter apart at all. All 6 are independent
# and run in parallel. A "whole" combo still excludes Meninges/Nerve
# bundle spots (region %in% c("GM","WM")), matching what the 4 original
# combos always excluded too -- the only difference from the original 4
# is whether GM and WM are pooled into one factorization or kept
# separate, not which spots are eligible in the first place.
#
# Cleaned up from a colleague's stale draft (hardcoded to motor
# cortex/white matter only, with a commented-out alternative
# CoLocatedGroupsSklearnNMF call left in place). Real bugs fixed, both
# confirmed directly against cell2location's own source (not guessed):
# - train_args={'sample_id': 'sample_id', ...} -- run_colocation() never
#   reads a key called "sample_id" at all; it reads "sample_name_col"
#   (confirmed in cell2location/run_colocation.py). As written, this key
#   was silently ignored, and run_colocation() would have fallen back to
#   collapsing every spot onto a single pseudo-sample name when
#   "sample_name_col" isn't given explicitly. Fixed to the real key,
#   pointing at "code" (this project's actual per-Visium-section id
#   column -- "sample_id" doesn't exist in this project's metadata at
#   all).
# - The active ArchetypalAnalysis call never passed model_kwargs, so
#   run_colocation()'s own sklearn-NMF-shaped defaults (init, alpha,
#   nmf_kwd_args) got forwarded into ArchetypalAnalysis.__init__(), which
#   doesn't accept any of them (confirmed directly against
#   cell2location/models/downstream/ArchetypalAnalysis.py's signature) --
#   this would have raised a TypeError before training ever started.
#   Fixed by passing model_kwargs={'random_state': 0} explicitly (the one
#   kwarg ArchetypalAnalysis actually takes that the old defaults also
#   set).
# - py_pcha (ArchetypalAnalysis.fit()'s lazily-imported dependency,
#   "please first run `pip install py_pcha`" per its own source) was
#   missing from envs/cell2location.yml entirely -- added there.
# - Read path fixed to match where 07b_c2l.py actually saves
#   adata_vis_after.h5ad (results/07_c2l/<tissue>/, no adata_objects/
#   subfolder) -- the draft read from a .../adata_objects/
#   adata_vis_after.h5ad path that 07b's own save call never wrote to,
#   confirmed by comparing the two scripts directly.
#
# Dropped the commented-out CoLocatedGroupsSklearnNMF alternative (dead
# code -- ArchetypalAnalysis is the one actually invoked, per the user)
# and the draft's unused flax/squidpy/anndata-adjacent imports this
# script never calls directly.
#
# LABEL_SET: 07a/07b now each run per (tissue, label set) combo (two
# scRNAseq annotation columns, "discrete" and "microglia_subtypes" --
# see 07a_c2l.py's header), so results/07_c2l/<tissue>/ no longer
# exists as a path at all -- it's results/07_c2l/<tissue>_<label_set>/.
# This script was originally only asked to split 4 ways by (tissue,
# region), now 6 with the whole-tissue combos above -- either way it
# needs ONE label set's cell-type abundances to run NMF on. Confirmed
# with the user: "microglia_subtypes" is the one this should actually
# run against (not just a placeholder default anymore).
#
# Directory rework, per explicit user request ("I hate the way Thomas
# saved data... I don't want many nested folders"): 07b_c2l.py's
# adata_vis_after.h5ad now lives under data/07_c2l/ (combo-prefixed),
# not results/07_c2l/<combo>/ -- see that script's header. This
# script's own data objects (AfterNMF_<combo>.h5ad, nmf_res_<combo>.pkl)
# move the same way, to data/08_nmf/. The results/08_nmf/<combo>/ folder
# passed below as export_args={'path': ...} keeps whatever internal
# subfolder structure (stability_plots/, predictive_accuracy/,
# factor_markers/, etc.) cell2location's own run_colocation() builds --
# confirmed directly against cell2location/run_colocation.py that this
# nesting is baked into the library's own export logic, not something
# this script controls, so it's left alone rather than fought.

import os
import sys
import pickle
import numpy as np
import scanpy as sc
from cell2location import run_colocation

sys.stdout.reconfigure(line_buffering=True)

os.chdir("/projects/b1169/boles/als_cns_visium")

# Which 07a/07b label-set combo this reads cell-type abundances from --
# confirmed with the user, see header note above.
LABEL_SET = "microglia_subtypes"

# (tissue, region) combos -- "mcx"/"sc" tissue codes and "GM"/"WM" region
# codes match 04_spot_annotation.R's own `region` column values exactly.
# "whole" isn't a real `region` value -- it's this script's own label for
# "don't split GM/WM apart", handled specially in the filter step below.
COMBOS = [
    ("mcx", "GM"), ("mcx", "WM"), ("sc", "GM"), ("sc", "WM"),
    ("mcx", "whole"), ("sc", "whole"),
]

task_id_str = os.environ.get("SLURM_ARRAY_TASK_ID")
if task_id_str is None:
    sys.exit("SLURM_ARRAY_TASK_ID is not set -- this script is meant to run "
              "as a SLURM job array (see jobs/08_nmf.sh), one task per "
              "(tissue, region) combo in COMBOS, not as a standalone "
              "invocation.")
task_id = int(task_id_str)
if task_id < 1 or task_id > len(COMBOS):
    sys.exit(f"SLURM_ARRAY_TASK_ID ({task_id}) is out of range for "
             f"{len(COMBOS)} combos -- check the --array range in "
             f"jobs/08_nmf.sh.")
tissue, region = COMBOS[task_id - 1]
combo = f"{tissue}_{region}"

print(f"INFO: Started! Combo: {combo} (task {task_id}/{len(COMBOS)}, label set: {LABEL_SET})")

input_dir = "data/07_c2l/"
data_dir = "data/08_nmf/"
results_folder = f"results/08_nmf/{combo}/"
os.makedirs(data_dir, exist_ok=True)
os.makedirs(results_folder, exist_ok=True)

adata_vis = sc.read_h5ad(f"{input_dir}{tissue}_{LABEL_SET}_adata_vis_after.h5ad")

if region == "whole":
    adata_vis = adata_vis[adata_vis.obs['region'].isin(["GM", "WM"])].copy()
else:
    adata_vis = adata_vis[adata_vis.obs['region'] == region].copy()

res_dict, adata_vis = run_colocation(
    adata_vis,
    model_name='ArchetypalAnalysis',
    verbose=True,
    train_args={
        'n_fact': np.arange(5, 20),  # IMPORTANT: use a wider range of the number of factors (5-30)
        'sample_name_col': 'code',   # per-Visium-section id -- see header note above
        'n_iter': 5000,
        'random_state': 0
    },
    # ArchetypalAnalysis only accepts random_state/pcha_kwd_args (not
    # run_colocation()'s own sklearn-NMF-shaped defaults) -- see header
    # note above.
    model_kwargs={'random_state': 0},
    export_args={'path': results_folder}
)

adata_vis.write_h5ad(f"{data_dir}{combo}_AfterNMF.h5ad")

with open(f"{data_dir}{combo}_nmf_res.pkl", "wb") as f:
    pickle.dump(res_dict, f)

print("INFO: Done")
