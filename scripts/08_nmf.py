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
#   ORIGINALLY "fixed" by passing model_kwargs={'random_state': 0}
#   explicitly -- this was NOT actually sufficient, and the real
#   production run of this script crashed with exactly the predicted
#   TypeError ("unexpected keyword argument 'init'") despite that fix
#   being in place. Root cause, confirmed directly against
#   cell2location/run_colocation.py's actual source: run_colocation()
#   builds its model_kwargs as `d_model_kwargs = {"init": "random",
#   "random_state": 0, "nmf_kwd_args": {...}, "alpha": 0.01}` and then
#   does `for k in model_kwargs.keys(): d_model_kwargs[k] =
#   model_kwargs[k]` -- a plain dict UPDATE, not a replacement. Passing
#   model_kwargs={'random_state': 0} only ever overwrote the
#   random_state entry; init/nmf_kwd_args/alpha were always still
#   present in the final dict forwarded to Model(**model_kwargs), no
#   matter what was passed in -- there's no "clear these keys" option in
#   run_colocation()'s own API. Actually fixed below by subclassing
#   ArchetypalAnalysis to absorb and discard those three kwargs
#   specifically, and passing the subclass itself (not the string
#   'ArchetypalAnalysis') as model_name -- run_colocation() accepts
#   either (`Model = getattr(models, model_name) if type(model_name) is
#   str else model_name`), so this bypasses the string lookup entirely
#   without touching the installed library.
# - py_pcha (ArchetypalAnalysis.fit()'s lazily-imported dependency,
#   "please first run `pip install py_pcha`" per its own source) was
#   missing from envs/cell2location.yml entirely -- added there.
#   REAL FAILURE after the _ArchetypalAnalysisCompat fix above got past
#   model creation: ModuleNotFoundError: No module named 'py_pcha', on
#   the actual cluster. This is NOT a gap in envs/cell2location.yml --
#   py_pcha==0.1.3 is already pinned there (confirmed by reading the
#   committed file directly) -- it means the conda env actually
#   installed at /projects/b1169/boles/als_cns_visium/envs/cell2location
#   predates that pin and was never rebuilt from the updated yml. Fix is
#   operational, not a code change: with that env active,
#   `pip install py_pcha==0.1.3` installs just the missing package
#   without a full env rebuild.
# - mod.plot_cell_type_loadings(**export_args["plot_cell_type_loadings_
#   kwargs"]) -- called unconditionally by run_colocation() (NOT wrapped
#   in a try/except, unlike its mod.plot_posterior_mu_vs_data() and
#   spatial-plotting calls) right after mod.print_gene_loadings(). This
#   method is NOT defined on ArchetypalAnalysis or its BaseModel parent
#   at all -- confirmed by grepping the entire installed cell2location
#   package -- it only exists on CoLocatedGroupsSklearnNMF, the model
#   run_colocation() was originally written around. Would have raised
#   AttributeError: 'ArchetypalAnalysis' object has no attribute
#   'plot_cell_type_loadings', likely after a long, expensive run
#   through fit()/evaluate_stability()/sample_posterior() for every
#   n_fact value in train_args['n_fact'] (5 through 19). Fixed by giving
#   _ArchetypalAnalysisCompat the same method CoLocatedGroupsSklearnNMF
#   defines (confirmed by reading its implementation directly) -- it's a
#   thin wrapper around self.plot_gene_loadings() (inherited from
#   BaseModel, which ArchetypalAnalysis also inherits) over
#   self.var_names_read/self.fact_filt/self.cell_type_fractions, all of
#   which ArchetypalAnalysis instances already have correctly populated
#   by the time run_colocation() reaches this call (self.cell_type_
#   fractions specifically is set by mod.sample2df(), called just
#   before) -- so the exact same implementation works unchanged.
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
#
# Real crash found after a full production run completed all 15 n_fact
# values (5 through 19) -- training/stability-eval/posterior-sampling/
# saving/plotting all succeeded for every single one, confirming every
# fix above this holds -- then run_colocation() itself crashed on its
# very last line, after the whole per-n_fact loop had already finished:
# `TypeError: Scalar datasets don't support chunk/filter options`,
# writing key 'SSE' to /uns/mod_coloc_n_fact5/post_sample_means, from
# `sp_data.write(filename=f"{save_path}sp.h5ad", compression="gzip")`
# (confirmed the only `compression=` write call anywhere in
# cell2location/run_colocation.py). Root cause: ArchetypalAnalysis.fit()
# stores PCHA's "SSE" value as a 0-dimensional (scalar) numpy array in
# self.samples["post_sample_means"] -- h5py simply does not support
# chunking/compression on scalar datasets (a well-known h5py limitation,
# not a cell2location-specific bug), and run_colocation() unconditionally
# passes compression="gzip" to this one internal write with no way to
# turn it off from the caller's side. Since this happens on the very
# last line of the function (after every n_fact value's real work,
# including each one's own already-pickled model_n_fact<n>.p, is already
# safely on disk), the only thing actually lost by this crash is
# run_colocation()'s own redundant combined results/.../anndata/sp.h5ad
# summary AND, critically, this script's own two final save calls below
# (data/08_nmf/<combo>_AfterNMF.h5ad, <combo>_nmf_res.pkl) -- since the
# exception propagates out of run_colocation() before it can `return
# res_dict, sp_data`, this script's own save calls never run either.
# Fixed with a narrow monkeypatch on h5py.Group.create_dataset (not on
# anything inside cell2location or anndata, so it's agnostic to both
# libraries' own internals) that drops compression/chunking kwargs
# specifically when the data being written is scalar, and behaves
# identically otherwise -- this is the standard workaround for this
# exact, well-known h5py limitation.

import os
import sys
import pickle
import numpy as np
import scanpy as sc
import h5py
from cell2location import run_colocation
from cell2location.models.downstream import ArchetypalAnalysis as _ArchetypalAnalysis

sys.stdout.reconfigure(line_buffering=True)

os.chdir("/projects/b1169/boles/als_cns_visium")

# See header note above. Scalar (0-d) datasets can't be chunked/
# compressed in h5py -- strip those kwargs only when they'd actually
# trigger that error, otherwise defer to h5py's normal behavior.
_orig_create_dataset = h5py.Group.create_dataset


def _create_dataset_allow_scalars(self, name, shape=None, dtype=None, data=None, **kwds):
    is_scalar = shape == () or (data is not None and np.asarray(data).shape == ())
    if is_scalar:
        for key in ("compression", "compression_opts", "chunks", "shuffle", "fletcher32", "scaleoffset"):
            kwds.pop(key, None)
    return _orig_create_dataset(self, name, shape=shape, dtype=dtype, data=data, **kwds)


h5py.Group.create_dataset = _create_dataset_allow_scalars


# Absorbs and discards run_colocation()'s own hardcoded sklearn-NMF-shaped
# model_kwargs defaults (init, nmf_kwd_args, alpha) that
# ArchetypalAnalysis.__init__() doesn't accept -- see header note above
# for why model_kwargs={'random_state': 0} alone isn't enough to avoid
# them.
#
# This class's __name__/__qualname__ were originally overridden to
# "ArchetypalAnalysis" purely so run_colocation()'s own run_name/
# path_name construction (which embeds mod.__class__.__name__) wouldn't
# change the exported folder names. REAL BUG this caused, found on the
# actual cluster once training/stability-eval/posterior-sampling/
# plot_cell_type_loadings (the fix above) all succeeded:
# _pickle.PicklingError: Can't pickle <class '__main__.ArchetypalAnalysis'>:
# attribute lookup ArchetypalAnalysis on __main__ failed, at
# run_colocation()'s own `pickle.dump(res_dict_1, ...)` call (res_dict_1
# includes the trained `mod` object itself). Root cause: Python's pickle
# module identifies a class by (module, __qualname__) and reconstructs
# it on load via getattr(sys.modules[module], qualname) -- it does NOT
# store the class object itself. Overriding __qualname__ to
# "ArchetypalAnalysis" made pickle look for an attribute literally named
# "ArchetypalAnalysis" on the __main__ module, which doesn't exist (this
# class is bound to the name _ArchetypalAnalysisCompat in this script's
# own namespace) -- so EVERY save of a trained model (including this
# script's own final nmf_res_<combo>.pkl, not just run_colocation()'s
# internal model_n_fact*.p files) would have failed to pickle. Fixed by
# dropping the __name__/__qualname__ override entirely: pickling now
# correctly resolves __main__._ArchetypalAnalysisCompat, and the only
# cost is a purely cosmetic one -- exported folder/file names and the
# model_name field inside saved results now read
# "_ArchetypalAnalysisCompat_..." instead of "ArchetypalAnalysis_...".
# Correctness (pickling actually working at all) matters far more than
# that cosmetic naming, so this was not worth re-attempting cleverly.
#
# A further real bug found on the same cluster once the pickling fix
# above let a run actually get through one full n_fact iteration (model
# creation, 5x training restart, stability eval, posterior sampling,
# saving, and even plot_cell_type_loadings/the spatial plots once
# plotnine was installed -- confirming every fix above this one holds):
# run_colocation() itself crashed with KeyError: 'run_time', confirmed
# NOT caused by anything in this script -- see the verbose=False
# comment at the run_colocation() call below for the exact root cause
# (a genuine typo in run_colocation()'s own source, present for ANY
# model_name whenever verbose=True) and why verbose=False is the fix.

class _ArchetypalAnalysisCompat(_ArchetypalAnalysis):
    def __init__(self, *args, init=None, nmf_kwd_args=None, alpha=None, **kwargs):
        super().__init__(*args, **kwargs)

    # Copied verbatim from CoLocatedGroupsSklearnNMF.plot_cell_type_loadings()
    # (cell2location/models/downstream/CoLocatedGroupsSklearnNMF.py) --
    # ArchetypalAnalysis doesn't define this method itself, but has every
    # attribute this implementation needs already in the right state by
    # the time run_colocation() calls it. See header note above.
    def plot_cell_type_loadings(self, selected_cell_types=None, **kwargs):
        if selected_cell_types is None:
            selected_cell_types = self.var_names_read

        self.plot_gene_loadings(
            sel_var_names=selected_cell_types,
            var_names=self.var_names_read,
            fact_filt=self.fact_filt,
            loadings_attr="cell_type_fractions",
            gene_fact_name="cell_type_fractions",
            fun_type="dotplot",
            cmap="RdPu",
            figsize=[5 + 0.12 * self.n_fact, 5 + 0.1 * self.n_var],
            **kwargs
        )


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
    model_name=_ArchetypalAnalysisCompat,
    # verbose=False, not True -- see header note above. This is a real,
    # confirmed bug in run_colocation() itself, unrelated to anything in
    # this script: with verbose=True, right after finishing each n_fact
    # value's full pipeline (training/stability-eval/posterior-sampling/
    # saving/plotting) it does `print("### Done ### - time " +
    # res_dict["run_time"])`, but res_dict (keyed by f"n_fact{n}") never
    # has a top-level "run_time" key at all -- only the per-n_fact
    # res_dict_1 dict stored at res_dict[f"n_fact{n}"] does. This raises
    # KeyError: 'run_time' unconditionally on the FIRST n_fact value
    # processed, for ANY model_name, crashing before the loop ever
    # reaches a second n_fact value -- which would have completely
    # defeated train_args['n_fact']'s whole 5-19 sweep below, regardless
    # of which model this script uses.
    verbose=False,
    train_args={
        'n_fact': np.arange(5, 20),  # IMPORTANT: use a wider range of the number of factors (5-30)
        'sample_name_col': 'code',   # per-Visium-section id -- see header note above
        'n_iter': 5000,
        'random_state': 0
    },
    # ArchetypalAnalysis only accepts random_state/pcha_kwd_args (not
    # run_colocation()'s own sklearn-NMF-shaped defaults) -- see header
    # note above for why model_kwargs alone can't filter those out, and
    # why _ArchetypalAnalysisCompat (passed as model_name above) is the
    # actual fix.
    model_kwargs={'random_state': 0},
    # plot_extension: png, not run_colocation()'s own "pdf" default, per
    # explicit user request -- covers every plot this call produces
    # (cell_type_factors/location_factors stability plots,
    # data_vs_posterior_mean, the cell_type_fractions_heatmap, and the
    # per-section spatial density plots).
    export_args={'path': results_folder, 'plot_extension': 'png'}
)

adata_vis.write_h5ad(f"{data_dir}{combo}_AfterNMF.h5ad")

with open(f"{data_dir}{combo}_nmf_res.pkl", "wb") as f:
    pickle.dump(res_dict, f)

print("INFO: Done")
