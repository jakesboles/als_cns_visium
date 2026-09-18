# NetRep module preservation between two hdWGCNA networks: this project's
# own Visium consensus network (wgcna_consensus.R) and a pre-built
# scRNAseq per-cell-type consensus network from the sister single-cell
# project (one .rds per cell type, already fit elsewhere -- not built by
# this repo). Driven by a SLURM array over
# jobs/wgcna_module_preservation_with_sc_params.txt (one row per
# comparison, 3 cell types x 2 directions = 6 rows); see that .txt and
# jobs/wgcna_module_preservation_with_sc.sh for how each task's 6
# arguments (obj1 path, obj2 path, name1, name2, type1, type2) are passed
# in. "ref"/reference is always argument 1, "test"/query is always
# argument 2 -- module preservation is assessed for the reference's
# modules as tested against the query dataset.
#
# Reworked from an earlier version of this script (kept, commented, only
# up to "stopped editing here" -- everything below that line is this
# rewrite) that predates this project's BPCells convention and loaded
# both sides as a single whole-object readRDS(). That's still correct for
# the scRNAseq side (`type == "sc"`): those per-cell-type .rds files are
# pre-built elsewhere, as whole Seurat objects, and aren't this repo's to
# change. It is NOT correct for the Visium side (`type == "spatial"`):
# wgcna_consensus.R never saves a whole Seurat object, only
# data/wgcna_consensus/wgcna_experiment.rds (just the hdWGCNA @misc
# entry -- see that script's own final comment). So for `type ==
# "spatial"`, `load_wgcna_obj()` below reconstructs the actual Seurat
# object wgcna_consensus.R ran consensus WGCNA on (04_spot_annotation.R's
# counts/metadata/images + 05_integration.R's CCA embedding, filtered to
# region %in% c("GM", "WM") -- the exact recipe documented in
# wgcna_consensus.R's own header and already reused in
# wgcna_consensus_viz.R) and reattaches the saved wgcna_experiment.rds
# into it, rather than treating the passed file path as a whole object.
# jobs/wgcna_module_preservation_with_sc_params.txt's `type == "spatial"`
# rows have been repointed at this actual save path accordingly (the old
# version pointed at a stale, pre-BPCells whole-object path that no
# longer corresponds to anything this repo produces).
#
# Other departures from the earlier version, beyond BPCells:
# - `group.by`/`group_name` for the spatial side: the earlier version
#   used a generic "all_cells_group"/"all_cells" convention (presumably
#   matching how that OLDER, now-superseded spatial WGCNA object was
#   itself built). This project's own wgcna_consensus.R instead built its
#   metacells with a constant placeholder column named `dummy`
#   (ident.group = "dummy", group_name = 1 throughout that script's own
#   SetMultiExpr()/ModuleConnectivity() calls) -- SetDatExpr(use_metacells
#   = T) subsets the METACELL object, which only carries whatever
#   grouping columns MetacellsByGroups() was actually built with, so a
#   brand new "all_cells_group" column added to `obj` after the fact
#   wouldn't propagate there. `dummy`/1 is reused instead, matching how
#   this object's metacells actually exist on disk.
# - The `names(obj@misc)[2] <- "ref"`/`"test"` rename (needed because
#   ModulePreservationNetRep() apparently wants both networks under
#   caller-chosen names, and the scRNAseq side's real internal wgcna_name
#   varies by cell type and isn't something this script controls) is
#   kept, since that's still the right strategy here, but is no longer
#   keyed to a hardcoded @misc list position -- it now renames whichever
#   entry `obj@misc$active_wgcna` actually names, which is robust to
#   @misc's element order (the earlier version's position-2 assumption
#   doesn't hold in general, and there's no indication it was ever
#   actually verified against a real object -- this script was still
#   mid-edit).
# - `SetActiveWGCNA()` (already this project's own established call, see
#   wgcna_consensus.R) is used to activate "wgcna_consensus" on the
#   reconstructed spatial object, rather than reaching into
#   obj@misc$active_wgcna directly.
# - `csv_dir`/`plots_dir` (referenced but never defined in the earlier
#   version -- a leftover from whatever directory layout the original,
#   pre-this-repo version of this script used) are dropped in favor of
#   this project's own flat data_dir/results_dir convention: everything
#   below writes directly into results_dir, matching every other script
#   in this repo (no further csvs/plots subdivision).
# - The two loading branches (spatial/sc) are factored into one
#   load_wgcna_obj() helper called once per object, instead of being
#   duplicated inline for obj1 and obj2 -- the branching logic is
#   identical either way, just parameterized by which file/name/type is
#   passed.
# - Dropped: a large commented-out, already-superseded block (a much
#   older per-array-index ProjectModules()-based design, never active in
#   this version) that was pure noise once the current design is in
#   place.
#
# NOTE: this session has no cluster access, real WGCNA objects, or R
# interpreter -- every line below is a careful read/reasoning-based
# rewrite of the earlier version plus wgcna_consensus.R's own documented
# save/reconstruction recipe, not something that's been executed. One
# real unknown: the scRNAseq side's actual assay name isn't verified here
# (SetDatExpr() is left to use that object's own DefaultAssay() for the
# `type == "sc"` branch, rather than guessing "RNA") -- if that's wrong,
# pass assay explicitly in that branch once you can check the object
# itself. Sanity-check the whole script on the cluster before treating it
# as final.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(BPCells)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(igraph)
  library(NetRep)
})

enableWGCNAThreads(nThreads = 16)
set.seed(256)

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/wgcna_cross_modality/"
dir.create(results_dir, showWarnings = F, recursive = T)

obj1_file <- commandArgs(trailingOnly = TRUE)[1]
obj2_file <- commandArgs(trailingOnly = TRUE)[2]
name1 <- commandArgs(trailingOnly = T)[3]
name2 <- commandArgs(trailingOnly = T)[4]
type1 <- commandArgs(trailingOnly = T)[5]
type2 <- commandArgs(trailingOnly = T)[6]

# Load + prepare one side of the comparison -----------------------------
# `type` is this argument's data modality ("spatial" = this project's own
# Visium consensus network, "sc" = a pre-built scRNAseq per-cell-type
# network from the sister project) -- NOT this project's own mcx/sc
# tissue code, despite the unfortunate naming overlap with "sc" =
# spinal cord elsewhere in this repo. That's jobs/
# wgcna_module_preservation_with_sc_params.txt's own established column
# values, kept as-is rather than renamed, since changing it would mean
# rewriting a params file that already has real, hand-verified cluster
# paths in it.

load_wgcna_obj <- function(file, type, name){

  if (type == "spatial"){

    # `file` here is data/wgcna_consensus/wgcna_experiment.rds -- just
    # the hdWGCNA @misc entry, not a whole Seurat object. Rebuild the
    # rest from this project's own known save locations (see header note
    # above) and graft the experiment back in.
    meta <- readRDS("data/04_spot_annotation/metadata.rds")

    raw_mat <- open_matrix_dir("data/04_spot_annotation/bpcells_data")
    raw_mat <- raw_mat[, rownames(meta)]

    images <- readRDS("data/04_spot_annotation/images.rds")

    obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta, assay = "Spatial")
    obj@images <- images

    cca <- readRDS("data/05_integration/cca.rds")
    cca@cell.embeddings <- cca@cell.embeddings[rownames(meta), ]
    obj[["cca"]] <- cca

    obj <- obj %>%
      subset(region %in% c("WM", "GM"))

    obj <- NormalizeData(obj)

    # SetDatExpr()/ModulePreservationNetRep() reach back into the "data"
    # layer for correlation-based steps that need a real CsparseMatrix --
    # BPCells' lazy matrix classes don't support that coercion (same
    # constraint documented in wgcna_consensus.R's own ModuleConnectivity()
    # call).
    obj[["Spatial"]]$data <- as(obj[["Spatial"]]$data, "dgCMatrix")

    obj@misc[["wgcna_consensus"]] <- readRDS(file)
    obj <- SetActiveWGCNA(obj, "wgcna_consensus")

    # Every spot is one group here -- there's no cell-type axis on the
    # spatial side, only the `dummy` constant wgcna_consensus.R's
    # metacells were actually built under (ident.group = "dummy"). See
    # header note above for why this is `dummy`/1, not a fresh
    # "all_cells_group"/"all_cells" column.
    obj$dummy <- 1

    obj <- SetDatExpr(obj,
                      group_name = 1,
                      group.by = "dummy",
                      use_metacells = T,
                      assay = "Spatial",
                      layer = "data")

  } else if (type == "sc"){

    meta <- readRDS("/projects/b1169/boles/als_cns_scrnaseq/data/18_full_integration/brain_sc/metadata.rds")
    meta_sub <- meta[meta$cell_type3 == "Microglia", ]
    
    raw_mat <- open_matrix_dir("/projects/b1169/boles/als_cns_scrnaseq/data/06_obj_reassembly/bpcells")
    raw_mat <- raw_mat[, rownames(meta_sub)]
    
    obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "RNA")
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj)
    
    harmony <- readRDS("/projects/b1169/boles/als_cns_scrnaseq/data/18_full_integration/brain_sc/harmony.rds")
    harmony@cell.embeddings <- harmony@cell.embeddings[rownames(meta_sub), ]
    obj[["harmony"]] <- harmony
    
    obj@misc[["wgcna_consensus"]] <- readRDS(file)
    obj <- SetActiveWGCNA(obj, "wgcna_consensus")

    obj <- SetDatExpr(obj,
                      group_name = name,
                      group.by = "cell_type3",
                      use_metacells = T)

  } else {
    stop(paste0("Unrecognized type \"", type, "\" -- expected \"spatial\" or \"sc\"."))
  }

  return(obj)
}

obj1 <- load_wgcna_obj(obj1_file, type1, name1)
obj2 <- load_wgcna_obj(obj2_file, type2, name2)

# Rename whichever @misc entry is actually active to a fixed "ref"/"test"
# name, rather than assuming a fixed list position (see header note
# above) -- ModulePreservationNetRep() needs both networks under
# caller-chosen names, and the scRNAseq side's real internal wgcna_name
# varies by cell type.

active1 <- obj1@misc$active_wgcna
names(obj1@misc)[names(obj1@misc) == active1] <- "ref"
obj1@misc$active_wgcna <- "ref"

active2 <- obj2@misc$active_wgcna
names(obj2@misc)[names(obj2@misc) == active2] <- "test"
obj2@misc$active_wgcna <- "test"

# 2 is the query object, 1 is the reference: module preservation is
# assessed for the reference's (obj1's) modules as tested against the
# query (obj2's) data.

obj2 <- ModulePreservationNetRep(
  obj2,
  seurat_ref = obj1,
  name = "test",
  n_permutations = 10000,
  TOM_use = "test",
  n_threads = 16,
  wgcna_name = "test",
  wgcna_name_ref = "ref"
)

results <- GetModulePreservation(obj2,
                                 "test",
                                 "test")

write.csv(results$p.values,
          file = paste0(results_dir, "pvals_", name1, "_into_", name2, ".csv"))

write.csv(results$observed,
          file = paste0(results_dir, "observed_", name1, "_into_", name2, ".csv"))

p <- PlotModulePreservationLollipop(
  obj2,
  name = "test",
  features = "average",
  wgcna_name = "test"
) +
  ggtitle(paste0(name1, " projected\nto ", name2)) +
  theme(axis.text = element_text(color = "black"))

ggsave(p,
       filename = paste0(results_dir, "preservation_", name1, "_into_", name2, ".png"),
       units = "in", dpi = 600,
       height = nrow(results$p.values) * 0.7,
       width = 5)
