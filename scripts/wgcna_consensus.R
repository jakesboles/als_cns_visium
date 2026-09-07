# Runs consensus hdWGCNA (co-expression modules found consistently across
# all 4 major anatomical compartments -- mcx GM, mcx WM, sc GM, sc WM --
# not fit on one pooled population) on 04_spot_annotation.R's spots,
# reattaching 05_integration.R's already-fit CCA embedding. Single
# run covering all 4 compartments at once (no SLURM array), matching this
# repo's own deseq2_by_compartment.R rather than the scRNAseq sibling
# repo's per-cell-type array-job scripts (wgcna_single.R,
# wgcna_consensus_cns.R) -- there's no separate "target" to loop over
# here the way there's a cell type there, since the goal is one joint
# consensus network across all 4 compartments in a single call.
#
# Modeled closely on als_cns_scrnaseq/r_scripts/wgcna_consensus_cns.R
# (consensus hdWGCNA across brain/spinal cord for shared cell types),
# adapted from "consensus across tissue, held constant per cell type" to
# "consensus across compartment, held constant overall":
# - That script's per-cell-type filtering (meta_sub <- meta[meta$cell_type3
#   == cell_type_target, ]) has no analog here -- there's no cell-type
#   axis to filter down to first, only the region filter
#   (region %in% c("GM", "WM")) needed to restrict to the 4 compartments
#   in the first place.
# - Its `cell_type3` column served two roles: MetacellsByGroups()'s
#   ident.group, and the constant single value passed to SetMultiExpr()/
#   ModuleConnectivity()'s group_name/group.by (there, the specific cell
#   type being analyzed; the tissue split lives entirely in
#   multi.group.by/multi_groups instead). Since there's no cell-type-like
#   axis here, `dummy` (a column that's just the constant 1) fills that
#   role -- same "constant placeholder column used as ident.group/
#   group_name" pattern that script already uses for `cell_type3` once
#   it's constant post-filtering, just with nothing to filter down from in
#   the first place.
# - `compartment` (tissue + region, e.g. "mcx_GM") plays the role
#   `tissue` plays there: the multi.group.by/multi_groups axis
#   SetMultiExpr()/TestSoftPowersConsensus() build separate per-group
#   networks across before finding their consensus.
# - Real raw counts and the CCA embedding come from this project's own
#   04_spot_annotation.R (bpcells_data, metadata.rds) and
#   05_integration.R (cca.rds) -- this project's actual pipeline
#   stages, not the scRNAseq repo's 06/17/18. 05_integration.R found CCA
#   integrated better than Harmony for this data, so that's what's
#   reattached here for metacell construction (MetacellsByGroups()'s
#   `reduction` below) -- but RunHarmonyMetacells() a bit further down is
#   NOT a leftover: hdWGCNA has no CCA equivalent for that specific
#   metacell-level correction step (every WGCNA script in the scRNAseq
#   sibling repo calls it unconditionally too, regardless of what
#   single-cell reduction was used upstream), and nothing downstream
#   (SetMultiExpr() pulls raw expression, not a reduction) actually
#   consumes its output, so it's effectively a required formality of
#   hdWGCNA's own pipeline rather than a real Harmony-vs-CCA choice.
# - This object's assay is "Spatial" throughout (this project's Visium
#   assay name), never "RNA" (the scRNAseq repo's assay name) -- watch
#   for this specifically when porting anything further from that repo.
# - Same ConstructNetwork() TOM.rda workaround as the scRNAseq repo's
#   wgcna scripts (unfixed hdWGCNA bug, smorabit/hdWGCNA#182: it writes a
#   temp .rda file to the working directory regardless of tom_outdir/
#   tom_name). Not strictly needed without a SLURM array here, but kept
#   since it's harmless and isolates the TOM output either way.
# - Added a per-compartment minimum-spot-count check before touching raw
#   counts, matching the abundance-safeguard convention already
#   established in this project's deseq2_by_compartment.R and the
#   scRNAseq repo's wgcna scripts -- a too-sparse compartment (spinal cord
#   WM is the likeliest candidate) would otherwise produce degenerate
#   metacells or fail deep inside TestSoftPowersConsensus()/
#   ConstructNetwork() with a much less clear error.

suppressMessages({
  library(hdWGCNA)
  library(Seurat)
  library(scCustomize)
  library(tidyverse)
  library(patchwork)
  library(UCell)
  library(cowplot)
  library(BPCells)
})

message2 <- function(text){
  v1 <- paste(rep("~", 15),
              collapse = "")
  message(paste0(v1, text, v1))
}

setwd("/projects/b1169/boles/als_cns_visium")

theme_set(theme_cowplot())
set.seed(256)
enableWGCNAThreads(nThreads = 16)

data_dir <- "data/wgcna_consensus/"
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- "results/wgcna_consensus/"
dir.create(results_dir, showWarnings = F, recursive = T)

# Filter to the 4 anatomical compartments before touching raw counts at all -
# See header note above.

message2("Reading in metadata")

meta <- readRDS("data/04_spot_annotation/metadata.rds")

message2("Reading in raw counts and CCA embedding")

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

obj$dummy <- 1
obj$compartment <- paste0(obj$tissue, "_", obj$region)

compartments <- unique(obj$compartment)

# Skip if any compartment is too sparse for stable metacell construction --
# see header note above. Checked per compartment, not just on the combined
# total, since SetMultiExpr() builds a separate metacell population and
# network per compartment before finding the consensus.

min_cells <- 200 # change as needed

cell_counts <- table(obj@meta.data$compartment)
if (any(cell_counts < min_cells)){
  stop(paste0("At least one compartment has too few spots for stable ",
              "metacell construction (min_cells = ", min_cells, "): ",
              paste(names(cell_counts), cell_counts, sep = " = ", collapse = ", ")))
}

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

obj <- ScaleData(obj)

# Identify genes expressed in at least 5% of these spots ---------------

message2("Selecting genes expressed in at least 5% of spots")

pe <- rowMeans(GetAssayData(obj, layer = "data", assay = "Spatial") > 0)
genes_keep <- names(pe)[pe > 0.05] # change this cutoff as needed

# Set up hdWGCNA -------------------------------------------------------
# obj is already filtered to just the 4 compartments -- see header note above.

message2("Setting up hdWGCNA")

obj <- SetupForWGCNA(obj,
                     gene_select = "custom",
                     features = genes_keep,
                     wgcna_name = "wgcna_consensus")

message2("Constructing metacells")

# group.by includes tissue (via compartment) so metacells never mix cells
# across tissue -- see header note above. dummy is kept as ident.group,
# same "constant placeholder column" pattern as wgcna_consensus_cns.R's
# use of cell_type3.
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("code", "compartment", "dummy"),
  reduction = "cca",
  k = 25, # change as needed
  max_shared = 10, # change as needed
  ident.group = "dummy"
)

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj, features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj, features = VariableFeatures(obj))

# Still Harmony despite the rest of this script being CCA -- not a
# leftover, hdWGCNA has no CCA option for this metacell-level correction
# step. See header note above.
obj <- RunHarmonyMetacells(obj, group.by.vars = "code")

obj <- SetMultiExpr(
  obj,
  group_name = 1,
  group.by = "dummy",
  multi.group.by = "compartment",
  multi_groups = compartments,
  assay = "Spatial",
  layer = "data",
  use_metacells = T
)

# Find soft power per compartment ------------------------------------------

message2("Testing soft powers")

obj <- TestSoftPowersConsensus(obj)

plot_list <- PlotSoftPowers(obj)

p_list <- lapply(seq_along(compartments), function(i){
  plot_list[[i]][[1]] +
    ggtitle(paste0("Compartment: ", compartments[i])) +
    theme(plot.title = element_text(hjust = 0.5))
})
p <- wrap_plots(p_list, ncol = 2)
ggsave(p,
       filename = paste0(results_dir, "soft_power.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)

power_table <- GetPowerTable(obj)
write.csv(power_table,
          file = paste0(results_dir, "soft_powers.csv"),
          row.names = F)

# Build consensus TOM and cluster genes into modules -------------------------
# Letting ConstructNetwork() pick the soft power automatically.
#
# ConstructNetwork() writes a temp .rda file to the working directory
# regardless of tom_outdir/tom_name (unfixed hdWGCNA bug,
# smorabit/hdWGCNA#182) -- see header note above. Isolating the working
# directory for just this call is defensive here (no SLURM array), but
# harmless either way.

message2("Constructing consensus network")

tom_dir <- paste0(data_dir, "tom/")
dir.create(tom_dir, showWarnings = F, recursive = T)

setwd(tom_dir)
tryCatch({
  obj <- ConstructNetwork(obj,
                          tom_name = "wgcna_consensus",
                          consensus = T,
                          overwrite_tom = T)
}, finally = {
  setwd("/projects/b1169/boles/als_cns_visium")
})

png(paste0(results_dir, "dendrogram.png"),
    height = 8, width = 8, units = "in", res = 600)
PlotDendrogram(obj, main = "Consensus dendrogram")
dev.off()

# Module stats ------------------------------------------------------------

message2("Computing module eigengenes and connectivity")

obj <- SetActiveWGCNA(obj, "wgcna_consensus")
obj <- ModuleEigengenes(obj, group.by.vars = "code")

# ModuleConnectivity() reaches back into the full single-cell "data" layer
# for its corSparse()-based correlation step, which needs a real
# CsparseMatrix -- BPCells' lazy matrix classes don't support that
# coercion. Deferred to right before this call, not earlier, so gene
# selection, SetupForWGCNA(), and metacell construction above keep the
# benefit of BPCells' lazy/streaming evaluation.
obj[["Spatial"]]$data <- as(obj[["Spatial"]]$data, "dgCMatrix")

obj <- ModuleConnectivity(obj, group_name = 1, group.by = "dummy")

p <- PlotKMEs(obj, ncol = 4, text_size = 4)
ggsave(p,
       filename = paste0(results_dir, "module_connectivity.png"),
       units = "in", dpi = 600,
       height = 12, width = 12)

mods <- obj@misc[["wgcna_consensus"]][["wgcna_modules"]]
write.csv(mods,
          file = paste0(results_dir, "modules.csv"),
          row.names = F)

# Module expression scores via UCell -----------------------------------------

message2("Scoring modules with UCell")

module_names <- setdiff(unique(mods$module), "grey")

gene_sets <- lapply(module_names, function(m){
  mods$gene_name[mods$module == m]
})
names(gene_sets) <- module_names

maxrank <- max(lengths(gene_sets))

obj <- AddModuleScore_UCell(obj, features = gene_sets, maxRank = maxrank)
obj <- SmoothKNN(obj,
                 signature.names = paste0(names(gene_sets), "_UCell"),
                 reduction = "cca")

scores <- obj@meta.data %>%
  dplyr::select(code, group, tissue, region, ptdp, pga, batch, age, sex, matches("_UCell_kNN$"))

write.csv(scores,
          file = paste0(results_dir, "module_scores_ucell.csv"),
          row.names = F)

# Harmonized module eigengenes -----------------------------------------------

message2("Saving module eigengenes")

hMEs <- GetMEs(obj, harmonized = T)
write.csv(hMEs,
          file = paste0(results_dir, "module_eigengenes.csv"))

png(paste0(results_dir, "module_eigengene_correlogram.png"),
    height = 8, width = 8, units = "in", res = 600)
ModuleCorrelogram(obj, features = "MEs")
dev.off()

# Save the hdWGCNA experiment for further downstream use ---------------------
# Just the hdWGCNA network/module state (@misc[[wgcna_name]]), not the
# whole Seurat object -- matches the scRNAseq repo's wgcna scripts (the
# object's own expression matrix is already saved separately by
# 04_spot_annotation.R).

message2("Saving hdWGCNA experiment object")

wgcna_experiment <- obj@misc[["wgcna_consensus"]]
saveRDS(wgcna_experiment,
        file = paste0(data_dir, "wgcna_experiment.rds"))
