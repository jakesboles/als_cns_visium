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

# Figure out which cell type this task handles ------------------------

data_dir <- "data/wgcna_consensus/"
dir.create(data_dir, showWarnings = F, recursive = T)

results_dir <- "results/wgcna_consensus/"
dir.create(results_dir, showWarnings = F, recursive = T)

# Filter to this cell type before touching raw counts at all ----------------
# See header note above.

message2("Reading in metadata and filtering to this cell type")

meta <- readRDS("data/04_spot_annotation/metadata.rds")
meta_sub <- meta %>% 
  filter(region %in% c("GM", "WM")) %>%
  mutate(compartment = paste0(tissue, "_", region),
         dummy = 1)

compartments <- unique(meta_sub$compartment)

message2("Reading in raw counts and Harmony embedding")

raw_mat <- open_matrix_dir("data/04_spot_annotation/bpcells_data")
raw_mat <- raw_mat[, rownames(meta_sub)]

images <- readRDS("data/04_spot_annotation/images.rds")

obj <- CreateSeuratObject(counts = raw_mat, meta.data = meta_sub, assay = "Spatial")
obj@images <- images

obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj)

harmony <- readRDS("data/05_integration/harmony.rds")
harmony@cell.embeddings <- harmony@cell.embeddings[rownames(meta_sub), ]
obj[["harmony"]] <- harmony

obj <- ScaleData(obj)

# Identify genes expressed in at least 5% of this cell type's cells ---------

message2("Selecting genes expressed in at least 5% of this cell type")

pe <- rowMeans(GetAssayData(obj, layer = "data", assay = "RNA") > 0)
genes_keep <- names(pe)[pe > 0.05] # change this cutoff as needed

# Set up hdWGCNA -------------------------------------------------------
# obj is already filtered to just this cell type -- see header note above.

message2("Setting up hdWGCNA")

obj <- SetupForWGCNA(obj,
                     gene_select = "custom",
                     features = genes_keep,
                     wgcna_name = "wgcna_consensus")

message2("Constructing metacells")

# group.by includes tissue (unlike wgcna_single.R) so metacells never mix
# cells across tissue -- see header note above. cell_type3 is dropped
# from group.by (constant post-filtering, same simplification as the
# single script) but kept as ident.group.
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("code", "compartment", "dummy"),
  reduction = "harmony",
  k = 25, # change as needed
  max_shared = 10, # change as needed
  ident.group = "dummy"
)

obj <- NormalizeMetacells(obj)
obj <- ScaleMetacells(obj, features = VariableFeatures(obj))
obj <- RunPCAMetacells(obj, features = VariableFeatures(obj))
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

# Find soft power per tissue -----------------------------------------------

message2("Testing soft powers")

obj <- TestSoftPowersConsensus(obj)

plot_list <- PlotSoftPowers(obj)

p_list <- lapply(seq_along(compartments), function(i){
  plot_list[[i]][[1]] +
    ggtitle(paste0("Tissue: ", compartments[i])) +
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
# Letting ConstructNetwork() pick the soft power automatically, matching
# the draft.
#
# Same TOM.rda SLURM-array collision as wgcna_single.R -- see that
# script's header for the full diagnosis (unfixed hdWGCNA bug,
# smorabit/hdWGCNA#182). Same workaround: isolate the working directory
# for just this call.

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

# ModuleConnectivity() reaches back into the full single-cell "data"
# layer for its corSparse()-based correlation step -- see wgcna_single.R's
# header for the full BPCells/CsparseMatrix diagnosis. Deferred to right
# before this call, not earlier, for the same reason.
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
# Same approach as wgcna_single.R -- see header note above.

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
                 reduction = "harmony")

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
# whole Seurat object -- see header note above.

message2("Saving hdWGCNA experiment object")

wgcna_experiment <- obj@misc[["wgcna_consensus"]]
saveRDS(wgcna_experiment,
        file = paste0(data_dir, "wgcna_experiment.rds"))
