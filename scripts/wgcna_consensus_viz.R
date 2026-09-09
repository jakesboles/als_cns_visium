# Plots of wgcna_consensus.R's output, modeled on
# als_cns_scrnaseq/r_scripts/wgcna_single_params.R-style per-module figures
# (kME bar chart, spot-level violin, pseudobulk quasirandom) pasted in from
# that repo as a starting point. Shares results/wgcna_consensus/ with
# wgcna_consensus.R and wgcna_consensus_analysis.R rather than getting its
# own folder, matching this project's established downstream-script
# convention.
#
# Adaptations from the pasted sibling-repo example:
# - No `celltype`/per-cell-type subfolder -- wgcna_consensus.R runs one
#   joint consensus network across all 4 anatomical compartments in a
#   single call (not one run per cell type), so there's no equivalent axis
#   to loop over or nest the results under here.
# - The sibling script's tissue-only facet (Motor cortex / Cervical spinal
#   cord) is replaced with `compartment` (mcx_GM, mcx_WM, sc_GM, sc_WM) --
#   this project's WGCNA run splits on compartment, not just tissue, so
#   that's the axis its module scores actually vary across. Built the same
#   way wgcna_consensus_analysis.R already builds it
#   (paste0(tissue, "_", region)), since module_scores_ucell.csv only
#   carries tissue/region separately.
# - Pseudobulk ("pb") is summarized per donor (`sample`), not per
#   `orig.ident` -- this project has no orig.ident-as-donor+tissue-key
#   column (see CLAUDE.md's assay/column-naming warnings); `sample` is
#   this project's actual donor id and is what every other pseudobulk
#   script here (deseq2_by_compartment.R, deseq2_by_pathology.R) groups
#   by, so module scores should be summarized the same way. Grouped by
#   (sample, compartment) together, not sample alone, since one donor can
#   contribute spots to more than one compartment.
# - Output file names: the sibling script's "_expression_sc.png" (spot/
#   cell-resolution violin, i.e. "single-cell") and "_expression_pb.png"
#   (pseudobulk) are renamed to "_expression_spot.png"/
#   "_expression_pseudobulk.png" -- "_sc" would be read as "spinal cord"
#   in this project's own vocabulary (tissue == "sc"), which is exactly
#   backwards from what it means in the sibling script.
# - Assay is "Spatial" throughout, never "RNA" -- this project's Visium
#   objects don't have an "RNA" assay (see CLAUDE.md's repeated warning
#   about this exact copy-paste mistake).
# - The Seurat object is reconstructed using wgcna_consensus.R's own
#   recipe (04_spot_annotation.R's counts/metadata/images +
#   05_integration.R's CCA embedding, then the same
#   region %in% c("GM", "WM") filter) -- this is the actual object
#   hdWGCNA/UCell were run on, not 04's full, unfiltered CNS object (which
#   still has Meninges/Nerve bundle spots that were never part of the
#   consensus network). The sibling script's FeaturePlot_scCustom() (UMAP-
#   embedding feature plot, scRNAseq) is replaced with
#   SpatialFeaturePlot_scCustom() -- this data has physical spot
#   coordinates, so a UMAP feature plot isn't the right visualization for
#   it; module scores are reattached from module_scores_ucell.csv (already
#   computed by wgcna_consensus.R on this exact object) rather than
#   recomputing them with a second AddModuleScore_UCell()/SmoothKNN() call.
# - The sibling script's commented-out, unfinished "ORA result plots"
#   section is dropped rather than duplicated -- wgcna_consensus_analysis.R
#   already produces real ORA dotplots per module
#   (results/wgcna_consensus/<module>_ora.png), so there's nothing left
#   for this script to do there.
#
# NOTE: this session has no cluster access, real data, or R interpreter --
# every line below is a careful read/reasoning-based adaptation of the
# pasted sibling script and wgcna_consensus.R's own save/reconstruction
# code, not something that's been executed. Sanity-check on the cluster,
# especially the SpatialFeaturePlot_scCustom()/scCustomize argument names
# (image.alpha, colors_use, multi-image `images =`) at the end, before
# treating any of this as final.

library(tidyverse)
library(ggplot2)
library(Seurat)
library(scCustomize)
library(paletteer)
library(BPCells)
library(ggbeeswarm)

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/wgcna_consensus/"
data_dir <- "data/wgcna_consensus/"

scores <- read.csv(paste0(results_dir, "module_scores_ucell.csv"))

scores <- scores %>%
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")),
         compartment = paste0(tissue, "_", region),
         compartment = factor(compartment, 
                              levels = sort(unique(compartment)),
                              labels = sort(unique(compartment)) %>% str_replace_all("_", " ") %>% str_to_upper()))

# One row per donor per compartment -- a donor with spots in more than one
# compartment (e.g. both mcx and sc) still contributes separately to each,
# matching how every other pseudobulk script in this project (
# deseq2_by_compartment.R, deseq2_by_pathology.R) treats compartment
# membership.
pb <- scores %>%
  group_by(sample, group, compartment) %>%
  mutate(across(where(is.numeric), median)) %>%
  distinct(sample, compartment, .keep_all = T)

modules <- read.csv(paste0(results_dir, "modules.csv"))

mois <- setdiff(unique(modules$color), "grey") # change to a specific subset once real module colors are known

# kME plots -----------------------------------------------------------------

for (i in seq_along(mois)){

modules %>%
  filter(color == mois[i]) %>%
  arrange(desc(!!sym(paste0("kME_", mois[i])))) %>%
  mutate(gene_name = fct_inorder(gene_name)) %>%
  slice_head(n = 30) %>%
  ggplot(aes(x = !!sym(paste0("kME_", mois[i])),
             y = gene_name)) +
  geom_col(color = "black",
           fill = mois[i]) +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(x = "kME") +
  theme_linedraw(base_size = 12) +
  theme(axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(results_dir, mois[i], "_kme_bars.png"),
       units = "in", dpi = 600,
       height = 6, width = 2.5)

# Expression plots (spot-level) ----------------------------------------------

scores %>%
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) +
  geom_violin(aes(fill = group)) +
  facet_wrap(. ~ compartment,
             nrow = 1) +
  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) +
  labs(y = paste0(str_to_title(mois[i]), " module score")) +
  theme_linedraw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(results_dir, mois[i], "_expression_spot.png"),
       units = "in", dpi = 600,
       height = 3, width = 7)

# Expression plots (pseudobulk) ----------------------------------------------

# fit <- lmer(green_UCell_kNN ~ group * compartment + (1|sample),
#           data = pb)
# 
# joint_tests(fit)

pb %>%
  ggplot(aes(x = group,
             y = !!sym(paste0(mois[i], "_UCell_kNN")))) +
  geom_quasirandom(aes(fill = group),
                   shape = 21,
                   size = 4,
                   alpha = 0.7) +
  stat_summary(fun = mean,
               geom = "crossbar") +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               linewidth = 1.2,
               width = 0.6) +
  facet_wrap(. ~ compartment,
             nrow = 1) +
  scale_fill_manual(values = c("#b8b0a8", "#0CAA00", "#CC00FF")) +
  labs(y = paste0(str_to_title(mois[i]), " module score")) +
  theme_linedraw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))
ggsave(filename = paste0(results_dir, mois[i], "_expression_pseudobulk.png"),
       units = "in", dpi = 600,
       height = 3, width = 7)
}

# Make the Seurat object consensus WGCNA was actually run on -----------------
# Exact recipe from wgcna_consensus.R itself: 04_spot_annotation.R's counts/
# metadata/images, 05_integration.R's CCA embedding subset to the same
# barcodes, then the same region %in% c("GM", "WM") filter -- NOT
# 04's full, unfiltered object, which still includes Meninges/Nerve bundle
# spots that were never part of this consensus network.

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

obj$compartment <- paste0(obj$tissue, "_", obj$region)

obj <- NormalizeData(obj)

obj <- AddMetaData(obj,
                   scores)

# Spatial FeaturePlots --------------------------------------------------------
# Two representative sections -- the same AN67-7 (motor cortex)/AN72-4
# (spinal cord) examples dataset_overview_figs.R already uses, reused here
# rather than picking a new pair so the same two sections recur across this
# project's example figures. ("." not "-" in the image slot name: Seurat
# sanitizes image names containing "-" on load.)
#
# One example module (the first entry in `mois`) is plotted spatially below
# -- both its UCell module score and its single strongest hub gene (the
# top row of the kME bar chart above, i.e. the highest-kME gene). Swap in
# any other module from `mois`/`modules` or any gene name to explore
# further; SpatialFeaturePlot_scCustom() is scCustomize's spatial analog of
# the sibling script's FeaturePlot_scCustom() (which plots onto a UMAP
# embedding -- not meaningful for spot data with real physical
# coordinates).

example_images <- c("JSB146.4", "AN72.8", "AN67.7", "AN72.4")

example_module <- "yellow"

SpatialFeaturePlot(obj,
                   images = example_images,
                   features = paste0(example_module, "_UCell_kNN"),
                   # colors_use = viridis_inferno_light_high,
                   image.alpha = 0,
                   pt.size = 3,
                   ncol = 2)
ggsave(filename = paste0(results_dir, "yellow_spatial_expression.png"),
       units = "in", dpi = 600,
       height = 8, width = 8)
