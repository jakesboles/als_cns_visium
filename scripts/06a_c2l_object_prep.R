# Prepares the raw counts + metadata inputs 06b_c2l_object_prep.py needs to
# assemble Cell2Location's .h5ad reference/query objects: per-tissue
# scRNAseq references (the sister project's cell-typed reference) and
# per-tissue Visium queries (this project's own annotated spots), written
# as .mtx counts + .csv metadata triplets rather than .h5ad directly --
# R-to-Python AnnData conversion (zellkonverter/reticulate) is a real
# source of version fragility on a shared cluster environment, so the
# handoff goes through plain matrix-market/CSV files instead, which 06b
# reassembles into .h5ad.
#
# Cleaned up from a colleague's now-stale draft (written against a
# different, pre-BPCells Visium pipeline, with absolute paths under
# /gpfs/projects/b1169/thomas/als_multitissue/...): consolidated 4
# near-identical copy-pasted blocks (each re-loading the same libraries
# and re-setwd()ing) into one script with shared setup; dropped
# zellkonverter/SingleCellExperiment (imported in the draft but never
# actually called anywhere -- the real R-to-Python handoff is the
# .mtx/.csv files, not either of those); repointed every path at this
# project's own data_dir and at the real upstream objects (see below).
#
# scRNAseq reference: data/20_subcluster_phenotyping/microglia/
# full_metadata.rds -- despite the "microglia" in the path, this IS the
# correct all-cell-type reference, not a microglia-only subset. Confirmed
# directly against als_cns_scrnaseq/r_scripts/20_subcluster_phenotyping.R
# (lines 254-271): that script starts from 18_full_integration's full,
# all-cell-type metadata, and only overwrites the `cell_type4` column for
# microglia specifically -- every other cell type keeps its coarser
# `cell_type3` label under `cell_type4` too. It's saved under a
# "microglia" subfolder only because that was the subclustering run that
# last touched it, not because the file itself is microglia-only. No
# other cell type has gone through this same per-cell-type phenotyping
# step yet (checked: no other als_cns_scrnaseq r_scripts/*.R reference a
# sibling data/20_subcluster_phenotyping/<other_celltype>/ path), so this
# is currently the single most up-to-date, all-cell-type annotated
# reference available.
#
# Per-tissue split: filters the scRNAseq reference and the Visium query
# to the SAME two tissues (motor cortex/"mcx", cervical spinal cord/"sc")
# so 07's two-task array (one task per tissue) can each train against its
# own tissue-matched reference, rather than one combined reference/query
# pair spanning both tissues.

suppressMessages({
  library(Seurat)
  library(BPCells)
  library(dplyr)
  library(Matrix)
})

setwd("/projects/b1169/boles/als_cns_visium")

data_dir <- "data/06_c2l_object_prep/"
dir.create(data_dir, showWarnings = F, recursive = T)

tissues <- data.frame(
  file = c("mcx", "sc"),
  scrna_tissue = c("Motor cortex", "Cervical spinal cord"),
  visium_tissue = c("mcx", "sc")
)

# scRNAseq references ------------------------------------------------------
# One per tissue, from the sister project's final cell-typed reference --
# see header note above.

message("Reading in scRNAseq reference metadata and counts")

scrna_meta <- readRDS("/projects/b1169/boles/als_cns_scrnaseq/data/20_subcluster_phenotyping/microglia/full_metadata.rds")
scrna_counts <- open_matrix_dir("/projects/b1169/boles/als_cns_scrnaseq/data/06_obj_reassembly/bpcells")

for (i in seq_len(nrow(tissues))){

  file <- tissues$file[i]

  message(paste0("scRNAseq reference: ", file))

  meta_sub <- scrna_meta %>%
    filter(tissue == tissues$scrna_tissue[i])

  counts_sub <- scrna_counts[, rownames(meta_sub)]

  obj <- CreateSeuratObject(counts = counts_sub, meta.data = meta_sub)

  counts_matrix <- as(GetAssayData(obj, layer = "counts"), "dgCMatrix")

  writeMM(counts_matrix, file = paste0(data_dir, file, "_sc_counts.mtx"))
  write.csv(obj@meta.data, file = paste0(data_dir, file, "_sc_cell_metadata.csv"), row.names = T)
  write.csv(data.frame(gene = rownames(counts_matrix)),
            file = paste0(data_dir, file, "_sc_gene_metadata.csv"), row.names = F)
}

# Visium queries ------------------------------------------------------------
# 04_spot_annotation.R's already-filtered/annotated object, split by
# tissue the same way -- see header note above.

message("Reading in Visium metadata, counts, and images")

visium_meta <- readRDS("data/04_spot_annotation/metadata.rds")
visium_counts <- open_matrix_dir("data/04_spot_annotation/bpcells_data")
visium_counts <- visium_counts[, rownames(visium_meta)]
visium_images <- readRDS("data/04_spot_annotation/images.rds")

visium_obj <- CreateSeuratObject(counts = visium_counts, meta.data = visium_meta, assay = "Spatial")
visium_obj@images <- visium_images

for (i in seq_len(nrow(tissues))){

  file <- tissues$file[i]

  message(paste0("Visium query: ", file))

  sub <- visium_obj %>%
    subset(tissue == tissues$visium_tissue[i])

  counts_matrix <- as(GetAssayData(sub, layer = "counts"), "dgCMatrix")

  writeMM(counts_matrix, file = paste0(data_dir, file, "_spatial_counts.mtx"))
  write.csv(sub@meta.data, file = paste0(data_dir, file, "_spatial_cell_metadata.csv"), row.names = T)
  write.csv(data.frame(gene = rownames(counts_matrix)),
            file = paste0(data_dir, file, "_spatial_gene_metadata.csv"), row.names = F)
}

message("Done")
