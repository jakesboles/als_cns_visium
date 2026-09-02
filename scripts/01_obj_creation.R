library(tidyverse)
library(Seurat)
library(scCustomize)
library(BPCells)

setwd("/projects/b1169/boles/als_motor_circuit_visium/")

data_dir <- "data/01_obj_creation/"
dir.create(data_dir,
           recursive = T,
           showWarnings = F)

# Each sample's counts matrix is written to its own on-disk BPCells directory
# right after it's loaded, so the merge below never holds every sample's
# in-memory sparse matrix at once.
bpcells_persample_dir <- paste0(data_dir, "bpcells_persample/")
dir.create(bpcells_persample_dir,
           recursive = T,
           showWarnings = F)

spaceranger <- list.dirs("/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/spaceranger",
                         recursive = F)

samples <- str_split_i(spaceranger, "/", i = 8)

# Excluded per current cohort definition
exclude <- c("137-1", "137-2", "AN16-1", "AN68-1", "AN68-2", "AN69-7", "AN69-8")
keep <- !(samples %in% exclude)
spaceranger <- spaceranger[keep]
samples <- samples[keep]

muscle_objs <- list()
cns_objs <- list()

muscle <- c("JSB147-7", "JSB147-8", paste0("AN69-", c(1:6)))

for (i in seq_along(spaceranger)){

  dir <- paste0(spaceranger[i], "/outs")

  print(samples[i])

  s <- Load10X_Spatial(dir,
                       slice = samples[i],
                       assay = "Spatial")

  spatial_df <- read.csv(paste0(dir, "/spatial/tissue_positions.csv"),
                         row.names = 1)
  spatial_df <- spatial_df[rownames(s@meta.data),]

  print(sum(rownames(spatial_df) != rownames(s@meta.data)))
  print(sum(is.na(spatial_df)))
  s <- AddMetaData(s, spatial_df)

  s <- RenameCells(object = s,
                   add.cell.id = samples[i])

  s$sample_id <- samples[i]

  # dgCMatrix always stores its @x slot as R type "double" regardless of the
  # values it holds, which leaves write_matrix_dir()'s compressed writer with
  # an ambiguously-typed matrix. Converting to an explicit integer type first
  # avoids that.
  counts <- convert_matrix_type(s[["Spatial"]]$counts, type = "uint32_t")
  bp_dir <- paste0(bpcells_persample_dir, samples[i])
  write_matrix_dir(mat = counts, dir = bp_dir)
  s[["Spatial"]]$counts <- open_matrix_dir(dir = bp_dir)

  if (samples[i] %in% muscle) {
    muscle_objs[i] <- s
  } else {
    cns_objs[i] <- s
  }

}

cns_objs <- compact(cns_objs)
muscle_objs <- compact(muscle_objs)

cns_obj <- Merge_Seurat_List(cns_objs)
muscle_obj <- Merge_Seurat_List(muscle_objs)

cns_obj <- JoinLayers(cns_obj)
muscle_obj <- JoinLayers(muscle_obj)

message("Saving CNS counts matrix as BPCells on-disk matrix")
counts_out <- convert_matrix_type(cns_obj[["Spatial"]]$counts, type = "uint32_t")
write_matrix_dir(mat = counts_out,
                 dir = paste0(data_dir, "bpcells_cns"))

message("Saving muscle counts matrix as BPCells on-disk matrix")
counts_out <- convert_matrix_type(muscle_obj[["Spatial"]]$counts, type = "uint32_t")
write_matrix_dir(mat = counts_out,
                 dir = paste0(data_dir, "bpcells_muscle"))

message("Saving metadata and images as RDS")
saveRDS(cns_obj@meta.data,
        file = paste0(data_dir, "cns_metadata.rds"))
saveRDS(muscle_obj@meta.data,
        file = paste0(data_dir, "muscle_metadata.rds"))

saveRDS(cns_obj@images,
        file = paste0(data_dir, "cns_images.rds"))
saveRDS(muscle_obj@images,
        file = paste0(data_dir, "muscle_images.rds"))

# Downstream scripts should reconstruct each object from these on-disk pieces:
#   counts <- open_matrix_dir(paste0(data_dir, "bpcells_cns"))
#   meta <- readRDS(paste0(data_dir, "cns_metadata.rds"))
#   images <- readRDS(paste0(data_dir, "cns_images.rds"))
#   obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
#   obj@images <- images
# (swap in the muscle_* files for the muscle object)
