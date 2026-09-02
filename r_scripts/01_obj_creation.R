library(tidyverse)
library(Seurat)
library(scCustomize)

setwd("/projects/b1169/boles/als_motor_circuit_visium/")

data_dir <- "data/01_obj_creation/"
dir.create(data_dir,
           recursive = T,
           showWarnings = F)

spaceranger <- list.dirs("/projects/b1042/Gate_Lab/boles/als_motor_circuit_visium/spaceranger",
                         recursive = F)
# spaceranger <- spaceranger[str_detect(spaceranger, "AN68|AN69-7|AN69-8|JSB14", negate = T)]

muscle_objs <- list()
cns_objs <- list()

samples <- str_split_i(spaceranger, "/", i = 8)
# samples <- samples[str_detect(samples, "AN68", negate = T)]
# samples <- samples[str_detect(samples, "AN69-7|AN69-8", negate = T)]
# samples <- samples[str_detect(samples, "JSB14", negate = T)]

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

saveRDS(cns_obj,
        file = paste0(data_dir, "cns_obj.rds"))
saveRDS(muscle_obj,
        file = paste0(data_dir, "muscle_obj.rds"))

