suppressMessages({
  library(tidyverse) 
  library(Seurat)
  library(readxl)
  library(scCustomize)
  library(BPCells)
})

setwd("/projects/b1169/boles/als_cns_visium")

in_dir <- "data/02_qc/"

results_dir <- "results/05_integration/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

data_dir <- "data/05_integration/"
dir.create(data_dir,
           showWarnings = F,
           recursive = T)

counts <- open_matrix_dir(paste0(in_dir, "bpcells_cns"))
meta <- readRDS("data/04_spot_annotation/metadata.rds")
images <- readRDS(paste0(in_dir, "cns_images.rds"))

counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
obj@images <- images

obj <- NormalizeData(obj) %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA()

p <- ElbowPlot(obj, ndims = 50)
ggsave(p,
       filename = paste0(results_dir, "pca_elbow.png"),
       units = "in", dpi = 600, bg = "white",
       height = 6, width = 6)

Iterate_PC_Loading_Plots(obj,
                         file_path = results_dir,
                         file_name = "pca_loadings")

obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj@meta.data$orig.ident)

obj <- IntegrateLayers(obj, 
                     method = HarmonyIntegration, 
                     assay = "Spatial", 
                     layers = "data", 
                     orig.reduction = "pca", 
                     new.reduction = "harmony",
                     dims = 1:15)

obj[["Spatial"]] <- JoinLayers(obj[["Spatial"]])

obj <- RunUMAP(obj,
               umap.method = "uwot",
               reduction = "harmony",
               dims = 1:15,
               # nn.name = "RNA.nn",
               metric = "euclidean",
               min.dist = 0.5,
               n.neighbors = 15L,
               # repulsion.strength = 0.5,
               # uwot.init = "random",
               reduction.name = "harmony_umap",
               return.model = F)

for (group in c("region", "tissue", "orig.ident", "group", "ptdp")){
  w <- if (group %in% c("orig.ident")) 15 else 11
  
  p <- DimPlot_scCustom(obj,
                        reduction = "harmony_umap",
                        group.by = group)
  ggsave(p,
         filename = paste0(results_dir, group, "_dimplot.png"),
         units = "in", dpi = 600,
         height = 8, width = w)
}

bpcells_data_dir <- paste0(data_dir, "bpcells_data")
if (dir.exists(bpcells_data_dir)){
  unlink(bpcells_data_dir, recursive = T)
}

write_matrix_dir(mat = obj[["RNA"]]$data,
                 dir = bpcells_data_dir)

saveRDS(obj[["harmony"]],
        file = paste0(data_dir, "harmony.rds"))

saveRDS(obj[["harmony_umap"]],
        file = paste0(data_dir, "harmony_umap.rds"))
