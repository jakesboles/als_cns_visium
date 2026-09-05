suppressMessages({
  library(tidyverse) 
  library(Seurat)
  library(readxl)
  library(scCustomize)
  library(BPCells)
})

setwd("/projects/b1169/boles/als_cns_visium")

in_dir <- "data/02_qc/"

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

obj[["Spatial"]] <- split(obj[["Spatial"]], f = obj@meta.data$orig.ident)

s <- IntegrateLayers(s, 
                     method = HarmonyIntegration, 
                     assay = "Spatial", 
                     layers = "data", 
                     orig.reduction = "pca", 
                     new.reduction = "harmony")

# saveRDS(s, "Harmony.rds")

s <- FindNeighbors(s, reduction = "harmony", dims = 1:30)

# s <- FindClusters(s, resolution = 0.5, algorithm = 4)

s <- RunUMAP(s, reduction = "harmony", reduction.name = "harmony_umap", dims = 1:30)

# pdf("Harmony_UMAP_CNS.pdf", width = 12, height = 12)
DimPlot(s, group.by = c("tissue"), reduction = "harmony_umap", raster = T, shuffle = T)
# dev.off()

FeaturePlot_scCustom(s,
                     features = "MOBP",
                     raster = T,
                     reduction = "harmony_umap")

saveRDS(s, "Harmony_CNS.rds")





meta <- s@meta.data %>%
  dplyr::select(c(sample_id, tissue)) %>%
  distinct()

meta$Counts <- meta$Features <- 0

for(sample in meta$sample_id){
  
  meta$Counts[meta$sample_id == sample] <- median(s@meta.data$nCount_Spatial[s@meta.data$sample_id == sample])
  
  meta$Features[meta$sample_id == sample] <- median(s@meta.data$nFeature_Spatial[s@meta.data$sample_id == sample])
  
}

meta <- meta %>% 
  pivot_longer(names_to = "Var", values_to = "Value", cols = c(3:4))

ggplot(meta, aes(x = sample_id, y = Value, fill = Var, color = tissue)) +
  geom_col(position = position_dodge2(1)) + scale_fill_manual(values = c("dodgerblue", "firebrick1")) +
  scale_color_manual(values = c("black", "green")) + 
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


meta_sorted <- meta %>%
  # 1. Get the 'Counts' value for each sample to use as our sorting key
  group_by(sample_id) %>%
  mutate(count_key = Value[Var == "Counts"]) %>%
  ungroup() %>%
  # 2. Arrange by tissue first, then by the Count key in descending order
  arrange(tissue, desc(count_key)) %>%
  # 3. Lock this specific order into the sample_id factor levels
  mutate(sample_id = fct_inorder(sample_id))

# Now use 'meta_sorted' in your ggplot code
ggplot(meta_sorted, aes(x = sample_id, y = Value, fill = Var, color = tissue)) +
  geom_col(position = position_dodge2(1)) + 
  scale_fill_manual(values = c("dodgerblue", "firebrick1")) +
  scale_color_manual(values = c("black", "green")) + 
  theme_cowplot() + ggtitle("Median") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))




