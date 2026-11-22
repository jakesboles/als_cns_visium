library(Seurat)
library(BPCells)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)
setwd("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects")
# Brain reference ---------------------------------------------------------

scrna_meta <- readRDS("/projects/b1169/boles/als_cns_scrnaseq/data/20_subcluster_phenotyping/microglia/full_metadata.rds")
scrna_counts <- open_matrix_dir("/projects/b1169/boles/als_cns_scrnaseq/data/06_obj_reassembly/bpcells")

brain_scrna_meta <- scrna_meta %>%
  filter(tissue == "Motor cortex")

brain_scrna_counts <- scrna_counts[, rownames(brain_scrna_meta)]

brain_scrna_obj <- CreateSeuratObject(counts = brain_scrna_counts,
                                      meta.data = brain_scrna_meta)

counts_matrix <- as(GetAssayData(brain_scrna_obj, layer = "counts"), "dgCMatrix")

writeMM(counts_matrix, file = "mcx_sc_counts.mtx")

write.csv(brain_scrna_obj@meta.data, file = "mcx_sc_cell_metadata.csv", row.names = TRUE)

write.csv(data.frame(gene = rownames(counts_matrix)), file = "mcx_sc_gene_metadata.csv", row.names = FALSE)

# DONE

# SC reference ------------------------------------------------------------

library(Seurat)
library(BPCells)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)
setwd("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects")

scrna_meta <- readRDS("/projects/b1169/boles/als_cns_scrnaseq/data/20_subcluster_phenotyping/microglia/full_metadata.rds")
scrna_counts <- open_matrix_dir("/projects/b1169/boles/als_cns_scrnaseq/data/06_obj_reassembly/bpcells")

sc_scrna_meta <- scrna_meta %>%
  filter(tissue == "Cervical spinal cord")

sc_scrna_counts <- scrna_counts[, rownames(sc_scrna_meta)]

sc_scrna_obj <- CreateSeuratObject(counts = sc_scrna_counts,
                                   meta.data = sc_scrna_meta)

counts_matrix <- as(GetAssayData(sc_scrna_obj, layer = "counts"), "dgCMatrix")

writeMM(counts_matrix, file = "sc_sc_counts.mtx")

write.csv(sc_scrna_obj@meta.data, file = "sc_sc_cell_metadata.csv", row.names = TRUE)

write.csv(data.frame(gene = rownames(counts_matrix)), file = "sc_sc_gene_metadata.csv", row.names = FALSE)

# DONE

# Brain Visium ------------------------------------------------------------

library(Seurat)
library(BPCells)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)
setwd("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects")

visium_counts <- open_matrix_dir("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/bpcells_data")
visium_meta <- readRDS("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/metadata.rds")
visium_images <- readRDS("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/images.rds")

visium_counts <- visium_counts[, rownames(visium_meta)]

visium_obj <- CreateSeuratObject(counts = visium_counts, meta.data = visium_meta, assay = "Spatial")
visium_obj@images <- visium_images

brain_visium_obj <- visium_obj %>%
  subset(tissue == "mcx")

counts_matrix <- as(GetAssayData(brain_visium_obj, layer = "counts"), "dgCMatrix")

writeMM(counts_matrix, file = "mcx_spatial_counts.mtx")

write.csv(brain_visium_obj@meta.data, file = "mcx_spatial_cell_metadata.csv", row.names = TRUE)

write.csv(data.frame(gene = rownames(counts_matrix)), file = "mcx_spatial_gene_metadata.csv", row.names = FALSE)

# DONE

# SC Visium ---------------------------------------------------------------

library(Seurat)
library(BPCells)
library(zellkonverter)
library(SingleCellExperiment)
library(dplyr)
library(Matrix)
setwd("/gpfs/projects/b1169/thomas/als_multitissue/Visium/Cell2Location2/CreateObjects")

visium_counts <- open_matrix_dir("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/bpcells_data")
visium_meta <- readRDS("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/metadata.rds")
visium_images <- readRDS("/projects/b1169/boles/als_cns_visium/data/04_spot_annotation/images.rds")

visium_counts <- visium_counts[, rownames(visium_meta)]

visium_obj <- CreateSeuratObject(counts = visium_counts, meta.data = visium_meta, assay = "Spatial")
visium_obj@images <- visium_images

sc_visium_obj <- visium_obj %>%
  subset(tissue == "sc")

counts_matrix <- as(GetAssayData(sc_visium_obj, layer = "counts"), "dgCMatrix")

writeMM(counts_matrix, file = "sc_spatial_counts.mtx")

write.csv(sc_visium_obj@meta.data, file = "sc_spatial_cell_metadata.csv", row.names = TRUE)

write.csv(data.frame(gene = rownames(counts_matrix)), file = "sc_spatial_gene_metadata.csv", row.names = FALSE)