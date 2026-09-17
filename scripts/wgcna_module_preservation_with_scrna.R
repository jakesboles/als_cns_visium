library(hdWGCNA)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)
library(igraph)
library(NetRep)

enableWGCNAThreads(16)

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- paste0("results/wgcna_cross_modality/")
dir.create(results_dir,
           showWarnings = F,
           recursive = T)
# stopped editing here

obj1_file <- commandArgs(trailingOnly = TRUE)[1]
obj2_file <- commandArgs(trailingOnly = TRUE)[2]
name1 <- commandArgs(trailingOnly = T)[3]
name2 <- commandArgs(trailingOnly = T)[4]
type1 <- commandArgs(trailingOnly = T)[5]
type2 <- commandArgs(trailingOnly = T)[6]

obj1 <- readRDS(obj1_file)
obj2 <- readRDS(obj2_file)

if (type1 == "sc") { 
  obj1 <- SetDatExpr(obj1,
                     group_name = name1,
                     group.by = "PredictedCellType",
                     use_metacells = T)
} else { 
    obj1 <- SetDatExpr(obj1,
                       group_name = "all_cells",
                       group.by = "all_cells_group",
                       use_metacells = T)
}

if (type2 == "sc") { 
  obj2 <- SetDatExpr(obj2,
                     group_name = name2,
                     group.by = "PredictedCellType",
                     use_metacells = T)
} else { 
  obj2 <- SetDatExpr(obj2,
                     group_name = "all_cells",
                     group.by = "all_cells_group",
                     use_metacells = T)
}

# cell1 <- str_split_i(param1, "_", i = 1)
# tissue1 <- str_split_i(param1, "_", i = 2)

# cell2 <- str_split_i(param2, "_", i = 1)
# tissue2 <- str_split_i(param2, "_", i = 2)

# 2 is the query object, 1 is the reference
# 1 --> 2 

# if (i == 1) { 
# obj1 <- readRDS(paste0("25_", tissue1, "/", cell1, "/wgcna_obj.rds"))
# } else {
# if(param1 != params$V1[i-1]) { 
# obj1 <- readRDS(paste0("25_", tissue1, "/", cell1, "/wgcna_obj.rds"))
# }
# }

# if (i == 1) { 
# obj2 <- readRDS(paste0("25_", tissue2, "/", cell2, "/wgcna_obj.rds"))
# } else {
# if(param1 != params$V2[i-1]) { 
# obj2 <- readRDS(paste0("25_", tissue2, "/", cell2, "/wgcna_obj.rds"))
# }
# }

# module_df <- read.csv(paste0("25_", tissue1, "/", cell1, "/csvs/Modules.csv"))

# obj2 <- ProjectModules(obj2,
#                        # modules = module_df,
#                        seurat_ref = obj1,
#                        wgcna = cell1,
#                        wgcna_name_proj = cell2,
#                        assay = "RNA")

# set expression matrix for reference dataset

# obj1 <- SetDatExpr(
#   obj1,
#   group_name = cell1,
#   group.by = "PredictedCellType",
#   use_metacells = T
# )

# set expression matrix for query dataset:
# obj2 <- SetDatExpr(
#   obj2,
#   group_name = cell2,
#   group.by = "PredictedCellType",
#   use_metacells = T
# )

names(obj1@misc)[2] <- "ref"
obj1@misc$active_wgcna <- "ref"

names(obj2@misc)[2] <- "test"
obj2@misc$active_wgcna <- "test"

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
results$p.values

write.csv(results$p.values,
          file = paste0(csv_dir, "pvals_", name1, "_into_", name2, ".csv"))

write.csv(results$observed,
          file = paste0(csv_dir, "observed_", name1, "_into_", name2, ".csv"))

p <- PlotModulePreservationLollipop(
  obj2, 
  name='test',
  features='average',
  wgcna_name="test"
) + 
  ggtitle(paste0(name1, " projected\nto ", name2)) + 
  theme(axis.text = element_text(color = "black"))

ggsave(p, 
       filename = paste0(plots_dir, "preservation_", name1, "_into_", name2, ".png"),
       units = "in", dpi = 600,
       height = nrow(results$p.values) * 0.7,
       width = 5)
# }