library(tidyverse)
library(Seurat)
library(scCustomize)

setwd("/projects/b1169/boles/als_motor_circuit_visium/")

data_dir <- "data/02_qc/"
results_dir <- "results/02_qc/"

dir.create(data_dir,
           F, T)
dir.create(results_dir,
           F, T)

obj <- readRDS("data/01_obj_creation/cns_obj.rds")

obj <- Add_Mito_Ribo(obj,
                     species = "Hs")
obj <- Add_Cell_Complexity(obj,
                           assay = "Spatial")

# Set cutoffs ---------------------------------------

samples <- unique(obj$sample_id)
n_samples <- length(samples)

thresh_df <- data.frame(sample_id = samples,
                        umi_med = c(rep(NA, n_samples)),
                        umi_mad = c(rep(NA, n_samples)),
                        feature_med = c(rep(NA, n_samples)),
                        feature_mad = c(rep(NA, n_samples)),
                        max_col_idx = c(rep(NA, n_samples)),
                        max_row_idx = c(rep(NA, n_samples)),
                        min_col_idx = c(rep(NA, n_samples)),
                        min_row_idx = c(rep(NA, n_samples)))

meta <- obj@meta.data

for (i in seq_along(samples)){
  message(paste0("Getting cutoffs for ", samples[i]))
  
  df <- meta %>%
    filter(sample_id == samples[i])
  
  thresh_df$umi_med[i] <- median(log10(df$nCount_Spatial))
  
  thresh_df$umi_mad[i] <- stats::mad(log10(df$nCount_Spatial))
  
  thresh_df$feature_med[i] <- median(log10(df$nFeature_Spatial))
  
  thresh_df$feature_mad[i] <- stats::mad(log10(df$nFeature_Spatial))
  
  thresh_df$max_col_idx[i] <- max(df$array_col)
  
  thresh_df$max_row_idx[i] <- max(df$array_row)
  
  thresh_df$min_col_idx[i] <- min(df$array_col)
  
  thresh_df$min_row_idx[i] <- min(df$array_row)
  
}

# Pretty "strict" cutoffs of 2 x MA
thresh_df <- thresh_df %>%
  mutate(umi_lower = umi_med - 3*umi_mad,
         feature_lower = feature_med - 3*feature_mad,
         mito_upper = 20,
         umi_upper = umi_med + 2*umi_mad,
         feature_upper = feature_med + 2*feature_mad) %>%
  mutate(feature_lower = case_when(10^(feature_lower) < 300 ~ log10(300),
                                   .default = feature_lower),
         umi_lower = case_when(10^(umi_lower) < 300 ~ log10(300),
                               .default = umi_lower))

meta <- meta %>%
  rownames_to_column(var = "cell") %>%
  left_join(thresh_df,
            by = "sample_id")

p <- meta %>%
  ggplot(aes(x = sample_id)) +
  geom_violin(aes(y = nFeature_Spatial)) +
  geom_point(aes(y = 10^(feature_lower))) +
  geom_point(aes(y = 10^(feature_upper))) +
  scale_y_log10() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p,
       filename = paste0(results_dir, "nfeature_violins.png"),
       units = "in", dpi = 600,
       height = 6, width = 15)

# meta %>%
#   group_by(sample_id) %>%
#   summarize(n = n()) %>%
#   ggplot(aes(x = sample_id,
#              y = n)) +
#   geom_col() +
#   theme_linedraw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p <- meta %>%
  ggplot(aes(x = sample_id)) +
  geom_violin(aes(y = nCount_Spatial)) +
  geom_point(aes(y = 10^(umi_lower))) +
  geom_point(aes(y = 10^(umi_upper))) +
  scale_y_log10() +
  theme_linedraw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p,
       filename = paste0(results_dir, "ncount_violins.png"),
       units = "in", dpi = 600,
       height = 6, width = 15)

p <- meta %>% 
  ggplot(aes(x = sample_id)) + 
  geom_violin(aes(y = percent_mito)) +
  geom_hline(yintercept = 20) + 
  theme_linedraw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave(p,
       filename = paste0(results_dir, "mito_violins.png"),
       units = "in", dpi = 600,
       height = 6, width = 15)

meta <- meta %>%
  mutate(discard = if_else(percent_mito > 20 | 
                             log10(nCount_Spatial) < umi_lower | 
                             log10(nCount_Spatial) > umi_upper | 
                             log10(nFeature_Spatial) < feature_lower | 
                             log10(nFeature_Spatial) > feature_upper |
                             array_col == min_col_idx | 
                             array_col == max_col_idx | 
                             array_row == max_row_idx | 
                             array_row == min_row_idx,
                           "discard", "keep"))

meta <- meta %>%
  column_to_rownames("cell")

meta %>% 
  filter(discard == "keep") %>% 
  rownames() %>%
  write.csv(file = paste0(results_dir, "spots_to_keep.csv"),
            row.names = F,
            quote = F)

pre_stats <- meta %>% 
  group_by(sample_id) %>% 
  summarize(raw_med_ncount = median(nCount_Spatial),
            raw_med_nfeature = median(nFeature_Spatial),
            raw_med_mito = median(percent_mito, na.rm = T),
            raw_spots = n())

post_stats <- meta %>% 
  dplyr::filter(discard == "keep") %>% 
  group_by(sample_id) %>% 
  summarize(filtered_med_ncount = median(nCount_Spatial),
            filtered_med_nfeature = median(nFeature_Spatial),
            filtered_med_mito = median(percent_mito, na.rm = T),
            filtered_spots = n())

stats <- full_join(pre_stats, post_stats,
                   by = "sample_id")

p <- stats %>% 
  pivot_longer(!sample_id) %>%
  mutate(stage = str_split_i(name, "_", i = 1),
         stat = str_remove_all(name, "filtered_|raw_")) %>% 
  mutate(stage = factor(stage,
                        levels = c("raw", "filtered"))) %>%
  ggplot(aes(x = stage,
             y = value)) + 
  geom_line(aes(group = sample_id)) + 
  geom_point(aes(fill = sample_id),
             show.legend = F,
             shape = 21,
             color = "black") + 
  facet_wrap(. ~ stat,
             nrow = 2,
             scales = "free_y") + 
  theme_linedraw()
ggsave(p,
       filename = paste0(results_dir, "qc_stats_change.png"),
       units = "in", dpi = 600,
       height = 6, width = 6)  

obj <- AddMetaData(obj,
                   meta)

# obj$log_ncount <- log(obj$nCount_Spatial)
# obj$log_nfeature <- log(obj$nFeature_Spatial)
# 
# SpatialDimPlot(obj,
#                group.by = "discard",
#                images = "X146.7")

filtered_obj <- subset(obj,
                       discard == "keep")

saveRDS(filtered_obj,
        file = paste0(data_dir, "filtered_cns_obj.rds"))

# 
# SpatialFeaturePlot(obj,
#                    features = "log_ncount",
#                    images = "X146.2")

# SpatialFeaturePlot(obj,
#                    features = "nCount_Spatial",
#                    images = "AN72.6")
