library(tidyverse)
library(ggplot2)
library(paletteer)

setwd("/projects/b1169/boles/als_cns_visium")

in_dir <- "results/deseq2_by_pathology"

ptdp_dirs <- list.dirs(in_dir, recursive = F)
ptdp_dirs <- ptdp_dirs[str_detect(ptdp_dirs, "ptdp_")]

for (i in ptdp_dirs){
  message(str_split_i(i, "/", i = 3))

  file <- list.files(i, pattern = ".csv")
  file <- file[str_detect(file, "shrunk|filtering", negate = T)]
  
  df <- read.csv(paste0(i, "/", file))
  
  print(table(df$padj < 0.05))
}

p_thresh <- 0.05
lfc_thresh <- 0
  

c9 <- read.csv(paste0(in_dir, "/ptdp_sc_C9orf72/C9orf72.csv")) %>%
  dplyr::select(X, log2FoldChange, padj)

sals <- read.csv(paste0(in_dir, "/ptdp_sc_sALS/sALS.csv")) %>% 
  dplyr::select(X, log2FoldChange, padj)

df <- inner_join(c9, sals, by = "X", suffix = c("_c9", "_sals"))

df <- df %>%
  mutate(sig_c9 = !is.na(padj_c9) & padj_c9 < p_thresh & abs(log2FoldChange_c9) > lfc_thresh,
         sig_sals = !is.na(padj_sals) & padj_sals < p_thresh & abs(log2FoldChange_sals) > lfc_thresh,
         sig_group = case_when(
           sig_c9 & sig_sals ~ "Both",
           sig_c9 & !sig_sals ~ "C9orf72-ALS",
           !sig_c9 & sig_sals ~ "sALS",
           TRUE ~ "Neither"
         )) %>%
  filter(sig_group != "Neither")

# Ranked by Euclidean distance from the origin, not by either axis alone
# -- that rewards a gene for being an outlier on EITHER comparison (or
# both), rather than always favoring whichever comparison happens to have
# larger fold changes overall.
label_df <- df %>%
  mutate(dist = sqrt(log2FoldChange_c9^2 + log2FoldChange_sals^2)) %>%
  slice_max(dist, n = 30, with_ties = F)

ggplot(df, aes(x = log2FoldChange_c9, y = log2FoldChange_sals, color = sig_group)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_label_repel(data = label_df,
                  aes(label = X, color = sig_group),
                  size = 2, show.legend = F, max.overlaps = Inf) +
  labs(x = "log2(fold change) in C9orf72-ALS",
       y = "log2(fold change) in sALS",
       color = "DEGs in pTDP+ vs\npTDP- spots in:") +
  scale_color_manual(values = c("darkslategrey", "#CC00FF", "#0CAA00")) +
  ggtitle("DEGs in spinal cord pTDP+ spots") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme_linedraw(base_size = 12) +
  theme(axis.text = element_text(color = "black"))
ggsave(filename = paste0(in_dir, "/sc_ptdp_multivolcano.png"),
       units = "in", dpi = 600,
       height = 4, width = 5)

# pGA volcano plot -------------------------------------------------------

pga <- read.csv(paste0(in_dir, "/pga_mcx_C9orf72/C9orf72.csv"))

table(pga$padj < 0.05)

pga <- pga %>% 
  mutate(sig = if_else(padj < 0.05, "sig", "ns") %>% 
           factor(levels = c("sig", "ns"),
                  labels = c("True", "False")))

label_df <- pga %>%
  mutate(dist = sqrt(log2FoldChange^2 + (-log10(padj))^2)) %>%
  slice_max(dist, n = 30, with_ties = F)


pga %>% 
  ggplot(aes(x = log2FoldChange,
             y = -log10(padj))) + 
  geom_point(aes(color = sig),
             size = 0.7) + 
  geom_label_repel(data = label_df,
                   aes(label = X, color = sig),
                   size = 2, show.legend = F, max.overlaps = Inf) + 
  scale_color_manual(values = c("firebrick1", "gray60")) +
  ggtitle("DEGs in pGA+ spots", subtitle = "C9orf72-ALS motor cortex") +
  labs(color = "Adjusted\np < 0.05") +
  theme_linedraw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.text = element_text(color = "black"))
ggsave(filename = paste0(in_dir, "/mcx_pga_multivolcano.png"),
       units = "in", dpi = 600,
       height = 4, width = 4)
