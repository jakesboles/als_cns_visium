library(tidyverse)
library(ggrepel)
library(scCustomize)
library(UpSetR)
library(ggplot2)
library(reshape2)
library(paletteer)

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/deseq2_by_compartment"

get_degs <- function(df,
                     pval_cutoff,
                     logfc_cutoff){
  genes <- df %>% 
    filter(padj < pval_cutoff & 
             abs(log2FoldChange) > logfc_cutoff) %>% 
    pull(X)
  
  return(genes)
}

p_thresh <- 0.05
lfc_thresh <- log2(1.5)

degs <- list()

dirs <- list.dirs(results_dir,
                  full.names = T,
                  recursive = F)

region <- str_split_i(dirs, "/", i = 3)

df <- matrix(nrow = length(dirs),
              ncol = 5)

for (i in seq_along(dirs)){
  
  c9 <- read.csv(paste0(dirs[i], "/C9orf72_vs_Control.csv"))
  sals <- read.csv(paste0(dirs[i], "/sALS_vs_Control.csv"))
  
  c9_degs <- get_degs(c9,
                      pval_cutoff = p_thresh,
                      logfc_cutoff = lfc_thresh)
  sals_degs <- get_degs(sals,
                        pval_cutoff = p_thresh,
                        logfc_cutoff = lfc_thresh)
  
  shared <- intersect(sals_degs,
                      c9_degs)
  
  c9_only <- setdiff(c9_degs,
                     sals_degs)
  
  sals_only <- setdiff(sals_degs,
                       c9_degs)
  
  df[i, ] <- c(region[i], length(shared), length(c9_only), length(sals_only), length(union(sals_degs, c9_degs)))
  
}

df <- as.data.frame(df)
colnames(df) <- c("region", "shared", "c9_only", "sals_only", "total")

df <- df %>% 
  mutate(tissue = str_split_i(region, "_", i = 1),
         compartment = str_split_i(region, "_", i = 2))

df %>%
  pivot_longer(c(shared, c9_only, sals_only)) %>%
  mutate(total = as.numeric(total),
         value = as.numeric(value)) %>%
  arrange(desc(total)) %>% 
  mutate(region = fct_inorder(region)) %>% 
  mutate(name = factor(name, 
                       levels = c("sals_only", "c9_only", "shared"),
                       labels = c("sALS only", "C9-ALS only", "Shared"))) %>%
  ggplot(aes(x = region,
             y = value)) +
  geom_bar(aes(fill = name),
           stat = "identity",
           color = "black",
           linewidth = 0.4) + 
  scale_fill_manual(values = c("#0CAA00", "#CC00FF", "darkslategrey")) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(
    labels = function(x) str_replace_all(x, "_", " ") %>% str_to_upper()
  ) +
  ylab("# DEGs") + 
  labs(fill = "DEGs versus\ncontrol in:") +
  theme_linedraw(base_size = 16) +
  theme(axis.title.x = element_blank(),
        # axis.title.y = element_text(face = "bold"),
        axis.text = element_text(color = "black", size = 14),
        # axis.text.x = element_text(angle = , hjust = 1, vjust = 1),
        # axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.title = element_text(size = 20, hjust = 0.5),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text = element_text(face = "bold", size = 20, color = "black"))
  