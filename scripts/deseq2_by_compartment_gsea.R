library(tidyverse)
library(ggrepel)
library(ggplot2)
library(paletteer)
library(msigdbr)
library(clusterProfiler)

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/deseq2_by_compartment/"

# Load GSEA pathways ------------------------------------------------------

m_t2g <- msigdbr(species = "Homo sapiens",
                 category = "C2")
m_t2g <- m_t2g %>%
  filter(gs_subcollection %in% c("CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:PID")) %>%
  dplyr::select(c(gs_name, gene_symbol))
# unique(m_t2g$gs_name)

go_t2g <- msigdbr(species = "Homo sapiens",
                  category = "C5")
go_t2g <- go_t2g %>%
  # filter(gs_subcollection %in% c("GO:BP", "GO:CC") %>%
  dplyr::select(c(gs_name, gene_symbol))

t2g <- rbind(go_t2g, m_t2g)

# Bring in LFC shrunk results ---------------------------------------------

get_ranked_genes <- function(data){
  rank <- data %>% 
    arrange(desc(log2FoldChange))
  
  genes <- rank$X
  rank <- rank$log2FoldChange
  
  names(rank) <- genes
  
  return(rank)
}

dirs <- list.dirs(results_dir,
                  full.names = F,
                  recursive = F)

for (j in dirs) {
  
  message(j)
  
  sals <- read.csv(paste0(results_dir, "/", j, "/sALS_vs_Control_lfc_shrunk.csv"))
  
  c9 <- read.csv(paste0(results_dir, "/", j, "/C9orf72_vs_Control_lfc_shrunk.csv"))
  
  sals_rank <- get_ranked_genes(sals)
  
  c9_rank <- get_ranked_genes(c9)
  
  suppressMessages({
    
    sals_res <- GSEA(sals_rank,
                     TERM2GENE = t2g)
    
    c9_res <- GSEA(c9_rank,
                   TERM2GENE = t2g)
    
  })
  
  if (sals_res@result$p.adjust %>% min() < 0.05){
    p <- dotplot(sals_res,
                 showCategory = 20,
                 x = "NES",
                 color = "p.adjust",
                 size = "setSize")
    ggsave(p,
           filename = paste0(results_dir, j, "/sALS_vs_Control_GSEA.png"),
           units = "in", dpi = 600,
           height = 10, width = 8)
  }
  
  if (c9_res@result$p.adjust %>% min() < 0.05){
    p <- dotplot(c9_res,
                 showCategory = 20,
                 x = "NES",
                 color = "p.adjust",
                 size = "setSize")
    ggsave(p,
           filename = paste0(results_dir, j, "/C9orf72_vs_Control_GSEA.png"),
           units = "in", dpi = 600,
           height = 10, width = 8)
  }
  
  write.csv(as.data.frame(sals_res@result),
            file = paste0(results_dir, j, "/sALS_vs_Control_GSEA.csv"))
  
  write.csv(as.data.frame(c9_res@result),
            file = paste0(results_dir, j, "/C9orf72_vs_Control_GSEA.csv"))
}
