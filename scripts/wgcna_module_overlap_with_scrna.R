library(tidyverse)

setwd("/projects/b1169/boles/als_cns_visium/")

st_modules <- read.csv("results/wgcna_consensus/modules.csv") %>%
  dplyr::select(c(gene_name, module))

scrna_modules <- read.csv("../als_cns_scrnaseq/results/wgcna_consensus/Microglia/modules.csv") %>%
  dplyr::select(c(gene_name, module))

tab <- full_join(st_modules,
                 scrna_modules,
                 by = "gene_name") %>%
  replace(is.na(.), "dropped") # %>% 
  # mutate(module.x = if_else(module.x == "yellow", module.x, "gold"),
  #        module.y = if_else(module.y %in% c("blue", "turquoise"), module.y, "gold"))

colors1 <- unique(tab[,2])
colors2 <- unique(tab[,3])

n1 <- length(colors1)
n2 <- length(colors2)

p_tab <- matrix(0, nrow = n1, ncol = n2)
count_tab <- matrix(0, nrow = n1, ncol = n2)

bonf <- n1 * n2

for (k in 1:n1){
  for(j in 1:n2){
    
    members1 <- tab[,2] == colors1[k]
    members2 <- tab[,3] == colors2[j]
    
    p_tab[k, j] = bonf * (fisher.test(members1, members2, alternative = "greater")$p.value)
    
    count_tab[k, j] = length(intersect(tab$gene[members1], tab$gene[members2]))
  }
}

p_tab_long <- p_tab %>%
  as.data.frame()

colnames(p_tab_long) <- colors2
rownames(p_tab_long) <- colors1

count_tab_long <- count_tab %>%
  as.data.frame() %>%
  rownames_to_column(var = "y") %>%
  mutate(y = as.numeric(y)) %>%
  mutate(y = (n1 + 1) - y) %>%
  pivot_longer(2:(ncol(count_tab) + 1),
               names_to = "x")

count_tab_long$x <- str_replace_all(count_tab_long$x, "V", "")
count_tab_long$x <- as.numeric(count_tab_long$x)

p_tab_long <- p_tab_long %>%
  rownames_to_column(var = "colors1") %>%
  pivot_longer(2:(ncol(p_tab_long) + 1),
               names_to = "colors2")
p_tab_long$colors1 <- factor(p_tab_long$colors1,
                             levels = rev(unique(p_tab_long$colors1)))
p_tab_long$colors2 <- factor(p_tab_long$colors2,
                             levels = unique(p_tab_long$colors2))

colors1_factor <- levels(p_tab_long$colors1)
colors2_factor <- levels(p_tab_long$colors2)

p_tab_long %>%
  mutate(value = if_else(value < 0.05, "P < 0.05", "P > 0.05")) %>%
  # dplyr::rename("consModules" = "consMod", "spModules" = "spMod") %>%
  # left_join(spGenes, by = "spModules") %>%
  # left_join(consGenes, by = "consModules") %>%
  # mutate(spModTotals = as.character(spModTotals)) %>%
  # mutate(consModTotals = as.character(consModTotals)) %>%
  # unite(col = "spModules", c(mbModules, spModTotals), sep = ": ") %>%
  # unite(col = "consModules", c(consModules, consModTotals), sep = ": ") %>%
  # mutate(spModules = factor(spModules,
  #                           levels = rev(unique(spModules)))) %>%
  # mutate(consModules = factor(consModules,
  #                             levels = unique(consModules))) %>%
  ggplot(aes(x = colors2, y = colors1)) + 
  geom_tile(aes(fill = value), color = "black") +
  scale_fill_manual(values = c("red", "white")) + 
  # labs(y = paste0(str_to_title(tissues[i]), "-specific modules"), x = "Consensus modules") +
  # ggtitle(paste0("Correspondence of ", tissues[i], "-specific modules and brain consensus modules")) +
  geom_text(data = count_tab_long, aes(x = x, y = y, label = value),
            size = 3) +
  theme_classic() +
  labs(x = "SC",
       y = "ST") + 
  # scale_y_discrete(limits = rev) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5),
    axis.text = element_text(color = "black")
  )
ggsave(filename = paste0("results/wgcna_cross_modality/spatial_vs_microglia_module_overlap.png"),
       units = "in", dpi = 600,
       height = 5, width = 8)
