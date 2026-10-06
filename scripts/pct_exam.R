library(Seurat)
library(tidyverse)
library(BPCells)
library(ggbeeswarm)

# Build Visium object -----------------------------------------------------

setwd("/projects/b1169/boles/als_cns_visium")

in_dir <- "data/04_spot_annotation/"

counts <- open_matrix_dir(paste0(in_dir, "bpcells_data"))
meta <- readRDS(paste0(in_dir, "metadata.rds"))
images <- readRDS(paste0(in_dir, "images.rds"))

counts <- counts[, rownames(meta)]

st <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
st@images <- images

# Build scRNAseq object ---------------------------------------------------

setwd("/projects/b1169/boles/als_cns_scrnaseq")

meta_list <- map(paste0(list.dirs("data/17_obj_reassembly",
                 recursive = F), "/metadata.rds"),
                 readRDS)

meta <- list_rbind(meta_list)

meta_sub <- meta %>%
  filter(cell_type3 == "Microglia")

raw_mat <- open_matrix_dir("data/06_obj_reassembly/bpcells")
raw_mat <- raw_mat[, rownames(meta_sub)]

sc <- CreateSeuratObject(counts = raw_mat, 
                         meta.data = meta_sub)

# Get WGCNA modules from scRNAseq -----------------------------------------

modules <- read.csv("results/wgcna_consensus/Microglia/modules.csv")

blue_genes <- modules %>% 
  filter(module == "blue") %>% 
  pull(gene_name)

turquoise_genes <- modules %>% 
  filter(module == "turquoise") %>% 
  pull(gene_name)

blue_kme <- modules %>% 
  filter(module == "blue") %>% 
  dplyr::select(c(gene_name, kME_blue))

turquoise_kme <- modules %>% 
  filter(module == "turquoise") %>% 
  dplyr::select(c(gene_name, kME_turquoise))

colnames(blue_kme) <- c("gene", "kme") -> colnames(turquoise_kme)

kme <- rbind(blue_kme, turquoise_kme)

exp_data <- FetchData(st, vars = blue_genes)
blue_st_pct <- colMeans(exp_data > 0) * 100

exp_data <- FetchData(st, vars = turquoise_genes)
turquoise_st_pct <- colMeans(exp_data > 0) * 100

exp_data <- FetchData(sc, vars = blue_genes)
blue_sc_pct <- colMeans(exp_data > 0) * 100

exp_data <- FetchData(sc, vars = turquoise_genes)
turquoise_sc_pct <- colMeans(exp_data > 0) * 100

print(blue_st_pct)

df <- tibble(
  pct = c(blue_st_pct, turquoise_st_pct, blue_sc_pct, turquoise_sc_pct),
  gene = c(names(blue_st_pct), names(turquoise_st_pct), names(blue_sc_pct), names(turquoise_sc_pct)),
  modality = c(rep("st", length(c(blue_st_pct, turquoise_st_pct))),
               rep("sc", length(c(blue_sc_pct, turquoise_sc_pct)))),
  module = c(rep("blue", length(blue_st_pct)), rep("turquoise", length(turquoise_st_pct)),
             rep("blue", length(blue_sc_pct)), rep("turquoise", length(turquoise_sc_pct))))

df <- df %>% 
  left_join(kme,
            by = "gene")

df$modality_factor <- factor(df$modality,
                      levels = c("sc", "st"),
                      labels = c("scRNAseq microglia", "ST spots"))

df$module <- factor(df$module,
                    levels = c("blue", "turquoise"),
                    labels = c("Microglia3 module", "Microglia4 module"))

df %>% 
  ggplot(aes(x = modality_factor,
             y = pct)) + 
  # geom_quasirandom() + 
  geom_violin(aes(fill = modality_factor),
              show.legend = F) +
  facet_wrap(. ~ module,
             ncol = 1) + 
  labs(y = "Percent expression") +
  theme_linedraw(base_size = 16) + 
  theme(axis.title.x = element_blank())

df %>%
  ggplot(aes(x = modality_factor,
             y = pct)) + 
  geom_line(aes(group = gene)) + 
  facet_wrap(. ~ module,
             ncol = 1) + 
  theme_linedraw()

blue_sc_pct %>% sort()

df %>% 
  dplyr::select(-modality_factor) %>%
  pivot_wider(names_from = modality,
              values_from = pct) %>%
  mutate(change = st/sc) %>%
  ggplot(aes(x = module,
             y = change)) + 
  # geom_quasirandom() +
  geom_violin(aes(fill = module),
              show.legend = F) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  scale_fill_manual(values = c("blue", "turquoise")) +
  labs(y = "ST pct / SC pct") +
  theme_linedraw(base_size = 16) + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

df %>% 
  pivot_wider(names_from = modality,
              values_from = pct) %>%
  filter(module == "blue") %>% 
  arrange(desc(sc)) %>%
  print(n = 100)

df %>% 
  pivot_wider(names_from = modality,
              values_from = pct) %>%
  ggplot(aes(x = sc,
             y = st)) + 
  geom_point() + 
  geom_abline(intercept = 0, slope = 1,
              color = "red") +
  facet_wrap(. ~ module) + 
  theme_linedraw()

df %>% 
  filter(modality == "st") %>% 
  ggplot(aes(x = kme,
             y = pct)) + 
  geom_point() + 
  labs(y = "Percent expression in ST spots",
       x = "kME in respective modules") +
  facet_wrap(. ~ module,
             ncol = 1,
             scales = "free_x") + 
  theme_linedraw(base_size = 16)
