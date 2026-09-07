library(Seurat)
library(scCustomize)
library(tidyverse)
library(BPCells)
library(readxl)

setwd("/projects/b1169/boles/als_cns_visium")

plots_dir <- "figures/"

# QC metrics --------------------------------------------------------------

meta <- readRDS("data/02_qc/prefilter_metadata.rds")

meta$code <- str_split_i(rownames(meta), "_", i = 2)
meta$batch <- str_split_i(meta$code, "-", i = 1)

# Get demographic data ----------------------------------------------------

key <- read_xlsx("tab_data/master.xlsx",
                 sheet = 1)

key <- key %>%
  dplyr::select(c(code, sample, tissue)) %>% 
  mutate(batch = str_split_i(code, "-", i = 1)) %>%
  mutate(code = if_else(str_detect(code, "AN|JSB"), code, paste0("JSB", code)))

demo <- read.csv("tab_data/target_als_demographics_compiled.csv")

demo <- demo %>% 
  mutate(Case.Number = str_remove_all(Case.Number, "-"),
         group = case_when(C9orf72.mutation == "Y" ~ "C9orf72",
                           Clinical.Diagnosis == "Control" ~ "Control",
                           .default = "sALS")) %>% 
  dplyr::rename("sample" = "Case.Number",
                "age" = "Age.at.Death",
                "sex" = "Sex") %>% 
  dplyr::select(c(sample, group, age, sex))

meta <- meta %>% 
  rownames_to_column(var = "barcode") %>%
  left_join(key, 
            by = "code") %>% 
  left_join(demo,
            by = "sample")

meta <- meta %>%
  mutate(group = factor(group,
                        levels = c("Control", "sALS", "C9orf72"),
                        labels = c("Control", "sALS", "C9orf72-ALS")))

meta <- meta %>% 
  mutate(tissue = factor(tissue,
                         levels = c("mcx", "sc"),
                         labels = c("Motor cortex", "Spinal cord")))

theme <- theme_linedraw(base_size = 12) + 
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "gray", color = "black"))

meta %>% 
  arrange(group) %>%
  mutate(sample = fct_inorder(sample)) %>%
  ggplot(aes(x = sample,
             y = nCount_Spatial)) + 
  facet_wrap(. ~ tissue,
             ncol = 1,
             scales = "free") + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "# UMIs per spot") +
  scale_y_log10() + 
  theme
ggsave(filename = paste0(plots_dir, "ncount.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)

meta %>% 
  arrange(group) %>%
  mutate(sample = fct_inorder(sample)) %>%
  ggplot(aes(x = sample,
             y = nFeature_Spatial)) + 
  facet_wrap(. ~ tissue,
             ncol = 1,
             scales = "free") + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "# genes per spot") +
  scale_y_log10() + 
  theme
ggsave(filename = paste0(plots_dir, "nfeature.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)

meta %>% 
  arrange(group) %>%
  mutate(sample = fct_inorder(sample)) %>%
  ggplot(aes(x = sample,
             y = percent_mito)) + 
  facet_wrap(. ~ tissue,
             ncol = 1,
             scales = "free") + 
  geom_violin(aes(fill = group)) +
  scale_fill_manual(values = c("#b8b0a8", "#CC00FF", "#0CAA00")) + 
  labs(y = "% mitochondrial genes per spot") +
  theme
ggsave(filename = paste0(plots_dir, "mito.png"),
       units = "in", dpi = 600,
       height = 6, width = 10)

# Show annotated spots ----------------------------------------------------

in_dir <- "data/04_spot_annotation/"

counts <- open_matrix_dir(paste0(in_dir, "bpcells_data"))
meta <- readRDS(paste0(in_dir, "metadata.rds"))
images <- readRDS(paste0(in_dir, "images.rds"))

counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")
obj@images <- images

obj$region <- factor(obj$region,
                     levels = c("GM", "WM", "Meninges", "Nerve bundle"))

codes <- key %>% 
  filter(sample == "GWF2156") %>% 
  pull(code)

p1 <- SpatialDimPlot(obj,
               images = "AN67.7",
               group.by = "region",
               image.alpha = 0,
               pt.size = 2) + 
  ggtitle("Motor cortex") +
  guides(fill = guide_legend(override.aes = list(size = 4)))

p2 <- SpatialDimPlot(obj,
               images = "AN72.4",
               group.by = "region",
               image.alpha = 0,
               pt.size = 2) +
  ggtitle("Spinal cord") +
  guides(fill = guide_legend(override.aes = list(size = 4)))

p1 + p2 + 
  plot_layout(ncol = 1) & 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(plots_dir, "spots_anatomy.png"),
       units = "in", dpi = 600,
       height = 6, width = 5)

p3 <- SpatialDimPlot(obj,
               images = "AN67.7",
               group.by = "ptdp",
               image.alpha = 0,
               pt.size = 2,
               cols = c("gray80", "midnightblue")) + 
  ggtitle("Motor cortex") +
  guides(fill = guide_legend(override.aes = list(size = 4)))

p4 <- SpatialDimPlot(obj,
               images = "AN72.4",
               group.by = "ptdp",
               image.alpha = 0,
               pt.size = 2,
               cols = c("gray80", "midnightblue")) + 
  ggtitle("Motor cortex") +
  guides(fill = guide_legend(override.aes = list(size = 4)))

p3 + p4 + 
  plot_layout(ncol = 1, guides = "collect") & 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(plots_dir, "spots_ptdp.png"),
       units = "in", dpi = 600,
       height = 6, width = 5)

p5 <- SpatialDimPlot(obj,
               images = "AN67.7",
               group.by = "pga",
               image.alpha = 0,
               pt.size = 2,
               cols = c("gray80", "firebrick")) + 
  ggtitle("Motor cortex") +
  guides(fill = guide_legend(override.aes = list(size = 4)))

p5 + plot_spacer() + 
  plot_layout(ncol = 1, guides = "collect") & 
  theme(legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5))
ggsave(filename = paste0(plots_dir, "spots_pga.png"),
       units = "in", dpi = 600,
       height = 6, width = 5)
