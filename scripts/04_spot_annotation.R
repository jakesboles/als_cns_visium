suppressMessages({
  library(tidyverse) 
  library(Seurat)
  library(readxl)
})

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/04_spot_annotation/"
dir.create(results_dir,
           showWarnings = F,
           recursive = T)

# Get meta data from 02 ---------------------------------------------------

meta <- readRDS("data/02_qc/cns_metadata.rds")

meta$code <- str_split_i(rownames(meta), "_", i = 2)

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

# Load anatomical annotations from 03b ------------------------------------

samples <- list.dirs("data/03b_make_halo_gdfs",
                     recursive = F,
                     full.names = F)

gdfs1 <- map(paste0("data/03b_make_halo_gdfs/", samples, "/results/in_roi.csv"),
             read.csv)

gdfs1 <- list_rbind(gdfs1)
gdfs1$barcode <- paste0("_", gdfs1$barcode)

# check that barcodes match
table(meta$barcode %in% gdfs1$barcode) # all good

gdfs1 <- gdfs1 %>%
  mutate(mcx_region = case_when(in_mcx_wm == "True" ~ "WM",
                                in_mcx_meninges == "True" ~ "Meninges",
                                in_mcx_wm == "True" & in_mcx_meninges == "True" ~ "Flag",
                                in_mcx_wm == "False" & 
                                  (in_mcx_meninges == "False" | is.na(in_mcx_meninges)) ~ "GM"),
         sc_region = case_when(in_sc_wm == "True" & in_sc_gm == "False" ~ "WM",
                               in_sc_gm == "True" ~ "GM",
                               in_sc_meninges == "True" ~ "Meninges",
                               in_sc_nerve_bundles == "True" ~ "Nerve bundle",
                               (in_sc_gm == "True" & in_sc_meninges == "True") | 
                                 (in_sc_gm == "True" & in_sc_nerve_bundles == "True") | 
                                 (in_sc_meninges == "True" & in_sc_nerve_bundles == "True") | 
                                 (in_sc_gm == "True" & in_sc_wm == "True") | 
                                 (in_sc_meninges == "True" & in_sc_wm == "True") | 
                                 (in_sc_nerve_bundles == "True" & in_sc_wm == "True") ~ "Flag",
                               in_sc_wm == "False" & in_sc_gm == "False" & 
                                 (in_sc_meninges == "False" | is.na(in_sc_meninges)) & 
                                 (in_sc_nerve_bundles == "False" | is.na(in_sc_nerve_bundles)) ~ "Remove"))

gdfs1 <- gdfs1 %>%
  mutate(region = coalesce(mcx_region, sc_region)) %>% 
  dplyr::select(c(barcode, region))

table(is.na(gdfs1$region))

meta <- meta %>% 
  left_join(gdfs1, 
            by = "barcode")

samples <- unique(meta$sample)

for (i in samples){
  
  n <- meta %>% 
    filter(sample == i) %>% 
    pull(tissue) %>%
    unique() %>% 
    length()
  
  p <- meta %>% 
    filter(sample == i) %>%
    ggplot(aes(x = array_col,
               y = array_row)) + 
    geom_point(aes(fill = region),
               shape = 21) +
    facet_wrap(. ~ tissue,
               ncol = n) +
    ggtitle(i) +
    theme_void(base_size = 12) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(p,
         filename = paste0(results_dir, i, "_spots_anatomy.png"),
         units = "in", dpi = 600,
         height = 6, width = 7*n)
}
# Load pathology/feature annotations from 03c -----------------------------

samples <- list.dirs("data/03c_make_halo_feature_gdfs",
                     recursive = F,
                     full.names = F)

gdfs2 <- map(paste0("data/03c_make_halo_feature_gdfs/", samples, "/results/in_roi.csv"),
             read.csv)

gdfs2 <- list_rbind(gdfs1)