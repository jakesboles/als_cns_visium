suppressMessages({
  library(tidyverse) 
  library(Seurat)
  library(readxl)
})

setwd("/projects/b1169/boles/als_cns_visium")

# Get meta data from 02 ---------------------------------------------------

meta <- readRDS("data/02_qc/cns_metadata.rds")

meta$code <- str_split_i(rownames(meta), "_", i = 2)

# Get demographic data ----------------------------------------------------

key <- read_xlsx("tab_data/master.xlsx",
                 sheet = 1)

key <- key %>%
  dplyr::select(c(code, sample, tissue)) %>% 
  mutate(batch = str_split_i(code, "-", i = 1)) %>%
  mutate(code = if_else(str_detect(code, "AN"), code, paste0("JSB", code)))

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

# Load pathology/feature annotations from 03c -----------------------------

samples <- list.dirs("data/03c_make_halo_feature_gdfs",
                     recursive = F,
                     full.names = F)

gdfs2 <- map(paste0("data/03c_make_halo_feature_gdfs/", samples, "/results/in_roi.csv"),
             read.csv)

gdfs2 <- list_rbind(gdfs1)


options(future.globals.maxSize=1048576000000)
load("/projects/b1169/projects/sea_ad_hypothalamus/results/preprocessing/qc/out_TW_05-04-2023/helperfunctions.RData")
setwd("/gpfs/projects/b1169/thomas/als_multitissue/Visium/HALO/attach")

s <- readRDS("/gpfs/projects/b1169/thomas/als_multitissue/Visium/CCA/CCA.rds")

s@meta.data$barcode <-  gsub("^([^_]*_[^_]*)_", "", rownames(s@meta.data))

s@meta.data$sample_id[startsWith(s@meta.data$sample_id, "14")] <- paste0("JSB", s@meta.data$sample_id[startsWith(s@meta.data$sample_id, "14")])

s@meta.data$sample_barcode <- paste0(s@meta.data$sample_id, "_", s@meta.data$barcode)

mngs <- c("JSB146-1", "JSB146-8", "AN67-3", "AN67-7", "AN68-1")

s@meta.data$Region <- NA

for(sample in unique(s@meta.data$sample_id)){
  
  wdat <- read.csv(paste0("/gpfs/projects/b1169/boles/als_motor_circuit_visium/halo_annotations/cortical_layers/", sample, "/results/in_roi_wm.csv"))
  
  wdat <- wdat[wdat$barcode %in% s@meta.data$sample_barcode,]
  
  rownames(wdat) <- wdat$barcode
  
  wdat <- wdat[s@meta.data$sample_barcode[s@meta.data$sample_id == sample],]
  
  all.equal(wdat$barcode, s@meta.data$sample_barcode[s@meta.data$sample_id == sample])
  
  wdat$WM <- ifelse(wdat$in_roi == "True", "White Matter", NA)
  
  if(sample %in% mngs){
    
    mdat <- read.csv(paste0("/gpfs/projects/b1169/boles/als_motor_circuit_visium/halo_annotations/cortical_layers/", sample, "/results/in_roi_mng.csv"))
    
    mdat <- mdat[mdat$barcode %in% s@meta.data$sample_barcode,]
    
    rownames(mdat) <- mdat$barcode
    
    mdat <- mdat[s@meta.data$sample_barcode[s@meta.data$sample_id == sample],]
    
    all.equal(mdat$barcode, s@meta.data$sample_barcode[s@meta.data$sample_id == sample])
    
    mdat$MNG <- ifelse(mdat$in_roi == "True", "Meninges", NA)
    
    dat <- cbind(wdat,"MNG" = mdat$MNG)
    
    dat$Tissue <- ifelse(is.na(dat$WM) & is.na(dat$MNG), "Gray Matter", NA)
    
    dat$Tissue[dat$WM == "White Matter"] <- "White Matter"
    
    dat$Tissue[dat$MNG == "Meninges"] <- "Meninges"
    
    s@meta.data$Region[s@meta.data$sample_id == sample] <- dat$Tissue
    
  }else{
    
    wdat$WM[is.na(wdat$WM)] <- "Gray Matter"
    
    s@meta.data$Region[s@meta.data$sample_id == sample] <- wdat$WM
    
  }
  
}

SpatialDimPlot(s, group.by = "Region", images = "X147.1")

DimPlot(s, group.by = "Region", reduction = "CCA_UMAP", raster = FALSE)

for(sample in names(s@images)){
  
  pdf(paste0("plots/", sample, ".pdf"), width = 8, height = 8)
  print(SpatialDimPlot(s, group.by = "Region", images = sample))
  dev.off()
  
}

saveRDS(s, "/gpfs/projects/b1169/boles/als_motor_circuit_visium/data/02_HALO_CCA/s_HALO_CCA.rds")

