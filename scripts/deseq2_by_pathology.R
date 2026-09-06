suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(DESeq2)
  library(apeglm)
  library(IHW)
})

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/deseq2_by_pathology/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_dir <- "data/deseq2_by_pathology/"
dir.create(data_dir, showWarnings = F, recursive = T)

# Load annotated metadata and raw counts -------------------------------------
# 04_spot_annotation.R writes out the region-filtered object's own pieces
# directly: metadata.rds's rownames are the real barcodes (not a separate
# "barcode" column), and bpcells_data is that same object's raw counts
# (pre-normalization), already subset to exactly the spots 04 kept -- no
# need to re-read 02_qc.R's unfiltered matrix and re-subset by hand here.

message("Reading in metadata and raw counts")

meta <- readRDS("data/04_spot_annotation/metadata.rds")

counts <- open_matrix_dir("data/04_spot_annotation/bpcells_data")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")

# Factorize grouping variable, Control as the reference level ---------------

obj$sex <- factor(obj$sex)
obj$ptdp <- factor(obj$ptdp, 
                   levels = c(F, T))
obj$pga <- factor(obj$pga,
                  levels = c(F, T))

obj <- obj %>% 
  subset(region == "GM")

# Define compartments ---------------------------------------------------
# Only gray and white matter are analyzed here -- Meninges/Nerve bundle
# spots exist in the annotated metadata but aren't part of this comparison.

guide <- tibble(
  tissues = c("mcx", "sc", "mcx", "sc", "mcx"),
  groups = c("sALS", "sALS", "C9orf72", "C9orf72", "C9orf72"),
  features = c("ptdp", "ptdp", "ptdp", "ptdp", "pga")
) %>%
  mutate(file = paste0(features, "_", tissues))

# A pseudobulk sample built from too few spots is mostly zero, which can
# make every gene contain a zero in some sample -- DESeq2's default
# median-of-ratios size factor estimation then fails outright ("every gene
# contains at least one zero, cannot compute log geometric means"). Below
# this per-donor spot count, that donor is dropped from the compartment's
# pseudobulk; if too few donors remain in any group after dropping, the
# whole compartment is skipped rather than run on an unreliable/unbalanced
# design.
min_spots_per_sample <- 10 # change as needed
min_samples_per_group <- 3 # change as needed

for (i in seq_len(nrow(guide))){
  
  message(paste0(guide$features[i], " in ", guide$tissues[i], " in ", guide$groups[i]))
  
  sub <- subset(obj,
                tissue == guide$tissues[i] & group == guide$groups[i])
  
  file <- guide$file[i]
  
  comp_results_dir <- paste0(results_dir, file, "/")
  dir.create(comp_results_dir, showWarnings = F, recursive = T)
  
  bulk <- AggregateExpression(sub,
                              assays = "Spatial",
                              return.seurat = F,
                              # layer = "counts",
                              group.by = c("sample", guide$features[i]))
  
  exp <- bulk$Spatial
  
  spot_counts <- sub@meta.data %>%
    dplyr::count(sample, !!sym(guide$features[i]), name = "n_spots")
  
  # Every donor present in this compartment, whether or not it survives the
  # min_spots_per_sample filter below -- saved so a skipped/thinned
  # compartment's sample composition can be checked later without rerunning
  # anything.
  sample_table <- sub@meta.data %>%
    dplyr::select(sample, guide$features[i]) %>%
    distinct() %>%
    left_join(spot_counts, by = c("sample", guide$features[i])) %>%
    mutate(retained = n_spots >= min_spots_per_sample) %>%
    arrange(sample)
  
  write.csv(sample_table,
            file = paste0(comp_results_dir, "sample_filtering.csv"),
            row.names = F)
  
  # if EITHER classification is below the threshold, dump the whole sample
  dump <- sample_table %>% 
    filter(retained == F) %>% 
    pull(sample) %>% 
    unique()
  
  meta_comp <- sub@meta.data %>% 
    dplyr::select(c(sample, guide$features[i], sex, age)) %>% 
    distinct() %>%
    filter(!(sample %in% dump)) %>% 
    left_join(spot_counts, by = c("sample", guide$features[i])) %>% 
    filter(n_spots > min_spots_per_sample)
  
  if (any(table(meta_comp$group) < min_samples_per_group)){
    message(paste0("Skipping ", compartments$title[i], " -- fewer than ",
                   min_samples_per_group, " samples per group have >= ",
                   min_spots_per_sample, " spots."))
    next
  }
  
  # Unlike als_cns_scrnaseq/deseq2.R's orig.ident (which pairs a donor id
  # with a tissue code using "_", e.g. "AU-066_b"), `sample` here is just
  # the donor id with hyphens already stripped in 04_spot_annotation.R's
  # demographics join -- it never contains an underscore, so there's
  # nothing for AggregateExpression()'s internal "_" -> "-" sanitization of
  # group.by values to touch, and colnames(exp) can be matched against it
  # directly.
  meta_comp <- meta_comp %>%
    mutate(sample2 = paste0(sample, "_", !!sym(guide$features[i]))) %>%
    dplyr::select(-n_spots) %>%
    mutate(age_scale = scale(age, center = T, scale = T)[,1]) %>% 
    dplyr::rename("feature" = guide$features[i])
  
  exp <- exp[, meta_comp$sample2, drop = F]
  
  idx <- match(colnames(exp), meta_comp$sample2)
  meta_comp <- meta_comp[idx, ]
  rownames(meta_comp) <- meta_comp$sample2
  
  # The abundance filter above catches the most common cause of DESeq2's
  # "every gene contains at least one zero" size factor error, but not
  # every case (e.g. a gene that's zero in every retained sample even
  # though each sample individually cleared min_spots_per_sample). Wrap the
  # whole DESeq2 pipeline so a failure on one compartment is logged and
  # skipped instead of killing the rest of the script.
  tryCatch({
    
    dds <- DESeqDataSetFromMatrix(countData = exp,
                                  colData = meta_comp,
                                  design = ~ sample + feature) # change this as needed
    
    keep <- rowSums(counts(dds) >= 10) >= 5 # change these cutoffs as needed
    
    dds <- dds[keep, ]
    
    dds <- DESeq(dds)
    
    saveRDS(dds,
            file = paste0(data_dir, file, "_dds.rds"))
    
    # resultsNames(dds)
    
    res <- results(dds,
                   contrast = c("feature", T, F),
                   filterFun = ihw,
                   independentFiltering = T)
    
    res <- as.data.frame(res)
    
    write.csv(res,
              file = paste0(comp_results_dir, guide$groups[i], ".csv"))
    
    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "ptdp_TRUE_vs_FALSE",
                              type = "apeglm")
    })
    
    write.csv(res_shrunk,
              file = paste0(comp_results_dir, guide$groups[i], "_lfc_shrunk.csv"))
    
  }, error = function(e){
    message(paste0("Skipping ", compartments$title[i], " -- DESeq2 pipeline failed: ",
                   conditionMessage(e)))
  })
  
}
