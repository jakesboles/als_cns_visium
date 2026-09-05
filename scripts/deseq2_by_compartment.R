# Pseudobulk DESeq2 differential expression within each anatomical
# compartment (motor cortex gray matter, motor cortex white matter, spinal
# cord gray matter, spinal cord white matter), comparing sALS vs Control and
# C9orf72 vs Control. Counts are pseudobulked per donor (sample), not per
# Visium section, so a donor with multiple sections of the same tissue
# still contributes one pseudobulk column. Modeled on
# als_cns_scrnaseq/r_scripts/deseq2.R's pseudobulk-per-cell-type design,
# adapted to pseudobulk-per-anatomical-region instead.
#
# Non-anatomical feature annotations (pTDP-43, pGA) are handled by a
# separate follow-up script once those annotations are complete -- this one
# only covers the region calls from 03b_make_halo_gdfs.py /
# 04_spot_annotation.R.

suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(BPCells)
  library(DESeq2)
  library(apeglm)
  library(IHW)
})

setwd("/projects/b1169/boles/als_cns_visium")

results_dir <- "results/deseq2_by_compartment/"
dir.create(results_dir, showWarnings = F, recursive = T)

data_dir <- "data/deseq2_by_compartment/"
dir.create(data_dir, showWarnings = F, recursive = T)

# Load annotated metadata and raw counts -------------------------------------
# 04_spot_annotation.R's metadata.rds carries region (GM/WM/Meninges/Nerve
# bundle), tissue, sample (donor), group, sex, and age for every spot it
# retained. Raw counts come from 02_qc.R's BPCells matrix (pre-normalization,
# unlike anything saved from 03_integration_harmony.R onward), subset to
# exactly the spots 04 kept.

message("Reading in metadata and raw counts")

meta <- readRDS("data/04_spot_annotation/metadata.rds")
rownames(meta) <- meta$barcode

counts <- open_matrix_dir("data/02_qc/bpcells_cns")
counts <- counts[, rownames(meta)]

obj <- CreateSeuratObject(counts = counts, meta.data = meta, assay = "Spatial")

# Factorize grouping variable, Control as the reference level ---------------

obj$group <- factor(obj$group, levels = c("Control", "sALS", "C9orf72"))
obj$sex <- factor(obj$sex)

# Define compartments ---------------------------------------------------
# Only gray and white matter are analyzed here -- Meninges/Nerve bundle
# spots exist in the annotated metadata but aren't part of this comparison.

compartments <- data.frame(
  tissue = c("mcx", "mcx", "sc", "sc"),
  region = c("GM", "WM", "GM", "WM"),
  title = c("Motor cortex gray matter", "Motor cortex white matter",
           "Spinal cord gray matter", "Spinal cord white matter"),
  file = c("mcx_gm", "mcx_wm", "sc_gm", "sc_wm")
)

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

for (i in seq_len(nrow(compartments))){

  message(paste0(compartments$title[i], " (", i, "/", nrow(compartments), ")"))

  sub <- subset(obj,
               tissue == compartments$tissue[i] & region == compartments$region[i])

  file <- compartments$file[i]

  comp_results_dir <- paste0(results_dir, file, "/")
  dir.create(comp_results_dir, showWarnings = F, recursive = T)

  bulk <- AggregateExpression(sub,
                              assays = "Spatial",
                              return.seurat = F,
                              # layer = "counts",
                              group.by = c("sample"))

  exp <- bulk$Spatial

  spot_counts <- sub@meta.data %>%
    dplyr::count(sample, name = "n_spots")

  # Every donor present in this compartment, whether or not it survives the
  # min_spots_per_sample filter below -- saved so a skipped/thinned
  # compartment's sample composition can be checked later without rerunning
  # anything.
  sample_table <- sub@meta.data %>%
    dplyr::select(sample, group) %>%
    distinct() %>%
    left_join(spot_counts, by = "sample") %>%
    mutate(retained = n_spots >= min_spots_per_sample) %>%
    arrange(group, sample)

  write.csv(sample_table,
            file = paste0(comp_results_dir, "sample_filtering.csv"),
            row.names = F)

  meta_comp <- sub@meta.data %>%
    dplyr::select(c(sample, group, sex, age)) %>%
    distinct() %>%
    left_join(spot_counts, by = "sample") %>%
    filter(n_spots >= min_spots_per_sample)

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
    dplyr::select(-n_spots) %>%
    mutate(age_scale = scale(age, center = T, scale = T)[,1])

  exp <- exp[, meta_comp$sample, drop = F]

  idx <- match(colnames(exp), meta_comp$sample)
  meta_comp <- meta_comp[idx, ]
  rownames(meta_comp) <- meta_comp$sample

  # The abundance filter above catches the most common cause of DESeq2's
  # "every gene contains at least one zero" size factor error, but not
  # every case (e.g. a gene that's zero in every retained sample even
  # though each sample individually cleared min_spots_per_sample). Wrap the
  # whole DESeq2 pipeline so a failure on one compartment is logged and
  # skipped instead of killing the rest of the script.
  tryCatch({

    dds <- DESeqDataSetFromMatrix(countData = exp,
                                  colData = meta_comp,
                                  design = ~ sex + age_scale + group) # change this as needed

    keep <- rowSums(counts(dds) >= 10) >= 10 # change these cutoffs as needed

    dds <- dds[keep, ]

    dds <- DESeq(dds)

    saveRDS(dds,
            file = paste0(data_dir, file, "_dds.rds"))

    # resultsNames(dds)

    res <- results(dds,
                   contrast = c("group", "sALS", "Control"),
                   filterFun = ihw,
                   independentFiltering = T)

    res <- as.data.frame(res)

    write.csv(res,
              file = paste0(comp_results_dir, "sALS_vs_Control.csv"))

    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "group_sALS_vs_Control",
                              type = "apeglm")
    })

    write.csv(res_shrunk,
              file = paste0(comp_results_dir, "sALS_vs_Control_lfc_shrunk.csv"))

    res <- results(dds,
                   contrast = c("group", "C9orf72", "Control"),
                   filterFun = ihw,
                   independentFiltering = T)

    res <- as.data.frame(res)

    write.csv(res,
              file = paste0(comp_results_dir, "C9orf72_vs_Control.csv"))

    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "group_C9orf72_vs_Control",
                              type = "apeglm")
    })

    write.csv(res_shrunk,
              file = paste0(comp_results_dir, "C9orf72_vs_Control_lfc_shrunk.csv"))

  }, error = function(e){
    message(paste0("Skipping ", compartments$title[i], " -- DESeq2 pipeline failed: ",
                   conditionMessage(e)))
  })

}
