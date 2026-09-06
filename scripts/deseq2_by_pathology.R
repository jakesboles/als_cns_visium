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

# Define comparisons ---------------------------------------------------
# NOTE: only GM is analyzed here (region == "GM" above) -- flagged
# separately, since this doesn't cover WM despite what a copy-pasted
# comment used to say.
#
# No mcx/sALS/pga row and no sc/pga row at all -- deliberate, not a gap:
# pGA (poly-GA) is a dipeptide repeat protein specific to the C9orf72
# hexanucleotide repeat expansion, so sALS donors (no C9orf72 mutation)
# aren't expected to have pGA pathology to compare, and there's no sc_pga
# Halo annotation category yet (see data/halo_annotations/features/).

guide <- tibble(
  tissues = c("mcx", "sc", "mcx", "sc", "mcx"),
  groups = c("sALS", "sALS", "C9orf72", "C9orf72", "C9orf72"),
  features = c("ptdp", "ptdp", "ptdp", "ptdp", "pga")
) %>%
  # Must include `groups`, not just `features`/`tissues` -- ptdp/mcx and
  # ptdp/sc each appear for two different groups above, and file is used
  # for both the results subdirectory AND the flat dds.rds path below. If
  # file collides across rows, the second row silently overwrites the
  # first row's dds.rds and sample_filtering.csv (the per-group results
  # CSVs were already fine, since those are separately named by
  # guide$groups[i] within comp_results_dir).
  mutate(file = paste0(features, "_", tissues, "_", groups))

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
  
  # tidyr::complete() fills in an explicit n_spots = 0 row for a donor
  # missing one feature status entirely (e.g. an sALS donor with zero
  # pTDP+ spots in GM) -- without it, a donor like that would never
  # appear as a row for the missing status at all, so it could never get
  # marked retained = FALSE / dumped below, and would sneak into meta_comp
  # with only one feature level. design = ~ sample + feature further down
  # needs every retained donor to contribute both levels (sample is used
  # as a paired/blocking term) -- an unpaired donor there makes that
  # donor's own sample indicator collinear with the feature indicator,
  # which DESeq2 will refuse to fit ("model matrix is not full rank").
  spot_counts <- sub@meta.data %>%
    dplyr::count(sample, !!sym(guide$features[i]), name = "n_spots") %>%
    tidyr::complete(sample, !!sym(guide$features[i]), fill = list(n_spots = 0))

  # Every donor x feature-status combination (including a completed-in
  # zero-spot one), whether or not it survives the min_spots_per_sample
  # filter below -- saved so a skipped/thinned comparison's sample
  # composition can be checked later without rerunning anything.
  sample_table <- spot_counts %>%
    mutate(retained = n_spots >= min_spots_per_sample) %>%
    arrange(sample)

  write.csv(sample_table,
            file = paste0(comp_results_dir, "sample_filtering.csv"),
            row.names = F)

  # If EITHER feature status is below the threshold for a donor (now
  # including a status with zero spots, thanks to complete() above), dump
  # the whole donor -- see the note above spot_counts for why.
  dump <- sample_table %>%
    filter(retained == F) %>%
    pull(sample) %>%
    unique()

  meta_comp <- sub@meta.data %>%
    dplyr::select(sample, all_of(guide$features[i]), sex, age) %>%
    distinct() %>%
    filter(!(sample %in% dump)) %>%
    left_join(spot_counts, by = c("sample", guide$features[i])) %>%
    filter(n_spots >= min_spots_per_sample)

  if (any(table(meta_comp[[guide$features[i]]]) < min_samples_per_group)){
    message(paste0("Skipping ", guide$features[i], " in ", guide$tissues[i],
                   " (", guide$groups[i], ") -- fewer than ",
                   min_samples_per_group, " samples per feature status have >= ",
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
  # sex/age are kept in colData for reference but deliberately not in the
  # design below -- sample is already a full per-donor blocking factor,
  # so a fixed per-donor covariate like sex or age would be perfectly
  # collinear with it (every one of a donor's pseudobulk rows shares the
  # same sex/age), which would make the model matrix rank-deficient.
  meta_comp <- meta_comp %>%
    mutate(sample2 = paste0(sample, "_", !!sym(guide$features[i]))) %>%
    dplyr::select(-n_spots) %>%
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
    
    keep <- rowSums(counts(dds) >= 10) >= 10 # change these cutoffs as needed
    
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
    
    # Coefficient name is "feature_TRUE_vs_FALSE" regardless of which
    # feature (ptdp/pga) this iteration is on, since the design column
    # was renamed to the generic "feature" above -- not "ptdp_TRUE_vs_FALSE"
    # (which never matches resultsNames(dds), even on a ptdp iteration).
    suppressMessages({
      res_shrunk <- lfcShrink(dds,
                              coef = "feature_TRUE_vs_FALSE",
                              type = "apeglm")
    })
    
    write.csv(res_shrunk,
              file = paste0(comp_results_dir, guide$groups[i], "_lfc_shrunk.csv"))
    
  }, error = function(e){
    message(paste0("Skipping ", guide$features[i], " in ", guide$tissues[i],
                   " (", guide$groups[i], ") -- DESeq2 pipeline failed: ",
                   conditionMessage(e)))
  })
  
}
