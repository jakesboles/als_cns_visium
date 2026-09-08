library(tidyverse)
library(ggrepel)
library(scCustomize)
library(ComplexUpset)
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
ggsave(filename = paste0(results_dir, "/degs_stacked_bars.png"),
       units = "in", dpi = 600,
       height = 5, width = 4)

# Fold-change scatter between any two DESeq2 comparisons --------------------
# Modeled on als_cns_scrnaseq/r_scripts/deseq_viz2.R's plot_fc_scatter() --
# a function rather than a fixed block, so any two (compartment, contrast)
# comparisons can be plotted against each other without editing the body
# each time. Uses the same raw (non-LFC-shrunk) results CSVs and "X" gene
# column as the bar chart above, and the same p_thresh/lfc_thresh already
# defined up top, so "DEG" means the same thing in every plot this script
# makes. Only genes significant in at least one of the two comparisons are
# plotted (same as the sibling repo's version), not the whole transcriptome
# background.
#
# Departure from the sibling repo's version: the top label_n genes (by
# distance from the origin -- see below) are labeled directly on the plot
# via ggrepel, colored to match their sig_group so a labeled gene's color
# still says which comparison(s) it's significant in.

plot_fc_scatter <- function(compartment_1, contrast_1, compartment_2, contrast_2,
                            label_1 = paste(compartment_1, contrast_1),
                            label_2 = paste(compartment_2, contrast_2),
                            label_n = 20){

  path_1 <- paste0(results_dir, "/", compartment_1, "/", contrast_1, ".csv")
  path_2 <- paste0(results_dir, "/", compartment_2, "/", contrast_2, ".csv")

  missing <- c(path_1, path_2)[!file.exists(c(path_1, path_2))]
  if (length(missing) > 0){
    stop(paste0("Missing DESeq2 results file(s): ",
                paste(missing, collapse = ", ")))
  }

  res_1 <- read.csv(path_1) %>%
    dplyr::select(X, log2FoldChange, padj)
  res_2 <- read.csv(path_2) %>%
    dplyr::select(X, log2FoldChange, padj)

  df <- inner_join(res_1, res_2, by = "X", suffix = c("_1", "_2"))

  df <- df %>%
    mutate(sig_1 = !is.na(padj_1) & padj_1 < p_thresh & abs(log2FoldChange_1) > lfc_thresh,
           sig_2 = !is.na(padj_2) & padj_2 < p_thresh & abs(log2FoldChange_2) > lfc_thresh,
           sig_group = case_when(
             sig_1 & sig_2 ~ "Both",
             sig_1 & !sig_2 ~ label_1,
             !sig_1 & sig_2 ~ label_2,
             TRUE ~ "Neither"
           )) %>%
    filter(sig_group != "Neither")

  # Ranked by Euclidean distance from the origin, not by either axis alone
  # -- that rewards a gene for being an outlier on EITHER comparison (or
  # both), rather than always favoring whichever comparison happens to have
  # larger fold changes overall.
  label_df <- df %>%
    mutate(dist = sqrt(log2FoldChange_1^2 + log2FoldChange_2^2)) %>%
    slice_max(dist, n = label_n, with_ties = F)

  ggplot(df, aes(x = log2FoldChange_1, y = log2FoldChange_2, color = sig_group)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_text_repel(data = label_df,
                    aes(label = X, color = sig_group),
                    size = 3, show.legend = F, max.overlaps = Inf) +
    labs(x = paste0("log2FC (", label_1, ")"),
        y = paste0("log2FC (", label_2, ")"),
        color = "Significant in") +
    theme_bw(base_size = 12) +
    theme(axis.text = element_text(color = "black"))
}

# Returns the ggplot object rather than saving it -- the specific pair of
# comparisons (and label_n) varies by call, so there's no single sensible
# fixed output filename the way the bar chart/upset plot have one. Example:
# p <- plot_fc_scatter("mcx_gm", "C9orf72_vs_Control", "sc_gm", "C9orf72_vs_Control",
#                      label_1 = "C9-ALS motor cortex GM", label_2 = "C9-ALS spinal cord GM",
#                      label_n = 20)
# p
# ggsave(p, filename = paste0(results_dir, "/fc_scatter_c9_mcxgm_vs_scgm.png"),
#        units = "in", dpi = 600, height = 6, width = 7)

# Upset plot of DEGs across all 8 sets ---------------------------------------
# 4 compartments (mcx GM, mcx WM, sc GM, sc WM) x 2 contrasts (sALS vs
# Control, C9orf72 vs Control) = 8 DEG sets, all pairwise/higher-order
# overlaps shown. Ported from als_cns_scrnaseq/r_scripts/deseq_viz2.R's
# ComplexUpset-based microglia upset plot (file_paths -> missing-file check
# -> per-set DEG lists -> intersect_df membership matrix -> upset() with a
# custom intersection_matrix/theme), generalized from that script's 4 sets
# (2 tissues x 2 groups) to this project's 8 (2 tissues x 2 anatomical
# compartments x 2 groups).
#
# The sibling script's version also colors specific "shared" intersections
# (a highlight_group fill + matching upset_query() calls on the matrix) --
# deliberately dropped here rather than mechanically extended, since that
# logic was a fully-specified AND/NOT condition across exactly 2 tissues x
# 2 groups (4 sets), and there's no single obvious way to generalize it to
# 4 compartments x 2 groups (8 sets): a literal extension (e.g. requiring
# significance in all 4 compartments for one contrast) is a much stricter
# condition than the original 2-tissue version and would likely highlight
# almost nothing. Everything else about the sibling script's structure and
# styling (labeller, matrix geom, theme) is kept as close to the original
# as possible.

file_paths <- c(
  sALS_mcx_gm = paste0(results_dir, "/mcx_gm/sALS_vs_Control.csv"),
  `C9-ALS_mcx_gm` = paste0(results_dir, "/mcx_gm/C9orf72_vs_Control.csv"),
  sALS_mcx_wm = paste0(results_dir, "/mcx_wm/sALS_vs_Control.csv"),
  `C9-ALS_mcx_wm` = paste0(results_dir, "/mcx_wm/C9orf72_vs_Control.csv"),
  sALS_sc_gm = paste0(results_dir, "/sc_gm/sALS_vs_Control.csv"),
  `C9-ALS_sc_gm` = paste0(results_dir, "/sc_gm/C9orf72_vs_Control.csv"),
  sALS_sc_wm = paste0(results_dir, "/sc_wm/sALS_vs_Control.csv"),
  `C9-ALS_sc_wm` = paste0(results_dir, "/sc_wm/C9orf72_vs_Control.csv")
)

missing <- file_paths[!file.exists(file_paths)]
if (length(missing) > 0){
  stop(paste0("Missing DESeq2 results file(s): ",
              paste(missing, collapse = ", "),
              " -- check whether deseq2_by_compartment.R's abundance filter ",
              "skipped this compartment."))
}

list_names <- names(file_paths)

list <- lapply(file_paths, read.csv)
names(list) <- list_names

# Same DEG definition (p_thresh/lfc_thresh, defined at the top of this
# script) as the bar chart and fold-change scatter above, not the sibling
# script's own hardcoded 0.05/log2(1.5) -- this script already
# parameterizes both, so every plot it makes should stay in sync if either
# threshold ever changes.
for (i in seq_along(list)){
  list[[i]] <- list[[i]] %>%
    filter(padj < p_thresh &
             abs(log2FoldChange) > lfc_thresh) %>%
    pull(X)
}

genes <- unique(list_c(list))

intersect_df <- data.frame(gene = genes)
for (i in seq_along(list)){

  col <- list_names[i]

  intersect_df <- intersect_df %>%
    mutate(x = if_else(gene %in% list[[i]], T, F)) %>%
    dplyr::rename(!!sym(col) := "x")
}

png(file = paste0(results_dir, "/degs_upset.png"),
    height = 7, width = 12,
    units = "in", res = 600)

upset(intersect_df, list_names,
     set_sizes = F,
     stripes = "white",
     wrap = T,
     encode_sets = F,
     height_ratio = 0.9,
     labeller = ggplot2::as_labeller(c(
       "sALS_mcx_gm" = "sALS\nMotor cortex GM",
       "C9-ALS_mcx_gm" = "C9orf72-ALS\nMotor cortex GM",
       "sALS_mcx_wm" = "sALS\nMotor cortex WM",
       "C9-ALS_mcx_wm" = "C9orf72-ALS\nMotor cortex WM",
       "sALS_sc_gm" = "sALS\nSpinal cord GM",
       "C9-ALS_sc_gm" = "C9orf72-ALS\nSpinal cord GM",
       "sALS_sc_wm" = "sALS\nSpinal cord WM",
       "C9-ALS_sc_wm" = "C9orf72-ALS\nSpinal cord WM"
     )),
     matrix = (intersection_matrix(
       geom = geom_point(shape = 19, size = 10),
       segment = geom_segment(linewidth = 1.5),
       outline_color = list(active = "white", inactive = "white")
     )),
     base_annotations = list(
       'Intersection size' = intersection_size(text = list(size = 0)) +
         scale_y_continuous(expand = c(0, 0)) +
         ylab("# DEGs")
     ),
     theme = upset_modify_themes(
       list(
         'Intersection size' = theme(axis.text = element_text(color = "black", size = 16),
                                     axis.title = element_text(size = 20),
                                     axis.ticks.y = element_line(),
                                     panel.border = element_rect(color = "black", fill = "transparent"),
                                     panel.grid = element_line(color = "gray70")),
         'intersections_matrix' = theme(axis.text = element_text(color = "black", size = 16),
                                        axis.title = element_blank())
       )
     )
) +
  ggtitle("DEGs by compartment and contrast") +
  theme(plot.title = element_text(hjust = 0.5, size = 20, face = "plain"))

dev.off()
