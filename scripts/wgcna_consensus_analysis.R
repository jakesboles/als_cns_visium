# Downstream analysis of wgcna_consensus.R's output: (1) a mixed model per
# module testing how UCell module expression differs by group (Control/
# sALS/C9orf72) within each of the 4 anatomical compartments, and (2) gene
# set overrepresentation analysis (ORA) on each module's gene membership.
# Reads directly from results/wgcna_consensus/ (module_scores_ucell.csv,
# modules.csv) and writes back into the same results_dir/data_dir, matching
# deseq2_by_compartment_viz1.R/deseq2_by_compartment_gsea.R's convention of
# a downstream analysis script sharing its parent step's output folder
# rather than getting its own.
#
# Adapted from a pasted draft of als_cns_scrnaseq's per-cell-type WGCNA
# analysis script. Design notes:
# - That draft's `tissue` (2 levels: Motor cortex/Cervical spinal cord) is
#   replaced with `compartment` (4 levels: mcx_GM, mcx_WM, sc_GM, sc_WM) --
#   this project's wgcna_consensus.R runs one consensus network jointly
#   across all 4 anatomical compartments (not per cell type, and not per
#   tissue alone), so "tissue/compartment-specific" here means testing the
#   group effect within each of those 4 compartments, the same way the
#   draft tested it within each of its 2 tissues. compartment isn't saved
#   directly in module_scores_ucell.csv (only its two components, tissue
#   and region, are), so it's rebuilt here the same way wgcna_consensus.R
#   built it in the first place.
# - The draft's `id` (donor) was derived by splitting `orig.ident` on "_".
#   This project already carries a real `sample` (donor id) column, and
#   `code` (Visium section id) plays the role orig.ident played there (the
#   finer-grained replicate nested under donor) -- so
#   (1 | id/orig.ident) becomes (1 | sample/code). NOTE: wgcna_consensus.R's
#   own `scores <- obj@meta.data %>% dplyr::select(...)` does not currently
#   include `sample` -- added there (see that script's diff) since this
#   script can't nest by donor without it. wgcna_consensus.R needs to be
#   rerun for module_scores_ucell.csv to actually carry that column before
#   this script will run.
# - Trimmed the library list to what's actually used -- the draft loaded
#   Seurat/scCustomize/UCell, but this script never touches a Seurat
#   object or calls a UCell function directly, only reads the CSVs
#   wgcna_consensus.R already produced from them.
# - Added write.csv() of both the joint (omnibus) test per module and the
#   full pairwise-comparison/compact-letter-display table -- the draft
#   only ever saved the summary plot, not the numbers behind it, but "for
#   interpretation later" (per the user) means the actual stats need to be
#   on disk too, not just a PNG.
# - Wrapped each module's model fit (task 1) and each module's enrichment
#   test (task 2) in tryCatch(), matching the per-iteration defensive
#   pattern already established in this project's deseq2_by_compartment.R/
#   deseq2_by_pathology.R -- one module failing to converge or enrich
#   shouldn't kill the whole run.
# - Renamed the ORA output files from "_gsea" to "_ora": this project
#   already has a genuinely rank-based GSEA script
#   (deseq2_by_compartment_gsea.R, using clusterProfiler::GSEA() on
#   continuous log2FC), and this task is overrepresentation analysis on a
#   discrete module gene list (clusterProfiler::enricher()) -- a
#   meaningfully different test, so it gets a different name rather than
#   reusing "_gsea" for both.
# - The msigdbr/t2g (MSigDB term-to-gene) setup block is unchanged from
#   both the draft and deseq2_by_compartment_gsea.R -- kept identical
#   rather than factored into a shared file, matching this project's
#   convention of self-contained scripts with no shared helper modules.

suppressMessages({
  library(tidyverse)
  library(lme4)
  library(emmeans)
  library(multcomp)
  library(msigdbr)
  library(clusterProfiler)
})

setwd("/projects/b1169/boles/als_cns_visium")

data_dir <- "data/wgcna_consensus/"
results_dir <- "results/wgcna_consensus/"

# Mixed model of module expression by group, per compartment ----------------

message("Reading in module UCell scores")

scores <- read.csv(paste0(results_dir, "module_scores_ucell.csv"))

scores <- scores %>%
  mutate(compartment = paste0(tissue, "_", region),
         compartment = factor(compartment, levels = sort(unique(compartment))),
         group = factor(group, levels = c("Control", "sALS", "C9orf72")))

cols <- colnames(scores)[str_detect(colnames(scores), "_UCell_kNN$")]
color <- str_remove_all(cols, "_UCell_kNN")

fit <- vector("list", length(cols))
emm <- vector("list", length(cols))
stats <- vector("list", length(cols))
names(fit) <- names(emm) <- names(stats) <- color

for (i in seq_along(cols)){

  message(color[i])

  df <- scores %>%
    dplyr::rename(active_col = all_of(cols[i]))

  # One module failing to converge (or erroring outright, e.g. a module
  # with near-constant scores in some compartment) shouldn't kill every
  # other module's fit.
  tryCatch({

    fit[[i]] <- lme4::lmer(active_col ~ group * compartment + (1 | sample/code),
                           data = df)

    suppressMessages({
      jt <- joint_tests(fit[[i]])
      print(jt)

      emm[[i]] <- emmeans(fit[[i]], pairwise ~ group | compartment)
    })

    write.csv(as.data.frame(jt),
              file = paste0(results_dir, color[i], "_lmer_joint_tests.csv"),
              row.names = F)

    stats[[i]] <- multcomp::cld(emm[[i]], Letters = letters) %>%
      mutate(.group = str_remove_all(.group, " ")) %>%
      mutate(module = color[i])

  }, error = function(e){
    message(paste0("Skipping ", color[i], " -- lmer/emmeans failed: ",
                   conditionMessage(e)))
  })

}

message("Saving model objects and comparison table")

saveRDS(fit,
        file = paste0(data_dir, "lmer_fits.rds"))
saveRDS(emm,
        file = paste0(data_dir, "lmer_emm.rds"))

stats_df <- list_rbind(compact(stats))

write.csv(stats_df,
          file = paste0(results_dir, "lmer_group_comparisons.csv"),
          row.names = F)

message("Plotting module expression by group and compartment")

p <- stats_df %>%
  ggplot(aes(x = compartment,
             y = emmean)) +
  geom_crossbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL,
                    fill = group),
                position = position_dodge(width = 1)) +
  geom_text(aes(label = .group,
                y = (asymp.UCL + emmean) / 2,
                group = group),
            position = position_dodge(width = 1)) +
  scale_fill_manual(values = c("dodgerblue1", "magenta1", "chartreuse1")) +
  facet_wrap(. ~ module,
             scales = "free_y") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        strip.text = element_text(color = "white", face = "bold"),
        strip.background = element_rect(fill = "black"))
ggsave(p,
       filename = paste0(results_dir, "lmer_stats_overview.png"),
       units = "in", dpi = 600,
       height = length(cols) * 0.9,
       width = length(cols) * 1.2)

# Gene set overrepresentation analysis on module membership ------------------

message("Loading MSigDB gene sets")

m_t2g <- msigdbr(species = "Homo sapiens",
                 category = "C2")
m_t2g <- m_t2g %>%
  filter(gs_subcollection %in% c("CP:BIOCARTA", "CP:KEGG_LEGACY", "CP:REACTOME", "CP:WIKIPATHWAYS", "CP:PID")) %>%
  dplyr::select(c(gs_name, gene_symbol))

go_t2g <- msigdbr(species = "Homo sapiens",
                  category = "C5")
go_t2g <- go_t2g %>%
  dplyr::select(c(gs_name, gene_symbol))

t2g <- rbind(go_t2g, m_t2g)

message("Reading in modules")

modules <- read.csv(paste0(results_dir, "modules.csv"))

background <- modules$gene_name %>% unique()

module_names <- unique(modules$module)
module_names <- module_names[module_names != "grey"]

for (j in seq_along(module_names)){

  message(paste0(str_to_title(module_names[j]), " module"))

  genes <- modules %>%
    filter(module == module_names[j]) %>%
    pull(gene_name)

  # A module could plausibly fail enricher() outright (e.g. too few genes
  # after intersecting with the term database) -- logged and skipped
  # rather than killing every other module's enrichment test.
  tryCatch({

    em <- enricher(genes,
                   TERM2GENE = t2g,
                   universe = background)

    write.csv(as.data.frame(em@result),
              file = paste0(results_dir, module_names[j], "_ora.csv"))

    if (nrow(em@result) > 0 && min(em@result$p.adjust) < 0.05){
      p <- dotplot(em,
                   showCategory = 20,
                   x = "FoldEnrichment",
                   color = "p.adjust",
                   size = "GeneRatio")
      ggsave(p,
             filename = paste0(results_dir, module_names[j], "_ora.png"),
             units = "in", dpi = 600,
             height = 10, width = 8)
    }

  }, error = function(e){
    message(paste0("Skipping ", module_names[j], " -- enricher() failed: ",
                   conditionMessage(e)))
  })

}

message("Done")
