# ==============================================================================
# Integrated Peak Annotation and DNA Methylation Analysis
# ==============================================================================
# This script creates a comprehensive visualization combining:
# 1. CCRE (cis-regulatory element) annotations
# 2. RepeatMasker (repetitive element) annotations
# 3. DNA methylation levels (M-values)
# across different histone modification and transcription factor combinations.
# ==============================================================================

# Load Required Libraries ======================================================
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(glue)
library(grid)

# Source Required Scripts ======================================================
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/dna_methylation/calculate_methylation_percentage.R")

# Parse Command Line Arguments =================================================
args <- commandArgs(trailingOnly = TRUE)
# peak_type <- args[1]
peak_dir <- args[1]
repeatmasker_rds_file <- args[2]
ccre_rds_file <- args[3]
out_dir <- args[4]
# Define Target Tags ===========================================================
# Select histone marks and transcription factors for analysis
tags <- c(
  "H3K4me1", "H3K27ac", "H3K9me3", 
  "H3K27me3", "H3K4me3", "POLR2AphosphoS2"
)

# Generate and filter target pairs
target_pair_list <- target_pair_generation(tags)
target_pair_list <- filter_target_pairs(0.25, target_pair_list)

# Define S4 Class for CCRE Annotations =========================================
setClass(
  "csCCREAnno",
  representation = representation(
    annoStat = "data.frame",
    peakNum = "numeric"
  )
)

# Define Biological Clusters ===================================================
# Group target pairs based on biological function and chromatin state
clusters <- list(
  # Group 1: Promoter/active transcription markers
  c(
    "H3K36me3-H3K4me3",
    "H3K4me3-YY1",
    "H3K4me3-Myc",
    "CTCF-H3K4me3",
    "H3K4me3-POLR2AphosphoS2",
    "H3K27ac-H3K4me3",
    "H3K4me3-NRF1",
    "H3K14ac-H3K4me3"
  ),
  
  # Group 2: Active enhancer/transcription markers
  c(
    "H3K27ac-H3K27ac",
    "H3K27ac-Myc",
    "H3K14ac-H3K9ac",
    "H3K9ac-Myc",
    "POLR2AphosphoS2-POLR2AphosphoS2",
    "H3K14ac-POLR2AphosphoS2",
    "H3K4me3-H3K4me3",
    "H3K27me3-H3K4me3"
  ),
  
  # Group 3: Polycomb repression markers
  c(
    "MLL1_KMT2A-MLL4_MLL2_KMT2B",
    "H3K36me3-H3S10ph",
    "CDK8-H2A_XS139ph",
    "PIM1-SETD2",
    "H3K27me3-cFos",
    "H3K27me3-NRF1",
    "H3K27me3-MSK1",
    "EHMT2-H3K27me3"
  ),
  
  # Group 4: Heterochromatin markers
  c(
    "H3K9me3-MSK2",
    "H3K9me3-Myc",
    "EHMT2-H3K9me3",
    "H2A_XS139ph-H3K9me3",
    "H3K9me3-H3K9me3",
    "H3K9me2-H3K9me3",
    "H3K9me2-MSK2",
    "H3K9me2-YY1"
  )
)

# Perform RepeatMasker-Based Clustering ========================================

# Load RepeatMasker annotations
# repeatmasker_rds_file <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/V/annotate.rds"
# )
repeatmasker_annotation_result <- readRDS(repeatmasker_rds_file)

# Cluster within each biological group
new_clusters <- list()

for (cluster_group in seq_along(clusters)) {
  cluster <- clusters[[cluster_group]]
  
  # Extract annotations for current cluster
  repeatmasker_annotation_result_i <- repeatmasker_annotation_result[cluster]
  anno <- lapply(repeatmasker_annotation_result_i, getAnnoStatCCRE)
  anno_mat <- getAnnotationMatrix(anno, repeatMaskerCategoriesOrder)
  
  # Perform hierarchical clustering
  h <- ComplexHeatmap::Heatmap(
    anno_mat,
    cluster_columns = FALSE,
    column_names_gp = grid::gpar(fontsize = 8),
    row_names_gp = grid::gpar(fontsize = 0.5)
  )
  ht <- draw(h)
  
  # Store ordered target pairs
  new_clusters[[cluster_group]] <- names(anno)[row_order(ht)]
}

# Combine all clusters in order
target_pair_list <- unlist(new_clusters)

# Create CCRE Annotation Plot ==================================================

# Load CCRE annotations
# ccre_rds_file <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/V/annotate.rds"
# )
ccre_annotation_result <- readRDS(ccre_rds_file)

# Order and prepare data
ccre_annotation_result <- ccre_annotation_result[target_pair_list]
anno <- lapply(ccre_annotation_result, getAnnoStatCCRE)
anno <- anno[target_pair_list]

# Convert to data frame
anno.df <- list_to_dataframe(anno)
anno.df$.id <- map_target_names(anno.df$.id, target_pair_mapping_df)
anno.df$.id <- factor(
  anno.df$.id,
  levels = rev(map_target_names(target_pair_list, target_pair_mapping_df))
)
anno.df$Feature <- factor(
  anno.df$Feature,
  levels = ChIPSeekerCCRECategoriesOrder
)

# Create CCRE plot
categoryColumn <- ".id"
p_ccre <- plotAnnoBar.data.frame(
  anno.df,
  categoryColumn = categoryColumn,
  colorOption = 1
) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 20)
  ) +
  ggtitle("Cis-regulatory Elements Annotation")

# Create RepeatMasker Annotation Plot ==========================================

# Order and prepare RepeatMasker data
repeatmasker_annotation_result <- repeatmasker_annotation_result[target_pair_list]
anno_repeat <- lapply(repeatmasker_annotation_result, getAnnoStatCCRE)

# Convert to data frame
anno.df <- list_to_dataframe(anno_repeat)
anno.df$Feature <- factor(
  anno.df$Feature,
  levels = c(repeatMaskerFeatures, "other")
)
anno.df$.id <- map_target_names(anno.df$.id, target_pair_mapping_df)
anno.df$.id <- factor(
  anno.df$.id,
  levels = rev(map_target_names(target_pair_list, target_pair_mapping_df))
)

# Create RepeatMasker plot
p_repeat <- plotAnnoBar.data.frame(
  anno.df,
  categoryColumn = categoryColumn,
  colorOption = 4
) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(size = 20),
    axis.text.y = element_blank()
  ) +
  ggtitle("Repetitive Elements Annotation")

# Calculate and Plot DNA Methylation Levels ====================================

# Initialize storage vectors
result_tags <- c()
result_methyl_levels <- c()

# Process each target pair
# peak_dir <- glue(
#   "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/peak/{peak_type}/"
# )

for (tag_pair in target_pair_list) {
  peak_file <- paste0(peak_dir, "V", "/", tag_pair, ".stringent.bed")
  peak <- readPeakFile(peak_file)
  
  # Calculate M-values (log2 ratio of methylated to unmethylated)
  result_grange <- calculate_M_values(peak, ignore_no_cpg_peaks = FALSE)
  result_grange$tag_pair <- tag_pair
  
  # Store results
  result_tags <- append(result_tags, result_grange$tag_pair)
  result_methyl_levels <- append(result_methyl_levels, result_grange$M_values)
}

# Create methylation data frame
plt_data <- data.frame(
  target_pair = result_tags,
  methylation = result_methyl_levels
)

# Save methylation data
# dir.create(
#   glue("/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/{peak_type}/"),
#   recursive = TRUE
# )
# saveRDS(
#   plt_data,
#   glue("/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/{peak_type}/plt_data_V_summary.rds")
# )
dir.create(
  out_dir,
  recursive = TRUE
)
saveRDS(
  plt_data,
  glue("{out_dir}}/plt_data_V_summary.rds")
)

# Filter and order methylation data
plt_data <- plt_data[plt_data$target_pair %in% target_pair_list, ]
plt_data$target_pair <- factor(plt_data$target_pair, levels = rev(target_pair_list))

# Create methylation violin plot
p_methyl <- ggplot(plt_data, aes(x = target_pair, y = methylation)) +
  geom_violin(fill = "#6CB0D6", alpha = 0.5) +
  ggtitle("Methylation Level Distribution") +
  theme(
    plot.title = element_text(size = 20),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_rect(fill = "white", colour = "black"),
    legend.position = "none",
    axis.title.y = element_blank()
  ) +
  coord_flip() +
  ylab("M value")

# Combine All Plots and Save ==================================================

# Create output directory
# dir.create(
#   glue("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/{peak_type}/"),
#   recursive = TRUE
# )

# Create combined plot
# pdf(
#   glue("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/{peak_type}/summary.pdf"),
#   width = 30,
#   height = 15
# )
pdf(
  glue("{out_dir}/summary.pdf"),
  width = 30,
  height = 15
)
plot_grid(
  p_ccre,
  p_repeat,
  p_methyl,
  align = "h",
  axis = "bt",
  ncol = 3,
  rel_widths = c(2.5, 2, 1)
)
dev.off()

cat("Analysis complete. Summary plot saved.\n")

# ==============================================================================
# End of Script
# ==============================================================================