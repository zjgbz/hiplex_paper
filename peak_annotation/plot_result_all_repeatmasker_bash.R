# ==============================================================================
# RepeatMasker Annotation Visualization with Clustering
# ==============================================================================
# This script visualizes RepeatMasker annotations for genomic peaks using a
# pre-computed clustering order from CCRE analysis. It creates ordered bar
# plots showing the distribution of repetitive elements across different
# histone modification and transcription factor combinations.
# ==============================================================================

# Load Required Libraries ======================================================
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(stats)
library(ggrepel)
library(glue)
library(dendextend)

# Source Required Scripts ======================================================
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/plotUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")

# Configuration ================================================================
filtered <- TRUE
# xV <- readRDS("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/data_peak_auc_0.05_extend_window_10000_binomial/V/annotate.rds")
# Parse Command Line Arguments =================================================
args <- commandArgs(trailingOnly = TRUE)
# peak_type <- args[1]
dataPeakAnnoDir <- args[1]
new_order_dir <- args[2]
out_dir <- args[3]
print(peak_type)

# Load RepeatMasker Annotation Data ============================================
# dataPeakAnnoDir <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/V/"
# )
rdsFileV <- paste0(dataPeakAnnoDir, "annotate.rds")
xV <- readRDS(rdsFileV)

# Create Output Directory ======================================================
# out_dir <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/plots/"
# )
dir.create(out_dir, recursive = TRUE)

# Extract Annotation Statistics ================================================
annoV <- lapply(xV, getAnnoStatCCRE)
cat("tag counts: ", length(annoV), "\n")

# Filter Targets if Required ===================================================
if (filtered) {
  filteredTags <- filter_target_pairs()
  annoV <- annoV[filteredTags]
  filterStr <- "filtered_"
  cat("after filtering: ", length(annoV), "\n")
} else {
  filterStr <- ""
}

# Load Pre-computed Clustering Order ===========================================
# Use the clustering order from CCRE analysis for consistency
# new_order_dir <- paste0(
#     glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/plots/"),
#     "CCRE_anno_V_order.rds"
#   )
new_order <- readRDS(
  new_order_dir
)

# Order Annotations by CCRE Clustering =========================================
annoVClusterOrdered <- annoV[new_order]

# Utility Functions ============================================================

#' Extract legend from ggplot object
#'
#' @param my_ggp ggplot object
#' @return grob containing the legend
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

#' Reverse heterotone tag name (swap tag order)
#'
#' @param tag Character string in format "tag1-tag2"
#' @return Character string in format "tag2-tag1"
getRevHeterotone <- function(tag) {
  tag <- as.character(tag)
  tags <- strsplit(tag, "-")[[1]]
  tag1 <- tags[1]
  tag2 <- tags[2]
  return(paste0(tag2, "-", tag1))
}

# Create RepeatMasker Annotation Plot ==========================================

# Convert to data frame
annoV.df <- list_to_dataframe(annoVClusterOrdered)
annoV.df$Feature <- factor(
  annoV.df$Feature,
  levels = repeatMaskerCategoriesOrder
)

# Map target names and set factor levels
annoV.df$.id <- map_target_names(annoV.df$.id, target_pair_mapping_df)
annoV.df$.id <- factor(
  annoV.df$.id,
  levels = rev(map_target_names(new_order, target_pair_mapping_df))
)

# Create bar plot
categoryColumn <- ".id"
gpV <- plotAnnoBar.data.frame(
  annoV.df,
  categoryColumn = categoryColumn,
  colorOption = 4
) +
  ggtitle("Vehicle") +
  scale_x_discrete(labels = TeX)

# Save Results =================================================================

# Save plot as PDF
pdf(
  paste0(out_dir, "repeatMasker_anno_V_clustered.pdf"),
  height = 140,
  width = 20
)
print(gpV)
dev.off()

# Save annotation data as CSV
write.csv(
  annoV.df,
  paste0(out_dir, "repeatMasker_anno_V_clustered.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Analysis complete. Results saved to:", out_dir, "\n")

# ==============================================================================
# End of Script
# ==============================================================================