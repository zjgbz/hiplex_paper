# ==============================================================================
# Peak Annotation Clustering and Visualization
# ==============================================================================
# This script performs hierarchical clustering on peak annotations and creates
# visualizations of CCRE (cis-regulatory element) distributions across different
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
library(latex2exp)
library(glue)
library(dendextend)

# Source Required Scripts ======================================================
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/plotUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")

# Configuration ================================================================
filtered <- TRUE
intermediatePlots <- FALSE

# Parse Command Line Arguments =================================================
args <- commandArgs(trailingOnly = TRUE)
# peak_type <- args[1]
CCREPeakAnnoDir <- args[1]
out_dir <- args[2]
# Load Annotation Data =========================================================
# CCREPeakAnnoDir <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/V/"
# )
CCRErdsFileV <- paste0(CCREPeakAnnoDir, "annotate.rds")
CCRExV <- readRDS(CCRErdsFileV)

# Extract Annotation Statistics ================================================
CCREannoV <- lapply(CCRExV, getAnnoStatCCRE)
cat("tag counts: ", length(CCREannoV), "\n")

# Filter Targets if Required ===================================================
if (filtered) {
  filteredTags <- filter_target_pairs()
  CCREannoV <- CCREannoV[filteredTags]
  filterStr <- "filtered_"
  cat("after filtering: ", length(CCREannoV), "\n")
} else {
  filterStr <- ""
}

# Create Output Directory ======================================================
# out_dir <- glue(
#   "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/plots/"
# )
dir.create(out_dir, recursive = TRUE)

# Perform Hierarchical Clustering ==============================================

# Create annotation matrix
resultMat <- getAnnotationMatrix(CCREannoV, ChIPSeekerCCRECategoriesOrder)

# Calculate distance matrix and perform hierarchical clustering
dist_matrix <- dist(resultMat, method = "euclidean")
hclust_result <- hclust(dist_matrix, method = "complete")

# Convert to dendrogram for visualization
dend <- as.dendrogram(hclust_result)

# Define Biological Clusters ===================================================
# Group 1: H3K4me3-related modifications (promoter/active transcription)
# Group 2: H3K27ac-related modifications (active enhancers/transcription)
# Group 3: H3K27me3-related modifications (polycomb repression)
# Group 4: H3K9me3/me2-related modifications (heterochromatin)

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

# Define colors for the clusters
cluster_colors <- c("red", "blue", "green", "orange")

# Color Dendrogram Labels by Cluster ===========================================

# Initialize all labels as black
label_colors <- rep("black", length(labels(dend)))
names(label_colors) <- labels(dend)

# Assign colors based on cluster membership
for (i in seq_along(clusters)) {
  label_colors[intersect(labels(dend), clusters[[i]])] <- cluster_colors[i]
}

# Apply colors to dendrogram labels
dend <- dendextend::set(dend, "labels_col", label_colors[labels(dend)])

# Get ordering from dendrogram
new_order <- rownames(resultMat)[order.dendrogram(as.dendrogram(hclust_result))]

# Save Dendrogram ==============================================================
pdf(glue("{out_dir}/dend.pdf"), width = 100, height = 10)
plot(
  dend,
  main = "Hierarchical Clustering Dendrogram",
  xlab = "",
  sub = "",
  cex = 0.8
)
dev.off()

# Order Annotations by Clustering ==============================================
CCREannoVClusterOrdered <- CCREannoV[new_order]

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

# Create Annotation Plot =======================================================

# Convert to data frame
CCREannoV.df <- list_to_dataframe(CCREannoVClusterOrdered)
CCREannoV.df$Feature <- factor(
  CCREannoV.df$Feature,
  levels = ChIPSeekerCCRECategoriesOrder
)

# Map target names and set factor levels
CCREannoV.df$.id <- map_target_names(CCREannoV.df$.id, target_pair_mapping_df)
CCREannoV.df$.id <- factor(
  CCREannoV.df$.id,
  levels = rev(map_target_names(new_order, target_pair_mapping_df))
)

# Create bar plot
categoryColumn <- ".id"
gpV <- plotAnnoBar.data.frame(
  CCREannoV.df,
  categoryColumn = categoryColumn,
  colorOption = 1
) +
  ggtitle("Vehicle") +
  scale_x_discrete(labels = TeX)

# Save Results =================================================================

# Save annotation data as CSV
write.csv(
  CCREannoV.df,
  paste0(out_dir, "CCRE_anno_V_clustered.csv"),
  row.names = FALSE,
  quote = FALSE
)

# Save plot as PDF
pdf(
  paste0(out_dir, "CCRE_anno_V_clustered.pdf"),
  height = 140,
  width = 20
)
print(gpV)
dev.off()

# Save ordering for future use
saveRDS(new_order, paste0(out_dir, "CCRE_anno_V_order.rds"))

cat("Analysis complete. Results saved to:", out_dir, "\n")

# ==============================================================================
# End of Script
# ==============================================================================