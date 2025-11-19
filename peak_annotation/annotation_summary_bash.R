library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(glue)

source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/dna_methylation/calculate_methylation_percentage.R")
args <- commandArgs(trailingOnly = TRUE)
peak_type <- args[1]

tags <- c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
          "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP",
          "EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1",
          "EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "POLR2AphosphoS2", "cJun", "cFos", "Max", "Myc",
          "USF1", "USF2", "NRF1", "YY1", "H3K9ac")
tags = c("H3K4me1", "H3K27ac", "H3K9me3", "H3K27me3", "H3K4me3", "POLR2AphosphoS2")

target_pair_list=target_pair_generation(tags)
target_pair_list = filter_target_pairs(0.25, target_pair_list)
setClass("csCCREAnno",
         representation=representation(
           annoStat="data.frame",
           peakNum="numeric"
         ))
cluster_start <- c(1, 30, 51, 133, 349, 393, length(target_pair_list))
clusters <- list(
  c("H3K36me3-H3K4me3",
      "H3K4me3-YY1",
      "H3K4me3-Myc",
      "CTCF-H3K4me3",
      "H3K4me3-POLR2AphosphoS2",
      "H3K27ac-H3K4me3",
      "H3K4me3-NRF1",
      "H3K14ac-H3K4me3"),
  c("H3K27ac-H3K27ac",
      "H3K27ac-Myc",
      "H3K14ac-H3K9ac",
      "H3K9ac-Myc",
      "POLR2AphosphoS2-POLR2AphosphoS2",
      "H3K14ac-POLR2AphosphoS2",
      "H3K4me3-H3K4me3",
      "H3K27me3-H3K4me3"
  ),
  c("MLL1_KMT2A-MLL4_MLL2_KMT2B", 
      "H3K36me3-H3S10ph",
      "CDK8-H2A_XS139ph",
      "PIM1-SETD2",
      "H3K27me3-cFos", 
      "H3K27me3-NRF1", 
      "H3K27me3-MSK1",
      "EHMT2-H3K27me3"
  ),
  c("H3K9me3-MSK2", 
      "H3K9me3-Myc",
      "EHMT2-H3K9me3", 
      "H2A_XS139ph-H3K9me3", 
      "H3K9me3-H3K9me3",
      "H3K9me2-H3K9me3", "H3K9me2-MSK2",
      "H3K9me2-YY1"
  )
)
## CCRE annotation
ccre_rds_file <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/V/annotate.rds")
ccre_annotation_result <- readRDS(ccre_rds_file)

ccre_annotation_result <- ccre_annotation_result[target_pair_list]
anno <- lapply(ccre_annotation_result, getAnnoStatCCRE)
anno_mat <- getAnnotationMatrix(anno, ChIPSeekerCCRECategoriesOrder)

h <- ComplexHeatmap::Heatmap(anno_mat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
                             row_names_gp = grid::gpar(fontsize = 0.5))
ht = draw(h)
target_pair_list <- names(anno)[row_order(ht)]
target_pair_list <- unlist(clusters)

repeatmasker_rds_file <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/V/annotate.rds")
repeatmasker_annotation_result <- readRDS(repeatmasker_rds_file)
## RepeatMasker annotation
new_clusters <- list()
for (cluster_group in seq_along(clusters)) {
  cluster <- clusters[[cluster_group]]
  repeatmasker_annotation_result_i <- repeatmasker_annotation_result[cluster]
  anno <- lapply(repeatmasker_annotation_result_i, getAnnoStatCCRE)
  anno_mat <- getAnnotationMatrix(anno, repeatMaskerCategoriesOrder)
  # mergeDf <- merge(mergeDfAnno, mergeDfCCRE, by = "Row.names", all = TRUE)
  
  h <- ComplexHeatmap::Heatmap(anno_mat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
                               row_names_gp = grid::gpar(fontsize = 0.5))
  ht = draw(h)
  new_clusters[[cluster_group]] <- names(anno)[row_order(ht)]
}
target_pair_list <- unlist(new_clusters)

# ccre_rds_file <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/V/annotate.rds"
ccre_annotation_result <- readRDS(ccre_rds_file)
ccre_annotation_result <- ccre_annotation_result[target_pair_list]
anno <- lapply(ccre_annotation_result, getAnnoStatCCRE)
anno <- anno[target_pair_list]
anno.df <- list_to_dataframe(anno)
anno.df$.id <- map_target_names(anno.df$.id, target_pair_mapping_df)
anno.df$.id <- factor(anno.df$.id, levels = rev(map_target_names(target_pair_list, target_pair_mapping_df)))
anno.df$Feature <- factor(anno.df$Feature, levels = c(ChIPSeekerCCRECategoriesOrder))
categoryColumn <- ".id"
p_ccre <- plotAnnoBar.data.frame(anno.df, categoryColumn=categoryColumn, colorOption = 1) + theme(legend.position="bottom", plot.title = element_text(size = 20,)) + ggtitle("Cis-regulatory Elements Annotation")


repeatmasker_annotation_result <- repeatmasker_annotation_result[target_pair_list]
anno_repeat <- lapply(repeatmasker_annotation_result, getAnnoStatCCRE)
anno.df <- list_to_dataframe(anno_repeat)
anno.df$Feature <- factor(anno.df$Feature, levels = c(repeatMaskerFeatures, "other"))
anno.df$.id <- map_target_names(anno.df$.id, target_pair_mapping_df)
anno.df$.id <- factor(anno.df$.id, levels = rev(map_target_names(target_pair_list, target_pair_mapping_df)))
categoryColumn <- ".id"

p_repeat <- plotAnnoBar.data.frame(anno.df, categoryColumn=categoryColumn, colorOption = 4) + theme(legend.position="bottom", plot.title = element_text(size = 20,), axis.text.y=element_blank()) + ggtitle("Repetitive Elements Annotation")

if (FALSE) {
  ## peak counts bar chart
  peak_counts <- c()
  for (tag_pair in names(ccre_annotation_result)) {
    peak_counts <- append(peak_counts, ccre_annotation_result[[tag_pair]]@peakNum)
  }
  names(peak_counts) <- names(ccre_annotation_result)
  
  peak_counts <- peak_counts[target_pair_list]
  
  peak_counts_df <- data.frame(tags=names(peak_counts), peakCounts=peak_counts, logPeakCounts=log2(peak_counts))
  peak_counts_df$tags <- factor(peak_counts_df$tags, levels = rev(peak_counts_df$tags))
  p_cnts <- ggplot(data=peak_counts_df, aes(x=tags, y=logPeakCounts, fill=tags)) +
    geom_bar(stat="identity") + coord_flip() + ggtitle("peak num")+ theme(plot.title = element_text(size = 20,), 
                                                                          axis.ticks.y=element_blank(),
                                                                          axis.text=element_text(size=15), 
                                                                          panel.border = element_blank(), 
                                                                          panel.grid.major = element_blank(),
                                                                          panel.grid.minor = element_blank(), 
                                                                          axis.line = element_line(colour = "black"), 
                                                                          panel.background = element_rect(fill = "white", 
                                                                                                          colour = "black"), 
                                                                          legend.position="none",
                                                                          axis.title.y=element_blank()) + geom_text(aes(label = peakCounts), hjust = 2, colour = "white") + ggtitle("Peak Counts") + ylab("Log(Peak Count)")
}




peak_dir <- glue("/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/peak/{peak_type}/")
result_tags <- c()
result_methyl_levels <- c()
for (tag_pair in target_pair_list) {
  peak_file <- paste0(peak_dir, "V", "/", tag_pair, ".stringent.bed")
  peak <- readPeakFile(peak_file)
  result_grange <- calculate_M_values(peak, ignore_no_cpg_peaks=FALSE)
  result_grange$tag_pair <- tag_pair
  result_tags <- append(result_tags, result_grange$tag_pair)
  result_methyl_levels <- append(result_methyl_levels, result_grange$M_values)
}
result_methyl_levels = result_methyl_levels 
plt_data <- data.frame(target_pair=result_tags, methylation=result_methyl_levels)
dir.create(glue("/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/{peak_type}/"), recursive=TRUE)
saveRDS(plt_data, glue("/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/{peak_type}/plt_data_V_summary.rds"))


# plt_data <- readRDS("/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/plt_data_V_summary.rds")
plt_data <- plt_data[plt_data$target_pair %in% target_pair_list,]
# saveRDS(plt_data, "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/avg_methylation/plt_data.rds")
plt_data$target_pair <- factor(plt_data$target_pair, levels = rev(target_pair_list))
p_methyl <- ggplot(plt_data, aes(x=target_pair, y=methylation)) + 
  geom_violin(fill="#6CB0D6", alpha=0.5) + ggtitle("Methylation Level Distribution") + theme(plot.title = element_text(size = 20), 
                                                        axis.ticks.y=element_blank(),
                                                        axis.text.x=element_text(size=15), 
                                                        axis.text.y=element_blank(),
                                                        panel.border = element_blank(), 
                                                        panel.grid.major = element_blank(),
                                                        panel.grid.minor = element_blank(), 
                                                        axis.line = element_line(colour = "black"), 
                                                        panel.background = element_rect(fill = "white", 
                                                                                        colour = "black"),
                                                        legend.position="none",
                                                        axis.title.y=element_blank()) + coord_flip() + ylab("M value")
dir.create(glue("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/{peak_type}/"), recursive=TRUE)
pdf(glue("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/{peak_type}/summary.pdf"), width = 30, height = 15)
plot_grid(p_ccre, p_repeat, p_methyl, align = "h", axis = "bt", ncol = 3, rel_widths = c(2.5,2,1))
dev.off()


