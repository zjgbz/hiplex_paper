library(ChIPseeker)
library(ComplexHeatmap)
library(glue)
library(latex2exp)

source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")


args <- commandArgs(trailingOnly = TRUE)
peak_dir <- args[1]
out_dir <- args[2]
samples <- c("V")

# tags <- c("H3K36me3", "H3K4me1", "H3K27ac", "H3K9ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
#           "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "POLR2AphosphoS2", "SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP",
#           "EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1",
#           "EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "cJun", "cFos", "Max", "Myc",
#           "USF1", "USF2", "NRF1", "YY1")
tags <- c("H3K4me1", "H3K4me3", "H3K9ac", "H3K9me2", "H3K9me3", "H3S10ph", "H3K14ac", "H3K27ac", "H3K27me3", "H3K36me3", "H3K79me3", "H2A_XS139ph",
          "POLR2AphosphoS2", "AURORA_Aurora_B", "CBP_CREBBP", "CDK8", "EHMT1", "EHMT2", "EP300", "EZH2", 
          "MLL1_KMT2A", "MLL4_MLL2_KMT2B", "MSK1", "MSK2", "PIM1", "SETD2", "SuVar39_SUV39H1",
          "CTCF", "Max", "Myc", "NRF1", "USF1", "USF2", "YY1", "cFos", "cJun")
splits <- c(rep(1, 12), 2, rep(3, 14), rep(4, 9))
tags_peak_num_mat <- matrix(0, nrow = length(tags), ncol = length(tags))
rownames(tags_peak_num_mat) <- tags
colnames(tags_peak_num_mat) <- tags



all_target_pair_list <- target_pair_generation(tags)
target_pair_list = filter_target_pairs(0.25)

filtered_target_pair_list <- all_target_pair_list[!all_target_pair_list %in% target_pair_list]


# peak_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/peak/data_peak_auc_0.05"
for (sample in samples) {
  for (target_pair in target_pair_list) {
    peak_file <- glue("{peak_dir}/{sample}/{target_pair}.stringent.bed")
    if (file.exists(peak_file) && length(peak_file) > 0 && file.size(peak_file) > 0) {
      peak <- ChIPseeker::readPeakFile(peak_file, as = "GRanges")   # Use ChIPseeker to read peak files
      peak_count <- length(peak)
    } else {
      peak_count <- 0
    }
    tag1 <- strsplit(target_pair, "-")[[1]][1]
    tag2 <- strsplit(target_pair, "-")[[1]][2]
    tags_peak_num_mat[tag1, tag2] <- length(peak)
    tags_peak_num_mat[tag2, tag1] <- length(peak)
  }
}
print(peak_file)
for (sample in samples) {
  for (target_pair in filtered_target_pair_list) {
    tag1 <- strsplit(target_pair, "-")[[1]][1]
    tag2 <- strsplit(target_pair, "-")[[1]][2]
    tags_peak_num_mat[tag1, tag2] <- NA
    tags_peak_num_mat[tag2, tag1] <- NA
  }
}
rownames(tags_peak_num_mat) <- map_target_names(rownames(tags_peak_num_mat), target_pair_mapping_df)
colnames(tags_peak_num_mat) <- map_target_names(colnames(tags_peak_num_mat), target_pair_mapping_df)

na_legend <- Legend(labels = "Filtered Targets", legend_gp = gpar(fill = "grey90"), title = "")
group_col <- c("#8ECFC9",
  "#FFBE7A",
  "#FA7F6F",
  "#82B0D2",
  "#BEB8DC")
row_annotation = rowAnnotation(
  foo = anno_block(gp = gpar(fill = group_col), 
                   labels = c("Histone Modification", "", "Writer", "Transcription Factor"),
                   labels_gp = gpar(col = "white", fontsize = 11))
)
col_annotation = HeatmapAnnotation(
  foo = anno_block(gp = gpar(fill = group_col), 
                   labels = c("Histone Modification", "", "Writer", "Transcription Factor"),
                   labels_gp = gpar(col = "white", fontsize = 11))
)
ht_opt$ROW_ANNO_PADDING =unit(0.4, "cm")
ht_opt$COLUMN_ANNO_PADDING =unit(0.4, "cm")
# rownames(tags_peak_num_mat)[12] <- "$\\gamma$H2AX"
pdf(glue("{out_dir}/V.pdf"), width = 12, height = 12)
ht <- Heatmap(log2(tags_peak_num_mat+1), rect_gp = gpar(type = "none"), na_col = "grey90", 
        cluster_rows = FALSE, cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
          }
        }, 
        row_labels = TeX(rownames(tags_peak_num_mat)),
        column_labels = TeX(colnames(tags_peak_num_mat)),
        row_names_side = "left",
        column_title = paste0('Target Pair Peak Count'),
        column_names_rot = 45, 
        col = circlize::colorRamp2(c(0, 7.5, 15), c("blue", "white", "red")),
        heatmap_legend_param = list(
          title = "log2(peak counts)",
          color_bar = "continuous",
          itle_gp = gpar(fontsize = 15)
        ),
        column_split = splits, 
        row_split = splits,
        row_title = rep("", length(unique(splits))),
        left_annotation = row_annotation,
        bottom_annotation = col_annotation,
        )
# legends <- packLegend(HeatmapLegend(ht), na_legend)

# Draw heatmap and add combined legend
draw(ht, annotation_legend_list = list(na_legend), merge_legend = TRUE)
dev.off()

