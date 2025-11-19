library(ggplot2)
library(cowplot)
library(GenomicRanges)
library(glue)
library(RJSONIO)
library(tidyverse)
library(rtracklayer)


options(error=function() { traceback(2); if(!interactive()) quit("no", status = 1, runLast = FALSE) })

source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")


options <- commandArgs(trailingOnly = TRUE)

if (length(options) != 0) {
  biclustering_result_file <- options[1]
  save_dir <- options[2]
}
print(biclustering_result_file)
print(save_dir)

biclustering_result <- read.table(biclustering_result_file, header = TRUE, sep = "\t")
biclustering_result <- biclustering_result[!biclustering_result$label %in% c("Epitope_Specific", "Background"), ]
pos_df <- do.call(rbind, (strsplit(biclustering_result$feature, "_")))
colnames(pos_df) <- c("seqnames", "start", "end")
biclustering_result <- cbind(pos_df, biclustering_result)
biclustering_gr <-
  makeGRangesFromDataFrame(
    biclustering_result,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = TRUE
  )

wg_overlap_res <- findOverlaps(biclustering_gr, wg_biclustering_gr)
wg_overlap_df <- data.frame(query=biclustering_gr$label[queryHits(wg_overlap_res)],
                             subject=wg_biclustering_gr$label[subjectHits(wg_overlap_res)])
wg_overlap_mat <- table(wg_overlap_df$query, wg_overlap_df$subject)
wg_overlap_mat <- wg_overlap_mat / rowSums(wg_overlap_mat)
wg_overlap_mat <- as.data.frame(wg_overlap_mat)
cluster_all_levels <-
  c(
    as.character(1:100),
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "Epitope_Specific",
    "non-target_specific",
    "background",
    "Background"
  )

unique_clusters <- sort(unique(biclustering_result$label), method = "radix")
chipseeker_ccre_annotations <- list()
chipseeker_ccre_celltype_agnostic_annotations <- list()
repeatmasker_annotations <- list()
chromhmm_short_annotations <- list()
chromhmm_full_annotations <- list()
chromhmm_group_annotations <- list()
for (cluster_id in unique_clusters) {
  biclustering_result_i <- biclustering_result[biclustering_result$label == cluster_id, ]
  biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  #####
  chipseeker_ccre_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(biclustering_grange_i, annotation, categories)
  chipseeker_ccre_celltype_agnostic_annotation <- annotatePeakByOverlappingChIPSeekerCCRE(biclustering_grange_i, annotation_celltype_agnostic, categories, featureColname="V6")
  repeatmasker_annotation <- annotatepeakByOverlappingRepeatMasker(biclustering_grange_i, annotationRepeatMasker, repeatMaskerFeatures)
  #####
  chromhmm_short_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, categoriesChromHMM, featureColname = "V4")
  chromhmm_full_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, unique(annotationChromHMM$full_anno), featureColname = "full_anno")
  chromhmm_group_annotation <- annotatepeakByOverlappingChromHMM(biclustering_grange_i, annotationChromHMM, unique(annotationChromHMM$group), featureColname = "group")
  chipseeker_ccre_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_annotation
  repeatmasker_annotations[[glue("{cluster_id}")]] <- repeatmasker_annotation
  chromhmm_short_annotations[[glue("{cluster_id}")]] <- chromhmm_short_annotation
  chromhmm_full_annotations[[glue("{cluster_id}")]] <- chromhmm_full_annotation
  chromhmm_group_annotations[[glue("{cluster_id}")]] <- chromhmm_group_annotation
  chipseeker_ccre_celltype_agnostic_annotations[[glue("{cluster_id}")]] <- chipseeker_ccre_celltype_agnostic_annotation
}


################## 
anno_ccre <- lapply(chipseeker_ccre_annotations, getAnnoStatCCRE)
anno_ccre.df <- list_to_dataframe(anno_ccre)
anno_ccre.df$Feature <- factor(anno_ccre.df$Feature, levels = c(ChIPSeekerCCRECategoriesOrder))
categoryColumn <- ".id"
cluster_levels <- cluster_all_levels[cluster_all_levels %in% unique(biclustering_result$label)]
anno_ccre.df$.id <- factor(anno_ccre.df$.id, levels = rev(cluster_levels))
print(anno_ccre.df$.id)
p_ccre <-
  plotAnnoBar.data.frame.one.target(anno_ccre.df,
                         categoryColumn = categoryColumn,
                         colorOption = 1,
                         features=ChIPSeekerCCRECategoriesOrder) + theme(
                           legend.position = "bottom",
                           plot.title = element_text(size = 20, ),
                           axis.text.x = element_text(size = 25, ),
                           axis.text.y = element_text(size = 25, ),
                           legend.text = element_text(size = 16),
                           legend.key.size = unit(1, 'cm'),
                           axis.title.x = element_text(size = 25)
                         ) + ggtitle("Cis-regulatory Elements K562-specific") + guides(fill = guide_legend(title = NULL, nrow =
                                                                                               5, reverse = TRUE))
################## plotting chipseeker + ccre
anno_ccre <- lapply(chipseeker_ccre_celltype_agnostic_annotations, getAnnoStatCCRE)
anno_ccre.df <- list_to_dataframe(anno_ccre)
anno_ccre.df$Feature <- factor(anno_ccre.df$Feature, levels = c(ChIPSeekerCCRECategoriesOrder))
categoryColumn <- ".id"

cluster_levels <- cluster_all_levels[cluster_all_levels %in% unique(biclustering_result$label)]
anno_ccre.df$.id <- factor(anno_ccre.df$.id, levels = rev(cluster_levels))
print(anno_ccre.df$.id)
p_ccre_agnostic <-
  plotAnnoBar.data.frame.one.target(anno_ccre.df,
                                    categoryColumn = categoryColumn,
                                    colorOption = 1,
                                    features=ChIPSeekerCCRECategoriesOrder) + theme(
                                      legend.position = "bottom",
                                      plot.title = element_text(size = 20, ),
                                      axis.text.x = element_text(size = 25, ),
                                      axis.text.y = element_text(size = 25, ),
                                      legend.text = element_text(size = 16),
                                      legend.key.size = unit(1, 'cm'),
                                      axis.title.x = element_text(size = 25)
                                    ) + ggtitle("Cis-regulatory Elements (cell type agnostic)") + guides(fill = guide_legend(title = NULL, nrow =
                                                                                                          5, reverse = TRUE))
##################

anno_repeatmasker <- lapply(repeatmasker_annotations, getAnnoStatCCRE)
anno_repeatmasker.df <- list_to_dataframe(anno_repeatmasker)
anno_repeatmasker.df$Feature <- factor(anno_repeatmasker.df$Feature, levels = c(repeatMaskerFeatures, "other"))
anno_repeatmasker.df$.id <- factor(anno_repeatmasker.df$.id, levels = rev(cluster_levels))
categoryColumn <- ".id"


p_repeat <-
  plotAnnoBar.data.frame.one.target(anno_repeatmasker.df,
                         categoryColumn = categoryColumn,
                         colorOption = 4, features=repeatMaskerCategoriesOrder) + theme(
                           legend.position = "bottom",
                           plot.title = element_text(size = 20, ),
                           axis.text.x = element_text(size = 25, ),
                           axis.text.y = element_blank(),
                           legend.text = element_text(size = 16),
                           legend.key.size = unit(1, 'cm'),
                           axis.title.x = element_text(size = 25)
                         ) + ggtitle("Repetitive Elements") + guides(fill = guide_legend(title = NULL, nrow = 5, reverse = TRUE))

################## plotting chromhmm
anno_chromhmm_short <- lapply(chromhmm_short_annotations, getAnnoStatCCRE)
anno_chromhmm_short.df <- list_to_dataframe(anno_chromhmm_short)
anno_chromhmm_short.df$Feature <- factor(anno_chromhmm_short.df$Feature, levels = c(unique(annotationChromHMM$V4), "other"))
categoryColumn <- ".id"
anno_chromhmm_short.df$.id <- factor(anno_chromhmm_short.df$.id, levels = rev(cluster_levels))
print(anno_chromhmm_short.df$.id)
p_chromhmm_short <-
  plotAnnoBar.data.frame.one.target(anno_chromhmm_short.df,
                                    categoryColumn = categoryColumn,
                                    colorOption = 6,
                                    features=c(unique(annotationChromHMM$V4), "other")) + theme(
                                      legend.position = "bottom",
                                      plot.title = element_text(size = 20, ),
                                      axis.text.x = element_text(size = 25, ),
                                      axis.text.y = element_text(size = 25, ),
                                      legend.text = element_text(size = 16),
                                      legend.key.size = unit(1, 'cm'),
                                      axis.title.x = element_text(size = 25)
                                    ) + ggtitle("ChromHMM Elements") + guides(fill = guide_legend(title = NULL, nrow =
                                                                                                          5, reverse = TRUE))
###################
anno_chromhmm_full <- lapply(chromhmm_full_annotations, getAnnoStatCCRE)
anno_chromhmm_full.df <- list_to_dataframe(anno_chromhmm_full)
anno_chromhmm_full.df$Feature <- factor(anno_chromhmm_full.df$Feature, levels = c(sort(unique(annotationChromHMM$full_anno)), "other"))
categoryColumn <- ".id"
anno_chromhmm_full.df$.id <- factor(anno_chromhmm_full.df$.id, levels = rev(cluster_levels))
print(anno_chromhmm_full.df$.id)
p_chromhmm_full <-
  plotAnnoBar.data.frame.one.target(anno_chromhmm_full.df,
                                    categoryColumn = categoryColumn,
                                    colorOption = 5,
                                    features=c(unique(annotationChromHMM$full_anno), "other")) + theme(
                                      legend.position = "bottom",
                                      plot.title = element_text(size = 20, ),
                                      axis.text.x = element_text(size = 25, ),
                                      axis.text.y = element_text(size = 25, ),
                                      legend.text = element_text(size = 16),
                                      legend.key.size = unit(1, 'cm'),
                                      axis.title.x = element_text(size = 25)
                                    ) + ggtitle("ChromHMM Elements (full ver.)") + guides(fill = guide_legend(title = NULL, nrow =
                                                                                                                   7, reverse = TRUE))

anno_chromhmm_group <- lapply(chromhmm_group_annotations, getAnnoStatCCRE)
anno_chromhmm_group.df <- list_to_dataframe(anno_chromhmm_group)
anno_chromhmm_group.df$Feature <- factor(anno_chromhmm_group.df$Feature, levels = c(unique(annotationChromHMM$group), "other"))
categoryColumn <- ".id"
anno_chromhmm_group.df$.id <- factor(anno_chromhmm_group.df$.id, levels = rev(cluster_levels))
print(anno_chromhmm_group.df$.id)
p_chromhmm_group <-
  plotAnnoBar.data.frame.one.target(anno_chromhmm_group.df,
                                    categoryColumn = categoryColumn,
                                    colorOption = 4,
                                    features=c(unique(annotationChromHMM$group), "other")) + theme(
                                      legend.position = "bottom",
                                      plot.title = element_text(size = 20, ),
                                      axis.text.x = element_text(size = 25, ),
                                      axis.text.y = element_text(size = 25, ),
                                      legend.text = element_text(size = 16),
                                      legend.key.size = unit(1, 'cm'),
                                      axis.title.x = element_text(size = 25)
                                    ) + ggtitle("ChromHMM Elements (group ver.)") + guides(fill = guide_legend(title = NULL, nrow =
                                                                                                                  5, reverse = TRUE))

promoter_range <- "0-0"
# promoter_range <- "5000-5000"

promoter_file <- glue("/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/ccre_region/promoter_-{promoter_range}.tsv")
promoter_table <- read.table(promoter_file, sep = "\t", header = TRUE)
promoter_table$seqnames <- paste0("chr", promoter_table$seqnames)
promoter_grange <- makeGRangesFromDataFrame(promoter_table, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand",keep.extra.columns = TRUE)





rnaseq_dir_filename ="/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/RNA_seq_TPM_all.csv"
# rnaseq_dir_filename = file.path(rnaseq_dir, rnaseq_filename)
rnaseq_raw = read.table(rnaseq_dir_filename, sep=",", check.names=FALSE, header=TRUE, row.names=1)
rnaseq_raw[, "V"] = (rnaseq_raw[, "veh.1.RNA.seq_S1"] + rnaseq_raw[, "veh.2.RNA.seq_S2"]) / 2
rnaseq_raw[, "T"] = (rnaseq_raw[, "NaB.1.RNA.seq_S3"] + rnaseq_raw[, "NaB.2.RNA.seq_S4"]) / 2
rnaseq <- rnaseq_raw[, c("gene_name", "V", "T")]
rownames(rnaseq) <- gsub("\\..*","",rownames(rnaseq))
rnaseq$V <- log10(rnaseq$V + 1)
rnaseq$T <- log10(rnaseq$T + 1)



diff_rnaseq_dir_filename ="/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/diff_expr_all_update.csv"
rnaseq_diff = read.table(diff_rnaseq_dir_filename, sep=",", check.names=FALSE, header=TRUE, row.names=1)
rnaseq_diff = rnaseq_diff[!is.na(rnaseq_diff$log2FoldChange),]
rnaseq_diff <- rnaseq_diff[which(abs(rnaseq_diff$log2FoldChange)>=1),]
rnaseq_diff <- rnaseq_diff[which(abs(rnaseq_diff$padj)<= 0.5) ,]
cluster_ids <- c()
rnaseq_vals <- c()

rnaseq_vals_T <- c()
gene_ids <- c()


overlapped_promoter_list <- list()
overlapped_gene_name_list <- list()
nearest_dist_cutoff <- 5000
for (cluster_id in unique_clusters) {
  biclustering_result_i <- biclustering_result[biclustering_result$label == cluster_id, ]
  biclustering_result_i[,"start"] <- round((as.numeric(biclustering_result_i[, "end"]) - as.numeric(biclustering_result_i[, "start"])) %/% 2) + as.numeric(biclustering_result_i[, "start"])
  biclustering_result_i[,"end"] <- biclustering_result_i[,"start"] 
  biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
  # overlap_res <- findOverlaps(biclustering_grange_i, promoter_grange)
  overlap_res <-
    distanceToNearest(biclustering_grange_i, promoter_grange)
  overlap_dist <- mcols(overlap_res)$distance
  overlap_res <- overlap_res[overlap_dist <= nearest_dist_cutoff]
  print(glue("num of bin: {length(biclustering_grange_i)}"))
  print(glue("num of overlapped genes: {length(unique(subjectHits(overlap_res)))}"))
  overlaped_promoters <- promoter_grange[unique(subjectHits(overlap_res))]
  overlapped_promoter_list[[glue("{cluster_id}")]] <- overlaped_promoters$gene_id
  overlapped_gene_name_list[[glue("{cluster_id}")]] <- overlaped_promoters$gene_name
  rnaseq_vals <- append(rnaseq_vals, rnaseq[overlaped_promoters$gene_id, "V"])

  rnaseq_vals_T <- append(rnaseq_vals_T, rnaseq[overlaped_promoters$gene_id, "T"])
  cluster_ids <- append(cluster_ids, rep(glue("{cluster_id}"), length(overlaped_promoters)))
  gene_ids <- append(gene_ids, overlaped_promoters$gene_id)
}
is_grouped_df_plot <- FALSE
if ("fold_change" %in% colnames(biclustering_result)) {
  if (length(unique(biclustering_result$fold_change))==2) {
    is_grouped_df_plot <- TRUE
  }
}
print(is_grouped_df_plot)
if (is_grouped_df_plot == TRUE) {
  log2fc_vals <- c()
  log2fc_cluster_ids <- c()
  fold_changes <- c()
  for (fold_change in unique(biclustering_result$fold_change)) {
    for (cluster_id in unique_clusters) {
      biclustering_result_i <-
        biclustering_result[biclustering_result$label == cluster_id,]
      biclustering_result_i <- biclustering_result_i[biclustering_result_i$fold_change==fold_change,]
      biclustering_result_i[,"start"] <- round((as.numeric(biclustering_result_i[, "end"]) - as.numeric(biclustering_result_i[, "start"])) %/% 2) + as.numeric(biclustering_result_i[, "start"])
      biclustering_result_i[,"end"] <- biclustering_result_i[,"start"] 
      biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
      # overlap_res <-
      #   findOverlaps(biclustering_grange_i, promoter_grange)
      overlap_res <-
        distanceToNearest(biclustering_grange_i, promoter_grange)
      overlap_dist <- mcols(overlap_res)$distance
      overlap_res <- overlap_res[overlap_dist <= nearest_dist_cutoff]
      print(glue("num of bin: {length(biclustering_grange_i)}"))
      print(glue(
        "num of overlapped genes: {length(unique(subjectHits(overlap_res)))}"
      ))
      overlaped_promoters <-
        promoter_grange[unique(subjectHits(overlap_res))]
      
      overlapped_genes <-
        overlaped_promoters$gene_id[overlaped_promoters$gene_id %in% rownames(rnaseq_diff)]
      
      log2fc_vals <-
        append(log2fc_vals, rnaseq_diff[overlapped_genes, "log2FoldChange"])
      log2fc_cluster_ids <-
        append(log2fc_cluster_ids, rep(glue("{cluster_id}"), length(rnaseq_diff[overlapped_genes, "log2FoldChange"])))
      fold_changes <- append(fold_changes, rep(fold_change, length(rnaseq_diff[overlapped_genes, "log2FoldChange"])))
    }
  }
  df_plt_data <- data.frame(cluster_ids=log2fc_cluster_ids, rnaseq_vals=log2fc_vals, fold_changes=fold_changes)
  df_plt_data$fold_changes <- factor(df_plt_data$fold_changes, levels = c("up", "down"))
  p_RNA_df <- ggplot(df_plt_data, aes(x = cluster_ids, y = rnaseq_vals, fill = fold_changes)) +
    geom_boxplot( alpha = 0.5) + ggtitle("mRNA log2FC") + 
    scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) + 
    theme(
      plot.title = element_text(size = 20, ),
      axis.title.x = element_text(size = 25),
      axis.ticks.y =
        element_blank(),
      axis.text.x =
        element_text(size = 25),
      axis.text.y =
        element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = "black"),
      legend.position =
        "bottom",
      axis.title.y =
        element_blank()
    )  + coord_flip() + ylab("Log2(FC)")
  
  p_RNA_df_ <- ggplot(df_plt_data, aes(x = cluster_ids, y = rnaseq_vals, fill = fold_changes)) +
    geom_boxplot( alpha = 0.5) + ggtitle("mRNA log2FC") + 
    scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) + 
    theme(
      plot.title = element_text(size = 20, ),
      axis.title.x = element_text(size = 25),
      axis.ticks.y =
        element_text(size = 25),
      axis.text.x =
        element_text(size = 25),
      axis.text.y =
        element_text(size = 25),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = "black"),
      legend.position =
        "bottom",
      axis.title.y =
        element_blank()
    )  +  xlab("Log2(FC)")
  
} else {
  log2fc_vals <- c()
  log2fc_cluster_ids <- c()
  for (cluster_id in unique_clusters) {
    biclustering_result_i <-
      biclustering_result[biclustering_result$label == cluster_id,]
    biclustering_result_i[,"start"] <- round((as.numeric(biclustering_result_i[, "end"]) - as.numeric(biclustering_result_i[, "start"])) %/% 2) + as.numeric(biclustering_result_i[, "start"])
    biclustering_result_i[,"end"] <- biclustering_result_i[,"start"] 
    biclustering_grange_i <- makeGRangesFromDataFrame(biclustering_result_i, seqnames.field = "seqnames", start.field = "start", end.field = "end", keep.extra.columns = TRUE)
    # overlap_res <-
    #   findOverlaps(biclustering_grange_i, promoter_grange)
    overlap_res <-
      distanceToNearest(biclustering_grange_i, promoter_grange)
    overlap_dist <- mcols(overlap_res)$distance
    overlap_res <- overlap_res[overlap_dist <= nearest_dist_cutoff]
    print(glue("num of bin: {length(biclustering_grange_i)}"))
    print(glue(
      "num of overlapped genes: {length(unique(subjectHits(overlap_res)))}"
    ))
    overlaped_promoters <-
      promoter_grange[unique(subjectHits(overlap_res))]
    
    overlapped_genes <-
      overlaped_promoters$gene_id[overlaped_promoters$gene_id %in% rownames(rnaseq_diff)]
    
    log2fc_vals <-
      append(log2fc_vals, rnaseq_diff[overlapped_genes, "log2FoldChange"])
    log2fc_cluster_ids <-
      append(log2fc_cluster_ids, rep(glue("{cluster_id}"), length(rnaseq_diff[overlapped_genes, "log2FoldChange"])))
  }
  
  
  df_plt_data <- data.frame(cluster_ids=log2fc_cluster_ids, rnaseq_vals=log2fc_vals)
  p_RNA_df <- ggplot(df_plt_data, aes(x = cluster_ids, y = rnaseq_vals)) +
    geom_boxplot(fill = "#6CB0D6", alpha = 0.5) + ggtitle("mRNA log2FC") + 
    scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) + 
    theme(
      plot.title = element_text(size = 20, ),
      axis.title.x = element_text(size = 25),
      axis.ticks.y =
        element_blank(),
      axis.text.x =
        element_text(size = 25),
      axis.text.y =
        element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = "black"),
      legend.position =
        "none",
      axis.title.y =
        element_blank()
    )  + coord_flip() + ylab("Log2(FC)")
  
  pdf(gsub(".pdf", "_dif_rnaseq.pdf", save_dir, fixed=TRUE), width = 38/2.4+5, height = 5)
  print(p_RNA_df)
  dev.off()
}


gene_select_names <- c("coding_cpg", "coding_non_cpg")

gene_select_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/dna_methylation/cpg_island"
gene_select_filename <- "gene_categories.json"
gene_select_dir_filename <- file.path(gene_select_dir, gene_select_filename)

gene_select_dict <- fromJSON(gene_select_dir_filename)


set1s <- c()
set2s <- c()
cluster_id1s <- c()
cluster_id2s <- c()
overlap_cts <- c()

for (cluster_id1 in names(overlapped_gene_name_list)) {
  for (cluster_id2 in names(overlapped_gene_name_list)) {
    if (cluster_id1 != cluster_id2) {
      set1 <- unique(overlapped_gene_name_list[[cluster_id1]])
      set2 <- unique(overlapped_gene_name_list[[cluster_id2]])
      total_gene <- unique(c(set1, set2))
      result_df <- as.data.frame(matrix(0, nrow=2, ncol=length(total_gene)))
      rownames(result_df) <- c(cluster_id1, cluster_id2)
      colnames(result_df) <- total_gene
      result_df[cluster_id1, set1] <- 1
      result_df[cluster_id2, set2] <- 1
      print(dim(result_df))
      ht <- Heatmap(result_df, name = "count", column_title = glue("{cluster_id1}: {length(set1)} vs {cluster_id2}: {length(set2)}\n overlap: {length(intersect(set1, set2))}"))
      set1s <- append(set1s, length(set1))
      set2s <- append(set2s, length(set2))
      cluster_id1s <- append(cluster_id1s, cluster_id1)
      cluster_id2s <- append(cluster_id2s, cluster_id2)
      overlap_cts <- append(overlap_cts, length(intersect(set1, set2)))
    }
  }
}
result_df <- data.frame(cluster_id1=cluster_id1s, cluster_id2=cluster_id2s, set1=set1s, set2=set2s, overlap=overlap_cts)
write.csv(result_df, gsub(".pdf", "_overlap_counts.csv", save_dir, fixed=TRUE), row.names = FALSE, quote = FALSE)

overlap_cnts <- result_df[, c("cluster_id1", "cluster_id2", "overlap")]
library(reshape2)
df_wide <- dcast(overlap_cnts, cluster_id1 ~ cluster_id2, value.var = "overlap")
rownames(df_wide) <- df_wide$cluster_id1
df_wide$cluster_id1 <- NULL
colname_df_wide <- colnames(df_wide)
colname_df_wide <- factor(colname_df_wide, levels=cluster_levels)
colname_df_wide <- as.vector(sort(colname_df_wide))
df_wide <- df_wide[colname_df_wide, colname_df_wide]


pdf(gsub(".pdf", "_overlap_counts.pdf", save_dir, fixed=TRUE))
print(Heatmap(df_wide, cluster_rows = FALSE, cluster_columns = FALSE, name = "shared_mRNA"))
dev.off()

plt_data <- data.frame(cluster_ids=cluster_ids, rnaseq_vals=rnaseq_vals, rnaseq_vals_T=rnaseq_vals_T ,gene_ids=gene_ids)
plt_data$gene_type <- NULL

for (gene_type in names(gene_select_dict)) {
  plt_data[plt_data$gene_ids %in% gene_select_dict[[gene_type]], "gene_type"] <- gene_type
}

gene_type_mat <- table(plt_data$cluster_ids, plt_data$gene_type)
gene_type_mat <- gene_type_mat / rowSums(gene_type_mat)
gene_type_mat <- as.data.frame(gene_type_mat)

gene_type_mat$Var1 <- factor(gene_type_mat$Var1, levels = rev(cluster_levels))

p_gene_percentage <- ggplot(gene_type_mat, aes(x = factor(Var1), y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Proportion", title = "Gene Percentage") +
  theme_minimal() + theme(
    plot.title = element_text(size = 20, ),
    axis.title.x = element_text(size = 25),
    axis.ticks.y =
      element_blank(),
    axis.text.x =
      element_text(size = 25),
    axis.text.y =
      element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_rect(fill = "white",
                                    colour = "black"),
    axis.title.y =
      element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, 'cm')
  ) + coord_flip() + guides(fill = guide_legend(title = NULL, nrow =
                                                  5, reverse = TRUE))



plt_data_list <- list()
for (cluster_id in unique(plt_data$cluster_ids)) {
  plt_data_i <- plt_data[plt_data$cluster_ids == glue("{cluster_id}"),]
  plt_data_i$shared <- NA
  plt_data_not_i <- plt_data[plt_data$cluster_ids != glue("{cluster_id}"),]
  plt_data_i[plt_data_i$gene_ids %in% plt_data_not_i$gene_ids, "shared"] <- "shared"
  plt_data_i[!plt_data_i$gene_ids %in% plt_data_not_i$gene_ids, "shared"] <- "unique"
  plt_data_list[[glue("{cluster_id}")]] <- plt_data_i
}

plt_data_extra <- do.call(rbind, plt_data_list)
plt_data_extra$cluster_ids <- factor(plt_data_extra$cluster_ids, levels = rev(cluster_levels))
plt_data$cluster_ids <- factor(plt_data$cluster_ids, levels = rev(cluster_levels))
showTV <- FALSE
if (showTV == TRUE) {
  plt_data_V <- plt_data[, c("cluster_ids", "rnaseq_vals")]
  plt_data_V$cond <- "V"
  plt_data_T <- plt_data[, c("cluster_ids", "rnaseq_vals_T")]
  colnames(plt_data_T)[2] <- "rnaseq_vals"
  plt_data_T$cond <- "T"
  plt_data <- rbind(plt_data_T, plt_data_V)
  plt_data$cond <- factor(plt_data$cond, levels = rev(c("V", "T")))
  p_RNA <- ggplot(plt_data, aes(x = cluster_ids, y = rnaseq_vals, fill = cond)) +
    geom_boxplot(alpha = 0.5) + ggtitle("mRNA level") + 
    scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) + 
    theme(
      plot.title = element_text(size = 20, ),
      axis.title.x = element_text(size = 25),
      axis.ticks.y =
        element_blank(),
      axis.text.x =
        element_text(size = 25),
      axis.text.y =
        element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = "black"),
      legend.position =
        "bottom",
      axis.title.y =
        element_blank()
    ) + coord_flip() + ylab("Log10(TPM)")
} else {
  p_RNA <- ggplot(plt_data, aes(x = cluster_ids, y = rnaseq_vals)) +
    geom_boxplot(fill = "#6CB0D6", alpha = 0.5) + ggtitle("mRNA level") + 
    scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) + 
    theme(
      plot.title = element_text(size = 20, ),
      axis.title.x = element_text(size = 25),
      axis.ticks.y =
        element_blank(),
      axis.text.x =
        element_text(size = 25),
      axis.text.y =
        element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(fill = "white",
                                      colour = "black"),
      legend.position =
        "none",
      axis.title.y =
        element_blank()
    ) + coord_flip() + ylab("Log10(TPM)")
}


p_RNA_shared <- ggplot(plt_data_extra, aes(x = cluster_ids, y = rnaseq_vals, fill = shared)) +
  geom_boxplot(position = position_dodge(0.8),  alpha = 0.5) + ggtitle("mRNA levels") + 
  scale_x_discrete("cluster_ids", breaks=factor(cluster_levels, levels=rev(cluster_levels)), drop=FALSE) +
  theme(
    plot.title = element_text(size = 20, ),
    axis.title.x = element_text(size = 25),
    axis.ticks.y =
      element_blank(),
    axis.text.x =
      element_text(size = 25),
    axis.text.y =
      element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    panel.background = element_rect(fill = "white",
                                    colour = "black"),
    legend.position =
      "bottom",
    axis.title.y =
      element_blank()
  ) + coord_flip() + ylab("Log10(TPM)")

get_legend2 <- function(plot, legend = NULL) {
  if (is.ggplot(plot)) {
    gt <- ggplotGrob(plot)
  } else {
    if (is.grob(plot)) {
      gt <- plot
    } else {
      stop("Plot object is neither a ggplot nor a grob.")
    }
  }
  pattern <- "guide-box"
  if (!is.null(legend)) {
    pattern <- paste0(pattern, "-", legend)
  }
  indices <- grep(pattern, gt$layout$name)
  not_empty <- !vapply(
    gt$grobs[indices], 
    inherits, what = "zeroGrob", 
    FUN.VALUE = logical(1)
  )
  indices <- indices[not_empty]
  if (length(indices) > 0) {
    return(gt$grobs[[indices[1]]])
  }
  return(NULL)
}

legend1 <- get_legend2(p_ccre, "bottom")
legend2 <- get_legend2(p_repeat, "bottom")
legend3 <- get_legend2(p_chromhmm_short, "bottom")
legend4 <- get_legend2(p_chromhmm_group, "bottom")
legend_gene_percentage <- get_legend2(p_gene_percentage, "bottom")
legend_wg_percentage <- get_legend2(p_wg_percentage, "bottom")
combined_plot <-
  plot_grid(
    p_ccre + theme(legend.position = "none"),
    p_chromhmm_short + theme(legend.position = "none", axis.text.y = element_blank()),
    p_RNA,
    align = "h",
    axis = "bt",
    ncol = 3,
    rel_widths = c(2.2,2, 1)
  )
final_plot <- plot_grid(
  combined_plot,
  plot_grid(legend1, legend3, nrow = 1), # Combine legends horizontally
  ncol = 1,
  rel_heights = c(1, 0.2) # Adjust height ratios
)

pdf(save_dir, width = 30/2.4+5, height = 70/2.4)
print(final_plot)
dev.off()


combined_plot <-
  plot_grid(
    p_ccre_agnostic + theme(legend.position = "none"),
    p_chromhmm_short + theme(legend.position = "none", axis.text.y = element_blank()),
    p_RNA,
    align = "h",
    axis = "bt",
    ncol = 3,
    rel_widths = c(2.2,2, 1)
  )
final_plot <- plot_grid(
  combined_plot,
  plot_grid(legend1, legend3, nrow = 1), # Combine legends horizontally
  ncol = 1,
  rel_heights = c(1, 0.2) # Adjust height ratios
)

pdf(gsub(".pdf", "_agnostic.pdf", save_dir, fixed=TRUE), width = 30/2.4+5, height = 70/2.4)
print(final_plot)
dev.off()


legend5 <- get_legend2(p_chromhmm_full, "bottom")
pdf(gsub(".pdf", "_chromhmm_full.pdf", save_dir, fixed=TRUE), width = 38/2.4+5, height = 40/2.4)
print(p_chromhmm_full)
dev.off()
