library(ggplot2)
library(ChIPseeker)
library("dplyr")
library(ComplexHeatmap)
library(grid)
library("gridExtra")
library(GenomicRanges)
library(stats)
library(ggrepel)
library(latex2exp)
library(glue)
library(dendextend)
# load functions for plotting
source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/plotUtils.R")
source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
filtered <- TRUE
tagSelection <- "EZH2"
tagSelection <- NA
chisq <- TRUE
args <- commandArgs(trailingOnly = TRUE)
peak_type <- args[1]


dataPeakAnnoDir <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/V/")

# rdsFileT <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/NaB_36x_CT_T/annotate.rds"
# rdsFileT <- paste0(dataPeakAnnoDir, "T/annotate.rds")
# xT <- readRDS(rdsFileT)
# 
# # rdsFileV <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/NaB_36x_CT_V/annotate.rds"
# rdsFileV <- paste0(dataPeakAnnoDir, "V/annotate.rds")
# xV <- readRDS(rdsFileV)

intermediatePlots <- FALSE


CCREPeakAnnoDir <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/V/")
# CCRErdsFileT <- paste0(CCREPeakAnnoDir, "T/annotate.rds")
# CCRExT <- readRDS(CCRErdsFileT)

CCRErdsFileV <- paste0(CCREPeakAnnoDir, "annotate.rds")
CCRExV <- readRDS(CCRErdsFileV)


allFeatures <- allFeatures
# annoT <- lapply(xT, getAnnoStat)
# annoV <- lapply(xV, getAnnoStat)

# CCREannoT <- lapply(CCRExT, getAnnoStatCCRE)
CCREannoV <- lapply(CCRExV, getAnnoStatCCRE)

cat("tag counts: ", length(CCREannoV), "\n")

if (filtered) {
  filteredTags <- filter_target_pairs()
  # filteredTags <- read.table("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/filteredTags.tsv", sep="\t")$x
  
  # filteredTags <- filteredTags[grepl(tagSelection, filteredTags)]
  
  # annoT <- annoT[filteredTags]
  # annoV <- annoV[filteredTags]
  # CCREannoT <- CCREannoT[filteredTags]
  CCREannoV <- CCREannoV[filteredTags]
  filterStr <- "filtered_"
  cat("after filtering: ", length(CCREannoV), "\n")
} else {
  filterStr <- ""
}


clusterBoth <- FALSE

out_dir <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/plots/")
dir.create(out_dir, recursive=TRUE)
if (clusterBoth) {
  # resultMatT <- getAnnotationMatrix(annoT, allFeatures)
  # resultMatV <- getAnnotationMatrix(annoV, allFeatures)
  
  CCREresultMatT <- getAnnotationMatrix(CCREannoT, ChIPSeekerCCRECategoriesOrder)
  CCREresultMatV <- getAnnotationMatrix(CCREannoV, ChIPSeekerCCRECategoriesOrder)
  
  # colnames(resultMatT) <- paste0(colnames(resultMatT), " T")
  # colnames(resultMatV) <- paste0(colnames(resultMatV), " V")
  # 
  # colnames(CCREresultMatT) <- paste0(colnames(CCREresultMatT), " T CCRE")
  # colnames(CCREresultMatV) <- paste0(colnames(CCREresultMatV), " V CCRE")
  
  # mergeDfAnno <- merge(resultMatT, resultMatV, by = "row.names", all = TRUE)
  mergeDfCCRE <- merge(CCREresultMatT, CCREresultMatV, by = "row.names", all = TRUE)
  # mergeDf <- merge(mergeDfAnno, mergeDfCCRE, by = "Row.names", all = TRUE)
  mergeMat <- as.matrix(mergeDfCCRE[,2:(dim(mergeDfCCRE)[2])])
  rownames(mergeMat) <- mergeDfCCRE$Row.names
  
  
  h <- ComplexHeatmap::Heatmap(mergeMat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
                               row_names_gp = grid::gpar(fontsize = 0.5))
  ht = draw(h)
  
  plotFilePrefix <- paste0("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/TV", "/", filterStr, "_addCCRE")
  if (!is.na(tagSelection)) {
    plotFilePrefix <- paste0(plotFilePrefix, tagSelection, "_")
  }
  plotNameClustered <- paste0(plotFilePrefix, "cluster_heatmap.pdf")
  pdf(plotNameClustered, height=120, width=8)
  print(ht)
  dev.off()
  # annoTClusterOrdered <- annoT[rownames(mergeMat)[row_order(ht)]]
  # annoVClusterOrdered <- annoV[rownames(mergeMat)[row_order(ht)]]
  CCREannoTClusterOrdered <- CCREannoT[rownames(mergeMat)[row_order(ht)]]
  CCREannoVClusterOrdered <- CCREannoV[rownames(mergeMat)[row_order(ht)]]
} else { # Cluster on V
  
  resultMat <- getAnnotationMatrix(CCREannoV, ChIPSeekerCCRECategoriesOrder)
  
  # h <- ComplexHeatmap::Heatmap(resultMat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
  #                              row_names_gp = grid::gpar(fontsize = 0.5))
  # ht = draw(h)
  # plotFilePrefix <- paste0("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/V_cluster", "/", filterStr)
  # plotNameClustered <- paste0(plotFilePrefix, "cluster_heatmap.pdf")
  # if (intermediatePlots) {
  #   pdf(plotNameClustered, height=18 * length(filteredTags) / 1444, width=8)
  #   print(ht)
  #   dev.off()
  # }

  # # annoTClusterOrdered <- annoT[rownames(mergeMat)[row_order(ht)]]
  # # annoVClusterOrdered <- annoV[rownames(mergeMat)[row_order(ht)]]
  # CCREannoTClusterOrdered <- CCREannoT[rownames(resultMat)[row_order(ht)]]
  # CCREannoVClusterOrdered <- CCREannoV[rownames(resultMat)[row_order(ht)]]
  
  dist_matrix <- dist(resultMat, method = "euclidean")
  
  # Step 2: Perform Hierarchical Clustering
  hclust_result <- hclust(dist_matrix, method = "complete") # 'method' can be "complete", "single", "average", etc.
  # Convert hclust object to dendrogram
  dend <- as.dendrogram(hclust_result)
  
  # Define clusters (your list)
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
      "H3K27me3-H3K4me3"),
    
    c("MLL1_KMT2A-MLL4_MLL2_KMT2B", 
      "H3K36me3-H3S10ph",
      "CDK8-H2A_XS139ph",
      "PIM1-SETD2",
      "H3K27me3-cFos", 
      "H3K27me3-NRF1", 
      "H3K27me3-MSK1",
      "EHMT2-H3K27me3"),
    
    c("H3K9me3-MSK2", 
      "H3K9me3-Myc",
      "EHMT2-H3K9me3", 
      "H2A_XS139ph-H3K9me3", 
      "H3K9me3-H3K9me3",
      "H3K9me2-H3K9me3", 
      "H3K9me2-MSK2",
      "H3K9me2-YY1")
  )

  # Define colors for the 4 clusters
  cluster_colors <- c("red", "blue", "green", "orange")

  # Initialize all labels as black
  label_colors <- rep("black", length(labels(dend)))
  names(label_colors) <- labels(dend)

  # Assign colors based on cluster membership
  for (i in seq_along(clusters)) {
    label_colors[intersect(labels(dend), clusters[[i]])] <- cluster_colors[i]
  }

  # Apply colors to dendrogram labels
  # dend <- dend %>%
  #   set("labels_col", label_colors[labels(dend)])
  dend <- dendextend::set(dend, "labels_col", label_colors[labels(dend)])
  new_order = rownames(resultMat)[order.dendrogram(as.dendrogram(hclust_result))]
  
  pdf(glue("{out_dir}/dend.pdf"), width=100, height=10)
  # plot(hclust_result, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "", cex = 0.8)
  plot(dend, main = "Hierarchical Clustering Dendrogram", xlab = "", sub = "", cex = 0.8)
  dev.off()
  # CCREannoTClusterOrdered <- CCREannoT[new_order]
  CCREannoVClusterOrdered <- CCREannoV[new_order]
}



extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}



# annoT.df <- list_to_dataframe(annoTClusterOrdered)
# annoT.df$Feature <- factor(annoT.df$Feature, levels = allFeatures)
# categoryColumn <- ".id"
# gpT <- plotAnnoBar.data.frame(annoT.df, categoryColumn=categoryColumn) + ggtitle("NaB_36x_CT_T")
# shared_legend <- extract_legend(gpT)
# gpTNoLegnd <- gpT + theme(legend.position = "none")
# 
# annoV.df <- list_to_dataframe(annoVClusterOrdered)
# annoV.df$Feature <- factor(annoV.df$Feature, levels = allFeatures)
# categoryColumn <- ".id"
# gpV <- plotAnnoBar.data.frame(annoV.df, categoryColumn=categoryColumn) + ggtitle("NaB_36x_CT_V")
# gpVNoLegnd <- gpV + theme(legend.position = "none")
# 
# resultPlot <- grid.arrange(arrangeGrob(gpVNoLegnd, gpTNoLegnd, ncol = 2),
#              shared_legend, ncol = 2, widths = c(8, 2))
# 
# 
# plotNameAnno <- paste0(plotFilePrefix, "anno_clustered.pdf")
# height <-18 * length(filteredTags) / 1444
# textSize = 0.5
# if (height < 10) {
#   height <- 10
#   textSize <- 3
# }
# textSize=16
# pdf(plotNameAnno, height=140, width=40)
# grid.arrange(arrangeGrob(gpVNoLegnd + theme(axis.text=element_text(size=textSize)), gpTNoLegnd + theme(axis.text=element_text(size=textSize)), ncol = 2),
#              shared_legend, ncol = 2, widths = c(1, 0.5))
# dev.off()

# CCREannoT.df <- list_to_dataframe(CCREannoTClusterOrdered)
# categoryColumn <- ".id"
# CCREannoT.df$.id <- map_target_names(CCREannoT.df$.id , target_pair_mapping_df)
# gpT <- plotAnnoBar.data.frame(CCREannoT.df, categoryColumn=categoryColumn, colorOption = 1) + ggtitle("Treated") + scale_x_discrete(labels = TeX)
# shared_legend <- extract_legend(gpT)
# gpTNoLegnd <- gpT + theme(legend.position = "none")

CCREannoV.df <- list_to_dataframe(CCREannoVClusterOrdered)
CCREannoV.df$Feature <- factor(CCREannoV.df$Feature, levels = ChIPSeekerCCRECategoriesOrder)
CCREannoV.df$.id <- map_target_names(CCREannoV.df$.id , target_pair_mapping_df)
CCREannoV.df$.id <- factor(CCREannoV.df$.id , levels=rev(map_target_names(new_order, target_pair_mapping_df)))
categoryColumn <- ".id"
gpV <- plotAnnoBar.data.frame(CCREannoV.df, categoryColumn=categoryColumn, colorOption = 1) + ggtitle("Viecle") + scale_x_discrete(labels = TeX)
write.csv(CCREannoV.df, paste0(out_dir, "CCRE_anno_V_clustered.csv"), row.names = FALSE, quote = FALSE)
pdf(paste0(out_dir, "CCRE_anno_V_clustered.pdf"), height=140, width=20)
print(gpV)
dev.off()
saveRDS(new_order, paste0(out_dir, "CCRE_anno_V_order.rds"))

# gpVNoLegnd <- gpV + theme(legend.position = "none")

# resultPlot <- grid.arrange(arrangeGrob(gpVNoLegnd, gpTNoLegnd, ncol = 2),
#                            shared_legend, ncol = 2, widths = c(8, 2))


# plotNameAnno <- paste0(plotFilePrefix, "CCRE_anno_clustered.pdf")
# height <-18 * length(filteredTags) / 1444
# textSize = 0.5
# if (height < 10) {
#   height <- 10
#   textSize <- 3
# }
# textSize=16
# pdf(plotNameAnno, height=140, width=40)
# grid.arrange(arrangeGrob(gpVNoLegnd + theme(axis.text=element_text(size=textSize)), gpTNoLegnd + theme(axis.text=element_text(size=textSize)), ncol = 2),
#              shared_legend, ncol = 2, widths = c(1, 0.5))
# dev.off()



# plts <- list()
# for (tag in unique(annoV.df$.id)) {
#   annoVTag <- annoV.df[annoV.df$.id == tag, ]
#   annoTTag <- annoT.df[annoT.df$.id == tag, ]
#   p1 <- plotAnnoBar.data.frame(annoVTag, categoryColumn=".id")+ theme(legend.position = "none") + theme(plot.title = element_text(size=1), plot.margin=unit(c(1,1,-0.5,1),  "cm"), panel.border = element_blank(), panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                         plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                         panel.grid.major = element_blank(),
#                                                                                                         panel.grid.minor = element_blank(),
#                                                                                                         legend.background = element_rect(fill='transparent'),
#                                                                                                         legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL)
#   p2 <- plotAnnoBar.data.frame(annoTTag, categoryColumn=".id") + geom_bar(colour = "red", stat="identity")+ theme(legend.position = "none")+ theme(plot.title = element_text(size=1), plot.margin=unit(c(-0.5,1,1,1), "cm"), panel.border = element_blank(),panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                                                                    plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                                                                    panel.grid.major = element_blank(),
#                                                                                                                                                    panel.grid.minor = element_blank(),
#                                                                                                                                                    legend.background = element_rect(fill='transparent'),
#                                                                                                                                                    legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL)
#   plts <- c(plts, list(p1))
#   plts <- c(plts, list(p2))
# }
# do.call("grid.arrange", c(plts[1:10], ncol=1))


# getRevHeterotone <- function(tag) {
#   tag <- as.character(tag)
#   tags <- strsplit(tag, "-")[[1]]
#   tag1 <- tags[1]
#   tag2 <- tags[2]
#   return(paste0(tag2, "-", tag1))
# }


# unqiueHeterotones <- c()
# for (tag in unique(annoV.df$.id)) {
#   if (substring(tag, 1, nchar("EZH2")) == "EZH2") {
#     unqiueHeterotones <- append(unqiueHeterotones, tag)
#   } else {
#     unqiueHeterotones <- append(unqiueHeterotones, getRevHeterotone(tag))
#   }
# }

# 
# unqiueHeterotones <- unique(unqiueHeterotones)
# 
# newAnnoStatV <- list()
# for (tag in unqiueHeterotones) {
#   annoStat <- xV[[tag]]@annoStat
#   peakNum <- xV[[tag]]@peakNum
#   annoStat$Frequency <- annoStat$Frequency * peakNum
#   revAnnoStat <- xV[[getRevHeterotone(tag)]]@annoStat
#   revPeakNum <- xV[[getRevHeterotone(tag)]]@peakNum
#   revAnnoStat$Frequency <- revAnnoStat$Frequency * revPeakNum
#   freq <- list()
#   for (f in allFeatures) {
#     freq[[f]] <- 0
#   }
#   for (f in annoStat$Feature) {
#     freq[[f]] <- freq[[f]] + annoStat$Frequency[annoStat$Feature == f]
#   }
#   for (f in revAnnoStat$Feature) {
#     freq[[f]] <- freq[[f]] + revAnnoStat$Frequency[revAnnoStat$Feature == f]
#   }
#   newannoStat <- data.frame(Feature=names(freq), Frequency=(as.numeric(freq)/sum(as.numeric(freq))*100))
#   newAnnoStatV[[tag]] <- newannoStat
# }
# 
# newAnnoStatV <- getAnnoStatForUniqueHeterotones(xV, allFeatures, unqiueHeterotones) 
# newAnnoStatT <- getAnnoStatForUniqueHeterotones(xT, allFeatures, unqiueHeterotones) 
# 
# 
# plotAnnoBarTV <- function(newAnnoStatV, newAnnoStatT) {
#   annoT.df <- list_to_dataframe(newAnnoStatT)
#   annoT.df$Feature <- factor(annoT.df$Feature, levels = allFeatures)
#   categoryColumn <- ".id"
# 
#   annoV.df <- list_to_dataframe(newAnnoStatV)
#   annoV.df$Feature <- factor(annoV.df$Feature, levels = allFeatures)
# 
#   plts <- list()
#   for (tag in unique(annoV.df$.id)) {
#     annoVTag <- annoV.df[annoV.df$.id == tag, ]
#     annoTTag <- annoT.df[annoT.df$.id == tag, ]
#     p1 <- plotAnnoBar.data.frame(annoVTag, categoryColumn=".id")+ theme(legend.position = "none") + theme(plot.title = element_blank(), plot.margin=unit(c(1,1,-0.5,1),  "cm"), panel.border = element_blank(), panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                           plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                           panel.grid.major = element_blank(),
#                                                                                                           panel.grid.minor = element_blank(),
#                                                                                                           legend.background = element_rect(fill='transparent'),
#                                                                                                           legend.box.background = element_rect(fill='transparent'), 
#                                                                                                           axis.text.x=element_blank(),  #remove y axis labels
#                                                                                                           axis.ticks.x=element_blank())+ labs(x=NULL, y=NULL)
#     p2 <- plotAnnoBar.data.frame(annoTTag, categoryColumn=".id") + geom_bar(colour = "red", stat="summary", fun = sum, fill = "transparent")+ theme(legend.position = "none")+ theme(plot.title = element_blank(), plot.margin=unit(c(-0.5,1,1,1), "cm"), panel.border = element_blank(),panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                                                                     plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                                                                     panel.grid.major = element_blank(),
#                                                                                                                                                     panel.grid.minor = element_blank(),
#                                                                                                                                                     legend.background = element_rect(fill='transparent'),
#                                                                                                                                                     legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL) + scale_alpha_manual(values = c(0.1,0.4))
#     # p2 <- plotAnnoBar.data.frame(annoTTag, categoryColumn=".id") + theme(legend.position = "none")+ theme(plot.title = element_blank(), plot.margin=unit(c(-0.5,1,1,1), "cm"), panel.border = element_blank(),panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#     #                                                                                                                                                            plot.background = element_rect(fill='transparent', color=NA),
#     #                                                                                                                                                            panel.grid.major = element_blank(),
#     #                                                                                                                                                            panel.grid.minor = element_blank(),
#     #                                                                                                                                                            legend.background = element_rect(fill='transparent'),
#     #                                                                                                                                                            legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL) + scale_alpha_manual(values = c(0.1,0.1))
#     plts <- c(plts, list(p1))
#     plts <- c(plts, list(p2))
#   }
#   # same width for plots
#   # https://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
#   grobs <- list()
#   widths <- list()
#   for (i in 1:length(plts)){
#     grobs[[i]] <- ggplotGrob(plts[[i]])
#     widths[[i]] <- grobs[[i]]$widths[2:5]
#   }
#   maxwidth <- do.call(grid::unit.pmax, widths)
#   for (i in 1:length(grobs)){
#     grobs[[i]]$widths[2:5] <- as.list(maxwidth)
#   }
#   return(grobs)
# }
# 
# annoT.df <- list_to_dataframe(newAnnoStatT)
# annoT.df$Feature <- factor(annoT.df$Feature, levels = allFeatures)
# categoryColumn <- ".id"
# 
# annoV.df <- list_to_dataframe(newAnnoStatV)
# annoV.df$Feature <- factor(annoV.df$Feature, levels = allFeatures)
# 
# 
# plts <- list()
# for (tag in unique(annoV.df$.id)) {
#   annoVTag <- annoV.df[annoV.df$.id == tag, ]
#   annoTTag <- annoT.df[annoT.df$.id == tag, ]
#   p1 <- plotAnnoBar.data.frame(annoVTag, categoryColumn=".id")+ theme(legend.position = "none") + theme(plot.title = element_blank(), plot.margin=unit(c(1,1,-0.5,1),  "cm"), panel.border = element_blank(), panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                         plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                         panel.grid.major = element_blank(),
#                                                                                                         panel.grid.minor = element_blank(),
#                                                                                                         legend.background = element_rect(fill='transparent'),
#                                                                                                         legend.box.background = element_rect(fill='transparent'), 
#                                                                                                         axis.text.x=element_blank(),  #remove y axis labels
#                                                                                                         axis.ticks.x=element_blank())+ labs(x=NULL, y=NULL)
#   p2 <- plotAnnoBar.data.frame(annoTTag, categoryColumn=".id") + geom_bar(colour = "red", stat="summary", fun = sum, fill = "transparent")+ theme(legend.position = "none")+ theme(plot.title = element_blank(), plot.margin=unit(c(-0.5,1,1,1), "cm"), panel.border = element_blank(),panel.background = element_rect(fill='transparent'), axis.text=element_text(size=8),
#                                                                                                                                                    plot.background = element_rect(fill='transparent', color=NA),
#                                                                                                                                                    panel.grid.major = element_blank(),
#                                                                                                                                                    panel.grid.minor = element_blank(),
#                                                                                                                                                    legend.background = element_rect(fill='transparent'),
#                                                                                                                                                    legend.box.background = element_rect(fill='transparent'))+ labs(x=NULL, y=NULL) + scale_alpha_manual(values = c(0.1,0.4))
#   plts <- c(plts, list(p1))
#   plts <- c(plts, list(p2))
# }
#  
# # same width for plots
# # https://stackoverflow.com/questions/13294952/left-align-two-graph-edges-ggplot
# grobs <- list()
# widths <- list()
# for (i in 1:length(plts)){
#   grobs[[i]] <- ggplotGrob(plts[[i]])
#   widths[[i]] <- grobs[[i]]$widths[2:5]
# }
# maxwidth <- do.call(grid::unit.pmax, widths)
# for (i in 1:length(grobs)){
#   grobs[[i]]$widths[2:5] <- as.list(maxwidth)
# }
# 
# # do.call("grid.arrange", c(grobs, ncol=1))
# p1 <- do.call("grid.arrange", c(grobs, ncol=1))
# 
# pdf("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/TV_cluster/ezh2.pdf", height=30, width=22)
# resultPlot <- grid.arrange(p1, shared_legend, ncol = 2, widths = c(7, 2))
# dev.off()

# test=as.data.frame(cbind(a=c(1,1,2,3), b=1:4, c=as.character(1:4)))
# ggplot(test) + geom_bar(aes(x=a, y=b, fill=c), colour="blue", stat="identity")
