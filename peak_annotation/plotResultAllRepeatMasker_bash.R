library(ggplot2)
library(ChIPseeker)
library("dplyr")
library(ComplexHeatmap)
library(grid)
library("gridExtra")
library(GenomicRanges)
library(stats)
library(ggrepel)
library(glue)
library(dendextend)

# load functions for plotting
source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/plotUtils.R")
source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/CCREUtils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
filtered <- TRUE
tagSelection <- "EZH2"
tagSelection <- NA
chisq <- TRUE

args <- commandArgs(trailingOnly = TRUE)
peak_type <- args[1]

print(peak_type)

dataPeakAnnoDir <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/V/")

# rdsFileT <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/NaB_36x_CT_T/annotate.rds"
# rdsFileT <- paste0(dataPeakAnnoDir, "T/annotate.rds")
# xT <- readRDS(rdsFileT)

# rdsFileV <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/NaB_36x_CT_V/annotate.rds"
rdsFileV <- paste0(dataPeakAnnoDir, "annotate.rds")
xV <- readRDS(rdsFileV)

intermediatePlots <- FALSE

out_dir <- glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/{peak_type}/plots/")
dir.create(out_dir, recursive=TRUE)

allFeatures <- allFeatures
# annoT <- lapply(xT, getAnnoStatCCRE)
annoV <- lapply(xV, getAnnoStatCCRE)


cat("tag counts: ", length(annoV), "\n")

if (filtered) {
  filteredTags <- filter_target_pairs()
  # filteredTags <- read.table("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/filteredTags.tsv", sep="\t")$x
  # if (!is.na(tagSelection)) {
  #   filteredTags <- filteredTags[grepl(tagSelection, filteredTags)]
  # }
  # annoT <- annoT[filteredTags]
  annoV <- annoV[filteredTags]
  filterStr <- "filtered_"
  cat("after filtering: ", length(annoV), "\n")
} else {
  filterStr <- ""
}


clusterBoth <- FALSE
if (clusterBoth) {
  resultMatT <- getAnnotationMatrix(annoT, repeatMaskerFeatures)
  resultMatV <- getAnnotationMatrix(annoV, repeatMaskerFeatures)
  
  # CCREresultMatT <- getAnnotationMatrix(CCREannoT, CCREFeatures)
  # CCREresultMatV <- getAnnotationMatrix(CCREannoV, CCREFeatures)
  
  # mergeDfAnno <- merge(resultMatT, resultMatV, by = "row.names", all = TRUE)
  mergeDfCCRE <- merge(resultMatT, resultMatV, by = "row.names", all = TRUE)
  # mergeDf <- merge(mergeDfAnno, mergeDfCCRE, by = "Row.names", all = TRUE)
  mergeMat <- as.matrix(mergeDfCCRE[,2:(dim(mergeDfCCRE)[2])])
  rownames(mergeMat) <- mergeDfCCRE$Row.names
  
  
  h <- ComplexHeatmap::Heatmap(mergeMat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
                               row_names_gp = grid::gpar(fontsize = 0.5))
  ht = draw(h)
  
  plotFilePrefix <- paste0("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/combineReplicates/homotone_heterotone/data_peak_all_auc_0.05/TV_cluster", "/", filterStr, "_addCCRE")
  if (!is.na(tagSelection)) {
    plotFilePrefix <- paste0(plotFilePrefix, tagSelection, "_")
  }
  plotNameClustered <- paste0(plotFilePrefix, "cluster_heatmap.pdf")
  pdf(plotNameClustered, height=18 * length(filteredTags) / 1444, width=16)
  print(ht)
  dev.off()
  annoTClusterOrdered <- annoT[rownames(mergeMat)[row_order(ht)]]
  annoVClusterOrdered <- annoV[rownames(mergeMat)[row_order(ht)]]
  # CCREannoTClusterOrdered <- annoT[rownames(mergeMat)[row_order(ht)]]
  # CCREannoVClusterOrdered <- CCREannoV[rownames(mergeMat)[row_order(ht)]]
} else { # Cluster on V
  
  new_order <- readRDS(paste0(glue("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/ChIPSeeker_CCRE/{peak_type}/plots/"), "CCRE_anno_V_order.rds"))
  resultMat <- getAnnotationMatrix(annoV, repeatMaskerFeatures)
  
  # h <- ComplexHeatmap::Heatmap(resultMat, cluster_columns = FALSE, column_names_gp = grid::gpar(fontsize = 8),
  #                              row_names_gp = grid::gpar(fontsize = 0.5))
  # ht = draw(h)
  # plotFilePrefix <- paste0("/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/repeatMasker/TV", "/", filterStr)
  # plotNameClustered <- paste0(plotFilePrefix, "cluster_heatmap.pdf")
  # if (intermediatePlots) {
  #   pdf(plotNameClustered, height=18 * length(filteredTags) / 1444, width=8)
  #   print(ht)
  #   dev.off()
  # }

  # annoTClusterOrdered <- annoT[rownames(mergeMat)[row_order(ht)]]
  # annoVClusterOrdered <- annoV[rownames(mergeMat)[row_order(ht)]]
  # annoTClusterOrdered <- annoT[new_order]
  annoVClusterOrdered <- annoV[new_order]
  # CCREannoTClusterOrdered <- CCREannoT[rownames(mergeMat)[row_order(ht)]]
  # CCREannoVClusterOrdered <- CCREannoV[rownames(mergeMat)[row_order(ht)]]
}



extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}



# annoT.df <- list_to_dataframe(annoTClusterOrdered)
# annoT.df$Feature <- factor(annoT.df$Feature, levels = repeatMaskerCategoriesOrder)
# annoT.df$.id <- map_target_names(annoT.df$.id , target_pair_mapping_df)
# categoryColumn <- ".id"
# gpT <- plotAnnoBar.data.frame(annoT.df, categoryColumn=categoryColumn, colorOption = 4) + ggtitle("Treated") + scale_x_discrete(labels = TeX)
# shared_legend <- extract_legend(gpT)
# gpTNoLegnd <- gpT + theme(legend.position = "none")

annoV.df <- list_to_dataframe(annoVClusterOrdered)
annoV.df$Feature <- factor(annoV.df$Feature, levels = repeatMaskerCategoriesOrder)
annoV.df$.id <- map_target_names(annoV.df$.id , target_pair_mapping_df)
annoV.df$.id <- factor(annoV.df$.id , levels=rev(map_target_names(new_order, target_pair_mapping_df)))
categoryColumn <- ".id"
gpV <- plotAnnoBar.data.frame(annoV.df, categoryColumn=categoryColumn, colorOption = 4) + ggtitle("Viecle") + scale_x_discrete(labels = TeX)
pdf(paste0(out_dir, "repeatMasker_anno_V_clustered.pdf"), height=140, width=20)
gpV
dev.off()
write.csv(annoV.df, paste0(out_dir, "repeatMasker_anno_V_clustered.csv"), row.names = FALSE, quote = FALSE)

# gpVNoLegnd <- gpV + theme(legend.position = "none")

# resultPlot <- grid.arrange(arrangeGrob(gpVNoLegnd, gpTNoLegnd, ncol = 2),
#              shared_legend, ncol = 2, widths = c(8, 2))


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
# CCREannoT.df$Feature <- factor(CCREannoT.df$Feature, levels = CCREFeatures)
# categoryColumn <- ".id"
# gpT <- plotAnnoBar.data.frame(CCREannoT.df, categoryColumn=categoryColumn, colorOption = 1) + ggtitle("NaB_36x_CT_T")
# shared_legend <- extract_legend(gpT)
# gpTNoLegnd <- gpT + theme(legend.position = "none")

# CCREannoV.df <- list_to_dataframe(CCREannoVClusterOrdered)
# CCREannoV.df$Feature <- factor(CCREannoV.df$Feature, levels = CCREFeatures)
# categoryColumn <- ".id"
# gpV <- plotAnnoBar.data.frame(CCREannoV.df, categoryColumn=categoryColumn, colorOption = 1) + ggtitle("NaB_36x_CT_V")
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


getRevHeterotone <- function(tag) {
  tag <- as.character(tag)
  tags <- strsplit(tag, "-")[[1]]
  tag1 <- tags[1]
  tag2 <- tags[2]
  return(paste0(tag2, "-", tag1))
}


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
