library(ggplot2)
library(ChIPseeker)
library("dplyr")
library(ComplexHeatmap)
library(grid)
library("gridExtra")
library(GenomicRanges)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tools)
library(stringr)

source("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/plotUtils.R")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

annotationFile <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/annotation_data/ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed" # https://downloads.wenglab.org/Registry-V3/Seven-Group/ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed
annotation_celltype_agnostic_file <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/annotation_data/GRCh38-cCREs.bed" # https://downloads.wenglab.org/V3/GRCh38-cCREs.bed
annotationFileChromHMM <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/chromHMM_annotation/annotation_data/hg38_genome_100_segments.bed" # https://github.com/ernstlab/full_stack_ChromHMM_annotations
annotationFileRepeatMasker <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/repeatMasker_annotation/annotation_data/hg38.fa.out" # https://www.repeatmasker.org/species/hg.html

annotationDf <- read.table(annotationFile,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
annotation <- makeGRangesFromDataFrame(annotationDf, 
                                       keep.extra.columns=TRUE,
                                       seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       strand.field="V6")
categories <- unique(annotation$V10)

annotation_celltype_agnostic_df <- read.table(annotation_celltype_agnostic_file,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
annotation_celltype_agnostic <- makeGRangesFromDataFrame(annotation_celltype_agnostic_df, 
                                       keep.extra.columns=TRUE,
                                       seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3")


annotationDfChromHMM <- read.table(annotationFileChromHMM,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
annotationChromHMM <- makeGRangesFromDataFrame(annotationDfChromHMM,
                                       keep.extra.columns=TRUE,
                                       seqnames.field="V1",
                                       start.field="V2",
                                       end.field="V3",
                                       strand.field="V6")
cleaned_strings <- gsub("\\d+_", "", annotationChromHMM$V4)

replace_strings <- function(strings, replacement_df) {
  # Match strings in 'original' column and replace with 'replacement'
  sapply(strings, function(x) {
    if (x %in% replacement_df$original) {
      replacement_df$replacement[match(x, replacement_df$original)]
    } else {
      x # Keep original if no match found
    }
  })
}

ChromHMM_state <- read.csv("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/utils/state_annotations_processed.csv")

replacement_df <- ChromHMM_state[, c("mneumonics", "Group")]
colnames(replacement_df) <- c("original", "replacement")
updated_strings <- replace_strings(cleaned_strings, replacement_df)
updated_strings <- unname(updated_strings)
annotationChromHMM$group <- updated_strings

annotationChromHMM$full_anno <- annotationChromHMM$V4
annotationChromHMM$full_anno <- gsub("^\\d+", "", annotationChromHMM$full_anno)
annotationChromHMM$full_anno <- gsub("_", "", annotationChromHMM$full_anno)

annotationChromHMM$V4 <- gsub("^\\d+|\\d+$", "", annotationChromHMM$V4)

annotationChromHMM$V4 <- gsub("_", "", annotationChromHMM$V4)
# annotationChromHMM$V4[annotationChromHMM$V4 %in% c("TxWk", "TxEx", "TxEnh")] <- "Tx"

# annotationChromHMM$V5 <- updated_strings

categoriesChromHMM <- unique(annotationChromHMM$V4)
# https://www.repeatmasker.org/species/hg.html
annotationRepeatMasker <- readRDS("/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/repeatMasker_annotation/annotation_data/repeatmasker.rds")
ignoreCTCFBound <- TRUE

setClass("csCCREAnno",
         representation=representation(
           annoStat="data.frame",
           peakNum="numeric",
           anno="GRanges"
         ))

getAnnoStatCCRE <- function(x) {
  if (!is(x, "csCCREAnno"))
    stop("not supported...")
  return(x@annoStat)
}

x <- new("csCCREAnno", annoStat = data.frame(a=c(1,2),b=c(3, 4)), peakNum=11)

if (ignoreCTCFBound) {
  newV10 <- sub(",.*", "", annotation$V10)
  annotation$V10 <- newV10
  annotation_celltype_agnostic$V6 <- sub(",.*", "", annotation_celltype_agnostic$V6)
  categories <- unique(annotation$V10)
  categoriesOrder <- c("dELS", "pELS", "PLS", "CTCF-only", "DNase-H3K4me3", "DNase-only", "Low-DNase")
  CCREFeatures <- c(categoriesOrder, "other")
  CCREChIPSeekerCategoriesOrder <- c(categoriesOrder, c("5' UTR", "3' UTR", "Exon", "Intron", "other"))
  ChIPSeekerCCRECategoriesOrder <- c("dELS", "pELS", "PLS", "5' UTR", "Exon", "Intron", "3' UTR", "DNase-H3K4me3", "Low-DNase", "DNase-only", "CTCF-only", "other")
  CCREChIPSeekerFeatures <- CCREChIPSeekerCategoriesOrder
  repeatMaskerFeatures <- unique(annotationRepeatMasker$X11)
  repeatMaskerCategoriesOrder <- c(repeatMaskerFeatures, "other")
}



annotatePeakByOverlapping <- function(peak, annotation, categories) {
  allPeakCts <- c()
  peakIndexSets <- c()
  for (category in categories) {
    categoryAnno <- annotation[annotation$V10 == category]
    overlapRes <- findOverlaps(peak, categoryAnno)
    overlapHits <- unique(queryHits(overlapRes))
    overlapHits <- overlapHits[!(overlapHits %in% peakIndexSets)]
    overlapCts <- length(overlapHits)
    allPeakCts <- append(allPeakCts, overlapCts)
    peakIndexSets <- unique(append(peakIndexSets, overlapHits))
  }
  allPeakCts <- append(allPeakCts, length(peak) - length(peakIndexSets))
  allPeakFreq <- allPeakCts / length(peak) * 100
  featureCol <- c(categories, "other")
  res <- data.frame(Feature=featureCol, Frequency = allPeakFreq)
  return(res)
}

annotatePeakByOverlappingCountRepeat <- function(peak, annotation, categories) {
  allPeakCts <- c()
  peakIndexSets <- c()
  for (category in categories) {
    categoryAnno <- annotation[annotation$V10 == category]
    overlapRes <- findOverlaps(peak, categoryAnno)
    overlapHits <- unique(queryHits(overlapRes))
    # overlapHits <- overlapHits[!(overlapHits %in% peakIndexSets)]
    overlapCts <- length(overlapHits)
    allPeakCts <- append(allPeakCts, overlapCts)
    peakIndexSets <- unique(append(peakIndexSets, overlapHits))
  }
  allPeakCts <- append(allPeakCts, length(peak) - length(peakIndexSets))
  allPeakFreq <- allPeakCts / sum(allPeakCts) * 100
  featureCol <- c(categories, "other")
  res <- data.frame(Feature=featureCol, Frequency = allPeakFreq)
  return(res)
}

findClosestGrange <- function(position, grs) {
  startsDistance <- abs(start(grs) - position)
  endsDistance <- abs(end(grs) - position)
  minDistance <- pmin(startsDistance, endsDistance)
  # distanceArray <- c(startsDistance, endsDistance)
  distanceArrayIndices <- 1:length(grs)
  closestGrangeIndex <- distanceArrayIndices[which.min(minDistance)]
  return(grs[closestGrangeIndex])
}

annotatePeakByOverlappingClosestFeatureNaive <- function(peak, annotation, categories) {
  # too slow dont run this, run annotatePeakByOverlappingClosestFeature
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation) 
  peakHits <- unique(queryHits(overlapRes))
  for (peakHit in peakHits) {
    overlapResHit <- overlapRes[queryHits(overlapRes) == peakHit]
    # print(peakHit)
    # print(peak)
    peakHit <- peak[peakHit]
    peakCenter <- ceiling((end(peakHit)-start(peakHit))/2 + start(peakHit))
    cloestAnno <- findClosestGrange(peakCenter, annotation[subjectHits(overlapResHit)])
    category <- cloestAnno$V10
    allPeakCtsList[[category]] <- allPeakCtsList[[category]] + 1
  }
  # print(allPeakCtsList)
  categories <- c(names(allPeakCtsList), "other")
  allPeakCts <- unname(unlist(allPeakCtsList))
  allPeakCts <- c(allPeakCts, (length(peak) - sum(allPeakCts)))
  allPeakFreq <- allPeakCts / sum(allPeakCts) * 100
  res <- data.frame(Feature=categories, Frequency = allPeakFreq)
  return(res)
}



annotatePeakByOverlappingClosestFeatureNaive2 <- function(peak, annotation, categories) {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation) 
  peakHitsUnique <- unique(queryHits(overlapRes))
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((end(peakHitsGR)-start(peakHitsGR))/2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)
  # cat("here\n")
  for (peakHit in peakHitsUnique) {
    overlapResHit <- overlapRes[peakHits == peakHit]
    annoHitsSelected <- annoHitsGR[peakHits == peakHit]
    minDistanceSelected <- minDistance[peakHits == peakHit]
    cloestAnno <- annoHitsSelected[which.min(minDistanceSelected)]
    category <- cloestAnno$V10
    allPeakCtsList[[category]] <- allPeakCtsList[[category]] + 1
    # if (category == "PLS") {
    #   if (sum(min(minDistanceSelected) == minDistanceSelected)>1) {
    #     print(peakHit)
    #     print(minDistanceSelected)
    #     print(annoHitsSelected$V10)
    #   }
    # }
  }
  # print(allPeakCtsList)
  categories <- c(names(allPeakCtsList), "other")
  allPeakCts <- unname(unlist(allPeakCtsList))
  allPeakCts <- c(allPeakCts, (length(peak) - sum(allPeakCts)))
  allPeakFreq <- allPeakCts / sum(allPeakCts) * 100
  res <- data.frame(Feature=categories, Frequency = allPeakFreq)
  return(res)
}

annotatePeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation) 
  peakHitsUnique <- unique(queryHits(overlapRes))
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((end(peakHitsGR)-start(peakHitsGR))/2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)
  
  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  # annoHitsDTGroup <- annoHitsDT[ , .SD[which.min(minDistance)], by = peakHits]
  annoHitsDT$queryStart <- start(peakHitsGR)
  annoHitsDT$queryEnd <- end(peakHitsGR)
  annoHitsDT$queryIndex <- peakHitsGR$index
  
  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]
  peakCategories <- annoHitsDTGroup[[featureColname]]
  
  peakCategoriesTable <- table(peakCategories)
  return(annoHitsDTGroup)
}

getAnnotateStatPeakByOverlappingClosestFeatureHelper <- function(peak, annotation, categories, featureColname = "V10") {
  allPeakCtsList <- list()
  for (category in categories) {
    allPeakCtsList[[category]] <- 0
  }
  overlapRes <- findOverlaps(peak, annotation) 
  peakHitsUnique <- unique(queryHits(overlapRes))
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  peakCenter <- ceiling((end(peakHitsGR)-start(peakHitsGR))/2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)
  
  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  # annoHitsDTGroup <- annoHitsDT[ , .SD[which.min(minDistance)], by = peakHits]
  
  annoHitsDTGroup <- annoHitsDT[annoHitsDT[ , .I[which.min(minDistance)], by = peakHits]$V1]
  
  peakCategories <- annoHitsDTGroup[[featureColname]]
  
  peakCategoriesTable <- table(peakCategories)
  return(peakCategoriesTable)
}
annotatePeakByOverlappingClosestFeature <- function(peak, annotation, categories, featureColname = "V10") {
  # almost 10 times faster than naive2. My take-away is : use less for-loop, use more vectorized operation and data.table
  peakCategoriesTable <- getAnnotateStatPeakByOverlappingClosestFeatureHelper(peak, annotation, categories, featureColname)
  otherLength <- length(peak)-sum(peakCategoriesTable)
  peakFreq <- c(unname(peakCategoriesTable), otherLength) / length(peak) * 100
  res <- data.frame(Feature=c(names(peakCategoriesTable), "other"), Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum=length(peak))
  return(x)
}




annotatePeakByOverlappingChIPSeekerCCRE <- function(peak, annotation, categories, featureColname="V10") {
  chipSeekerAnno <- annotatePeak(peak, TxDb=txdb, tssRegion=c(-1000, 1000), verbose=FALSE)
  chipSeekerAnnoFreq <-  c()
  chipSeekerAnnoFreq[allFeatures] <- 0
  chipSeekerAnnoFreq[as.character(chipSeekerAnno@annoStat$Feature)] <- chipSeekerAnno@annoStat$Frequency
  condensedFeatures <- c("5' UTR", "3' UTR", "Exon", "Intron")
  chipSeekerAnnoFreqCondensed <- c()
  chipSeekerAnnoFreqCondensed["5' UTR"] <- chipSeekerAnnoFreq["5' UTR"]
  chipSeekerAnnoFreqCondensed["3' UTR"] <- chipSeekerAnnoFreq["3' UTR"]
  chipSeekerAnnoFreqCondensed["Exon"] <- chipSeekerAnnoFreq["1st Exon"] + chipSeekerAnnoFreq["Other Exon"]
  chipSeekerAnnoFreqCondensed["Intron"] <- chipSeekerAnnoFreq["1st Intron"] + chipSeekerAnnoFreq["Other Intron"]
  chipSeekerAnnoFreqCondensed <- chipSeekerAnnoFreqCondensed * chipSeekerAnno@peakNum / 100
  
  chipSeekerAnno@anno$annotation <- gsub("\\s*\\([^\\)]+\\)","",chipSeekerAnno@anno$annotation)
  finalAnno <- rep("other", length(peak))
  finalAnno[chipSeekerAnno@anno$annotation %in% condensedFeatures] <- chipSeekerAnno@anno$annotation[chipSeekerAnno@anno$annotation %in% condensedFeatures]
  # print(finalAnno)
  chipSeekerAnno@anno$index <- 1: length(peak)
  disposeFeatures <- c("Promoter", "Downstream (<=300)", "Distal Intergenic")
  unmappedPeakGR <- chipSeekerAnno@anno[chipSeekerAnno@anno$annotation %in% disposeFeatures]
  if (length(unmappedPeakGR) > 0) {
    annoResult <- annotatePeakByOverlappingClosestFeatureHelper(unmappedPeakGR, annotation, categories, featureColname)
    # print(annoResult)
    peakCategories <- annoResult[[featureColname]]
    peakCategoriesTable <- table(peakCategories)
    finalAnno[annoResult$queryIndex] <- annoResult[[featureColname]]
  }
  
  if (length(unmappedPeakGR)==length(peak)) {
    otherCt <- length(peak) - sum(peakCategoriesTable)
    peakFeature <- c(names(peakCategoriesTable), "other")
    peakFreq <- c(unname(peakCategoriesTable), otherCt) / length(peak) * 100
  } else if (length(unmappedPeakGR)==0) {
    otherCt <- length(peak) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  } else {
    otherCt <- length(peak) - sum(peakCategoriesTable) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(peakCategoriesTable), names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(peakCategoriesTable), unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  }
  peak$annotation <- finalAnno
  res <- data.frame(Feature=peakFeature, Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum=length(peak), anno=peak)
  return(x)
}





annotatepeakByOverlappingRepeatMasker <- function(peak, annotationRepeatMasker, repeatMaskerFeatures) {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationRepeatMasker, repeatMaskerFeatures, featureColname = "X11"))
}

annotatepeakByOverlappingChromHMM <- function(peak, annotationChromHMM, categoriesChromHMM, featureColname = "V4") {
  return(annotatePeakByOverlappingClosestFeature(peak, annotationChromHMM, categoriesChromHMM, featureColname = featureColname))
}

annotatePeakFromDir <- function(peakDir, saveDir, annotation, categories, annotatePeakFUN = annotatePeakByOverlappingClosestFeature) {
  # peakDir <- c("/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/") # input dir 
  # saveDir <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/wgc/800/" # output dir
  fnames <- list.files(path = peakDir, pattern = "\\.bed$")
  peakAnnos <- list()
  for (fname in fnames) {
    peakFile <- paste0(peakDir, fname)
    peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")
    # annotatePeak(peak, annotation, categories)
    peakAnno <- annotatePeakFUN(peak, annotation, categories)
    peakAnnos[[fname]] <- peakAnno
  }
  anno <- lapply(peakAnnos, getAnnoStatCCRE)
  annoSorted <- anno[sort(names(anno))]
  anno.df <- list_to_dataframe(annoSorted)
  # anno.df <- list_to_dataframe(anno)
  anno.df$Feature <- factor(anno.df$Feature, levels = CCREChIPSeekerCategoriesOrder)
  categoryColumn <- ".id"
  p <- plotAnnoBar.data.frame(anno.df, categoryColumn=categoryColumn, colorOption = 1) + theme(panel.background = element_rect(fill='transparent'),
                                                                                               plot.background = element_rect(fill='transparent', color=NA),
                                                                                               panel.grid.major = element_blank(),
                                                                                               panel.grid.minor = element_blank(),
                                                                                               legend.background = element_rect(fill='transparent'),
                                                                                               legend.box.background = element_rect(fill='transparent'))
  return(p)
}

annotatePeakForLoopFromDir <- function(peakDir, saveDir, annotation, categories, annotatePeakFUN = annotatePeakByOverlappingClosestFeature) {
  # peakDir <- c("/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/") # input dir 
  # saveDir <- "/dcl02/hongkai/data/kyu/multitag_scripts/data_peakAnnotate/wgc/800/" # output dir
  fnames <- list.files(path = peakDir, pattern = "\\.bed$")
  peakAnnos <- list()
  for (fname in fnames) {
    peakFile <- paste0(peakDir, fname)
    peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")
    # annotatePeak(peak, annotation, categories)
    peakAnno <- annotatePeakFUN(peak, annotation, categories)
    # annotatepeakByOverlappingRepeatMasker(peak, annotationRepeatMasker, repeatMaskerFeatures)
    tagName <- gsub(".bed", "", fname)
    peakAnnos[[tagName]] <- peakAnno
  }
  anno <- lapply(peakAnnos, getAnnoStatCCRE)
  annoSorted <- anno[sort(names(anno))]
  for (annoName in names(anno)) {
    annoSingle <- annoSorted[annoName]
    anno.df <- list_to_dataframe(annoSingle)
    anno.df$Feature <- factor(anno.df$Feature, levels = c(categories, "other"))
    categoryColumn <- ".id"
    p <- plotAnnoBar.data.frame.one.target(anno.df,
                                      xlab="",
                                      ylab="Percentage(%)",
                                      title="Feature Distribution",
                                      categoryColumn = ".id", 
                                      colorOption = 1, 
                                      c(categories, "other"))
    saveFile <- paste0(saveDir, annoName, ".pdf")
    pdf(saveFile, width = 20, height = 5)
    print(p)
    dev.off()
  }
  return(NULL) # might need to modify this 
}


testing <- FALSE
if (testing) {
  peakDir <- "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/peak/data_peak_auc_0.05/"
  scen <- "V"
  peakFile <- paste0(peakDir, scen, "/", "H3K27me3-H3K4me3", ".stringent.bed")
  peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")
  start.time <- Sys.time()
  print(annotatePeakByOverlappingCCREChIPSeeker(peak, annotation, categories))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("annotatePeakByOverlappingClosestFeature ", time.taken, "\n")
  start.time <- Sys.time()
  print(annotatePeakByOverlappingClosestFeatureNew(peak, annotation, categories))
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("annotatePeakByOverlappingClosestFeatureTest ", time.taken, "\n")
}



