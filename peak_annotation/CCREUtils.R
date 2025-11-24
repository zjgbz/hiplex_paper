# ==============================================================================
# Peak Annotation Script for Hi-Plex CUT&Tag Data
# ==============================================================================
# This script provides functions for annotating genomic peaks using various
# reference annotations including CCRE, ChromHMM, and RepeatMasker databases.
# ==============================================================================

# Load Required Libraries ======================================================
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(GenomicRanges)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(tools)
library(stringr)

# Source utilities and load reference database
source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/plotUtils.R")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotation File Paths ========================================================
annotationFile <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/annotation_data/ENCFF414OGC_ENCFF806YEZ_ENCFF849TDM_ENCFF736UDR.7group.bed"
annotation_celltype_agnostic_file <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/ccre_annotation/annotation_data/GRCh38-cCREs.bed"
annotationFileChromHMM <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/chromHMM_annotation/annotation_data/hg38_genome_100_segments.bed"
annotationFileRepeatMasker <- "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/repeatMasker_annotation/annotation_data/hg38.fa.out"

# Load and Process CCRE Annotations ===========================================

# Cell-type specific CCRE annotations
annotationDf <- read.table(
  annotationFile,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

annotation <- makeGRangesFromDataFrame(
  annotationDf,
  keep.extra.columns = TRUE,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  strand.field = "V6"
)

# Cell-type agnostic CCRE annotations
annotation_celltype_agnostic_df <- read.table(
  annotation_celltype_agnostic_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

annotation_celltype_agnostic <- makeGRangesFromDataFrame(
  annotation_celltype_agnostic_df,
  keep.extra.columns = TRUE,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3"
)

# Load and Process ChromHMM Annotations ========================================

annotationDfChromHMM <- read.table(
  annotationFileChromHMM,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

annotationChromHMM <- makeGRangesFromDataFrame(
  annotationDfChromHMM,
  keep.extra.columns = TRUE,
  seqnames.field = "V1",
  start.field = "V2",
  end.field = "V3",
  strand.field = "V6"
)

# Clean ChromHMM state names
cleaned_strings <- gsub("\\d+_", "", annotationChromHMM$V4)

# Helper function to replace strings based on a lookup dataframe
replace_strings <- function(strings, replacement_df) {
  sapply(strings, function(x) {
    if (x %in% replacement_df$original) {
      replacement_df$replacement[match(x, replacement_df$original)]
    } else {
      x
    }
  })
}

# Load ChromHMM state annotations and create replacement mapping
ChromHMM_state <- read.csv(
  "/dcs05/hongkai/data/next_cutntag/script/peak_annotation/utils/state_annotations_processed.csv"
)

replacement_df <- ChromHMM_state[, c("mneumonics", "Group")]
colnames(replacement_df) <- c("original", "replacement")

# Apply replacements and clean up annotation fields
updated_strings <- replace_strings(cleaned_strings, replacement_df)
updated_strings <- unname(updated_strings)
annotationChromHMM$group <- updated_strings

# Store full annotation and create cleaned version
annotationChromHMM$full_anno <- gsub("^\\d+", "", annotationChromHMM$V4)
annotationChromHMM$full_anno <- gsub("_", "", annotationChromHMM$full_anno)
annotationChromHMM$V4 <- gsub("^\\d+|\\d+$", "", annotationChromHMM$V4)
annotationChromHMM$V4 <- gsub("_", "", annotationChromHMM$V4)

categoriesChromHMM <- unique(annotationChromHMM$V4)

# Load RepeatMasker Annotations ================================================
annotationRepeatMasker <- readRDS(
  "/dcl02/hongkai/data/kyu/multitag_scripts/scripts_peakAnnotation/repeatMasker_annotation/annotation_data/repeatmasker.rds"
)

# Configuration ================================================================
ignoreCTCFBound <- TRUE

# Define Custom S4 Class =======================================================
setClass(
  "csCCREAnno",
  representation = representation(
    annoStat = "data.frame",
    peakNum = "numeric",
    anno = "GRanges"
  )
)

# Accessor function for annotation statistics
getAnnoStatCCRE <- function(x) {
  if (!is(x, "csCCREAnno")) {
    stop("Input must be a csCCREAnno object")
  }
  return(x@annoStat)
}

# Process CTCF-bound categories if needed ======================================
if (ignoreCTCFBound) {
  # Remove CTCF-bound suffix from categories
  annotation$V10 <- sub(",.*", "", annotation$V10)
  annotation_celltype_agnostic$V6 <- sub(",.*", "", annotation_celltype_agnostic$V6)
  
  categories <- unique(annotation$V10)
  
  # Define category orders for plotting
  categoriesOrder <- c(
    "dELS", "pELS", "PLS", "CTCF-only",
    "DNase-H3K4me3", "DNase-only", "Low-DNase"
  )
  
  CCREFeatures <- c(categoriesOrder, "other")
  
  CCREChIPSeekerCategoriesOrder <- c(
    categoriesOrder,
    "5' UTR", "3' UTR", "Exon", "Intron", "other"
  )
  
  ChIPSeekerCCRECategoriesOrder <- c(
    "dELS", "pELS", "PLS",
    "5' UTR", "Exon", "Intron", "3' UTR",
    "DNase-H3K4me3", "Low-DNase", "DNase-only",
    "CTCF-only", "other"
  )
  
  CCREChIPSeekerFeatures <- CCREChIPSeekerCategoriesOrder
  repeatMaskerFeatures <- unique(annotationRepeatMasker$X11)
  repeatMaskerCategoriesOrder <- c(repeatMaskerFeatures, "other")
}

# Core Annotation Functions ====================================================

#' Annotate peaks by overlapping features (non-exclusive)
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @return Data frame with Feature and Frequency columns
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
  
  res <- data.frame(Feature = featureCol, Frequency = allPeakFreq)
  return(res)
}

#' Annotate peaks by overlapping features (allowing repeat counts)
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @return Data frame with Feature and Frequency columns
annotatePeakByOverlappingCountRepeat <- function(peak, annotation, categories) {
  allPeakCts <- c()
  peakIndexSets <- c()
  
  for (category in categories) {
    categoryAnno <- annotation[annotation$V10 == category]
    overlapRes <- findOverlaps(peak, categoryAnno)
    overlapHits <- unique(queryHits(overlapRes))
    overlapCts <- length(overlapHits)
    allPeakCts <- append(allPeakCts, overlapCts)
    peakIndexSets <- unique(append(peakIndexSets, overlapHits))
  }
  
  allPeakCts <- append(allPeakCts, length(peak) - length(peakIndexSets))
  allPeakFreq <- allPeakCts / sum(allPeakCts) * 100
  featureCol <- c(categories, "other")
  
  res <- data.frame(Feature = featureCol, Frequency = allPeakFreq)
  return(res)
}

#' Find the closest genomic range to a given position
#'
#' @param position Integer genomic position
#' @param grs GRanges object to search
#' @return GRanges object with the closest range
findClosestGrange <- function(position, grs) {
  startsDistance <- abs(start(grs) - position)
  endsDistance <- abs(end(grs) - position)
  minDistance <- pmin(startsDistance, endsDistance)
  distanceArrayIndices <- 1:length(grs)
  closestGrangeIndex <- distanceArrayIndices[which.min(minDistance)]
  return(grs[closestGrangeIndex])
}

#' Get annotation statistics using closest feature approach (optimized)
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param featureColname Column name containing feature labels
#' @return Table of peak category counts
getAnnotateStatPeakByOverlappingClosestFeatureHelper <- function(
    peak,
    annotation,
    categories,
    featureColname = "V10"
) {
  # Find all overlaps
  overlapRes <- findOverlaps(peak, annotation)
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  
  # Calculate distances from peak centers to annotation boundaries
  peakCenter <- ceiling((end(peakHitsGR) - start(peakHitsGR)) / 2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)
  
  # Convert to data.table for efficient grouping
  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  
  # Group by peak and select closest annotation
  annoHitsDTGroup <- annoHitsDT[
    annoHitsDT[, .I[which.min(minDistance)], by = peakHits]$V1
  ]
  
  # Extract categories and create table
  peakCategories <- annoHitsDTGroup[[featureColname]]
  peakCategoriesTable <- table(peakCategories)
  
  return(peakCategoriesTable)
}

#' Annotate detailed peak information using closest feature approach
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param featureColname Column name containing feature labels
#' @return Data.table with detailed annotation information
annotatePeakByOverlappingClosestFeatureHelper <- function(
    peak,
    annotation,
    categories,
    featureColname = "V10"
) {
  # Find all overlaps
  overlapRes <- findOverlaps(peak, annotation)
  peakHits <- queryHits(overlapRes)
  peakHitsGR <- peak[peakHits]
  annoHits <- subjectHits(overlapRes)
  annoHitsGR <- annotation[annoHits]
  
  # Calculate distances from peak centers to annotation boundaries
  peakCenter <- ceiling((end(peakHitsGR) - start(peakHitsGR)) / 2 + start(peakHitsGR))
  annoStarts <- start(annoHitsGR)
  annoEnds <- end(annoHitsGR)
  startsDistance <- abs(peakCenter - annoStarts)
  endsDistance <- abs(peakCenter - annoEnds)
  minDistance <- pmin(startsDistance, endsDistance)
  
  # Convert to data.table and add query information
  annoHitsDT <- as.data.table(annoHitsGR)
  annoHitsDT$peakHits <- peakHits
  annoHitsDT$minDistance <- minDistance
  annoHitsDT$queryStart <- start(peakHitsGR)
  annoHitsDT$queryEnd <- end(peakHitsGR)
  annoHitsDT$queryIndex <- peakHitsGR$index
  
  # Group by peak and select closest annotation
  annoHitsDTGroup <- annoHitsDT[
    annoHitsDT[, .I[which.min(minDistance)], by = peakHits]$V1
  ]
  
  return(annoHitsDTGroup)
}

#' Annotate peaks using closest overlapping feature
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param featureColname Column name containing feature labels
#' @return csCCREAnno object with annotation statistics
annotatePeakByOverlappingClosestFeature <- function(
    peak,
    annotation,
    categories,
    featureColname = "V10"
) {
  # Get category counts using optimized helper function
  peakCategoriesTable <- getAnnotateStatPeakByOverlappingClosestFeatureHelper(
    peak,
    annotation,
    categories,
    featureColname
  )
  
  # Calculate frequencies
  otherLength <- length(peak) - sum(peakCategoriesTable)
  peakFreq <- c(unname(peakCategoriesTable), otherLength) / length(peak) * 100
  
  # Create result data frame
  res <- data.frame(
    Feature = c(names(peakCategoriesTable), "other"),
    Frequency = peakFreq
  )
  
  # Return as custom S4 object
  x <- new("csCCREAnno", annoStat = res, peakNum = length(peak))
  return(x)
}



#' Annotate peaks using closest overlapping feature
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param featureColname Column name containing feature labels
#' @return csCCREAnno object with annotation statistics
annotatePeakByOverlappingClosestFeatureV2 <- function(
    peak,
    annotation,
    categories,
    featureColname = "V10"
) {
  # Get category counts using optimized helper function
  annoHitsDTGroup <- annotatePeakByOverlappingClosestFeatureHelper(
    peak,
    annotation,
    categories,
    featureColname
  )
  print(annoHitsDTGroup)
  # Extract categories and create table
  peakCategories <- annoHitsDTGroup[[featureColname]]
  peakCategoriesTable <- table(peakCategories)
  result_gr <- makeGRangesFromDataFrame(
    annoHitsDTGroup,
    keep.extra.columns = TRUE,  # This keeps all metadata!
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )

  # peakCategoriesTable <- getAnnotateStatPeakByOverlappingClosestFeatureHelper(
  #   peak,
  #   annotation,
  #   categories,
  #   featureColname
  # )
  
  # Calculate frequencies
  otherLength <- length(peak) - sum(peakCategoriesTable)
  peakFreq <- c(unname(peakCategoriesTable), otherLength) / length(peak) * 100
  
  # Create result data frame
  res <- data.frame(
    Feature = c(names(peakCategoriesTable), "other"),
    Frequency = peakFreq
  )
  
  # Return as custom S4 object
  x <- new("csCCREAnno", annoStat = res, peakNum = length(peak), anno=result_gr)
  return(x)
}

#' Annotate peaks using combined ChIPSeeker and CCRE annotations
#'
#' @param peak GRanges object containing peak regions
#' @param annotation GRanges object containing CCRE annotation regions
#' @param categories Character vector of CCRE categories
#' @param featureColname Column name containing feature labels
#' @return csCCREAnno object with combined annotation statistics
annotatePeakByOverlappingChIPSeekerCCRE <- function(
    peak,
    annotation,
    categories,
    featureColname = "V10"
) {
  # Run ChIPSeeker annotation
  chipSeekerAnno <- annotatePeak(
    peak,
    TxDb = txdb,
    tssRegion = c(-1000, 1000),
    verbose = FALSE
  )
  
  # Condense ChIPSeeker features
  condensedFeatures <- c("5' UTR", "3' UTR", "Exon", "Intron")
  
  chipSeekerAnnoFreqCondensed <- c()
  chipSeekerAnnoFreqCondensed["5' UTR"] <- 
    chipSeekerAnno@annoStat$Frequency[chipSeekerAnno@annoStat$Feature == "5' UTR"]
  chipSeekerAnnoFreqCondensed["3' UTR"] <- 
    chipSeekerAnno@annoStat$Frequency[chipSeekerAnno@annoStat$Feature == "3' UTR"]
  
  # Combine exon and intron annotations
  exonFreq <- sum(
    chipSeekerAnno@annoStat$Frequency[
      chipSeekerAnno@annoStat$Feature %in% c("1st Exon", "Other Exon")
    ]
  )
  intronFreq <- sum(
    chipSeekerAnno@annoStat$Frequency[
      chipSeekerAnno@annoStat$Feature %in% c("1st Intron", "Other Intron")
    ]
  )
  
  chipSeekerAnnoFreqCondensed["Exon"] <- exonFreq
  chipSeekerAnnoFreqCondensed["Intron"] <- intronFreq
  chipSeekerAnnoFreqCondensed <- chipSeekerAnnoFreqCondensed * chipSeekerAnno@peakNum / 100
  
  # Clean annotation labels
  chipSeekerAnno@anno$annotation <- gsub(
    "\\s*\\([^\\)]+\\)",
    "",
    chipSeekerAnno@anno$annotation
  )
  
  # Initialize final annotation
  finalAnno <- rep("other", length(peak))
  finalAnno[chipSeekerAnno@anno$annotation %in% condensedFeatures] <- 
    chipSeekerAnno@anno$annotation[chipSeekerAnno@anno$annotation %in% condensedFeatures]
  
  chipSeekerAnno@anno$index <- 1:length(peak)
  
  # Identify peaks to be annotated with CCRE
  disposeFeatures <- c("Promoter", "Downstream (<=300)", "Distal Intergenic")
  unmappedPeakGR <- chipSeekerAnno@anno[chipSeekerAnno@anno$annotation %in% disposeFeatures]
  
  # Annotate unmapped peaks with CCRE
  if (length(unmappedPeakGR) > 0) {
    annoResult <- annotatePeakByOverlappingClosestFeatureHelper(
      unmappedPeakGR,
      annotation,
      categories,
      featureColname
    )
    peakCategories <- annoResult[[featureColname]]
    peakCategoriesTable <- table(peakCategories)
    finalAnno[annoResult$queryIndex] <- annoResult[[featureColname]]
  }
  
  # Calculate final frequencies
  if (length(unmappedPeakGR) == length(peak)) {
    # All peaks annotated with CCRE
    otherCt <- length(peak) - sum(peakCategoriesTable)
    peakFeature <- c(names(peakCategoriesTable), "other")
    peakFreq <- c(unname(peakCategoriesTable), otherCt) / length(peak) * 100
  } else if (length(unmappedPeakGR) == 0) {
    # All peaks annotated with ChIPSeeker
    otherCt <- length(peak) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(names(chipSeekerAnnoFreqCondensed), "other")
    peakFreq <- c(unname(chipSeekerAnnoFreqCondensed), otherCt) / length(peak) * 100
  } else {
    # Mixed annotation
    otherCt <- length(peak) - sum(peakCategoriesTable) - sum(chipSeekerAnnoFreqCondensed)
    peakFeature <- c(
      names(peakCategoriesTable),
      names(chipSeekerAnnoFreqCondensed),
      "other"
    )
    peakFreq <- c(
      unname(peakCategoriesTable),
      unname(chipSeekerAnnoFreqCondensed),
      otherCt
    ) / length(peak) * 100
  }
  
  # Add annotation to peak object
  peak$annotation <- finalAnno
  
  # Create result
  res <- data.frame(Feature = peakFeature, Frequency = peakFreq)
  x <- new("csCCREAnno", annoStat = res, peakNum = length(peak), anno = peak)
  
  return(x)
}

#' Annotate peaks using RepeatMasker annotations
#'
#' @param peak GRanges object containing peak regions
#' @param annotationRepeatMasker GRanges object with RepeatMasker annotations
#' @param repeatMaskerFeatures Character vector of repeat types
#' @return csCCREAnno object with RepeatMasker annotation statistics
annotatepeakByOverlappingRepeatMasker <- function(
    peak,
    annotationRepeatMasker,
    repeatMaskerFeatures
) {
  return(
    annotatePeakByOverlappingClosestFeature(
      peak,
      annotationRepeatMasker,
      repeatMaskerFeatures,
      featureColname = "X11"
    )
  )
}

#' Annotate peaks using ChromHMM annotations
#'
#' @param peak GRanges object containing peak regions
#' @param annotationChromHMM GRanges object with ChromHMM annotations
#' @param categoriesChromHMM Character vector of ChromHMM states
#' @param featureColname Column name containing feature labels
#' @return csCCREAnno object with ChromHMM annotation statistics
annotatepeakByOverlappingChromHMM <- function(
    peak,
    annotationChromHMM,
    categoriesChromHMM,
    featureColname = "V4"
) {
  return(
    annotatePeakByOverlappingClosestFeature(
      peak,
      annotationChromHMM,
      categoriesChromHMM,
      featureColname = featureColname
    )
  )
}

# Batch Processing Functions ===================================================

#' Annotate all peak files in a directory and create summary plot
#'
#' @param peakDir Directory containing .bed peak files
#' @param saveDir Directory to save output plots
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param annotatePeakFUN Function to use for annotation
#' @return ggplot object with annotation summary
annotatePeakFromDir <- function(
    peakDir,
    saveDir,
    annotation,
    categories,
    annotatePeakFUN = annotatePeakByOverlappingClosestFeature
) {
  # Get all .bed files in directory
  fnames <- list.files(path = peakDir, pattern = "\\.bed$")
  peakAnnos <- list()
  
  # Process each peak file
  for (fname in fnames) {
    peakFile <- paste0(peakDir, fname)
    peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")
    peakAnno <- annotatePeakFUN(peak, annotation, categories)
    peakAnnos[[fname]] <- peakAnno
  }
  
  # Extract annotation statistics and sort
  anno <- lapply(peakAnnos, getAnnoStatCCRE)
  annoSorted <- anno[sort(names(anno))]
  
  # Convert to data frame for plotting
  anno.df <- list_to_dataframe(annoSorted)
  anno.df$Feature <- factor(anno.df$Feature, levels = CCREChIPSeekerCategoriesOrder)
  
  # Create bar plot
  categoryColumn <- ".id"
  p <- plotAnnoBar.data.frame(anno.df, categoryColumn = categoryColumn, colorOption = 1) +
    theme(
      panel.background = element_rect(fill = "transparent"),
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.background = element_rect(fill = "transparent"),
      legend.box.background = element_rect(fill = "transparent")
    )
  
  return(p)
}

#' Annotate all peak files in a directory and save individual plots
#'
#' @param peakDir Directory containing .bed peak files
#' @param saveDir Directory to save output plots
#' @param annotation GRanges object containing annotation regions
#' @param categories Character vector of annotation categories
#' @param annotatePeakFUN Function to use for annotation
#' @return NULL (saves plots to disk)
annotatePeakForLoopFromDir <- function(
    peakDir,
    saveDir,
    annotation,
    categories,
    annotatePeakFUN = annotatePeakByOverlappingClosestFeature
) {
  # Get all .bed files in directory
  fnames <- list.files(path = peakDir, pattern = "\\.bed$")
  peakAnnos <- list()
  
  # Process each peak file
  for (fname in fnames) {
    peakFile <- paste0(peakDir, fname)
    peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")
    peakAnno <- annotatePeakFUN(peak, annotation, categories)
    tagName <- gsub(".bed", "", fname)
    peakAnnos[[tagName]] <- peakAnno
  }
  
  # Extract annotation statistics and sort
  anno <- lapply(peakAnnos, getAnnoStatCCRE)
  annoSorted <- anno[sort(names(anno))]
  
  # Create and save individual plots
  for (annoName in names(anno)) {
    annoSingle <- annoSorted[annoName]
    anno.df <- list_to_dataframe(annoSingle)
    anno.df$Feature <- factor(anno.df$Feature, levels = c(categories, "other"))
    
    categoryColumn <- ".id"
    p <- plotAnnoBar.data.frame.one.target(
      anno.df,
      xlab = "",
      ylab = "Percentage(%)",
      title = "Feature Distribution",
      categoryColumn = ".id",
      colorOption = 1,
      c(categories, "other")
    )
    
    saveFile <- paste0(saveDir, annoName, ".pdf")
    pdf(saveFile, width = 20, height = 5)
    print(p)
    dev.off()
  }
  
  return(NULL)
}

# ==============================================================================
# End of Script
# ==============================================================================