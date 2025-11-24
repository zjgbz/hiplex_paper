source("/dcs05/hongkai/data/next_cutntag/script/peak_annotation/scripts_for_manuscript_v2/CCREUtils.R")
library(glue)

args <- commandArgs(trailingOnly = TRUE)
# peakType <- args[1]
# print(peakType)
# peakDir <- glue("/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/peak/{peakType}/")

peakDir <- args[1]
outDir <- args[2]
print(peakDir)
print(outDir)

tags <- c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
          "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP",
          "EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1",
          "EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "POLR2AphosphoS2", "cJun", "cFos", "Max", "Myc",
          "USF1", "USF2", "NRF1", "YY1", "H3K9ac")    #igG removed
merged <- FALSE

if (!merged) {
  tmpTags <- c()
  for (tagA in tags) {
    for (tagB in tags) {
      sortedTagATagB <- sort(c(tagA, tagB), method = "radix")
      tmpTags <- append(tmpTags, paste0(sortedTagATagB[1], "-", sortedTagATagB[2]))
    }
  }
  tags <- unique(tmpTags)
} else {
  tags <- tags
}

for (scen in scens) {
  peakAnnos <- list()
  for (tag in tags) {
    peakFile <- paste0(peakDir, scen, "/", tag, ".stringent.bed")
    print(peakFile)
    if (file.exists(peakFile)) {
      if (length(peakFile) > 0 && !is.na(file.size(peakFile)) && file.size(peakFile) > 0) {
        peak <- ChIPseeker::readPeakFile(peakFile, as = "GRanges")   # Use ChIPseeker to read peak files
        peakAnno <- annotatePeakByOverlappingClosestFeatureV2(peak, annotationRepeatMasker, repeatMaskerFeatures, featureColname = "X11")
        peakAnnos[[tag]] <- peakAnno
      } 
    } 
  }
  dir.create(outDir, recursive = TRUE)
  RDSFile <- glue("{outDir}/annotate.rds")
  # dir.create(paste0("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/test/", peakType, "/", scen), recursive = TRUE)
  # RDSFile <- paste0("/dcs05/hongkai/data/next_cutntag/bulk/peak_annotation/test/", peakType, "/", scen, "/annotate.rds")
  saveRDS(peakAnnos, RDSFile)

}
