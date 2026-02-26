# /dcs05/hongkai/data/next_cutntag/script/dna_methylation/calculate_methylation_percentage.Rrm(list=ls())
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Bioconductor packages
    bioc_packages <- c("GenomicRanges", "plyranges", "rtracklayer")
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

packages <- c("arrow", "ggplot2", "GenomicRanges", "plyranges", 
              "rtracklayer", "data.table")
lapply(packages, install_if_missing)

suppressPackageStartupMessages({
  library(arrow)
  library(ggplot2)
  library(GenomicRanges)
  library(plyranges)
  library(rtracklayer)
  library(data.table)
})

methyl_bed_file <- "../data/peak_annotation/ENCFF660IHA.bed" 
methyl_bed_df <- read.table(methyl_bed_file, sep="\t")
# methyl_bed_df <- as.data.frame(data.table::fread(methyl_bed_file))
methyl_granges <- makeGRangesFromDataFrame(methyl_bed_df, start.field = "V2", end.field = "V3", seqnames.field = "V1", strand.field = "V6", keep.extra.columns = TRUE)
methyl_granges$unmethylated <- methyl_granges$V10 - methyl_granges$V10 * methyl_granges$V11 / 100
methyl_granges$methylated <- methyl_granges$V10 * methyl_granges$V11 / 100
calculate_methylation_percentage <- function(search_window_gr, ignore_no_cpg_peaks = TRUE) {
  res <- findOverlaps(search_window_gr, methyl_granges, ignore.strand = TRUE)
  overlapped_methyl_granges <- methyl_granges[subjectHits(res)]
  overlapped_methyl_granges$search_window_index <- queryHits(res)
  overlapped_methyl_dt <- as.data.table(overlapped_methyl_granges)
  avg_percentage_methyl <- aggregate(overlapped_methyl_dt[, "V11"], list(overlapped_methyl_dt$search_window_index), mean)
  if (ignore_no_cpg_peaks) {
    search_window_gr <- search_window_gr[avg_percentage_methyl$Group.1]
    search_window_gr$methyl_level <- (avg_percentage_methyl$V11/100)
  } else {
    methy_array <- rep(0, length(search_window_gr))
    methy_array[avg_percentage_methyl$Group.1] <- avg_percentage_methyl[["V11"]]
    methy_array <- methy_array/100
    search_window_gr$methyl_level <- methy_array
  }
  return(search_window_gr)
}

calculate_M_values <- function(search_window_gr, ignore_no_cpg_peaks = TRUE) {
  res <- findOverlaps(search_window_gr, methyl_granges, ignore.strand = TRUE)
  overlapped_methyl_granges <- methyl_granges[subjectHits(res)]
  overlapped_methyl_granges$search_window_index <- queryHits(res)
  overlapped_methyl_dt <- as.data.table(overlapped_methyl_granges)
  avg_percentage_methyl <- aggregate(overlapped_methyl_dt[, "V11"], list(overlapped_methyl_dt$search_window_index), mean)
  # sum_methyl <- aggregate(overlapped_methyl_dt[, "methylated"], list(overlapped_methyl_dt$search_window_index), sum)
  # sum_unmethyl <- aggregate(overlapped_methyl_dt[, "unmethylated"], list(overlapped_methyl_dt$search_window_index), sum)
  sum_methyl <- aggregate(overlapped_methyl_dt[, c("methylated", "unmethylated")], list(overlapped_methyl_dt$search_window_index), sum)
  if (ignore_no_cpg_peaks) {
    search_window_gr <- search_window_gr[avg_percentage_methyl$Group.1]
    search_window_gr$M_values <- log2((sum_methyl$methylated+1)/(sum_methyl$unmethylated+1))
    
  } else {
    methy_array <- rep(0, length(search_window_gr))
    methy_array[avg_percentage_methyl$Group.1] <- log2((sum_methyl$methylated+1)/(sum_methyl$unmethylated+1))
    search_window_gr$M_values <- methy_array
  }
  return(search_window_gr)
}



calculate_methylation_percentage_from_promoter_file <- function(promoter_file) {
  promoter_df <- read.table(promoter_file, header = TRUE, sep = "\t")
  promoter_df$seqnames <- paste0("chr", promoter_df$seqnames)
  search_window_gr <- makeGRangesFromDataFrame(promoter_df, start.field = "start", end.field = "end", seqnames.field = "seqnames", strand.field = "strand", keep.extra.columns = TRUE)
  return(calculate_methylation_percentage(search_window_gr))
}