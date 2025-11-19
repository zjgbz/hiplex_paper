#usage:#Rscript code_row.R "/Users/kyu/Documents/thesis/data/bigwig_normalized/V" "/Users/kyu/Documents/thesis/data/peak/data_peak_auc_0.05/V" "H3K27me3" "H3K4me3"
args <- commandArgs(trailingOnly = TRUE)

if (FALSE) {
  args <- c(
    "/Users/kyu/Documents/thesis/data/bigwig_normalized/V",
    "/Users/kyu/Documents/thesis/data/peak/data_peak_auc_0.05_extend_window_10000_binomial/V/",
    "H3K27me3", "H3K27ac", "#97BE5A", "#f5ac8e", "#2878B5", "1", 
    "/Users/kyu/Documents/thesis/scripts/enriched_heatmap/result/data_peak_auc_0.05_ctrl-singletone_stringent/"
  )
}

dir_bw <- args[1]
dir_peak <- args[2]
target1_name <- args[3]
target2_name <- args[4]
colt_t12 <- args[5]
col_t1 <- args[6]
col_t2 <- args[7]
maxlimit<- as.numeric(args[8])
out_dir <- args[9]
subsample_t1 <- args[10]
subsample_t2 <- args[11]

suppressPackageStartupMessages({
  library(GenomicAlignments)
  #library(ChIPseeker)
  library(Biostrings)
  library(EnrichedHeatmap)
  library(rtracklayer)
  library(data.table)
  library(circlize)
  library(glue)
})

hide_top_annotation <- TRUE

##### normed
peak1_files <- c(paste0(dir_peak,"/",target1_name,"-",target1_name, ".stringent.bed"),
                 paste0(dir_peak,"/",target1_name,"-",target1_name, ".relaxed.bed"))
peak2_files <- c(paste0(dir_peak,"/",target2_name,"-",target2_name, ".stringent.bed"),
                 paste0(dir_peak,"/",target2_name,"-",target2_name, ".relaxed.bed"))

print(peak1_files)
print(peak2_files)

peak1_file <- peak1_files[file.exists(peak1_files)]
peak2_file <- peak2_files[file.exists(peak2_files)]

target1_file <- paste0(dir_bw,"/",target1_name,"-",target1_name, ".bw")
target2_file <- paste0(dir_bw,"/",target2_name,"-",target2_name,".bw")

# Create targets and other necessary objects
peaks1 <- fread(peak1_file, header = FALSE, data.table = FALSE)
peaks2 <- fread(peak2_file, header = FALSE, data.table = FALSE)

peak1 <- makeGRangesFromDataFrame(
  df = peaks1,
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

peak2 <- makeGRangesFromDataFrame(
  df = peaks2,
  seqnames.field = "V1", start.field = "V2", end.field = "V3")

target1_signal <- rtracklayer::import(target1_file,
                                      format = "BigWig")
target2_signal <- rtracklayer::import(target2_file,
                                      format = "BigWig")

#log transformation, exp base -> log(x+1)
target1_signal$score <- log10(target1_signal$score+1)
target2_signal$score <- log10(target2_signal$score+1)

#peaks
overlap <- findOverlaps(peak1, peak2)

if (FALSE) {
  InterPeaks <- append(peak1[queryHits(overlap)],peak2[subjectHits(overlap)])
  InterPeaks <- GenomicRanges::reduce(InterPeaks)
} else {
  InterPeaks <- pintersect(peak1[queryHits(overlap)], peak2[subjectHits(overlap)])
}

peak1_no_overlap <- setdiff(peak1, InterPeaks)
peak2_no_overlap <- setdiff(peak2, InterPeaks)
peak1_no_overlap_length <- length(peak1_no_overlap)
peak2_no_overlap_length <- length(peak2_no_overlap)
interpeak_length <- length(InterPeaks)
partition_percentage <- c(interpeak_length, peak1_no_overlap_length, peak2_no_overlap_length)
partition_percentage <- partition_percentage/sum(partition_percentage)*100
names(partition_percentage) <- c("Intersection",glue("{target1_name} only"), glue("{target2_name} only"))

# Subsample peak1_no_overlap and peak2_no_overlap if subsample parameters are provided
if (!is.na(subsample_t1) && subsample_t1 != "") {
  subsample_ratio_t1 <- as.numeric(subsample_t1)
  if (subsample_ratio_t1 > 0 && subsample_ratio_t1 < 1) {
    n_subsample_t1 <- round(length(peak1_no_overlap) * subsample_ratio_t1)
    set.seed(42)  # For reproducibility
    subsample_idx_t1 <- sample(1:length(peak1_no_overlap), n_subsample_t1)
    peak1_no_overlap <- peak1_no_overlap[subsample_idx_t1]
    print(paste0("Subsampled ", target1_name, " only peaks to ", subsample_ratio_t1*100, "% (", n_subsample_t1, " peaks)"))
  }
}

if (!is.na(subsample_t2) && subsample_t2 != "") {
  subsample_ratio_t2 <- as.numeric(subsample_t2)
  if (subsample_ratio_t2 > 0 && subsample_ratio_t2 < 1) {
    n_subsample_t2 <- round(length(peak2_no_overlap) * subsample_ratio_t2)
    set.seed(43)  # Different seed for independence
    subsample_idx_t2 <- sample(1:length(peak2_no_overlap), n_subsample_t2)
    peak2_no_overlap <- peak2_no_overlap[subsample_idx_t2]
    print(paste0("Subsampled ", target2_name, " only peaks to ", subsample_ratio_t2*100, "% (", n_subsample_t2, " peaks)"))
  }
}

# Check if there is no overlap
has_overlap <- length(InterPeaks) > 0

if (has_overlap) {
  # Original code path: with overlap
  peaks_df <- rbind(peaks1, peaks2)
  peaks <- makeGRangesFromDataFrame(
    df = peaks_df,
    seqnames.field = "V1", start.field = "V2", end.field = "V3")
  
  peaksTri <- c(InterPeaks, peak1_no_overlap, peak2_no_overlap)
  
  #mat for plotting
  mat_target1_signal <- normalizeToMatrix(signal = target1_signal,
                                          target = resize(peaksTri, fix = "center", width = 1),
                                          background = 0,
                                          mean_mode = "w0",
                                          value_column = "score",
                                          extend = 2000,
                                          smooth = FALSE)
  
  mat_target2_signal <- normalizeToMatrix(signal =  target2_signal,
                                          target = resize(peaksTri, fix = "center", width = 1),
                                          background = 0,
                                          mean_mode = "w0",
                                          value_column = "score",
                                          extend = 2000,
                                          smooth = FALSE)
  
  #prepare order for peak 1 and peak 2
  axis_name <- c("-2000", "center", "2000")
  mat_pre1 <- normalizeToMatrix(signal = target1_signal,
                                target = resize(peak1_no_overlap, fix = "center", width = 1),
                                background = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = 2000,
                                smooth = FALSE)
  
  mat_pre2 <- normalizeToMatrix(signal =  target2_signal,
                                target = resize(peak2_no_overlap, fix = "center", width = 1),
                                background = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = 2000,
                                smooth = FALSE)
  
  draw1 <-  EnrichedHeatmap(mat_pre1,
                            column_title = target1_name,
                            axis_name = axis_name,
                            name = target1_name)
  
  draw2 <-  EnrichedHeatmap(mat_pre2,
                            column_title = target2_name,
                            axis_name = axis_name,
                            name = target2_name)
  
  ht1 = draw(draw1); ht2 = draw(draw2);   
  
  ##intersect order
  mat_pre1Inter <- normalizeToMatrix(signal = target1_signal,
                                     target = resize(InterPeaks, fix = "center", width = 1),
                                     background = 0,
                                     mean_mode = "w0",
                                     value_column = "score",
                                     extend = 2000,
                                     smooth = FALSE)
  
  mat_pre2Inter <- normalizeToMatrix(signal =  target2_signal,
                                     target = resize(InterPeaks, fix = "center", width = 1),
                                     background = 0,
                                     mean_mode = "w0",
                                     value_column = "score",
                                     extend = 2000,
                                     smooth = FALSE)
  
  draw1Inter <-  EnrichedHeatmap(mat_pre1Inter,
                                 column_title = target1_name,
                                 axis_name = axis_name,
                                 name = target1_name)
  
  draw2Inter <-  EnrichedHeatmap(mat_pre2Inter,
                                 column_title = target2_name,
                                 axis_name = axis_name,
                                 name = target2_name)
  
  ht1Inter = draw(draw1Inter); ht2Inter = draw(draw1Inter); 
  
  order1Inter <- row_order(ht1Inter)
  order2Inter <- row_order(ht2Inter)
  
  average_ranksInter <- numeric(length(order1Inter)) #use average rank as order
  for (i in 1:length(order1Inter)) {
    rank_1 <- which(order1Inter == i)
    rank_2 <- which(order2Inter == i)
    average_ranksInter[i] <- mean(c(rank_1, rank_2))
  }
  ordersInter <- order(average_ranksInter)
  
  #combine orders
  orders <- c(ordersInter, length(ordersInter)+ row_order(ht1),
              length(peak1_no_overlap)+length(ordersInter)+ row_order(ht2))
  
  partition = rep(NA, length(orders))
  
  if (target1_name=="POLR2AphosphoS2"){
    target1_name="RNAPII"
  }
  if (target2_name=="POLR2AphosphoS2"){
    target2_name="RNAPII"
  }
  
  partition[ordersInter] <- "Intersection"
  partition[length(ordersInter)+ row_order(ht1)] <- glue("{target1_name} only")
  partition[length(peak1_no_overlap)+length(ordersInter)+ row_order(ht2)] <- glue("{target2_name} only")

  # partition_percentage <- table(partition)/sum(table(partition)) * 100
  
  name_list <- c("Intersection",glue("{target1_name} only"), glue("{target2_name} only"))
  
  for (partition_name in name_list) {
    partition[partition == partition_name] <- glue("{partition_name}\n{round(partition_percentage[partition_name], 2)}%")
  }
  
  partition <- factor(partition, levels = glue("{name_list}\n{round(partition_percentage[name_list], 2)}%"))
  
  cluster_colors <- c(colt_t12, col_t1, col_t2)
  
} else {
  # No overlap case: only target1 and target2 regions
  print("No overlap detected between targets. Plotting only target1 and target2 specific regions.")
  
  peaksBi <- c(peak1_no_overlap, peak2_no_overlap)
  
  #mat for plotting
  mat_target1_signal <- normalizeToMatrix(signal = target1_signal,
                                          target = resize(peaksBi, fix = "center", width = 1),
                                          background = 0,
                                          mean_mode = "w0",
                                          value_column = "score",
                                          extend = 2000,
                                          smooth = FALSE)
  
  mat_target2_signal <- normalizeToMatrix(signal =  target2_signal,
                                          target = resize(peaksBi, fix = "center", width = 1),
                                          background = 0,
                                          mean_mode = "w0",
                                          value_column = "score",
                                          extend = 2000,
                                          smooth = FALSE)
  
  #prepare order for peak 1 and peak 2
  axis_name <- c("-2000", "center", "2000")
  mat_pre1 <- normalizeToMatrix(signal = target1_signal,
                                target = resize(peak1_no_overlap, fix = "center", width = 1),
                                background = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = 2000,
                                smooth = FALSE)
  
  mat_pre2 <- normalizeToMatrix(signal =  target2_signal,
                                target = resize(peak2_no_overlap, fix = "center", width = 1),
                                background = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = 2000,
                                smooth = FALSE)
  
  draw1 <-  EnrichedHeatmap(mat_pre1,
                            column_title = target1_name,
                            axis_name = axis_name,
                            name = target1_name)
  
  draw2 <-  EnrichedHeatmap(mat_pre2,
                            column_title = target2_name,
                            axis_name = axis_name,
                            name = target2_name)
  
  ht1 = draw(draw1); ht2 = draw(draw2);
  
  #combine orders (no intersection orders)
  orders <- c(row_order(ht1), length(peak1_no_overlap) + row_order(ht2))
  
  partition = rep(NA, length(orders))
  
  if (target1_name=="POLR2AphosphoS2"){
    target1_name="RNAPII"
  }
  if (target2_name=="POLR2AphosphoS2"){
    target2_name="RNAPII"
  }
  
  partition[row_order(ht1)] <- glue("{target1_name} only")
  partition[length(peak1_no_overlap) + row_order(ht2)] <- glue("{target2_name} only")
  
  partition_percentage <- table(partition)/sum(table(partition)) * 100
  
  name_list <- c(glue("{target1_name} only"), glue("{target2_name} only"))
  
  for (partition_name in name_list) {
    partition[partition == partition_name] <- glue("{partition_name}\n{round(partition_percentage[partition_name], 2)}%")
  }
  
  partition <- factor(partition, levels = glue("{name_list}\n{round(partition_percentage[name_list], 2)}%"))
  
  cluster_colors <- c(col_t1, col_t2)
}

#############
print(paste0(out_dir, "/", target1_name, "_", target2_name, ".pdf"))
pdf(paste0(out_dir, "/", target1_name, "_", target2_name, ".pdf"), width = 6, height = 10)

# Define annotation range
annotation_min <- 0
annotation_max <- maxlimit

line_cluster_colors <- cluster_colors
if (has_overlap) {
  line_cluster_colors[1] <- "#ffffff00"
}

if (!hide_top_annotation) {
  top_anno1 <-
    HeatmapAnnotation(lines = anno_enriched(
      ylim = c(0, annotation_max),
      gp = gpar(col = line_cluster_colors),
      axis_param = list(side = "left")
    ))
  
  top_anno2 <- HeatmapAnnotation(lines = anno_enriched(
    ylim = c(0, annotation_max),
    gp = gpar(col = line_cluster_colors),
    axis_param = list(side = "right")
  ))
} else {
  top_anno1 <- NULL
  top_anno2 <- NULL
}

ht_list = Heatmap(partition,
                  row_title_gp = gpar(col = cluster_colors),
                  col = cluster_colors, show_heatmap_legend = FALSE, show_row_names = FALSE,   
                  show_column_names = FALSE, width = unit(3, "mm")) +
  EnrichedHeatmap(mat_target1_signal,
                  column_title = target1_name,
                  axis_name = axis_name,
                  row_title = NULL,
                  show_row_names = F,
                  cluster_rows = FALSE, show_row_dend = FALSE,
                  top_annotation = top_anno1,
                  heatmap_legend_param = list(legend_width = unit(6, "cm"), title = NULL,
                                              color_bar = "continuous", legend_direction = "horizontal",
                                              legend_width = unit(5, "cm")),
                  row_order = orders) +
  EnrichedHeatmap(mat_target2_signal,
                  column_title = target2_name,
                  axis_name = axis_name,
                  row_title = NULL,
                  show_row_names = F,
                  cluster_rows = FALSE, show_row_dend = FALSE,
                  top_annotation = top_anno2,
                  heatmap_legend_param = list(legend_width = unit(6, "cm"), title = NULL,
                                              color_bar = "continuous", legend_direction = "horizontal",
                                              legend_width = unit(5, "cm")),
                  row_order = orders) 

draw(ht_list, heatmap_legend_side = "bottom", annotation_legend_side = "right",
     split = partition)

dev.off()