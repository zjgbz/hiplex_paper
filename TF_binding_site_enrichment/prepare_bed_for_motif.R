library(rtracklayer)
library(GenomicRanges)
library(glue)

# input_file <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_row_clusters_gw_manually_reorganized_sparsity-10_cor_quantile-0.75_w_sparse.tsv"
# save_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/bed/reorganized_wgc_final/"
options <- commandArgs(trailingOnly = TRUE)
if (length(options) != 0) {
  input_file <- options[1]
  save_dir <- dirname(input_file)
  filename <- tools::file_path_sans_ext(basename(input_file))
  save_dir <- glue("{save_dir}/{filename}")

} else {
  filename <- tools::file_path_sans_ext(basename(input_file))
}
print(filename)

cluster_result <- read.table(input_file, sep = "\t", header = TRUE)
cluster_labels <- unique(cluster_result$label)


dir.create(save_dir, recursive = TRUE)
for (cluster_label in cluster_labels) {
  cluster_result_i <- cluster_result[cluster_result$label==cluster_label,]
  cluster_result_i <- as.data.frame(do.call(rbind, strsplit(cluster_result_i$feature, "[_]")))
  colnames(cluster_result_i) <- c("seqnames", "start", "end")
  cluster_result_i_gr <- makeGRangesFromDataFrame(cluster_result_i)
  export.bed(cluster_result_i_gr,glue("{save_dir}/{cluster_label}.bed"))
}
 