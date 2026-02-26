install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Bioconductor packages
    bioc_packages <- c("GenomicRanges", "ComplexHeatmap")
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

# 所有需要的包（确保与下面library()一致）
packages <- c(
  # CRAN packages
  "arrow", "glue", "tibble", "ggpointdensity", "viridis",
  "ggplot2", "ggpubr", "ggrastr", "data.table", "dplyr", 
  "circlize", "svglite", "tidyr", "latex2exp",
  # Bioconductor packages
  "GenomicRanges", "ComplexHeatmap"
)

lapply(packages, install_if_missing)

suppressPackageStartupMessages({
    library(arrow)
    library(GenomicRanges)
    library(glue)
    library(tibble)
    library(ggpointdensity)
    library(viridis)
    library(ComplexHeatmap)
    library(ggplot2)
    library(ggpubr)
    library(ggrastr)
    library(data.table)
    library(dplyr)
    library(circlize)
    library(svglite)
    library(tidyr)
    library(latex2exp)
})

source("../hiplex_paper/utils.R")

options <- commandArgs(trailingOnly = TRUE)
if (length(options) != 0) {
  wgc_col_cluster_dir <- options[1]
  out_dir <- options[2]
}

method = "limma"
mean_per_thres = 0.25
fdr_thres = "0.25"
cluster_idx_list = c(1:8, 10:15)
# cluster_idx_list = c(1:8, 10, 11, 13, 14)
width_base = 8
col_size_coef = 20
colnames_fontsize = 10
random_seed = 42
font_size = 50
l2fc_thres = 0.5
show_heatmap_legend = "on"
show_colnames = "off"
if (show_heatmap_legend == "on") {
    show_heatmap_legend_bool = TRUE
} else if (show_heatmap_legend == "off") {
    show_heatmap_legend_bool = FALSE
}

wgc_dir = paste0(out_dir, "/column_cluster_wgc")
fig_dir = paste0(out_dir, "/column_cluster_fig")

col_cluster_full_filename = "bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_column_clusters_manually_reorganized.tsv"
col_cluster_full_dir_filename = file.path(wgc_col_cluster_dir, col_cluster_full_filename)
col_cluster_full = read.table(col_cluster_full_dir_filename, header=TRUE, sep="\t", row.names=NULL)

for (cluster_idx in cluster_idx_list) {
    # Plot heatmaps for log2 data
    V1_log2_wgc_filename = glue("V1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
    V2_log2_wgc_filename = glue("V2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
    T1_log2_wgc_filename = glue("T1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
    T2_log2_wgc_filename = glue("T2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
    V1_log2_wgc_dir_filename = file.path(wgc_dir, V1_log2_wgc_filename)
    V2_log2_wgc_dir_filename = file.path(wgc_dir, V2_log2_wgc_filename)
    T1_log2_wgc_dir_filename = file.path(wgc_dir, T1_log2_wgc_filename)
    T2_log2_wgc_dir_filename = file.path(wgc_dir, T2_log2_wgc_filename)
    V1_log2_wgc = read_feather(V1_log2_wgc_dir_filename)
    V2_log2_wgc = read_feather(V2_log2_wgc_dir_filename)
    T1_log2_wgc = read_feather(T1_log2_wgc_dir_filename)
    T2_log2_wgc = read_feather(T2_log2_wgc_dir_filename)
    V1_log2_wgc = column_to_rownames(V1_log2_wgc, var="pos")
    V2_log2_wgc = column_to_rownames(V2_log2_wgc, var="pos")
    T1_log2_wgc = column_to_rownames(T1_log2_wgc, var="pos")
    T2_log2_wgc = column_to_rownames(T2_log2_wgc, var="pos")
    colnames(V1_log2_wgc) <- map_target_names(colnames(V1_log2_wgc), target_pair_mapping_df)
    colnames(V2_log2_wgc) <- map_target_names(colnames(V2_log2_wgc), target_pair_mapping_df)
    colnames(T1_log2_wgc) <- map_target_names(colnames(T1_log2_wgc), target_pair_mapping_df)
    colnames(T2_log2_wgc) <- map_target_names(colnames(T2_log2_wgc), target_pair_mapping_df)

    df_V1 <- get_cluster_df(V1_log2_wgc, "V1")
    df_V2 <- get_cluster_df(V2_log2_wgc, "V2")
    df_T1 <- get_cluster_df(T1_log2_wgc, "T1")
    df_T2 <- get_cluster_df(T2_log2_wgc, "T2")

    col_cluster_df <- bind_rows(df_V1, df_V2, df_T1, df_T2)
    col_cluster_df <- col_cluster_df %>%
        mutate(label = factor(label, levels = c("V1", "V2", "T1", "T2")))
    col_cluster_df <- col_cluster_df %>%
        mutate(order = match(nonprefix, col_cluster_full$feature))
    col_cluster_df <- col_cluster_df %>% arrange(label, order)

    col_order = col_cluster_df$feature
    col_split = col_cluster_df$label

    colnames(V1_log2_wgc) = paste0("V1:", colnames(V1_log2_wgc))
    colnames(V2_log2_wgc) = paste0("V2:", colnames(V2_log2_wgc))
    colnames(T1_log2_wgc) = paste0("T1:", colnames(T1_log2_wgc))
    colnames(T2_log2_wgc) = paste0("T2:", colnames(T2_log2_wgc))

    wgc_log2_cbind = cbind(V1_log2_wgc, V2_log2_wgc, T1_log2_wgc, T2_log2_wgc)
    wgc_log2_cbind = wgc_log2_cbind[, col_order]

    # col_fun = colorRamp2(c(min(wgc_log2_cbind), 0.9, 1.8), c("blue", "white", "red"))
    # col_fun = colorRamp2(c(min(wgc_log2_cbind), 0.9, 1.8), c("#3155C3", "white", "#AF0525"))

    if (show_colnames == "on") {
        show_colnames_bool = TRUE
    } else if (show_colnames == "off") {
        show_colnames_bool = FALSE
    }
    heatmap_prefix = glue("diff_heatmap_col_cluster-{cluster_idx}_colname-{show_colnames}_col-reorder_size-{col_size_coef}_limmanorm-manually")

    col_num = ncol(wgc_log2_cbind)
    # if (col_num > 80) {
    #     col_size_coef = 80
    # } else if (col_num < 60 && col_num > 10) {
    #     col_size_coef = 60
    # } else if (col_num < 10) {
    #     col_size_coef = 25
    # } else {
    #     col_size_coef = ncol(wgc_log2_cbind)
    # }

    ht = Heatmap(as.matrix(wgc_log2_cbind), name="log2",
        # col = col_fun,
    	show_row_names=FALSE,
    	show_column_names=show_colnames_bool,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        column_order = col_order,
        column_split = col_split,
    	width=col_size_coef * unit(width_base, "mm"),
    	# height=50 * unit(10, "mm"),
        height= (3750 / 89) * unit(10, "mm"),
        column_title = c("C1", "C2", "T1", "T2"),
        column_title_gp=gpar(col = c("#3155C3", "#3155C3", "#AF0525", "#AF0525"), fontsize=90),
        column_gap=unit(8, "mm"),
        column_names_gp=gpar(fontsize=colnames_fontsize),
    	show_row_dend=FALSE,
        show_heatmap_legend=show_heatmap_legend_bool,
        heatmap_legend_param = list(
        		title = "log2",
                grid_width = 5*unit(5, "mm"),
                legend_height = 77 * 1.3*unit(5, "mm"),
        	    title_gp = gpar(fontsize = 10),
        	    labels_gp = gpar(fontsize = 40),
                legend_direction = "vertical"),
        use_raster=TRUE
    	)

    set.seed(random_seed)
    size = calc_ht_size(ht, unit = "inch", show_annotation_legend = FALSE)

    pdf_h_heatmap_filename = glue("{heatmap_prefix}.pdf")
    pdf_h_heatmap_dir_filename = file.path(fig_dir, pdf_h_heatmap_filename)
    pdf(pdf_h_heatmap_dir_filename, width = 1.01 * size[1], height = 1.01 * size[2])
    set.seed(random_seed)
    draw(ht, background = "transparent", show_annotation_legend = FALSE)
    dev.off()
    print(glue("cluster_idx: {cluster_idx}; {nrow(wgc_log2_cbind)}"))
}