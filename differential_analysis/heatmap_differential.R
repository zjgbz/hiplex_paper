rm(list=ls())
library(arrow)
library(tibble)
library(glue)
library(svglite)
library(ComplexHeatmap)
library(circlize)
library(tidyr)
library(dplyr)

ht_opt$COLUMN_ANNO_PADDING <- unit(6, "mm")

calc_ht_size = function(ht, unit = "inch", show_annotation_legend=FALSE, column_title=NULL, column_title_fontsize=14) {
	pdf(NULL)
    if (show_annotation_legend) {
        ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize), merge_legend = TRUE, annotation_legend_side = "top")
    } else {
        ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize), show_annotation_legend = FALSE)
    }
	ht = draw(ht, background = "transparent", column_title=column_title, column_title_gp = gpar(fontsize=column_title_fontsize))
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly = TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly = TRUE)
	dev.off()

	c(w, h)
}

region_target_heatmap <- function(matrix, name, col_order, col_split, heatmap_prefix, fig_dir, font_size, col_fun=NULL, show_heatmap_legend=FALSE, show_colnames=FALSE, random_seed=42) {
    set.seed(random_seed)

    if (is.null(col_fun)) {
        ht = Heatmap(matrix, name=name,
        	show_row_names=FALSE,
        	show_column_names=show_colnames,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_order = col_order,
            column_split = col_split,
        	width=ncol(matrix) * unit(5, "mm"),
        	height=77 * 1.3 * unit(5, "mm"),
        	# height=nrow(matrix)*unit(5, "mm"),
            column_title_gp=gpar(fontsize=0),
            column_gap=unit(8, "mm"),
        	show_row_dend=FALSE,
            show_heatmap_legend=show_heatmap_legend,
            heatmap_legend_param = list(
            		title = "log2",
                    grid_width = 5*unit(5, "mm"),
                    legend_height = 77 * 1.3*unit(5, "mm"),
            	    # title_position = "topcenter",
            	    title_gp = gpar(fontsize = 40),
            	    labels_gp = gpar(fontsize = 40),
                    legend_direction = "vertical")
        	)
    } else {
        ht = Heatmap(matrix, name=name,
            col = col_fun,
        	show_row_names=FALSE,
        	show_column_names=show_colnames,
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            column_order = col_order,
            column_split = col_split,
        	width=ncol(matrix) * unit(5, "mm"),
        	height=77 * 1.3 * unit(5, "mm"),
        	# height=nrow(matrix)*unit(5, "mm"),
            column_title_gp=gpar(fontsize=0),
            column_gap=unit(8, "mm"),
        	show_row_dend=FALSE,
            show_heatmap_legend=show_heatmap_legend,
            heatmap_legend_param = list(
            		title = "log2",
                    grid_width = 5*unit(5, "mm"),
                    legend_height = 77 * 1.3*unit(5, "mm"),
            	    # title_position = "topcenter",
            	    title_gp = gpar(fontsize = 40),
            	    labels_gp = gpar(fontsize = 40),
                    legend_direction = "vertical")
        	)
    }
    set.seed(random_seed)
    size = calc_ht_size(ht, unit = "inch")

    pdf_h_heatmap_filename = glue("{heatmap_prefix}.pdf")
    pdf_h_heatmap_dir_filename = file.path(fig_dir, pdf_h_heatmap_filename)
    pdf(pdf_h_heatmap_dir_filename, width = size[1], height = size[2])
    set.seed(random_seed)
    ht_draw = draw(ht, background="transparent")
    dev.off()
}

method_sig_list = list("limma"=list("0.25"=c(1:8, 10:15)))
# method_sig_list = list("limma"=list("0.25"=c(15)))
method_filter_list = list("limma"=0.25)
wgc_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_wgc"
fig_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_fig"
l2fc_thres = 0.5
show_heatmap_legend_switch = "on"
show_colnames_switch = "on"

method_list = names(method_sig_list)
for (method in method_list) {
    mean_per_thres = method_filter_list[[method]]
    fdr_thres_list = names(method_sig_list[[method]])
    for (fdr_thres in fdr_thres_list) {
        cluster_idx_list = method_sig_list[[method]][[fdr_thres]]
        for (cluster_idx in cluster_idx_list) {
            # Plot heatmaps for log2 data
            V1_log2_wgc_filename = glue("V1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
            V2_log2_wgc_filename = glue("V2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
            T1_log2_wgc_filename = glue("T1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
            T2_log2_wgc_filename = glue("T2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
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

            colnames(V1_log2_wgc) = paste0("V1:", colnames(V1_log2_wgc))
            colnames(V2_log2_wgc) = paste0("V2:", colnames(V2_log2_wgc))
            colnames(T1_log2_wgc) = paste0("T1:", colnames(T1_log2_wgc))
            colnames(T2_log2_wgc) = paste0("T2:", colnames(T2_log2_wgc))

            wgc_log2_cbind = cbind(V1_log2_wgc, V2_log2_wgc, T1_log2_wgc, T2_log2_wgc)

            col_cluster_df = as.data.frame(colnames(wgc_log2_cbind))
            colnames(col_cluster_df) = "feature"
            col_cluster_df <- col_cluster_df %>% separate(col = feature, into = c("label", "target_pair"), sep = ":", remove = FALSE)
            col_cluster_df <- col_cluster_df %>% arrange(label, target_pair)

            col_order = col_cluster_df$feature
            col_split = col_cluster_df$label
            wgc_log2_cbind = wgc_log2_cbind[, col_order]

            # col_fun = colorRamp2(c(min(wgc_log2_cbind), 0.9, 1.8), c("blue", "white", "red"))
            if (show_heatmap_legend_switch == "on") {
                show_heatmap_legend = TRUE
            } else if (show_heatmap_legend_switch == "off") {
                show_heatmap_legend = FALSE
            }
            if (show_colnames_switch == "on") {
                show_colnames = TRUE
            } else if (show_colnames_switch == "off") {
                show_colnames = FALSE
            }
            heatmap_prefix = glue("publish_heatmap_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{cluster_idx}_{method}_FDR-{fdr_thres}_logFC-{l2fc_thres}_legend-{show_heatmap_legend_switch}_colnames-{show_colnames_switch}")
            # region_target_heatmap(as.matrix(wgc_log2_cbind), glue("log2"), col_order, col_split, heatmap_prefix, fig_dir, font_size = 50, col_fun=col_fun, show_heatmap_legend=show_heatmap_legend, show_colnames=show_colnames)
            region_target_heatmap(as.matrix(wgc_log2_cbind), glue("log2"), col_order, col_split, heatmap_prefix, fig_dir, font_size = 50, show_heatmap_legend=show_heatmap_legend, show_colnames=show_colnames)

            print(glue("method: {method}, fdr_thres: {fdr_thres}, cluster_idx: {cluster_idx} completed!"))
        }
    }
}