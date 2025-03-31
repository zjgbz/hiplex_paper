rm(list=ls())
library(ComplexHeatmap)
library(circlize)
library(glue)
library(dplyr)
library(tidyr)
library(latex2exp)
source("/projects/foundation_model_for_single_cell_multiomics_data/hiplex/script/utils/map_target_pair_names.R")

calc_ht_size <- function(ht, heatmap_legend_side=NULL, title=NULL, fontsize=NULL, unit = "inch") {
	pdf(NULL)
	ht = ComplexHeatmap::draw(ht, heatmap_legend_side=heatmap_legend_side, column_title=title, column_title_gp=gpar(fontsize=fontsize), background="transparent")
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly=TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly=TRUE)
	dev.off()

	c(w, h)
}

draw_heatmap <- function(df, cell_reads_thres, perm_num, cooccur_scen, cluster_i, heatmap_name, enrich_pseudocount, col_fun, fig_dir) {
    set.seed(123)
    df = as.matrix(df)
    df_ht = Heatmap(df,
        name=glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_{heatmap_name}_pseudocount-{enrich_pseudocount}"),
        col=col_fun,
    	cluster_columns=FALSE, cluster_rows=FALSE,
        show_row_names=TRUE, show_column_names=TRUE,
        row_names_side = "left",
        rect_gp = gpar(col = "white", lwd = 4, type = "none"),
        na_col = "gray74",
    	width=ncol(df) * unit(8, "mm"),
    	height=nrow(df) * unit(8, "mm"),
        row_names_gp = gpar(fontsize = 14),
        column_names_gp = gpar(fontsize = 14),
        # column_names_rot = 45,
        cell_fun = function(j, i, x, y, w, h, fill) {
            if(i >= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "white", lwd=4))
        }},
    	show_heatmap_legend=TRUE,
	    heatmap_legend_param = list(
		        title = heatmap_name,
    	        legend_height = 10*unit(5, "mm"),
		        grid_width = 1.3*unit(5, "mm"),
                # legend_width = 15*unit(5, "mm"),
                # grid_height = 1*unit(5, "mm"),
    	        title_position = "topleft",
		        title_gp = gpar(fontsize = 12),
		        labels_gp = gpar(fontsize = 10),
                legend_direction = "vertical"),
        column_labels=TeX(colnames(df)),
        row_labels=TeX(rownames(df))
    	)
    set.seed(123)
    df_ht_size = calc_ht_size(df_ht, heatmap_legend_side="right", title=heatmap_name, fontsize=17)
    
    df_ht_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_{heatmap_name}_pseudocount-{enrich_pseudocount}.pdf")
    df_ht_dir_filename = file.path(fig_dir, df_ht_filename)
    pdf(df_ht_dir_filename, width=df_ht_size[1], height=df_ht_size[2])
    set.seed(123)
    df_ht_draw = ComplexHeatmap::draw(df_ht, heatmap_legend_side="right", column_title=heatmap_name, column_title_gp=gpar(fontsize=17), background="transparent")
    dev.off()
}

obs_thres = 5
FDR_thres = 0.05
perm_num = 10000
tissue_name = "K562"
frag_type = "mixed"
cooccur_scen = "region"
# cooccur_scen = "cell"
# cluster_list = c("all", LETTERS[1:15])
cluster_list = c("all")
cell_reads_thres_list = c(0, 1)
enrich_pseudocount_list = c(0)
sc_root_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
wgc_root_dir = file.path(sc_root_dir, "wgc")
cooccur_root_dir = file.path(wgc_root_dir, tissue_name, frag_type, "cooccur")
test_result_root_dir = file.path(cooccur_root_dir, "test_result")
test_result_scen_root_dir = file.path(test_result_root_dir, cooccur_scen)

coexp_num_thres = 1
observation_root_dir = file.path(cooccur_root_dir, "observation")
target_pair_avail_root_dir = file.path(observation_root_dir, "target_pair_avail")
target_pair_avail_dir = file.path(target_pair_avail_root_dir, cooccur_scen)
fig_root_dir = file.path(cooccur_root_dir, "figure")
fig_dir = file.path(fig_root_dir, cooccur_scen)

for (cell_reads_thres in cell_reads_thres_list) {
    for (enrich_pseudocount in enrich_pseudocount_list) {
        target_pair_avail_final_list = c()
        for (cluster_i in cluster_list) {
            target_pair_avail_filename = glue("target_pair_avail_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{cluster_i}_coexp_num_thres-{coexp_num_thres}.tsv")
            target_pair_avail_dir_filename = file.path(target_pair_avail_dir, target_pair_avail_filename)
            target_pair_avail_df = read.table(target_pair_avail_dir_filename, header=TRUE, sep="\t", row.names=NULL)
            target_pair_avail_list = c(target_pair_avail_df$target_pair_1, target_pair_avail_df$target_pair_2)
            target_pair_avail_final_list = c(target_pair_avail_final_list, target_pair_avail_list)
        }

        target_pair_avail_final_list = unique(target_pair_avail_final_list)
        target_pair_avail_final_sort_list = sort(target_pair_avail_final_list, method="radix")

        for (cluster_i in cluster_list) {
            cooccur_log10p = data.frame(matrix(NaN, nrow=length(target_pair_avail_final_sort_list), ncol=length(target_pair_avail_final_sort_list)))
            rownames(cooccur_log10p) = target_pair_avail_final_sort_list
            colnames(cooccur_log10p) = target_pair_avail_final_sort_list
            cooccur_beta = cooccur_log10p
            cooccur_ob = cooccur_log10p
            cooccur_sim = cooccur_log10p

            target_pair_avail_filename = glue("target_pair_avail_cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{cluster_i}_coexp_num_thres-{coexp_num_thres}.tsv")
            target_pair_avail_dir_filename = file.path(target_pair_avail_dir, target_pair_avail_filename)
            target_pair_avail_df = read.table(target_pair_avail_dir_filename, header=TRUE, sep="\t", row.names=NULL)

            case_num = nrow(target_pair_avail_df)
            test_result_df_list = list()
            for (case_idx in 1:case_num) {
                target_pair_1 = target_pair_avail_df[case_idx, "target_pair_1"]
                target_pair_2 = target_pair_avail_df[case_idx, "target_pair_2"]

                target_pair_12 = c(target_pair_1, target_pair_2)
                target_pair_12_sort = sort(target_pair_12, method="radix")
                target_pair_sort_1 = target_pair_12_sort[1]
                target_pair_sort_2 = target_pair_12_sort[2]

                test_result_dir = file.path(test_result_scen_root_dir, cluster_i)
                tmp_test_result_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{target_pair_1}_{target_pair_2}_{cluster_i}_{cooccur_scen}.tsv")
                tmp_test_result_dir_filename = file.path(test_result_dir, tmp_test_result_filename)
                tmp_test_result_df = read.table(tmp_test_result_dir_filename, header=TRUE, sep="\t", row.names=NULL, check.names=FALSE)

                log10p_value = tmp_test_result_df[1, "-log10(pval)"]
                if (is.infinite(log10p_value)) {
                    log10p_value = -log10(1/(perm_num + 1))
                }
                tmp_test_result_df[1, "-log10(pval)"] = log10p_value
                test_result_df_list[[case_idx]] = tmp_test_result_df
            }
            test_result_df_full = do.call(rbind, test_result_df_list)
            test_result_df = test_result_df_full[test_result_df_full[, "1|1_region"] > obs_thres, ]
            test_result_df[, "p"] = 10^(-test_result_df[, "-log10(pval)"])
            test_result_df[, "FDR"] = p.adjust(test_result_df[, "p"], method="BH")
            test_result_df[, "log10_FDR"] = -log10(test_result_df[, "FDR"])
            test_result_df[, "enrich"] = (test_result_df[, glue("1|1_{cooccur_scen}")] + enrich_pseudocount) / (test_result_df[, glue("1|1_{cooccur_scen}_sim")] + enrich_pseudocount)
            test_result_df[, "log2_enrich"] = log2(test_result_df[, "enrich"])
            test_result_df_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_test_result_pseudocount-{enrich_pseudocount}_FDR.tsv")
            test_result_df_dir_filename = file.path(test_result_scen_root_dir, test_result_df_filename)
            write.table(test_result_df, test_result_df_dir_filename, sep="\t", row.names=FALSE, quote=FALSE)

            test_result_df_sig = test_result_df[test_result_df$FDR < FDR_thres, ]
            test_result_df_sig_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_test_result_pseudocount-{enrich_pseudocount}_FDR_sig_enrich.tsv")
            test_result_df_sig_dir_filename = file.path(test_result_scen_root_dir, test_result_df_sig_filename)
            write.table(test_result_df_sig, test_result_df_sig_dir_filename, sep="\t", row.names=FALSE, quote=FALSE)

            target_pair_name_sig_list = unique(c(test_result_df_sig$target_pair_1, test_result_df_sig$target_pair_2))
            log10_FDR_sym = data.frame(matrix(NaN, nrow=length(target_pair_name_sig_list), ncol=length(target_pair_name_sig_list)))
            rownames(log10_FDR_sym) = target_pair_name_sig_list
            colnames(log10_FDR_sym) = target_pair_name_sig_list
            log2_enrich_sym = log10_FDR_sym
            # enrich_sym = log10_FDR_sym
            sig_case_num = nrow(test_result_df_sig)
            for (sig_case_idx in c(1:sig_case_num)) {
                target_pair_1 = test_result_df_sig[sig_case_idx, "target_pair_1"]
                target_pair_2 = test_result_df_sig[sig_case_idx, "target_pair_2"]
                log10_FDR_sym[target_pair_1, target_pair_2] = test_result_df_sig[sig_case_idx, "log10_FDR"]
                log10_FDR_sym[target_pair_2, target_pair_1] = test_result_df_sig[sig_case_idx, "log10_FDR"]
                log2_enrich_sym[target_pair_1, target_pair_2] = test_result_df_sig[sig_case_idx, "log2_enrich"]
                log2_enrich_sym[target_pair_2, target_pair_1] = test_result_df_sig[sig_case_idx, "log2_enrich"]
                # enrich_sym[target_pair_1, target_pair_2] = test_result_df_sig[sig_case_idx, "enrich"]
                # enrich_sym[target_pair_2, target_pair_1] = test_result_df_sig[sig_case_idx, "enrich"]
            }

            log10_FDR_sym_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_pseudocount-{enrich_pseudocount}_log10_FDR_sym.tsv")
            log10_FDR_sym_dir_filename = file.path(test_result_scen_root_dir, log10_FDR_sym_filename)
            write.table(log10_FDR_sym, log10_FDR_sym_dir_filename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

            log2_enrich_sym_filename = glue("cell_reads_thres-{cell_reads_thres}_perm-{perm_num}_{cooccur_scen}_{cluster_i}_pseudocount-{enrich_pseudocount}_log2_enrich_sym.tsv")
            log2_enrich_sym_dir_filename = file.path(test_result_scen_root_dir, log2_enrich_sym_filename)
            write.table(log2_enrich_sym, log2_enrich_sym_dir_filename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

            # enrich_sym_filename = glue("cell_reads_thres-{cell_reads_thres}_{cooccur_scen}_{cluster_i}_enrich_sym.tsv")
            # enrich_sym_dir_filename = file.path(test_result_scen_root_dir, enrich_sym_filename)
            # write.table(enrich_sym, enrich_sym_dir_filename, sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)

            colnames(log10_FDR_sym) <- map_target_names(colnames(log10_FDR_sym), target_pair_mapping_df)
            colnames(log2_enrich_sym) <- map_target_names(colnames(log2_enrich_sym), target_pair_mapping_df)
            # colnames(enrich_sym) <- map_target_names(colnames(enrich_sym), target_pair_mapping_df)
            rownames(log10_FDR_sym) <- map_target_names(rownames(log10_FDR_sym), target_pair_mapping_df)
            rownames(log2_enrich_sym) <- map_target_names(rownames(log2_enrich_sym), target_pair_mapping_df)
            # rownames(enrich_sym) <- map_target_names(rownames(enrich_sym), target_pair_mapping_df)

            col_fun_log10_FDR = colorRamp2(range(log10_FDR_sym, na.rm=TRUE), hcl_palette = "Purp", reverse = TRUE)
            col_fun_log2_enrich = colorRamp2(range(log2_enrich_sym, na.rm=TRUE), hcl_palette = "PinkYl", reverse = TRUE)
            draw_heatmap(log10_FDR_sym, cell_reads_thres, perm_num, cooccur_scen, cluster_i, "-log10(FDR)", enrich_pseudocount, col_fun_log10_FDR, fig_dir)
            draw_heatmap(log2_enrich_sym, cell_reads_thres, perm_num, cooccur_scen, cluster_i, "log2(Enrichment)", enrich_pseudocount, col_fun_log2_enrich, fig_dir)
            # draw_heatmap(enrich_sym, cell_reads_thres, cooccur_scen, cluster_i, "Enrichment", fig_dir)

            print(glue("{cooccur_scen} {cluster_i} {cell_reads_thres} {enrich_pseudocount}"))
        }
    }
}