rm(list=ls())
library(arrow)
library(tibble)
library(glue)
library(latex2exp)
library(edgeR)
library(matrixStats)

frag_len = 800
frag_type = "mixed"
target_qc_type = "all-qc"
process_type = "orig"
pseudocount = 0.5
normalization_factor = 1E6
pseudocount_for_log = 1
lowess_span = 0.5
l2fc_thres = 0.5

coldata = data.frame("condition"=c("untreated", "untreated", "treated", "treated"), row.names = c("V1", "V2", "T1", "T2"))
group = factor(coldata$condition)
group1 = c("V1", "V2")
group2 = c("T1", "T2")

mm <- model.matrix(~0 + group)

wgc_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800"
sig_result_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_result"
sig_wgc_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_wgc"
sig_plot_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_fig"
col_cluster_filename = "bicluster_V_mixed_all-qc_kmeans_euclidean_row_num-15_column_num-16_heatmap_column_clusters_raw_colnames_manually_reorganized.tsv"
col_cluster_dir_filename = file.path(wgc_dir, col_cluster_filename)
col_cluster = read.table(col_cluster_dir_filename, header=TRUE, sep="\t", row.names=NULL)
col_label_list = unique(col_cluster$label)

V1_wgc_orig_filename = glue("V1_{frag_type}_{frag_len}_colQC-{target_qc_type}_orig.feather")
V2_wgc_orig_filename = glue("V2_{frag_type}_{frag_len}_colQC-{target_qc_type}_orig.feather")
T1_wgc_orig_filename = glue("T1_{frag_type}_{frag_len}_colQC-{target_qc_type}_orig.feather")
T2_wgc_orig_filename = glue("T2_{frag_type}_{frag_len}_colQC-{target_qc_type}_orig.feather")
V1_wgc_orig_dir_filename = file.path(wgc_dir, V1_wgc_orig_filename)
V2_wgc_orig_dir_filename = file.path(wgc_dir, V2_wgc_orig_filename)
T1_wgc_orig_dir_filename = file.path(wgc_dir, T1_wgc_orig_filename)
T2_wgc_orig_dir_filename = file.path(wgc_dir, T2_wgc_orig_filename)
V1_wgc_orig = read_feather(V1_wgc_orig_dir_filename)
V2_wgc_orig = read_feather(V2_wgc_orig_dir_filename)
T1_wgc_orig = read_feather(T1_wgc_orig_dir_filename)
T2_wgc_orig = read_feather(T2_wgc_orig_dir_filename)
V1_wgc_orig = column_to_rownames(V1_wgc_orig, var="pos")
V2_wgc_orig = column_to_rownames(V2_wgc_orig, var="pos")
T1_wgc_orig = column_to_rownames(T1_wgc_orig, var="pos")
T2_wgc_orig = column_to_rownames(T2_wgc_orig, var="pos")

V1_wgc_transformed_filename = glue("V1_{frag_type}_{frag_len}_colQC-{target_qc_type}_libnorm_log2.feather")
V2_wgc_transformed_filename = glue("V2_{frag_type}_{frag_len}_colQC-{target_qc_type}_libnorm_log2.feather")
T1_wgc_transformed_filename = glue("T1_{frag_type}_{frag_len}_colQC-{target_qc_type}_libnorm_log2.feather")
T2_wgc_transformed_filename = glue("T2_{frag_type}_{frag_len}_colQC-{target_qc_type}_libnorm_log2.feather")
V1_wgc_transformed_dir_filename = file.path(wgc_dir, V1_wgc_transformed_filename)
V2_wgc_transformed_dir_filename = file.path(wgc_dir, V2_wgc_transformed_filename)
T1_wgc_transformed_dir_filename = file.path(wgc_dir, T1_wgc_transformed_filename)
T2_wgc_transformed_dir_filename = file.path(wgc_dir, T2_wgc_transformed_filename)
V1_wgc_transformed = read_feather(V1_wgc_transformed_dir_filename)
V2_wgc_transformed = read_feather(V2_wgc_transformed_dir_filename)
T1_wgc_transformed = read_feather(T1_wgc_transformed_dir_filename)
T2_wgc_transformed = read_feather(T2_wgc_transformed_dir_filename)
V1_wgc_transformed = column_to_rownames(V1_wgc_transformed, var="pos")
V2_wgc_transformed = column_to_rownames(V2_wgc_transformed, var="pos")
T1_wgc_transformed = column_to_rownames(T1_wgc_transformed, var="pos")
T2_wgc_transformed = column_to_rownames(T2_wgc_transformed, var="pos")

# mean_per_thres_dict = list(c(0.85, 0.91, 0.8, 0.81, 0.91, 0.7, 0.8, 0.8, 0.7, 0.92, 0.7, 0.8, 0.85, 0.8, 0.87, 0.8),
#                             c(0.85, 0.75, 0.7, 0.87, 0.9, 0.8, 0.7, 0.85, 0.9, 0.85, 0.68, 0.8, 0.85, 0.85, 0.85, 0.85))
# var_per_thres = 0.1

# mean_per_thres_list = mean_per_thres_dict[[min_total_count_cutoff]]
# mean_per_thres_list = c(0.14, 0.25)
mean_per_thres_list = c(0.25)
fdr_thres_list = c(0.25, 0.1)
for (col_label in col_label_list) {
    target_pair_select_list = col_cluster[col_cluster$label == col_label, "feature"]

    V1_wgc_orig_select_pre = V1_wgc_orig[, target_pair_select_list]
    V2_wgc_orig_select_pre = V2_wgc_orig[, target_pair_select_list]
    T1_wgc_orig_select_pre = T1_wgc_orig[, target_pair_select_list]
    T2_wgc_orig_select_pre = T2_wgc_orig[, target_pair_select_list]
    tmp_combine_orig = data.frame("V1"=rowSums(V1_wgc_orig_select_pre),
                                  "V2"=rowSums(V2_wgc_orig_select_pre),
                                  "T1"=rowSums(T1_wgc_orig_select_pre),
                                  "T2"=rowSums(T2_wgc_orig_select_pre)
                                 )
    
    group1_w_zero_boolean <- rowSums(tmp_combine_orig[, group1] == 0) == 0
    group2_w_zero_boolean <- rowSums(tmp_combine_orig[, group2] == 0) == 0
    tmp_combine = tmp_combine_orig[group1_w_zero_boolean | group2_w_zero_boolean, ]

    d0 <- DGEList(tmp_combine)
    d <- calcNormFactors(d0)
    # extract information
    norm_info = d$samples
    lib.size = norm_info$lib.size
    norm.factors = norm_info$norm.factors
    effect_libsize = norm_info$lib.size * norm_info$norm.factors
    # mean-variance plot
    voom_plot_filename = glue("voom_plot_lowess_span-{lowess_span}_post-filter-one_condition_nonzero-2_column_cluster-{col_label}.pdf")
    voom_plot_dir_filename = file.path(sig_plot_dir, voom_plot_filename)
    pdf(voom_plot_dir_filename)
    y <- voom(d, mm, span=lowess_span, plot = TRUE)
    dev.off()

    mean_log2_cpm <- rowMeans(y$E) + log2(exp(mean(log(lib.size + 1)))) - log2(normalization_factor)

    fit <- lmFit(y, mm)
    fitted_values <- fit$coefficients %*% t(mm)
    residuals <- y$E - fitted_values
    std_dev_residuals <- apply(residuals, 1, sd)
    sqrt_res_std <- sqrt(std_dev_residuals)

    # Step 2: Filter rows in tmp_combines based on the 0.9 percentile
    # tmp_combine_filtered <- tmp_combine[(mean_log2_cpm >= quantile(mean_log2_cpm, mean_per_thres)) & (sqrt_res_std > quantile(sqrt_res_std, var_per_thres)), ]
    for (mean_per_thres in mean_per_thres_list) {
        tmp_combine_filtered <- tmp_combine[mean_log2_cpm >= quantile(mean_log2_cpm, mean_per_thres), ]
        
        # Step 3: Re-process the filtered data
        d0_filtered <- DGEList(tmp_combine_filtered)
        d_filtered <- calcNormFactors(d0_filtered)
        norm_info_filtered = d_filtered$samples
        lib_size_filtered = norm_info_filtered$lib.size
        norm_factors_filtered = norm_info_filtered$norm.factors
        effect_libsize_filtered = lib_size_filtered * norm_factors_filtered
        
        # Perform voom transformation on the filtered data
        filtered_voom_plot_filename = glue("voom_plot_lowess_span-{lowess_span}_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}.pdf")
        filtered_voom_plot_dir_filename = file.path(sig_plot_dir, filtered_voom_plot_filename)
        pdf(filtered_voom_plot_dir_filename)
        y_filtered <- voom(d_filtered, mm, span = lowess_span, plot = TRUE)
        dev.off()
        E_value_filtered = y_filtered$E
        
        fit <- lmFit(y_filtered, mm)
        contr <- makeContrasts(grouptreated - groupuntreated, levels = colnames(coef(fit)))
        tmp <- contrasts.fit(fit, contr)
        tmp <- eBayes(tmp)
        print(glue("column cluster: {col_label}, df: {unique(tmp$df.total)}, original: {nrow(tmp_combine_orig)}, filtered: {nrow(tmp_combine)}, second filtered: {nrow(tmp_combine_filtered)}"))
        top.table <- topTable(tmp, sort.by = "P", n = Inf)
        
        # save the variance histogram
        E_value_filtered_var = rowVars(E_value_filtered)
        var_hist_filename = glue("hist_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_for_limma.pdf")
        var_hist_dir_filename = file.path(sig_plot_dir, var_hist_filename)
        pdf(var_hist_dir_filename)
        hist(E_value_filtered_var, breaks=100, main=glue("filter: {mean_per_thres}, col cluster: {col_label}\nregion #: {nrow(E_value_filtered)}, dof: {unique(tmp$df.total)}"), xlab="Variance", col="blue", border=FALSE)
        dev.off()
        
        significance_plot_filename = glue("limma_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_significance_plot.pdf")
        significance_plot_dir_filename = file.path(sig_plot_dir, significance_plot_filename)
        pdf(significance_plot_dir_filename)
        # Plot two histograms
        top.table$P.Value[top.table$P.Value == 0] = 1E-300
        top.table$adj.P.Val[top.table$adj.P.Val == 0] = 1E-300
        p_list = -log10(top.table$P.Value)
        fdr_list = -log10(top.table$adj.P.Val)
        
        top.table_for_save = rownames_to_column(top.table, var = "pos")
        top.table_for_save_filename = glue("result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_limma.feather")
        top.table_for_save_dir_filename = file.path(sig_result_dir, top.table_for_save_filename)
        write_feather(top.table_for_save, top.table_for_save_dir_filename)
        
        breaks <- seq(min(c(p_list, fdr_list)), max(c(p_list, fdr_list)), length.out = 100)
        
        # Create first histogram
        color1 = "#454597"
        color2 = "#6a936a"
        color1_trans <- adjustcolor(color1, alpha.f = 0.5)
        color2_trans <- adjustcolor(color2, alpha.f = 0.5)
        # title = paste(target_pair_select_list, collapse = ",")
        # wrapped_title <- paste(strwrap(title, width = 20), collapse = "\n")
        hist_obj = hist(p_list, breaks = breaks, col = color1_trans, border = NA, 
             xlab = "", main = "")
        
        # Overlay second histogram using `add = TRUE`
        hist(fdr_list, breaks = breaks, col = color2_trans, border = NA, add = TRUE)
        
        # Add vertical dashed red lines at specified x values
        vlines <- c(-log10(0.05), -log10(0.1), -log10(0.25))
        vline_labels <- c("0.05", "0.1", "0.25")  # Labels for vertical lines
        abline(v = vlines, col = "#b40e0e", lty = 2, lwd = 2)
        title(main = glue("limma, filter: {mean_per_thres}, col cluster: {col_label}"), cex.main = 0.7)
        y_max_hist = max(hist_obj$counts)
        text(x = vlines + 0.1, y = c(y_max_hist * 0.45, y_max_hist * 0.65, y_max_hist * 0.85), labels = vline_labels, col = "#5500ff", pos = 3, cex = 0.8)
        legend("topright", legend = c("-log10(p-value)", "-log10(FDR)"), 
                fill = c("#454597", "#6a936a"), border = c("#454597", "#6a936a"),
                bty = "n", cex = 0.8)
        dev.off()
        
        for (fdr_thres in fdr_thres_list) {
            top.table_sig = top.table[(top.table$adj.P.Val < fdr_thres) & (abs(top.table$logFC) > l2fc_thres), ]
            print(glue("column cluster: {col_label}, sig regions: {nrow(top.table_sig)}"))
            sig_region_num = nrow(top.table_sig)
        
            if (sig_region_num == 0) {
                print(glue("Will not save for cluster {col_label}!"))
            } else if (sig_region_num > 0) {
                top.table_sig_filename = glue("result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.tsv")
                top.table_sig_dir_filename = file.path(sig_result_dir, top.table_sig_filename)
                top.table_sig = top.table_sig[order(top.table_sig$logFC), ]
                sig_region_list = rownames(top.table_sig)
                write.table(top.table_sig, top.table_sig_dir_filename, sep="\t", quote=FALSE, row.names=TRUE, col.names=NA)
        
                V1_wgc_orig_select = V1_wgc_orig_select_pre[sig_region_list, target_pair_select_list]
                V2_wgc_orig_select = V2_wgc_orig_select_pre[sig_region_list, target_pair_select_list]
                T1_wgc_orig_select = T1_wgc_orig_select_pre[sig_region_list, target_pair_select_list]
                T2_wgc_orig_select = T2_wgc_orig_select_pre[sig_region_list, target_pair_select_list]
        
                V1_wgc_log2_select = log2(V1_wgc_orig_select + pseudocount_for_log)
                V2_wgc_log2_select = log2(V2_wgc_orig_select + pseudocount_for_log)
                T1_wgc_log2_select = log2(T1_wgc_orig_select + pseudocount_for_log)
                T2_wgc_log2_select = log2(T2_wgc_orig_select + pseudocount_for_log)
        
                V1_wgc_log2_select = rownames_to_column(V1_wgc_log2_select, var = "pos")
                V2_wgc_log2_select = rownames_to_column(V2_wgc_log2_select, var = "pos")
                T1_wgc_log2_select = rownames_to_column(T1_wgc_log2_select, var = "pos")
                T2_wgc_log2_select = rownames_to_column(T2_wgc_log2_select, var = "pos")
        
                V1_wgc_log2_select_filename = glue("V1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                V2_wgc_log2_select_filename = glue("V2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                T1_wgc_log2_select_filename = glue("T1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                T2_wgc_log2_select_filename = glue("T2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_log2_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
        
                V1_wgc_log2_select_dir_filename = file.path(sig_wgc_dir, V1_wgc_log2_select_filename)
                V2_wgc_log2_select_dir_filename = file.path(sig_wgc_dir, V2_wgc_log2_select_filename)
                T1_wgc_log2_select_dir_filename = file.path(sig_wgc_dir, T1_wgc_log2_select_filename)
                T2_wgc_log2_select_dir_filename = file.path(sig_wgc_dir, T2_wgc_log2_select_filename)
        
                write_feather(V1_wgc_log2_select, V1_wgc_log2_select_dir_filename)
                write_feather(V2_wgc_log2_select, V2_wgc_log2_select_dir_filename)
                write_feather(T1_wgc_log2_select, T1_wgc_log2_select_dir_filename)
                write_feather(T2_wgc_log2_select, T2_wgc_log2_select_dir_filename)
        
                # sum_wgc_orig = tmp_combine[sig_region_list, ]
                # sum_wgc_orig = rownames_to_column(sum_wgc_orig, var = "pos")
                # sum_wgc_orig_filename = glue("sum_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_orig_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                # sum_wgc_orig_dir_filename = file.path(sig_wgc_dir, sum_wgc_orig_filename)
                # write_feather(sum_wgc_orig, sum_wgc_orig_dir_filename)
        
                V1_wgc_limmanorm_select = log2((V1_wgc_orig_select + 0.5) * 1E6 / effect_libsize[1])
                V2_wgc_limmanorm_select = log2((V2_wgc_orig_select + 0.5) * 1E6 / effect_libsize[2])
                T1_wgc_limmanorm_select = log2((T1_wgc_orig_select + 0.5) * 1E6 / effect_libsize[3])
                T2_wgc_limmanorm_select = log2((T2_wgc_orig_select + 0.5) * 1E6 / effect_libsize[4])
        
                V1_wgc_limmanorm_select = rownames_to_column(V1_wgc_limmanorm_select, var = "pos")
                V2_wgc_limmanorm_select = rownames_to_column(V2_wgc_limmanorm_select, var = "pos")
                T1_wgc_limmanorm_select = rownames_to_column(T1_wgc_limmanorm_select, var = "pos")
                T2_wgc_limmanorm_select = rownames_to_column(T2_wgc_limmanorm_select, var = "pos")
        
                V1_wgc_limmanorm_select_filename = glue("V1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                V2_wgc_limmanorm_select_filename = glue("V2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                T1_wgc_limmanorm_select_filename = glue("T1_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                T2_wgc_limmanorm_select_filename = glue("T2_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm-manually_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
        
                V1_wgc_limmanorm_select_dir_filename = file.path(sig_wgc_dir, V1_wgc_limmanorm_select_filename)
                V2_wgc_limmanorm_select_dir_filename = file.path(sig_wgc_dir, V2_wgc_limmanorm_select_filename)
                T1_wgc_limmanorm_select_dir_filename = file.path(sig_wgc_dir, T1_wgc_limmanorm_select_filename)
                T2_wgc_limmanorm_select_dir_filename = file.path(sig_wgc_dir, T2_wgc_limmanorm_select_filename)
        
                write_feather(V1_wgc_limmanorm_select, V1_wgc_limmanorm_select_dir_filename)
                write_feather(V2_wgc_limmanorm_select, V2_wgc_limmanorm_select_dir_filename)
                write_feather(T1_wgc_limmanorm_select, T1_wgc_limmanorm_select_dir_filename)
                write_feather(T2_wgc_limmanorm_select, T2_wgc_limmanorm_select_dir_filename)
        
                # sum_wgc_norm = as.data.frame(E_value_filtered[sig_region_list, ])
                # sum_wgc_norm = rownames_to_column(sum_wgc_norm, var = "pos")
                # sum_wgc_norm_filename = glue("sum_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-{mean_per_thres}_limmanorm_column_cluster-{col_label}_limma_FDR-{fdr_thres}_logFC-{l2fc_thres}.feather")
                # sum_wgc_norm_dir_filename = file.path(sig_wgc_dir, sum_wgc_norm_filename)
                # write_feather(sum_wgc_norm, sum_wgc_norm_dir_filename)
            }
        }
    }
}