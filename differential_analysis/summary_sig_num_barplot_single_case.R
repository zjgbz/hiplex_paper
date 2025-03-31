rm(list=ls())
library(glue)
library(tidyr)
library(ggplot2)

# fdr_thres_list = c(0.25, 0.1)
# mean_filter_per_list = c(0.14, 0.25)
fdr_thres_list = c(0.25)
mean_filter_per_list = c(0.25)
l2fc_thres = 0.5

fig_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_fig"
result_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_result"

for (mean_filter_per in mean_filter_per_list) {
    for (fdr_thres in fdr_thres_list) {
        result_table_filename = glue("summary_sig_num_{mean_filter_per}_{fdr_thres}_l2fc-{l2fc_thres}.tsv")
        result_table_dir_filename = file.path(result_dir, result_table_filename)
        result_table = read.table(result_table_dir_filename, header=TRUE, sep="\t")
        colnames(result_table) = c("column_cluster", "sig_num")

        result_figure = ggplot(result_table, aes(y=sig_num, x=factor(column_cluster, levels=c(1:16)))) + 
                                geom_bar(position="dodge", stat="identity") +
                                geom_text(aes(label = sig_num), 
                                          position = position_dodge(width = 0.9), 
                                          vjust = -0.5, 
                                          size = 3) +
                                scale_x_discrete(expand = c(0, 0)) +
                                labs(
                                    title = glue("FDR < {fdr_thres}, abs(log2FC) > {l2fc_thres}"),
                                    x = "Column Clusters",
                                    y = "# Significant Regions",
                                    )
        result_figure_filename = glue("summary_sig_num_{mean_filter_per}_{fdr_thres}_l2fc-{l2fc_thres}.pdf")
        result_figure_dir_filename = file.path(fig_dir, result_figure_filename)
        ggsave(result_figure_dir_filename, result_figure, width=12, height=5)
    }
}