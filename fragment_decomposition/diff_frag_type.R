# Input: Fragment count and percentage data from ../data/frag_decomposition/
#        - frag_split_fastq-demux_num_{V,T}_peak-valley-notrimer.tsv (fragment counts)
#        - frag_split_fastq-demux_per_{V,T}_peak-valley-notrimer.tsv (fragment percentages)
# Output: Differential analysis results and volcano plots saved to ../results/Figure2/
#        - diff_{subnucleo,monomer,dimer}.tsv (Fisher's exact test results)
#        - volcano_{per_l2fc,log_OR}_{zoom,full}_{subnucleo,monomer,dimer}.pdf (volcano plots)

rm(list=ls())

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

packages <- c("ggplot2", "ggrepel", "glue")
lapply(packages, install_if_missing)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(glue)
})


args <- commandArgs(trailingOnly = TRUE)
frag_len_dir <- args[1]
out_dir <- args[2]


dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

before_num_filename = "frag_split_fastq-demux_num_V_peak-valley-notrimer.tsv"
after_num_filename = "frag_split_fastq-demux_num_T_peak-valley-notrimer.tsv"
before_num_dir_filename = file.path(frag_len_dir, before_num_filename)
after_num_dir_filename = file.path(frag_len_dir, after_num_filename)
before_num = read.table(before_num_dir_filename, sep="\t", check.names=FALSE, header=TRUE, row.names=1)
after_num = read.table(after_num_dir_filename, sep="\t", check.names=FALSE, header=TRUE, row.names=1)

before_per_filename = "frag_split_fastq-demux_per_V_peak-valley-notrimer.tsv"
after_per_filename = "frag_split_fastq-demux_per_T_peak-valley-notrimer.tsv"
before_per_dir_filename = file.path(frag_len_dir, before_per_filename)
after_per_dir_filename = file.path(frag_len_dir, after_per_filename)
before_per = read.table(before_per_dir_filename, sep="\t", check.names=FALSE, header=TRUE, row.names=1)
after_per = read.table(after_per_dir_filename, sep="\t", check.names=FALSE, header=TRUE, row.names=1)

frag_type_list = c("subnucleo", "monomer", "dimer")
target_pair_list = rownames(before_num)

l2fc_thres = 0.1
log_OR_thres = 0.25

for (frag_type_i in frag_type_list) {
    frag_diff_i = data.frame()
    for (target_pair_i in target_pair_list) {
        before_num_in = before_num[target_pair_i, frag_type_i]
        after_num_in = after_num[target_pair_i, frag_type_i]
        before_num_out = sum(before_num[target_pair_i, ]) - before_num_in
        after_num_out = sum(after_num[target_pair_i, ]) - after_num_in
        tmp_contengency = matrix(c(after_num_in, after_num_out, before_num_in, before_num_out), nrow=2, byrow=TRUE)
        tmp_test = fisher.test(tmp_contengency)
        or = tmp_test$estimate
        log_or = log2(or)
        tmp_num_log2fc = log2((after_num_in + 1) / (before_num_in + 1))

        before_per_in = before_per[target_pair_i, frag_type_i]
        after_per_in = after_per[target_pair_i, frag_type_i]
        tmp_per_log2fc = log2((after_per_in) / (before_per_in))
        tmp_diff_i = data.frame("target_pair"=target_pair_i, "pvalue"=tmp_test$p.value, "odds_ratio"=tmp_test$estimate, "log_OR"=log_or, "num_l2fc"=tmp_num_log2fc, "per_l2fc"=tmp_per_log2fc)
        frag_diff_i = rbind(frag_diff_i, tmp_diff_i)
    }

    frag_diff_i_num_l2fc_max = max(abs(frag_diff_i$num_l2fc))
    frag_diff_i_per_l2fc_max = max(abs(frag_diff_i$per_l2fc))
    frag_diff_i_num_l2fc_thres = frag_diff_i_num_l2fc_max * 1.01
    frag_diff_i_per_l2fc_thres = frag_diff_i_per_l2fc_max * 1.01

    frag_diff_i_log_OR_max = max(abs(frag_diff_i$log_OR))
    frag_diff_i_log_OR_thres = frag_diff_i_log_OR_max * 1.01
    
    frag_diff_i[, "pvalue_adj"] = p.adjust(frag_diff_i$pvalue, method = "BH")
    frag_diff_i[, "log10_pvalue_adj"] = -log10(frag_diff_i$pvalue_adj)

    # 保存到 out_dir
    frag_diff_i_filename = glue("diff_{frag_type_i}.tsv")
    frag_diff_i_dir_filename = file.path(out_dir, frag_diff_i_filename)
    write.table(frag_diff_i, frag_diff_i_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)

    # l2fc per; zoom in y axis
    frag_diff_i$pvalue[frag_diff_i$pvalue == 0] = 1E-300
    frag_diff_i[, "log10_pvalue"] = -log10(frag_diff_i$pvalue)
    frag_diff_i[, "diffexpressed_per"] = "NO"
    frag_diff_i[frag_diff_i$per_l2fc > l2fc_thres & frag_diff_i$pvalue < 0.05, "diffexpressed_per"] = "UP"
    frag_diff_i[frag_diff_i$per_l2fc < -l2fc_thres & frag_diff_i$pvalue < 0.05, "diffexpressed_per"] = "DOWN"
    diff_per_zoom_plot = ggplot(data = frag_diff_i, aes(x = per_l2fc, y = log10_pvalue, col = diffexpressed_per)) +
                            geom_vline(xintercept = c(-l2fc_thres, l2fc_thres), col = "gray", linetype = 'dashed') +
                            geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
                            geom_point(size = 1) + ylim(-0.01, 10) + xlim(-frag_diff_i_per_l2fc_thres, frag_diff_i_per_l2fc_thres)
    if (frag_type_i != "subnucleo") {
        diff_per_zoom_plot = diff_per_zoom_plot + scale_color_manual(values = c("#00AFBB", "grey", "#BB0E08"),
                                                                        labels = c("V", "Not significant", "T"))
    } else if (frag_type_i == "subnucleo") {
        diff_per_zoom_plot = diff_per_zoom_plot + scale_color_manual(values = c("grey", "#BB0E08"),
                                                                        labels = c("Not significant", "T"))
    }
    diff_per_zoom_plot_filename = paste0("volcano_per_l2fc_zoom_", frag_type_i, ".pdf")
    diff_per_zoom_plot_dir_filename = file.path(out_dir, diff_per_zoom_plot_filename)
    ggsave(diff_per_zoom_plot_dir_filename, diff_per_zoom_plot, width = 5, height = 5)

    # l2fc per; full y axis
    diff_per_full_plot = ggplot(data = frag_diff_i, aes(x = per_l2fc, y = log10_pvalue, col = diffexpressed_per)) +
                            geom_vline(xintercept = c(-l2fc_thres, l2fc_thres), col = "gray", linetype = 'dashed') +
                            geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
                            geom_point(size = 1) + xlim(-frag_diff_i_per_l2fc_thres, frag_diff_i_per_l2fc_thres)
    if (frag_type_i != "subnucleo") {
        diff_per_full_plot = diff_per_full_plot + scale_color_manual(values = c("#00AFBB", "grey", "#BB0E08"),
                                                                        labels = c("V", "Not significant", "T"))
    } else if (frag_type_i == "subnucleo") {
        diff_per_full_plot = diff_per_full_plot + scale_color_manual(values = c("grey", "#BB0E08"),
                                                                        labels = c("Not significant", "T"))
    }
    diff_per_full_plot_filename = paste0("volcano_per_l2fc_full_", frag_type_i, ".pdf")
    diff_per_full_plot_dir_filename = file.path(out_dir, diff_per_full_plot_filename)
    ggsave(diff_per_full_plot_dir_filename, diff_per_full_plot, width = 5, height = 5)

    # log OR; zoom in y axis
    frag_diff_i[, "diffexpressed_log_OR"] = "NO"
    frag_diff_i[frag_diff_i$log_OR > log_OR_thres & frag_diff_i$pvalue_adj < 0.05, "diffexpressed_log_OR"] = "UP"
    frag_diff_i[frag_diff_i$log_OR < -log_OR_thres & frag_diff_i$pvalue_adj < 0.05, "diffexpressed_log_OR"] = "DOWN"
    frag_diff_i_select = frag_diff_i[frag_diff_i$diffexpressed_log_OR == "UP" | frag_diff_i$diffexpressed_log_OR == "DOWN", ]

    diff_log_OR_zoom_plot = ggplot(data = frag_diff_i, aes(x = log_OR, y = log10_pvalue_adj, col = diffexpressed_log_OR)) +
                            geom_vline(xintercept = c(-log_OR_thres, log_OR_thres), col = "gray", linetype = 'dashed') +
                            geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
                            geom_text_repel(data=frag_diff_i_select, aes(x=log_OR, y=log10_pvalue_adj,label=target_pair), size = 1.5) +
                            geom_point(size = 1) + ylim(-0.01, 10) + xlim(-frag_diff_i_log_OR_thres, frag_diff_i_log_OR_thres)
    if (frag_type_i != "subnucleo") {
        diff_log_OR_zoom_plot = diff_log_OR_zoom_plot + scale_color_manual(values = c("#00AFBB", "grey", "#BB0E08"),
                                                                        labels = c("V", "Not significant", "T"))
    } else if (frag_type_i == "subnucleo") {
        diff_log_OR_zoom_plot = diff_log_OR_zoom_plot + scale_color_manual(values = c("grey", "#BB0E08"),
                                                                        labels = c("Not significant", "T"))
    }
    diff_log_OR_zoom_plot_filename = paste0("volcano_log_OR_zoom_", frag_type_i, ".pdf")
    diff_log_OR_zoom_plot_dir_filename = file.path(out_dir, diff_log_OR_zoom_plot_filename)
    ggsave(diff_log_OR_zoom_plot_dir_filename, diff_log_OR_zoom_plot, width = 5, height = 5)

    # log OR; full y axis
    diff_log_OR_full_plot = ggplot(data = frag_diff_i, aes(x = log_OR, y = log10_pvalue_adj, col = diffexpressed_log_OR)) +
                            geom_vline(xintercept = c(-log_OR_thres, log_OR_thres), col = "gray", linetype = 'dashed') +
                            geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
                            geom_text_repel(data=frag_diff_i_select, aes(x=log_OR, y=log10_pvalue_adj,label=target_pair), size = 1.5) +
                            geom_point(size = 1) + xlim(-frag_diff_i_log_OR_thres, frag_diff_i_log_OR_thres)
    if (frag_type_i != "subnucleo") {
        diff_log_OR_full_plot = diff_log_OR_full_plot + scale_color_manual(values = c("#00AFBB", "grey", "#BB0E08"),
                                                                        labels = c("V", "Not significant", "T"))
    } else if (frag_type_i == "subnucleo") {
        diff_log_OR_full_plot = diff_log_OR_full_plot + scale_color_manual(values = c("grey", "#BB0E08"),
                                                                        labels = c("Not significant", "T"))
    }
    diff_log_OR_full_plot_filename = paste0("volcano_log_OR_full_", frag_type_i, ".pdf")
    diff_log_OR_full_plot_dir_filename = file.path(out_dir, diff_log_OR_full_plot_filename)
    ggsave(diff_log_OR_full_plot_dir_filename, diff_log_OR_full_plot, width = 5, height = 5)

    print(c(frag_type_i, max(frag_diff_i$log_OR), min(frag_diff_i$log_OR), max(frag_diff_i$per_l2fc), min(frag_diff_i$per_l2fc)))
}