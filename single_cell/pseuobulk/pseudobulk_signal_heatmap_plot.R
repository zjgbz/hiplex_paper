rm(list=ls())
library(glue)
library(arrow)
library(tibble)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

calc_ht_size = function(ht, unit = "inch") {
	pdf(NULL)
	ht = draw(ht, background = "transparent")
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly = TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly = TRUE)
	dev.off()

	c(w, h)
}

# signal_visual_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell/wgc/K562/mixed"
signal_visual_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/40000"
frag_len = 40000
process = "orig"
target_pair_selected_list = c("H3K27me3-H3K27me3", "H3K27me3-H3K4me3", "H3K4me3-H3K4me3")
width = 240
height = 30
cell_number_kept = 100

final_region_chr = "chr1"
final_region_start = 108280001 #108276801
final_region_end = 114120000 #114093600

log10_coverage_list = c()
for (target_pair_selected in target_pair_selected_list) {
	signal_visual_df_filename = glue("signal_visual_{target_pair_selected}_{final_region_chr}_{final_region_start}_{final_region_end}_{frag_len}_{process}.feather")
	signal_visual_df_dir_filename = file.path(signal_visual_dir, signal_visual_df_filename)
	signal_visual_df = read_feather(signal_visual_df_dir_filename)
	signal_visual_df = column_to_rownames(signal_visual_df, var="cell_name")
	signal_visual_df[, "log10(Coverage)"] = rowSums(signal_visual_df)
	signal_visual_df = log10(signal_visual_df + 1)
	signal_visual_df = signal_visual_df[order(signal_visual_df[, "log10(Coverage)"], decreasing=TRUE),]
	cell_num = nrow(signal_visual_df)
	if (cell_num > cell_number_kept) {
		signal_visual_df = signal_visual_df[1:cell_number_kept, ]
	}
	signal_visual_sort = select(signal_visual_df, -c("log10(Coverage)"))
	log10_coverage_list = c(log10_coverage_list, max(signal_visual_sort))
}

col_fun = colorRamp2(c(0, max(log10_coverage_list)), c("white", "blue"))

for (target_pair_selected in target_pair_selected_list) {
	signal_visual_df_filename = glue("signal_visual_{target_pair_selected}_{final_region_chr}_{final_region_start}_{final_region_end}_{frag_len}_{process}.feather")
	signal_visual_df_dir_filename = file.path(signal_visual_dir, signal_visual_df_filename)
	signal_visual_df = read_feather(signal_visual_df_dir_filename)
	signal_visual_df = column_to_rownames(signal_visual_df, var="cell_name")
	signal_visual_df[, "log10(Coverage)"] = rowSums(signal_visual_df)
	signal_visual_df = log10(signal_visual_df + 1)
	signal_visual_df = signal_visual_df[order(signal_visual_df[, "log10(Coverage)"], decreasing=TRUE),]
	cell_num = nrow(signal_visual_df)
	if (cell_num > cell_number_kept) {
		signal_visual_df = signal_visual_df[1:cell_number_kept, ]
	}
	cell_num_plot = nrow(signal_visual_df)
	# remove column 'Coverage' from data frame signal_visual_df
	signal_visual_sort = select(signal_visual_df, -c("log10(Coverage)"))
	signal_visual_sort = as.matrix(signal_visual_sort)
	coverage = signal_visual_df[, "log10(Coverage)"]
	print(glue("{target_pair_selected}, {max(signal_visual_sort)}, cell number: {nrow(signal_visual_sort)}"))

	ht = Heatmap(signal_visual_sort, name="log10(Coverage)",
	    show_row_names=FALSE, show_column_names=FALSE,
		cluster_rows=FALSE, cluster_columns=FALSE,
	    right_annotation=rowAnnotation(Coverage=anno_barplot(coverage)),
	    border_gp=gpar(col="black", lwd=1), 
		col=col_fun,
	    width=unit(width, "mm"),
		height=unit(height, "mm"),
		show_row_dend=FALSE,
		show_column_dend=FALSE
		)
	size = calc_ht_size(ht)

	ht_filename = glue("signal_visual_{target_pair_selected}_{final_region_chr}_{final_region_start}_{final_region_end}_{frag_len}_{process}_{cell_number_kept}_height-{height}_width-{width}.pdf")
	ht_dir_filename = file.path(signal_visual_dir, ht_filename)
	pdf(ht_dir_filename, width=size[1], height=size[2])
	ht = draw(ht, background="transparent", heatmap_legend_side="left")
	dev.off()
}