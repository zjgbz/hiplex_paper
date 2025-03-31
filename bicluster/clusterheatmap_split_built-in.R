rm(list=ls())
library(ComplexHeatmap)
library(arrow)
library(rtracklayer)
library(tibble)
library(latex2exp)
library(glue)
source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")

calc_ht_size = function(ht, heatmap_legend_side=NULL, unit = "inch") {
	pdf(NULL)
	ht = draw(ht, heatmap_legend_side=heatmap_legend_side, background = "transparent")
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly = TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly = TRUE)
	dev.off()

	c(w, h)
}

col_km = 16
row_km = 15
qc_type = "all-qc"
random_seed = 42
frag_type = "mixed"
clustering_method = "kmeans"
# distance_method_list = c("pearson", "euclidean")
distance_method_list = c("euclidean")
show_dend = "on"

wgc_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800"
heatmap_dir = wgc_dir
wgc_filename = "V1V2_mixed_800_colQC-all-qc_libnorm_noAllZero_log2_qnorm_gam-40_per-0.01_mean-0.99.feather"
scen = "V"
wgc_dir_filename = file.path(wgc_dir, wgc_filename)
wgc = read_feather(wgc_dir_filename)
wgc = column_to_rownames(wgc, var="pos")
colnames(wgc) <- map_target_names(colnames(wgc), target_pair_mapping_df)

if (show_dend == "on") {
	show_dend_boolean = TRUE
} else if (show_dend == "off") {
	show_dend_boolean = FALSE
}

for (distance_method in distance_method_list) {
	heatmap_object_prefix = glue("bicluster_{scen}_{frag_type}_{qc_type}_{clustering_method}_{distance_method}_row_num-{row_km}_column_num-{col_km}")
	
	set.seed(random_seed)
	ht = Heatmap(as.matrix(wgc), name=scen, border=TRUE,
		row_km=row_km, column_km=col_km,
		show_row_names=FALSE, show_column_names=TRUE,
		row_title=LETTERS[1:row_km],
		row_title_gp=gpar(fontsize=120),
		column_title_gp=gpar(fontsize=60),
		width=ncol(wgc)*unit(5, "mm"),
		height=ncol(wgc)*unit(5, "mm"),
		row_gap=unit(12, "mm"),
		column_gap=unit(12, "mm"),
		clustering_distance_rows=distance_method,
		clustering_distance_columns=distance_method,
		row_title_rot=0,
		show_row_dend=show_dend_boolean,
		show_column_dend=show_dend_boolean,
		show_heatmap_legend=TRUE,
		heatmap_legend_param = list(
				title = "Normalized Read Counts",
				legend_width = ncol(wgc)/5*unit(5, "mm"),
				grid_height = 7*unit(5, "mm"),
				title_position = "topcenter",
				title_gp = gpar(fontsize = 100),
				labels_gp = gpar(fontsize = 80),
				legend_direction = "horizontal"),
		column_labels=TeX(colnames(wgc)),
		use_raster=TRUE
		)
	set.seed(random_seed)
	size = calc_ht_size(ht, heatmap_legend_side="bottom")
	
	heatmap_filename = glue("{heatmap_object_prefix}_heatmap_dend-{show_dend}.pdf")
	heatmap_dir_filename = file.path(heatmap_dir, heatmap_filename)
	pdf(heatmap_dir_filename, width = size[1], height = size[2])
	set.seed(random_seed)
	ht_draw = draw(ht, heatmap_legend_side="bottom", background = "transparent")
	dev.off()
	
	# ht_filename = glue("{heatmap_object_prefix}.RData")
	# ht_dir_filename = file.path(heatmap_dir, ht_filename)
	# save(ht, file=ht_dir_filename)
	
	row_order_list = row_order(ht_draw)
	row_order_df = stack(row_order_list)
	colnames(row_order_df) = c("feature", "label")
	row_order_df$feature = rownames(wgc)[row_order_df$feature]
	# row_order_df$label = as.integer(as.character(row_order_df$label))
	row_order_df$label = LETTERS[row_order_df$label]
	row_order_df_filename = glue("{heatmap_object_prefix}_heatmap_row_clusters.tsv")
	row_order_df_dir_filename = file.path(heatmap_dir, row_order_df_filename)
	write.table(row_order_df, row_order_df_dir_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

	col_order_list = column_order(ht_draw)
	col_order_df = stack(col_order_list)
	colnames(col_order_df) = c("feature", "label")
	col_order_df$feature = colnames(wgc)[col_order_df$feature]
	col_order_df$label = as.integer(as.character(col_order_df$label))
	col_order_df_filename = glue("{heatmap_object_prefix}_heatmap_column_clusters.tsv")
	col_order_df_dir_filename = file.path(heatmap_dir, col_order_df_filename)
	write.table(col_order_df, col_order_df_dir_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

	row_cluster_filename = glue("{heatmap_object_prefix}_heatmap_row_clusters.tsv")
	row_cluster_dir_filename = file.path(heatmap_dir, row_cluster_filename)
	row_cluster = read.table(row_cluster_dir_filename, header=TRUE, sep="\t")
	row_cluster$label <- gsub("E", "P", row_cluster$label)
	row_cluster$label <- gsub("F", "Q", row_cluster$label)
	row_cluster$label <- gsub("D", "F", row_cluster$label)

	row_cluster$label <- gsub("P", "D", row_cluster$label)
	row_cluster$label <- gsub("Q", "E", row_cluster$label)

	row_cluster$label <- gsub("J", "R", row_cluster$label)
	row_cluster$label <- gsub("K", "S", row_cluster$label)
	row_cluster$label <- gsub("L", "T", row_cluster$label)
	row_cluster$label <- gsub("M", "U", row_cluster$label)
	row_cluster$label <- gsub("N", "V", row_cluster$label)
	row_cluster$label <- gsub("O", "W", row_cluster$label)
	row_cluster$label <- gsub("I", "O", row_cluster$label)

	row_cluster$label <- gsub("R", "I", row_cluster$label)
	row_cluster$label <- gsub("S", "J", row_cluster$label)
	row_cluster$label <- gsub("T", "K", row_cluster$label)
	row_cluster$label <- gsub("U", "L", row_cluster$label)
	row_cluster$label <- gsub("V", "M", row_cluster$label)
	row_cluster$label <- gsub("W", "N", row_cluster$label)
	row_cluster = row_cluster[order(row_cluster$label), ]
	row_cluster_reorg_filename = glue("{heatmap_object_prefix}_heatmap_row_clusters_manually_reorganized.tsv")
	row_cluster_reorg_dir_filename = file.path(heatmap_dir, row_cluster_reorg_filename)
	write.table(row_cluster, row_cluster_reorg_dir_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	row_order = row_cluster$feature
	row_split = factor(row_cluster$label, levels=unique(row_cluster$label))

	col_cluster_filename = glue("{heatmap_object_prefix}_heatmap_column_clusters.tsv")
	col_cluster_dir_filename = file.path(heatmap_dir, col_cluster_filename)
	col_cluster = read.table(col_cluster_dir_filename, header=TRUE, sep="\t")
	col_cluster$label <- gsub(16, "A", col_cluster$label)
	col_cluster$label <- gsub(15, "B", col_cluster$label)
	col_cluster$label <- gsub(10, "C", col_cluster$label)
	col_cluster$label <- gsub(11, "D", col_cluster$label)
	col_cluster$label <- gsub(12, "E", col_cluster$label)
	col_cluster$label <- gsub(13, "F", col_cluster$label)
	col_cluster$label <- gsub(14, "G", col_cluster$label)
	col_cluster$label <- gsub(9, "H", col_cluster$label)
	col_cluster$label <- gsub(8, "I", col_cluster$label)
	col_cluster$label <- gsub(6, "J", col_cluster$label)
	col_cluster$label <- gsub(7, "K", col_cluster$label)
	col_cluster$label <- gsub(3, "L", col_cluster$label)
	col_cluster$label <- gsub(4, "M", col_cluster$label)
	col_cluster$label <- gsub(5, "N", col_cluster$label)
	col_cluster$label <- gsub(1, "O", col_cluster$label)
	col_cluster$label <- gsub(2, "P", col_cluster$label)

	col_cluster$label <- gsub("A", 1, col_cluster$label)
	col_cluster$label <- gsub("B", 9, col_cluster$label)
	col_cluster$label <- gsub("C", 5, col_cluster$label)
	col_cluster$label <- gsub("D", 3, col_cluster$label)
	col_cluster$label <- gsub("E", 2, col_cluster$label)
	col_cluster$label <- gsub("F", 8, col_cluster$label)
	col_cluster$label <- gsub("G", 7, col_cluster$label)
	col_cluster$label <- gsub("H", 4, col_cluster$label)
	col_cluster$label <- gsub("I", 6, col_cluster$label)
	col_cluster$label <- gsub("J", 10, col_cluster$label)
	col_cluster$label <- gsub("K", 11, col_cluster$label)
	col_cluster$label <- gsub("L", 14, col_cluster$label)
	col_cluster$label <- gsub("M", 12, col_cluster$label)
	col_cluster$label <- gsub("N", 13, col_cluster$label)
	col_cluster$label <- gsub("O", 16, col_cluster$label)
	col_cluster$label <- gsub("P", 15, col_cluster$label)

	# convert col_cluster$label to integer
	col_cluster$label = as.integer(as.character(col_cluster$label))
	col_cluster = col_cluster[order(col_cluster$label), ]
	col_cluster_reorg_filename = glue("{heatmap_object_prefix}_heatmap_column_clusters_manually_reorganized.tsv")
	col_cluster_reorg_dir_filename = file.path(heatmap_dir, col_cluster_reorg_filename)
	write.table(col_cluster, col_cluster_reorg_dir_filename, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
	col_order = col_cluster$feature
	col_split = col_cluster$label

	set.seed(random_seed)
	wgc_reorder = wgc[row_order, col_order]
	ht = Heatmap(as.matrix(wgc_reorder), name=scen, border=TRUE,
		row_order=row_order, row_split=row_split,
		column_order=col_order, column_split=col_split,
		cluster_rows=FALSE, cluster_columns=FALSE,
		show_row_names=FALSE, show_column_names=TRUE,
		row_title_gp=gpar(fontsize=120),
		column_title_gp=gpar(fontsize=60),
		width=ncol(wgc_reorder)*unit(5, "mm"),
		height=ncol(wgc_reorder)*unit(5, "mm"),
		row_gap=unit(12, "mm"),
		column_gap=unit(12, "mm"),
		clustering_distance_rows=distance_method,
		clustering_distance_columns=distance_method,
		row_title_rot=0,
		show_heatmap_legend=TRUE,
		heatmap_legend_param = list(
				title = "Normalized Read Counts",
				legend_width = ncol(wgc_reorder)/5*unit(5, "mm"),
				grid_height = 7*unit(5, "mm"),
				title_position = "topcenter",
				title_gp = gpar(fontsize = 100),
				labels_gp = gpar(fontsize = 80),
				legend_direction = "horizontal"),
		column_labels=TeX(colnames(wgc_reorder)),
		use_raster=TRUE
		)
	set.seed(random_seed)
	size = calc_ht_size(ht, heatmap_legend_side="bottom")
	
	current_time = format(Sys.time(), "%Y-%m-%d_%H-%M-%S-%Z")
	heatmap_filename = glue("{heatmap_object_prefix}_heatmap_row-LETTERS_dend-off_manually_reorganized_{current_time}.pdf")
	heatmap_dir_filename = file.path(heatmap_dir, heatmap_filename)
	pdf(heatmap_dir_filename, width = size[1], height = size[2])
	set.seed(random_seed)
	ht_draw = draw(ht, heatmap_legend_side="bottom", background = "transparent")
	dev.off()
}