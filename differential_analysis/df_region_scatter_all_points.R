library(arrow)
library(GenomicRanges)
library(glue)
library(ggpointdensity)
library(viridis)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(ggrastr)
library(data.table)
library(dplyr)
library(circlize)

# read rnaseq dif result
rnaseq_dir_filename ="/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/diff_expr_all_update.csv"
rnaseq_raw = read.table(rnaseq_dir_filename, sep=",", check.names=FALSE, header=TRUE, row.names=1)
rnaseq_raw = rnaseq_raw[!is.na(rnaseq_raw$log2FoldChange),]
rnaseq_raw[, "log2FoldChange_TPM"] = rnaseq_raw$log2_T - rnaseq_raw$log2_V
# get tss annotation
promoter_range <- "0-0"
promoter_file <- glue("/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/ccre_region/promoter_-{promoter_range}.tsv")
promoter_table <- read.table(promoter_file, sep = "\t", header = TRUE)
promoter_table$seqnames <- paste0("chr", promoter_table$seqnames)
promoter_grange <- makeGRangesFromDataFrame(promoter_table, seqnames.field = "seqnames", start.field = "start", end.field = "end", strand.field = "strand",keep.extra.columns = TRUE)


save_dir <- '/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_fig'
dir.create(save_dir, recursive = TRUE)
plts <- list()
ylim_top <- 7
ylim_bottom <- -4
xlim_top <- 3
xlim_bottom <- -3

# cluster_id_list = c(1:8, 10:15)
cluster_id_list = c(1, 12)
load_dir = "/dcs05/hongkai/data/next_cutntag/bulk/df_analysis/800/column_cluster_result"

l2fc_thres = 0.5

for (cluster_id in cluster_id_list) {
	file_name = glue("result_post-limmanorm_post-filter-one_condition_nonzero-2_rowmean-0.25_column_cluster-{cluster_id}_limma.feather")

	file_path = file.path(load_dir, file_name)
	pc <- read_feather(file_path)

	colnames(pc)[1] <- "pos"
	pc <- pc[, c("pos", "adj.P.Val", "logFC")]
	pc <- cbind(pc, do.call(rbind, strsplit(pc$pos, "[_]")))
	df_wgc <- pc[, c("logFC", "adj.P.Val", "1", "2", "3", "pos")]
	colnames(df_wgc)[ncol(df_wgc)] <- "genomic_bin_pos"
	df_wgc[,"mid"] <- round((as.numeric(df_wgc[, "3"]) - as.numeric(df_wgc[, "2"])) %/% 2) + as.numeric(df_wgc[, "2"])
	gr_wgc <- makeGRangesFromDataFrame(df_wgc, seqnames.field = "1", start.field = "mid", end.field = "mid", keep.extra.columns = TRUE)

	nearest_dist_cutoff = 5000

	overlap_res <- distanceToNearest(gr_wgc, promoter_grange)
	overlap_dist <- mcols(overlap_res)$distance
	overlap_res <- overlap_res[overlap_dist <= nearest_dist_cutoff]
	overlap_df <- cbind(mcols(gr_wgc[queryHits(overlap_res),]), mcols(promoter_grange[subjectHits(overlap_res),]))
	overlap_df <- cbind(overlap_df,rnaseq_raw[overlap_df$gene_id, c("log2FoldChange_TPM", "padj", "NaB.1.RNA.seq_S3", "NaB.2.RNA.seq_S4", "veh.1.RNA.seq_S1", "veh.2.RNA.seq_S2")])
	# Load necessary library
	overlap_df <- overlap_df[which(overlap_df$padj <= 0.05),]
	# overlap_df <- overlap_df[which(abs(overlap_df$log2FoldChange_TPM) >= 1),]
	overlap_df <-
		overlap_df[, c(
			"logFC",
			"adj.P.Val",
			"genomic_bin_pos",
			"pos",
			"gene_id",
			"gene_name",
			"log2FoldChange_TPM",
			"padj",
			"NaB.1.RNA.seq_S3",
			"NaB.2.RNA.seq_S4",
			"veh.1.RNA.seq_S1",
			"veh.2.RNA.seq_S2"
		)]

	df_to_draw <- overlap_df
	df_to_draw <- df_to_draw[order(df_to_draw$logFC),]
	df_to_draw_raw <- as.data.frame(df_to_draw)
	df_summarized <- df_to_draw_raw %>%
		group_by(gene_id) %>%
		summarise(
			gene_name = first(gene_name),
			pos = first(pos),
			genomic_bin_pos = paste(unique(genomic_bin_pos), collapse = ","),
			adj.P.Val = min(adj.P.Val, na.rm = TRUE),
			logFC = mean(logFC, na.rm = TRUE),
			padj = first(padj),
			log2FoldChange_TPM = first(log2FoldChange_TPM),
			.groups = "drop"
			)

	df_to_draw <- as.data.frame(df_summarized) %>%
		mutate(
			quadrant = case_when(
				logFC >= l2fc_thres & log2FoldChange_TPM >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 1",
				logFC < -l2fc_thres & log2FoldChange_TPM >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 2",
				logFC < -l2fc_thres & log2FoldChange_TPM < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 3",
				logFC >= l2fc_thres & log2FoldChange_TPM < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 4"
			)
		)

	df_to_draw_clean <- df_to_draw[!is.na(df_to_draw$quadrant), ]
	df_to_draw_clean_filename = glue("result_merge_bins_column_cluster-{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.tsv")
	df_to_draw_clean_dir_filename = file.path(load_dir, df_to_draw_clean_filename)
	write.table(df_to_draw_clean, df_to_draw_clean_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)

	# Create all expressed genes
	df_to_nonfilter <- as.data.frame(df_summarized) %>%
		mutate(
			quadrant = case_when(
				logFC >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 1",
				logFC < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 2",
				logFC < -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 3",
				logFC >= l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 4"
			)
		)

	df_to_nonfilter_clean <- df_to_nonfilter[!is.na(df_to_nonfilter$quadrant), ]
	df_to_nonfilter_clean_filename = glue("result_gene-non-filter_merge_bins_column_cluster-{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.tsv")
	df_to_nonfilter_clean_dir_filename = file.path(load_dir, df_to_nonfilter_clean_filename)
	write.table(df_to_nonfilter_clean, df_to_nonfilter_clean_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)

	# Create non filtered expressed genes
	df_to_nondiff <- as.data.frame(df_summarized) %>%
		mutate(
			quadrant = case_when(
				logFC >= l2fc_thres & log2FoldChange_TPM < l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 1",
				logFC < -l2fc_thres & log2FoldChange_TPM < l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 2",
				logFC < -l2fc_thres & log2FoldChange_TPM > -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 3",
				logFC >= l2fc_thres & log2FoldChange_TPM > -l2fc_thres & adj.P.Val <= 0.25 ~ "Quadrant 4"
			)
		)

	df_to_nondiff_clean <- df_to_nondiff[!is.na(df_to_nondiff$quadrant), ]
	df_to_nondiff_clean_filename = glue("result_gene-non-diff_merge_bins_column_cluster-{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.tsv")
	df_to_nondiff_clean_dir_filename = file.path(load_dir, df_to_nondiff_clean_filename)
	write.table(df_to_nondiff_clean, df_to_nondiff_clean_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)

	# Create line chart
	df_to_draw[is.na(df_to_draw$quadrant), "quadrant"] <- "NONE"
	quadrant_counts_tmp <- df_to_draw %>%
		group_by(quadrant) %>%
		summarise(count = n())
	quadrant_counts_tmp <- as.data.frame(quadrant_counts_tmp)
	rownames(quadrant_counts_tmp) <- quadrant_counts_tmp$quadrant
	all_quadrants <- c("Quadrant 1", "Quadrant 2", "Quadrant 3", "Quadrant 4")
	quadrant_counts <- data.frame(quadrant=all_quadrants, count=0)
	rownames(quadrant_counts) <- quadrant_counts$quadrant
	quadrant_counts[all_quadrants,"count"] <- quadrant_counts_tmp[all_quadrants, "count"]
	quadrant_counts[is.na(quadrant_counts$count), "count"] <- 0
	test_tb <-
		matrix(
			c(
				quadrant_counts$count[2],
				quadrant_counts$count[1],
				quadrant_counts$count[3],
				quadrant_counts$count[4]
			),
			nrow=2, 
			byrow = TRUE
		)
	test_tb <- test_tb + 1
	fisher_test_result <- fisher.test(test_tb)
	odds_ratio <- 1/fisher_test_result$estimate
	pval <- fisher_test_result$p.value
	if (pval < 0.01) {
		annotation_text <- paste0("OR=", round(odds_ratio, 2), 
	                            ", P<0.01 ")
	} else {
		annotation_text <- paste0("Odds Ratio=", round(odds_ratio, 2), 
	                            ", P=", round(pval, 4))
	}
	annotation_text <- glue("Cluster {cluster_id}: {annotation_text}")
	# max_x <- max(abs(df_to_draw$logFC))+0.4
	# max_y <- max(abs(df_to_draw$log2FoldChange_TPM))+0.4
	# min_y <- min((df_to_draw$log2FoldChange_TPM))-0.4
	max_x <- xlim_top-0.4
	min_x <- xlim_bottom+0.4
	max_y <- ylim_top-0.4
	min_y <- ylim_bottom+0.4
	quad_positions <- data.frame(
		quadrant = c("Quadrant 1", "Quadrant 2", "Quadrant 3", "Quadrant 4"),
		x = c(max_x, min_x, min_x, max_x),  # X position of the text
		y = c(max_y, max_y, min_y, min_y)   # Y position of the text
	)
	overlap_df <- as.data.frame(overlap_df)
	overlap_df[, c("NaB.1.RNA.seq_S3", "NaB.2.RNA.seq_S4", "veh.1.RNA.seq_S1", "veh.2.RNA.seq_S2")] <- log10(overlap_df[, c("NaB.1.RNA.seq_S3", "NaB.2.RNA.seq_S4", "veh.1.RNA.seq_S1", "veh.2.RNA.seq_S2")]+1)
	quad_positions <- merge(quad_positions, quadrant_counts, by = "quadrant")
	df_to_draw_subset <- df_to_draw[df_to_draw$quadrant!="NONE",]
	p <- ggplot(df_to_draw, aes(x = logFC, y = log2FoldChange_TPM)) +
		ggrastr::rasterise(geom_point(size = 1.5, color="#e9ecef")) +  # Add points for better visualization
		ggrastr::rasterise(geom_pointdensity(data = df_to_draw_subset, aes(logFC, log2FoldChange_TPM), size = 1.5)) +  # Add points for better visualization
		geom_text(data = quad_positions, aes(x = x, y = y, label = paste(count)), size = 6.5, color = "#367DB0", fontface = "bold") +
		scale_color_viridis()+
		labs(title = "",
			x = "Hi-Plex log2FC",
			y = "Gene Expression log2FC") +
		theme_minimal() + # Use a clean theme
		geom_hline(yintercept = 0, color = "#adb5bd", linewidth = 1, linetype = "dashed") +  # Add y = 0 line
		geom_vline(xintercept = 0, color = "#adb5bd", linewidth = 1, linetype = "dashed") +  # Add x = 0 line
		theme(legend.position = "none",
			axis.text.x = element_text(size = 16),
			axis.title.x = element_text(size = 18),
			axis.title.y = element_text(size = 18),
			axis.text.y = element_text(size = 16),
			panel.grid = element_blank(),  # Remove all grid lines
			panel.grid.major = element_blank(),  # Remove major grid lines
			panel.grid.minor = element_blank(),  # Remove minor grid lines
			panel.background = element_rect(fill = "white", color = "black"),
			panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
			plot.title = element_text(size = 18)
		) + 
		geom_smooth(data = df_to_draw_subset, aes(logFC, log2FoldChange_TPM), method = "lm", color = "red", se = TRUE) + 
		ggtitle(annotation_text) + xlim(xlim_bottom, xlim_top) + ylim(ylim_bottom, ylim_top)
	# annotate("text", x = -max_x+1, y = min_y + 2, 
	#          label = annotation_text, size = 4.75, fontface = "italic") + xlim(-max_x-0.5, max_x+0.5)
	scatter_plot_single_filename = glue("scatter_{cluster_id}_{nearest_dist_cutoff}_l2fc-{l2fc_thres}.pdf")
	scatter_plot_single_dir_filename = file.path(save_dir, scatter_plot_single_filename)
	pdf(scatter_plot_single_dir_filename, width=8/1.5, height=8/1.5)
	print(p)
	dev.off()
}