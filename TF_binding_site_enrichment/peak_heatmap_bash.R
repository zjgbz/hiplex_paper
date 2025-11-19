library(ComplexHeatmap)
library(plyranges)
library(rtracklayer)
library(matrixTests)
library(latex2exp)
library(glue)
library(circlize)

# Rscript /dcs05/hongkai/data/next_cutntag/script/motif_analysis/peak_heatma_bash.R /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw /dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/hm/gw/peak_analysis/
# options <- c("/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/gw", "/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/hm/gw/peak_analysis/")
# motif_result_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/motif_analysis/df_analysis/limma_fdr_0.1_down"
# folders <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O")
options <- commandArgs(trailingOnly = TRUE)
if (length(options) != 0) {
  motif_result_dir <- options[1]
  out_dir <- options[2]
}
print(motif_result_dir)
motif_results_dirs <- Sys.glob(glue("{motif_result_dir}/*/peak_result/result.tsv"))
# FDR_cutoff <- 0.8
# odds_ratio_cutoff <- 1.2

# read motif result when there is only one motif file for each cluster
motif_results <- list()
folders <- c()
for (motif_result_file in motif_results_dirs) {
  motif_result <- read.table(motif_result_file, sep="\t")
  folder <- gsub(glue("{motif_result_dir}/"), "", motif_result_file)
  folder <- gsub(glue("/peak_result/result.tsv"), "", folder)
  motif_results[[folder]] <- motif_result
  folders <- append(folders, folder)
}
if (is.na(as.numeric(folders[[1]]))) {
  folders <- sort(folders)
} else {
  folders <- as.character(sort(as.numeric(folders)))
}

motif_results <- motif_results[folders]

motif_odds_ratio <- list()
motif_result_file_i <- motif_results_dirs[1]
motif_result <- read.table(motif_result_file_i, sep="\t")
motif_id_order <- rownames(motif_result)
for (folder in folders) {
  motif_file_results <- list()
  motif_result_file_i <- glue("{motif_result_dir}/{folder}/peak_result/result.tsv")
  # motif_result_file_i <- paste0(motif_result_dir, folder, "/peak_result/result.tsv")
  motif_result <- read.table(motif_result_file_i, sep="\t")

  result_odds_ratio <- motif_result[motif_id_order, "odds_ratio"]
  motif_odds_ratio[[folder]] <- result_odds_ratio
}



motid_id_order <- rownames(motif_results[[folders[1]]])
motif_odds_ratio <- list()
FDRs <- list()
for (folder in folders) {
  motif_odds_ratio[[folder]] <- motif_results[[folder]][motid_id_order,]$odds_ratio
  FDRs[[folder]] <- motif_results[[folder]][motid_id_order,]$FDR
}
motif_odds_ratio_df <- as.data.frame(motif_odds_ratio, check.names=FALSE)
FDR_df <- as.data.frame(FDRs, check.names=FALSE)
rownames(motif_odds_ratio_df) <- motid_id_order
rownames(FDR_df) <- motid_id_order
log10_fdr_df <- -log10(FDR_df)


target_hits <- list()
for (folder in folders) {
  target_hits[[folder]] <- motif_results[[folder]][motid_id_order,]$target_hit
}
target_hits_df <- as.data.frame(target_hits)
rownames(target_hits_df) <- motid_id_order
target_hits_thresh <- unname(quantile(as.vector(as.matrix(target_hits_df)), 0.1))


motif_odds_ratio_df <- motif_odds_ratio_df[rowSums(motif_odds_ratio_df >= 2) >=1,]
# motif_odds_ratio_df <- motif_odds_ratio_df[rowSums(motif_odds_ratio_df >= 2) >=1,] 
FDR_df <- FDR_df[rowSums(FDR_df <= 0.05) >= 1,]
target_hits_df <- target_hits_df[rowSums(target_hits_df >= target_hits_thresh) >= 1,]
shared_sets <- intersect(rownames(motif_odds_ratio_df), rownames(log10_fdr_df))
shared_sets <- intersect(shared_sets, rownames(target_hits_df))
motif_odds_ratio_df <- motif_odds_ratio_df[shared_sets,]
FDR_df <- FDR_df[shared_sets,]

treatments <- do.call(rbind, strsplit(rownames(motif_odds_ratio_df), ")", fixed=TRUE))[,1]
targets <- do.call(rbind, strsplit(rownames(motif_odds_ratio_df), ")", fixed=TRUE))[,2]
treatments[!grepl("[(]", treatments)] <- "no_treatment"

target_treat_df <- data.frame(
  targets=targets,
  treatments=treatments
)
motif_odds_combined <- cbind(target_treat_df, motif_odds_ratio_df)

no_treat_targets <- target_treat_df$targets[target_treat_df$treatments!="no_treatment"]
treat_targets <- target_treat_df$targets[target_treat_df$treatments=="no_treatment"]
target_treat_df_2_cond <- no_treat_targets[no_treat_targets %in% treat_targets]
target_treat_df_2_cond_df <- target_treat_df[target_treat_df$targets %in% target_treat_df_2_cond,]

result_list <- list() 
splits <- c()
result_df <- data.frame(matrix(ncol = ncol(motif_odds_ratio_df), nrow = 0))
colnames(result_df) <- colnames(motif_odds_ratio_df)
unique_targets <- unique(target_treat_df_2_cond_df$targets)
for (i in seq_along(unique_targets)) {
  target <- unique_targets[i]
  motif_odds_ratio_df_subset <- motif_odds_combined[motif_odds_combined$targets==target, !colnames(motif_odds_combined) %in% c("targets", "treatments")]
  motif_odds_ratio_df_subset <- motif_odds_ratio_df_subset[order(rownames(motif_odds_ratio_df_subset), decreasing = TRUE),]
  result_list[[target]] <- motif_odds_ratio_df_subset
  print(motif_odds_ratio_df_subset)
  print(nrow(motif_odds_ratio_df_subset))
  splits <- append(splits, rep(i, nrow(motif_odds_ratio_df_subset)))
  print(rep(i, nrow(motif_odds_ratio_df_subset)))
  result_df <- rbind(result_df, motif_odds_ratio_df_subset)
}


result_df_no_crispr <- data.frame(matrix(ncol = ncol(motif_odds_ratio_df), nrow = 0))
colnames(result_df_no_crispr) <- colnames(motif_odds_ratio_df)
unique_targets <- unique(target_treat_df_2_cond_df$targets)
splits <- c()
j <- 1
for (i in seq_along(unique_targets)) {
  target <- unique_targets[i]
  motif_odds_ratio_df_subset <- motif_odds_combined[motif_odds_combined$targets==target, !colnames(motif_odds_combined) %in% c("targets", "treatments")]
  motif_odds_ratio_df_subset <- motif_odds_ratio_df_subset[order(rownames(motif_odds_ratio_df_subset), decreasing = TRUE),]
  motif_odds_ratio_df_subset <- motif_odds_ratio_df_subset[!grepl("CRISPR", rownames(motif_odds_ratio_df_subset)), ]
  if (nrow(motif_odds_ratio_df_subset) > 1) {
    print(motif_odds_ratio_df_subset)
    print(nrow(motif_odds_ratio_df_subset))
    splits <- append(splits, rep(i, nrow(motif_odds_ratio_df_subset)))
    print(rep(j, nrow(motif_odds_ratio_df_subset)))
    j <- j + 1
    result_df_no_crispr <- rbind(result_df_no_crispr, motif_odds_ratio_df_subset)
  }
}
rownames(result_df_no_crispr) <- gsub(" N-terminal eGFP-tagged", "", rownames(result_df_no_crispr), fixed=TRUE)
rownames(result_df_no_crispr)  <- gsub(".*\\((.*?)\\).*", "\\1", rownames(result_df_no_crispr))


selected_tfs <- c("NR2C2", "JUND", "ZNF24", "DDX20", "BACH1", "GATA2", "PYGO2")
selected_exp_tfs <- glue("stably expressing {selected_tfs}")
up_arrow <- "$\\uparrow$"
result_df_no_crispr_selected <- result_df_no_crispr[c(rbind(selected_tfs, selected_exp_tfs)),]
rownames(result_df_no_crispr_selected) <- gsub("stably expressing ", up_arrow, rownames(result_df_no_crispr_selected), fixed = TRUE)
splits = c(rbind(1:(
  nrow(result_df_no_crispr_selected) / 2
), 1:(
  nrow(result_df_no_crispr_selected) / 2
)))

motif_odds_ratio_df_no_perturb <- motif_odds_ratio_df[!grepl("(", rownames(motif_odds_ratio_df), fixed = TRUE),]
FDR_df <- FDR_df[rownames(motif_odds_ratio_df_no_perturb), colnames(motif_odds_ratio_df_no_perturb)]

rownames(motif_odds_ratio_df_no_perturb) <- gsub("POLR2AphosphoS2", "POLR2A pS2", rownames(motif_odds_ratio_df_no_perturb))
rownames(motif_odds_ratio_df_no_perturb) <- gsub("POLR2AphosphoS5", "POLR2A pS5", rownames(motif_odds_ratio_df_no_perturb))
rownames(FDR_df) <- gsub("POLR2AphosphoS2", "POLR2A pS2", rownames(FDR_df))
rownames(FDR_df) <- gsub("POLR2AphosphoS5", "POLR2A pS5", rownames(FDR_df))


log2_df <- log2(motif_odds_ratio_df_no_perturb)

col_fun = colorRamp2(c(-max(abs(log2_df)), 0, max(abs(log2_df))), c("#3155C3", "white", "#AF0525"))


pdf(glue("{out_dir}/peak_no_perturb.pdf"), width = 8, height=45)

h <- Heatmap(log2(motif_odds_ratio_df_no_perturb), name = "log(Odds Ratio)", col = col_fun, cluster_columns = FALSE)
ht <- draw(h)
dev.off()
row_order <- rownames(motif_odds_ratio_df_no_perturb)[row_order(ht)]

write.csv(log2_df[row_order,], glue("{out_dir}/peak_no_perturb.csv"), quote = FALSE)
write.csv(FDR_df[row_order,], glue("{out_dir}/peak_no_perturb_fdr.csv"), quote = FALSE)


top_ns <- 3:10
selected_folders <- folders
heatmap_dir <- out_dir
dir.create(heatmap_dir, recursive = TRUE)


i = 5

selected_motif_list <- list()
selected_motifs <- c()
for (folder in selected_folders) {
  motif_odds_ratio_mat_sub <- motif_odds_ratio_df_no_perturb
  odds_ratio_i <- motif_odds_ratio_mat_sub[, folder, drop=FALSE]
  odds_ratio_i <- odds_ratio_i[order(odds_ratio_i[[folder]], decreasing = TRUE),,drop=FALSE]
  selected_motifs <- append(selected_motifs, rownames(odds_ratio_i)[1:i])
  selected_motif_list[[folder]] <- rownames(odds_ratio_i)[1:i]

}
selected_motifs <- unique(selected_motifs)


motif_odds_ratio_mat_to_draw <- motif_odds_ratio_df_no_perturb[selected_motifs,]

motif_odds_ratio_mat_to_draw_transpose <- t(motif_odds_ratio_mat_to_draw)
log2_motif_odds_ratio_mat_to_draw <- log2(motif_odds_ratio_mat_to_draw_transpose)
h <- ComplexHeatmap::Heatmap(log2_motif_odds_ratio_mat_to_draw, 
                              column_names_gp = grid::gpar(fontsize = 25),
                              row_names_gp = grid::gpar(fontsize = 25),
                              cluster_columns = F, 
                              cluster_rows = F, 
                              column_names_rot = 45,
                              row_names_side = "left", 
                              col = col_fun,
                              heatmap_legend_param = list(
                                title = "log(Odds Ratio)",
                                color_bar = "continuous",
                                title_gp = gpar(fontsize = 24),
                                legend_direction = "horizontal",
                                labels_gp = gpar(fontsize = 20) ,
                                legend_width = unit(6, "cm")
                              ))
pdf(glue("{heatmap_dir}/reorganized_heatmap_filtered_rna_subtraction_{i}_t.pdf"), width = 36, height = 19)
draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
dev.off()

scaled_df <- as.data.frame(scale(log2_motif_odds_ratio_mat_to_draw))
col_fun_z = colorRamp2(c(min(scaled_df), 0, max(scaled_df)), c("#3155C3", "white", "#AF0525"))

h <- ComplexHeatmap::Heatmap(scaled_df, 
                              column_names_gp = grid::gpar(fontsize = 25),
                              row_names_gp = grid::gpar(fontsize = 25),
                              cluster_columns = F, 
                              cluster_rows = F, 
                              column_names_rot = 45,
                              row_names_side = "left", 
                              col = col_fun_z,
                              heatmap_legend_param = list(
                                title = "Z-score",
                                color_bar = "continuous",
                                title_gp = gpar(fontsize = 24),
                                legend_direction = "horizontal",
                                labels_gp = gpar(fontsize = 20) ,
                                legend_width = unit(6, "cm")
                              ))
pdf(glue("{heatmap_dir}/reorganized_heatmap_filtered_rna_subtraction_{i}_t_zscore.pdf"), width = 36, height = 19*2/3)
draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
dev.off()

h <- ComplexHeatmap::Heatmap(scaled_df, 
                              column_names_gp = grid::gpar(fontsize = 25),
                              row_names_gp = grid::gpar(fontsize = 25),
                              cluster_columns = TRUE, 
                              cluster_rows = F, 
                              column_names_rot = 45,
                              row_names_side = "left", 
                              col = col_fun_z,
                              heatmap_legend_param = list(
                                title = "Z-score",
                                color_bar = "continuous",
                                title_gp = gpar(fontsize = 24),
                                legend_direction = "horizontal",
                                labels_gp = gpar(fontsize = 20) ,
                                legend_width = unit(6, "cm")
                              ))
ht_drawn <- draw(h)

# Extract column order
col_order <- column_order(ht_drawn)

# Get the actual column names in that order
ordered_colnames <- rev(colnames(scaled_df)[col_order])
h <- ComplexHeatmap::Heatmap(scaled_df[,ordered_colnames], 
                              column_names_gp = grid::gpar(fontsize = 25),
                              row_names_gp = grid::gpar(fontsize = 25),
                              cluster_columns = F, 
                              cluster_rows = F, 
                              column_names_rot = 45,
                              row_names_side = "left", 
                              col = col_fun_z,
                              heatmap_legend_param = list(
                                title = "Z-score",
                                color_bar = "continuous",
                                title_gp = gpar(fontsize = 24),
                                legend_direction = "horizontal",
                                labels_gp = gpar(fontsize = 20) ,
                                legend_width = unit(6, "cm")
                              ))
pdf(glue("{heatmap_dir}/reorganized_heatmap_filtered_rna_subtraction_{i}_t_zscore_clustered.pdf"), width = 36, height = 19*2/3)
draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
dev.off()

h <- ComplexHeatmap::Heatmap(log2(motif_odds_ratio_mat_to_draw_transpose), 
                              column_names_gp = grid::gpar(fontsize = 25),
                              row_names_gp = grid::gpar(fontsize = 25),
                              cluster_columns = T, 
                              cluster_rows = F, 
                              column_names_rot = 45,
                              row_names_side = "left", 
                              col = col_fun,
                              heatmap_legend_param = list(
                                title = "log(Odds Ratio)",
                                color_bar = "continuous",
                                title_gp = gpar(fontsize = 24),
                                legend_direction = "horizontal",
                                labels_gp = gpar(fontsize = 20) ,
                                legend_width = unit(6, "cm")
                              ))
ht <- draw(h, heatmap_legend_side="bottom", padding = unit(c(8, 8, 8, 8), "mm"))
pdf(glue("{heatmap_dir}/reorganized_heatmap_filtered_rna_subtraction_{i}_t_cluster.pdf"), width = 36, height = 19)
print(ht)
dev.off()

