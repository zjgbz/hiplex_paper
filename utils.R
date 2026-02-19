library(ComplexHeatmap)

# calculate the size of a complextheatmap object,
# so that we can suave it to a file with the correct size
calc_ht_size = function(ht, unit = "inch") {
	pdf(NULL)
	ht = draw(ht)
	w = ComplexHeatmap:::width(ht)
	w = convertX(w, unit, valueOnly = TRUE)
	h = ComplexHeatmap:::height(ht)
	h = convertY(h, unit, valueOnly = TRUE)
	dev.off()
	c(w, h)
}

# generate a list of target pairs based on the vector
# of target names, or two vectors of target names, with
# considering the order of the target names () # to do: what is this kind of order?
single_target_pair_sort_generation <- function(target_1, target_2) {
	sort_target_1_target_2 = sort(c(target_1, target_2), method="radix")
	target_pair = paste(sort_target_1_target_2, collapse='-')
	return (target_pair)
}

target_pair_generation <- function(target_list) {
	if (class(target_list) == "character") {
		target_list_1 = target_list
		target_list_2 = target_list
	} else if (class(target_list) == "list") {
		target_type_list = names(target_list)
		target_list_1 = target_list[[target_type_list[1]]]
		target_list_2 = target_list[[target_type_list[2]]]
	}
	target_pair_list = c()
	for (target_1 in target_list_1) {
		for (target_2 in target_list_2) {
			target_pair_i = single_target_pair_sort_generation(target_1, target_2)
			target_pair_list = c(target_pair_list, target_pair_i)
		}
	}
	target_pair_list = unique(target_pair_list)
	return (target_pair_list)
}

# # convert the rownames of the dataframe to the first column
# rownames2first_col <- function(df, colname) {
# 	rownames_df = as.data.frame(rownames(df))
# 	colnames(rownames_df) = colname
# 	df = cbind(rownames_df, df)
# 	rownames(df) = NULL
# 	return (df)
# }

# # read the feather dataframe and convert the first column
# # to rownames
# read_wgc <- function(wgc_dir, wgc_filename, remove_all_zero_rows=TRUE) {
#     wgc_dir_filename = file.path(wgc_dir, wgc_filename)
#     wgc = read_feather(wgc_dir_filename)
#     pos_list = wgc$pos
#     wgc = wgc[, -1]
#     rownames(wgc) = pos_list
#     if (remove_all_zero_rows) {
#         wgc = wgc[rowSums(wgc) != 0, ]
#         pos_list = rownames(wgc)
#     }
#     target_pair_list = colnames(wgc)

#     out = list(target_pair_list=target_pair_list, wgc=wgc, pos_list=pos_list)
#     return(out)
# }

k_means_optimal <- function(df, k_max=20, iter.max=1000) {
	# Use map_dbl to run many models with varying value of k (centers)
	tot_withinss <- map_dbl(1:k_max, function(k){
		set.seed(123)
		# model <- kmeans(x=wgc_standard_row, iter.max=1000, centers=k)
		model <- kmeans(x=df, iter.max=1000, centers=k)
		model$tot.withinss
	})
	# Generate a data frame containing both k and tot_withinss
	elbow_df <- data.frame(
		k = 1:k_max,
		tot_withinss = tot_withinss
	)

	# Angle function
	wss_list = elbow_df$tot_withinss
	cluster_list = elbow_df$k
	cluster_num = length(wss_list)

	angle_list = c(360)
	end_idx = cluster_num - 1
	start_idx = 2
	for (cluster_idx in start_idx :end_idx) {
		upper_point = c(cluster_list[cluster_idx - 1], wss_list[cluster_idx - 1])
		middle_point = c(cluster_list[cluster_idx], wss_list[cluster_idx])
		lower_point = c(cluster_list[cluster_idx + 1], wss_list[cluster_idx + 1])
		angle_i = Angle(upper_point, middle_point, lower_point)
		angle_list = c(angle_list, angle_i)
	}

	optimal_k = cluster_list[which(angle_list == min(angle_list))]

	# Plot the elbow plot
	elbow_plot = ggplot(elbow_df, aes(x = k, y = tot_withinss)) +
						geom_line() + geom_point() +
						geom_vline(xintercept = optimal_k, linetype = 2) +
						scale_x_continuous(breaks = 1:k_max) + 
						theme(axis.title.y = element_text(angle = 90))

	output = list("optimal_k"=optimal_k, "elbow_df"=elbow_df, "elbow_plot"=elbow_plot)
	return (output)
}

kmeans_result2label <- function(df, optimal_k, to_letters, label_colname="Cluster", iter.max=1000, seed=123) {
	set.seed(seed)
	model_opt = kmeans(df, iter.max=iter.max, centers=optimal_k, nstart=50)
	cluster_label = as.data.frame(model_opt$cluster)
	row_order = rownames(cluster_label)
	if (to_letters == TRUE) {
		row_split = LETTERS[cluster_label[, 1]]
	} else if (to_letters == FALSE) {
		row_split = cluster_label[, 1]
	}
	colnames(cluster_label) = label_colname
	output = list("cluster_label"=cluster_label, "row_order"=row_order, "row_split"=row_split)
	return(output)
}

filter_target_pairs <- function(percentage_cutoff = 0.25, target_pairs=NULL, frag_len_num_file="../data/utils/frag_split_fastq-demux_num_V_peak-valley.tsv") {
  frag_len_num <- read.table(frag_len_num_file, sep = "\t", header = TRUE, row.names = 1)
  frag_len_num$total_frag <- rowSums(frag_len_num)
  percentile_res <- quantile(frag_len_num$total_frag, type=3, probs = c(percentage_cutoff))
  # print(percentile_res)
  frag_len_num_cutoff <- unname(percentile_res)
  frag_len_num_filtered <- frag_len_num[frag_len_num$total_frag >= frag_len_num_cutoff, ]
  filtered_target_pairs <- rownames(frag_len_num_filtered)
  if (!is.null(target_pairs)) {
    filtered_target_pairs <- filtered_target_pairs[filtered_target_pairs %in% target_pairs]
  }
  return(filtered_target_pairs)
}



target_pair_mapping_file_36x <- "../data/utils/target_pair_short_hand.csv"
target_pair_mapping_df <- read.table(target_pair_mapping_file_36x, header = TRUE, check.names = FALSE)

map_target_names <- function(target_pair_list, target_pair_mapping_df, from = "targets", to = "shorthand" ) {
  cur_names <- target_pair_mapping_df[[from]]
  new_names <- target_pair_mapping_df[[to]]
  result <- target_pair_list
  for (i in seq_along(cur_names)) {
    cur_name <- cur_names[i]
    new_name <- new_names[i]
    result <- gsub(cur_name, new_name, result, fixed = TRUE)
  }
  return(result)
}

gamma_utf8 <- "\u03FF"