library(arrow)
library(matrixStats)
library(mgcv)
library(ggplot2)
library(cowplot)
library(ggExtra)
library(ggpointdensity)

args <- commandArgs(trailingOnly = TRUE)
fitting_model = args[1]
robust_check = as.numeric(args[2])
if (fitting_model == "gam") {
	k = as.numeric(args[3])
} else if (fitting_model == "loess") {
	span = as.numeric(args[3])
}
seed = as.numeric(args[4])
nrow_sample = as.numeric(args[5])
frag_type = args[6]
bin_size = args[7]
scen = args[8]
qc_type = args[9]
post_process = args[10]

font_size = 10

wgc_root_dir = "/dcs05/hongkai/data/next_cutntag/bulk/wgc/"
wgc_dir = file.path(wgc_root_dir, frag_type, bin_size)
# wgc_qnorm_dir = "/data/hji7/minzhi/multitag/wgc/pair_end/cobinding/mixed/800"
# wgc_qnorm_filename_prefix = paste0(scen, "_", frag_type, "_", bin_size, "_colQC-", qc_type, "_", post_process, "_qnorm")
wgc_qnorm_filename_prefix = paste0(scen, "_", frag_type, "_", bin_size, "_colQC-", qc_type, "_", post_process)
wgc_qnorm_filename = paste0(wgc_qnorm_filename_prefix, ".feather")
wgc_qnorm_dir_filename = file.path(wgc_dir, wgc_qnorm_filename)

# wgc_qnorm_raw = read_parquet(wgc_qnorm_dir_filename)
# wgc_qnorm = wgc_qnorm_raw[sample(nrow(wgc_qnorm_raw), 50000), ]

# wgc_qnorm = read_parquet(wgc_qnorm_dir_filename)
wgc_qnorm = read_feather(wgc_qnorm_dir_filename)

pos_list <- wgc_qnorm[, "pos"]
wgc_qnorm = wgc_qnorm[, -1]

data <- as.matrix(wgc_qnorm)

clip_max <- sqrt(ncol(data))
print(paste0("using max clipping value: ",clip_max))

# wgc_qnorm_MAV_stat_filename_prefix_prefix = paste0(wgc_qnorm_filename_prefix, "_MAV_stat_", "seed-", robust_check, "_noabs-logmean-logvar_", fitting_model)
wgc_qnorm_MAV_stat_filename_prefix_prefix = paste0(wgc_qnorm_filename_prefix, "_", fitting_model)
if (fitting_model == "gam") {
	# wgc_qnorm_MAV_stat_filename_prefix = paste0(wgc_qnorm_MAV_stat_filename_prefix_prefix, "_k-", k)
	wgc_qnorm_MAV_stat_filename_prefix = paste0(wgc_qnorm_MAV_stat_filename_prefix_prefix, "-", k)
} else if (fitting_model == "loess") {
	# wgc_qnorm_MAV_stat_filename_prefix = paste0(wgc_qnorm_MAV_stat_filename_prefix_prefix, "_span-", format(span, nsmall=3))
	wgc_qnorm_MAV_stat_filename_prefix = paste0(wgc_qnorm_MAV_stat_filename_prefix_prefix, "-", format(span, nsmall=3))
}

gene_mean_all <- rowMeans(data)
gene_var_all <- rowVars(data)

data_filter <- data
# data_filter[gene_mean_all < 0, ] = -data[gene_mean_all < 0, ]

gene_mean <- rowMeans(data_filter)
gene_var <- rowVars(data_filter)

data_fit <- data.frame(X=gene_mean, Y=gene_var)

fit_model_filename = paste0(wgc_qnorm_MAV_stat_filename_prefix, ".RData")
fit_model_dir_filename = file.path(wgc_dir, fit_model_filename)

if (file.exists(fit_model_dir_filename)) {
	load(fit_model_dir_filename)
	print("Model exists and has been loaded!")
} else {
	set.seed(robust_check)
	if (fitting_model == "gam") {
		fit_model <- gam(formula = log2(x=Y) ~ s(log2(x=X), k=k), data = data_fit)
		# fit_model <- gam(formula = log2(x=Y) ~ s(X, k=k), data = data_fit)
	} else if (fitting_model == "loess") {
		fit_model <- loess(formula = log2(x=Y) ~ log2(x=X), data=data_fit, span=span)
		# fit_model <- loess(formula = log2(x=Y) ~ X, data=data_fit, span=span)
	}
	save(fit_model, file=fit_model_dir_filename)
	print("Model has been calculated and saved!")
}

gene_var_expect <- 2^(fit_model$fitted)
gene_sd_expect <- sqrt(gene_var_expect)

gene_var_norm <- (data_filter - gene_mean)/gene_sd_expect
gene_var_norm_mean = rowMeans(gene_var_norm)
gene_var_norm_var = rowVars(gene_var_norm)
gene_hyper_var <- rowSums(gene_var_norm^2)/(ncol(data_filter) - 1)
gene_var_norm_clip = gene_var_norm

# record the clipped rows
clip_or_not = matrix(0, nrow = dim(gene_var_norm)[1], ncol = dim(gene_var_norm)[2])
clip_or_not[which(gene_var_norm > clip_max)] <- 1
clip_or_not_vec = rowSums(clip_or_not)
gene_var_norm_clip[which(gene_var_norm > clip_max)] <- clip_max

gene_var_norm_mean_clip = rowMeans(gene_var_norm_clip)
gene_var_norm_var_clip = rowVars(gene_var_norm_clip)
gene_hyper_var_clip <- rowSums(gene_var_norm_clip^2)/(ncol(data_filter) - 1)

result <- data.frame(pos=pos_list, clip=clip_or_not_vec, mean_orig=gene_mean_all, mean=gene_mean, var=gene_var, norm_mean=gene_var_norm_mean, norm_var=gene_var_norm_var,
	var_expect=gene_var_expect, hypervar=gene_hyper_var, norm_mean_clip=gene_var_norm_mean_clip, norm_var_clip=gene_var_norm_var_clip, hypervar_clip=gene_hyper_var_clip)

wgc_qnorm_MAV_stat_filename = paste0(wgc_qnorm_MAV_stat_filename_prefix, ".feather")
wgc_qnorm_MAV_stat_dir_filename = file.path(wgc_dir, wgc_qnorm_MAV_stat_filename)
# write.table(result, wgc_qnorm_MAV_stat_dir_filename, row.names=FALSE, col.names=TRUE, sep='\t')
# write_parquet(result, wgc_qnorm_MAV_stat_dir_filename)
write_feather(result, wgc_qnorm_MAV_stat_dir_filename)

p1 <- ggplot(result, aes(log2(mean), log2(var))) + geom_point(alpha = 1/20) + 
				geom_point(data=result,aes(log2(mean),log2(var_expect)), color="red", size=0.1) +
				theme_bw() +
				theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

p2 <- ggplot(result, aes(log2(mean), hypervar)) + geom_point(alpha = 1/20) +
				theme_bw() +
				theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))

combined_plot <- plot_grid(p1, p2, labels = c('A', 'B'))
print(combined_plot)

fig_filename = paste0(wgc_qnorm_MAV_stat_filename_prefix, ".png")
fig_dir_filename = file.path(wgc_dir, fig_filename)
ggsave(fig_dir_filename, plot=combined_plot)

set.seed(seed)
nrow_sample = min(nrow(result), nrow_sample)
result_sample = result[sample(nrow(result), nrow_sample), ]

p3 <- ggplot(result_sample, aes(log2(mean), log2(var))) +
				geom_pointdensity() + scale_color_viridis_c() + 
				geom_point(data=result,aes(log2(mean),log2(var_expect)), color="red", size=0.1) +
				theme_bw() +
				theme(axis.text = element_text(size=font_size), axis.title = element_text(size=font_size))
				print(p3)

density_fig_filename = paste0(wgc_qnorm_MAV_stat_filename_prefix, "_smooth-", format(nrow_sample, scientific=FALSE), "_seed-", format(seed, scientific=FALSE), ".png")
density_fig_dir_filename = file.path(wgc_dir, density_fig_filename)
ggsave(density_fig_dir_filename, plot=p3)