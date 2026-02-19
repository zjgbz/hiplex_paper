# Input: ../data/frag_decomposition/{scen}-premerge_all_frag_lens.RData (fragment length data: bamWidths)
# Output: ../results/Figure2/{scen}_all_frag_hist_density.pdf (histogram with density overlay)

rm(list=ls())
# library(GenomicAlignments)
# library(mixtools)
# library(mclust)

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

packages <- c("ggplot2")
lapply(packages, install_if_missing)

suppressPackageStartupMessages({
  library(ggplot2)
})

scen_list = c("V")
# scen_list = c("T")
args <- commandArgs(trailingOnly = TRUE)
frag_dir <- args[1]
save_dir <- args[2]

for (scen in scen_list) {
	bam_overplot_filename = paste0(scen, "_all_frag_hist_density.pdf")
	bam_overplot_dir_filename = file.path(save_dir, bam_overplot_filename)

	bam_bundle_filename = paste0(scen, "-premerge_all_frag_lens.RData")
	bam_bundle_dir_filename = file.path(frag_dir, bam_bundle_filename)
	load(bam_bundle_dir_filename)
	bamWidths_df_tmp = data.frame("fragment_length"=bamWidths)

	dens_frq_plot = ggplot(bamWidths_df_tmp, aes(x = fragment_length)) + 
		geom_histogram(aes(y = ..density..), colour = 1, fill = "white", bins=60, linewidth=0.1) + 
		# geom_histogram(aes(y = ..count..), colour = 1, fill = "white", bins=60) + 
		# geom_density(colour = 4, fill = 4, alpha = 0.25, kernel="gaussian", linewidth=0.5, bw=25) +
		geom_density(colour = 4, fill = 4, alpha = 0.25, kernel="gaussian", linewidth=0.5) +
		ylab("Density") + xlab("Fragment Length") +
		theme(
			panel.background = element_rect(fill='transparent'),
			plot.background = element_rect(fill='transparent', color=NA),
			panel.grid.major = element_blank(),
			panel.grid.minor = element_blank(),
			legend.background = element_rect(fill='transparent'),
			legend.box.background = element_rect(fill='transparent'),
			axis.title.x = element_text(size = 24, face = "plain"),
			axis.title.y = element_text(size = 24, face = "plain")
			)

	ggsave(bam_overplot_dir_filename, plot=dens_frq_plot)
}