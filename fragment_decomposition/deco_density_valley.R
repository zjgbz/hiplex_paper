rm(list=ls())
library(sfsmisc)

derivative <- function(x, y) { # create a function with the name my_function
	x_0 = x[1:(length(x) - 1)]
	x_1 = x[2:(length(x))]
	y_0 = y[1:(length(y) - 1)]
	y_1 = y[2:(length(y))]
	dif_x = x_1 - x_0
	dif_y = y_1 - y_0
	out = list(x_dif=x[1:(length(x) - 1)], y_dif=dif_y / dif_x)
	return (out)
}

# scen_list = c("T", "V")
scen_list = c("all")
frag_dir = "/dcs05/hongkai/data/next_cutntag/bulk/frag_len"
save_dir = "/dcs05/hongkai/data/next_cutntag/bulk/frag_len"
dens_reso = 2^15
density_kernel = "gaussian"

c1 = rgb(173, 216, 230, max=255, names="lt.blue")
n_breaks = 200

lb_00 = 40
ub_00 = 60
hist_y_00_abs_min_x = 48.1065

lb_1 = 115
ub_1 = 135
hist_y_1_abs_min_x = 120.6945
# rough_x_1 = 128.2624
# rough_x_1 = 128.2097
# rough_x_1 = 126.2733
# rough_x_1 = 127.8093
# rough_x_1 = 127.8368
# rough_x_1 = 126.2516
rough_x_1 = 128.2625

lb_10 = 185
ub_10 = 205
hist_y_10_abs_min_x = 197.0955
hist_y_10_abs_min_x_mod = 209.7355

# lb_2 = 292
# ub_2 = 304
lb_2 = 280
ub_2 = 320
hist_y_2_abs_min_x = 297.2205
# rough_x_2 = 295.2996
# rough_x_2 = 295.2979
# rough_x_2 = 296.5215
# rough_x_2 = 295.1505
# rough_x_2 = 295.1761
# rough_x_2 = 296.5226
rough_x_2 = 295.2996

lb_20 = 365
ub_20 = 375
hist_y_20_abs_min_x = 371.3645

lb_3 = 420
ub_3 = 550
rough_x_3 = 462.3368
hist_x_3_max = 519.0855
# hist_x_3m_max = 473.22
hist_x_3m_max = 452.179
hist_x_3o_max = 445
hist_x_3d2_max = 444.6345

for (scen in scen_list) {
	# bam_bundle_filename = paste0(scen, "-premerge_all-qc_frag_lens.RData")
	# bam_bundle_filename = paste0(scen, "-premerge_all_frag_lens.RData")
    bam_bundle_filename = paste0(scen, "_all_frag_lens.RData") # if scen = "all" without any QC
	bam_bundle_dir_filename = file.path(frag_dir, bam_bundle_filename)
	load(bam_bundle_dir_filename)
	# bam_bundle_density = density(bamWidths, kernel=density_kernel, n=dens_reso, width=50)
	bam_bundle_density = density(bamWidths, kernel=density_kernel, n=dens_reso)
	# bam_bundle_density = density(bamWidths, kernel=density_kernel, n=dens_reso, bw=20)
	hist_x = bam_bundle_density$x
	hist_count = bam_bundle_density$y
	f1 = D1ss(x=hist_x, y=hist_count, spar.off=0.0)
	f2_raw <- D1D2(x=hist_x, y=hist_count, deriv=2, xout=hist_x)
	f2 = f2_raw$D2

	hist_x_00 = hist_x[hist_x < ub_00 & hist_x > lb_00]
	hist_count_00 = hist_count[hist_x < ub_00 & hist_x > lb_00]
	f1_00 = f1[hist_x < ub_00 & hist_x > lb_00]
	f2_00 = f2[hist_x < ub_00 & hist_x > lb_00]

	hist_x_1 = hist_x[hist_x < ub_1 & hist_x > lb_1]
	hist_count_1 = hist_count[hist_x < ub_1 & hist_x > lb_1]
	f1_1 = f1[hist_x < ub_1 & hist_x > lb_1]
	f2_1 = f2[hist_x < ub_1 & hist_x > lb_1]

	hist_x_10 = hist_x[hist_x < ub_10 & hist_x > lb_10]
	hist_count_10 = hist_count[hist_x < ub_10 & hist_x > lb_10]
	f1_10 = f1[hist_x < ub_10 & hist_x > lb_10]

	hist_x_2 = hist_x[hist_x < ub_2 & hist_x > lb_2]
	hist_count_2 = hist_count[hist_x < ub_2 & hist_x > lb_2]
	f1_2 = f1[hist_x < ub_2 & hist_x > lb_2]
	f2_2 = f2[hist_x < ub_2 & hist_x > lb_2]

	hist_x_20 = hist_x[hist_x < ub_20 & hist_x > lb_20]
	hist_count_20 = hist_count[hist_x < ub_20 & hist_x > lb_20]
	f1_20 = f1[hist_x < ub_20 & hist_x > lb_20]
	f2_20 = f2[hist_x < ub_20 & hist_x > lb_20]

	hist_x_3 = hist_x[hist_x < ub_3 & hist_x > lb_3]
	hist_count_3 = hist_count[hist_x < ub_3 & hist_x > lb_3]
	f1_3 = f1[hist_x < ub_3 & hist_x > lb_3]
	f2_3 = f2[hist_x < ub_3 & hist_x > lb_3]

	hist(bamWidths, breaks=n_breaks, col=c1, border=c1)
	par(new=TRUE)
	# # plot(hist_x, hist_count, pch=20, cex=0.1)
	# # par(new=TRUE)
	# # plot(hist_x, f1, pch=20, cex=0.1, col="red")
	# # par(new=TRUE)
	# # plot(hist_x, f2, pch=20, cex=0.1, col="blue")
	plot(hist_x, hist_count, pch=20, cex=0.1)
	par(new=TRUE)
	abline(v=rough_x_1, col="darkorchid1", lty=2)
	par(new=TRUE)
	abline(v=rough_x_2, col="darkorchid1", lty=2)
	# par(new=TRUE)
	# abline(v=rough_x_3, col="darkorchid1", lty=2)

	# plot(hist_x_1, hist_count_1, pch=20, cex=0.1)
	# rough_x_1 = hist_x_1[hist_count_1 == min(hist_count_1)]
	# par(new=TRUE)
	# abline(v=rough_x_1, col="darkorchid1", lty=2)
	# par(new=TRUE)
	# plot(hist_x_1, f1_1, pch=20, cex=0.1, col="red")
	# # par(new=TRUE)
	# # plot(hist_x_1, f2_1, pch=20, cex=0.1, col="blue")

	# plot(hist_x_2, hist_count_2, pch=20, cex=0.1)
	# rough_x_2 = hist_x_2[hist_count_2 == min(hist_count_2)]
	# par(new=TRUE)
	# abline(v=rough_x_2, col="darkorchid1", lty=2)

	# # plot(x_3, y_3, pch=20, cex=0.1)
	# # par(new=TRUE)
	# # plot(dif_x_3, dif_y_3, pch=20, cex=0.1, col="red")
	# # par(new=TRUE)
	# # plot(diff_x_3, diff_y_3, pch=20, cex=0.1, col="blue")
}
# png(bam_overplot_dir_filename)
# hist(bamWidths, breaks = 160)
# dev.off()

# ggplot(bamWidths_df, aes(fragment_length, fill = Scenario)) + 
# 	geom_histogram(alpha = 0.5, bins = 160, position = 'identity')

# ggplot(bamWidths_df, aes(fragment_length, fill = Scenario)) + geom_density(alpha = 0.2)