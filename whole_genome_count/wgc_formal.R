start_time <- Sys.time()

list.of.packages <- c("data.table", "arrow") # libraries from CRAN
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

listOfBiocPackages = c("GenomicAlignments", 
	"GenomicRanges", 
	"Biostrings", 
	"BSgenome.Hsapiens.NCBI.GRCh38",
	"plyranges") # libraries from bioconductor
notInstalled <- which(!listOfBiocPackages %in% rownames(installed.packages()))

if( length(notInstalled) ) {
	BiocManager::install(listOfBiocPackages[notInstalled])
}

suppressPackageStartupMessages({
	library(R.utils)
	library(GenomicAlignments)
	library(GenomicRanges)
	library(Biostrings)
	library(BSgenome.Hsapiens.NCBI.GRCh38)
	library(plyranges)
	library(arrow)
	library(preprocessCore)
	library(tibble)
})
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")
target_pairs_remained = filter_target_pairs(percentage_cutoff = 0.25)

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

rownames2first_col <- function(df, colname) {
	rownames_df = as.data.frame(rownames(df))
	colnames(rownames_df) = colname
	df = cbind(rownames_df, df)
	rownames(df) = NULL
	return (df)
}

args <- commandArgs(trailingOnly = TRUE)
scen <- args[1]
frag_type <- args[2]
BINSIZE <- strtoi(args[3])
target_qc_type <- args[4]
libnorm_type <- args[5]
save_dir <- args[6]
align_dir <- args[7]

# scen <- "V"
# frag_type <- "mixed"
# BINSIZE <- 800
# target_qc_type <- "hm-only"

# save_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800"
# align_dir <- "/dcs05/hongkai/data/next_cutntag/bulk/homotone_heterotone_merged/data_align"

if (libnorm_type == "libnorm-mean") {
	norm_scale_factor = 51633
} else if (libnorm_type == "libnorm-median") {
	norm_scale_factor = 5704
} else if (libnorm_type == "libnorm") {
	norm_scale_factor = 1E6
}

pos_colname = "pos"

target_list = c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
				"H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP",
				"EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1",
				"EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "POLR2AphosphoS2", "cJun", "cFos", "Max", "Myc",
				"USF1", "USF2", "NRF1", "YY1", "H3K9ac")

# hm_list = c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2", "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac")
# hm_writer_list = c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3",
# 					"H3K9me2", "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac", "SETD2",
# 					"MLL4_MLL2_KMT2B", "CBP_CREBBP", "EP300", "MSK1", "MSK2", "PIM1", "CDK8",
# 					"AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1", "EHMT1", "EZH2", "MLL1_KMT2A",
# 					"CTCF", "POLR2AphosphoS2")

# target_list = c('cJun', 'H3K9ac', 'H3K27me3', 'USF2', 'YY1', 'H3K4me3', 'cFos', 'POLR2AphosphoS2', 'NRF1',
# 				'H3K27ac', 'H3K14ac', 'Myc', 'H3K9me3', 'USF1', 'Max', 'CTCF')

# hm_list = c('H3K9ac', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K14ac', 'H3K9me3')
# hm_writer_list = c('H3K9ac', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K14ac', 'H3K9me3', 'CTCF', 'POLR2AphosphoS2')

if (target_qc_type == "all") {
	tags = target_pair_generation(target_list)
} else if (target_qc_type == "hm-only") {
	tags = target_pair_generation(hm_list)
} else if (target_qc_type == "hm-writer") {
	tags = target_pair_generation(hm_writer_list)
} else if (target_qc_type == "all-qc") {
	tags = intersect(target_pair_generation(target_list), target_pairs_remained)
}

bamFiles = c()

for (tag_i in tags) {
	target_pair_filename = paste0(tag_i, ".bam")
	target_pair_dir_filename = file.path(align_dir, scen, target_pair_filename)
	bamFiles = c(bamFiles, target_pair_dir_filename)
}

# dir.create(save_dir, recursive = TRUE)

refGenome <- BSgenome.Hsapiens.NCBI.GRCh38
chrSizes <- seqlengths(refGenome)[1:24]

chr_list = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
# chr_list = c("22", "Y")
chr_num = length(chr_list)
binChriDataframe_list = list()
# binChriDataframe_libnorm_list = list()
# binChriDataframe_libnorm_log2_list = list()
for (chr_idx in 1:chr_num) {
	chr_i = chr_list[chr_idx]
	print(paste0("compute at chr ", chr_i))
	
	chrSizei <- chrSizes[chr_i]
	print(chrSizei)
	bin <- tileGenome(chrSizei, tilewidth=BINSIZE, cut.last.tile.in.chrom=T)
	binChriDataframe <- as.data.frame(bin)
	binChriDataframe <- binChriDataframe[c("start", "end")]
	seqlevels(bin) <- paste0("chr", seqlevels(bin))
	
	grpName <- names(chrSizei)
	# dir.create(paste0(save_dir, grpName), recursive = TRUE)
	
	empty_target_pair_list = c()
	for (k in seq_along(bamFiles)) {
		bamFile <- bamFiles[k]
		temp <- readGAlignmentPairs(bamFile)
		locus <- data.frame(first_start = start(temp@first), 
							first_end = end(temp@first),
							last_start = start(temp@last), 
							last_end = end(temp@last))
		strand <- '*'

		if (dim(locus)[1] == 0) {
			empty_bamfile_target_pair = sub('\\.bam$', '', basename(bamFile))
			empty_target_pair_list = c(empty_target_pair_list, empty_bamfile_target_pair)
			bamContent <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
			print(paste0(empty_bamfile_target_pair, " is empty at chr", chr_i, "!"))
			# make empty grange
		} else {
			start <- rowMin(as.matrix(locus))
			end <- rowMax(as.matrix(locus))
			mid = ceiling((start + end) / 2)
			seqnames <- as.vector(seqnames(temp))
			bamContent <- makeGRangesFromDataFrame(data.frame(seqnames = seqnames, strand = strand, start = mid, end = mid))
			# bamContent %>% mutate(orig_start = start, orig_end = end)
		}

		overlapCount <- as.data.frame(GenomicRanges::countOverlaps(bin, bamContent))

		# if (exists("bamContent")) {
		# 	overlapCount <- as.data.frame(GenomicRanges::countOverlaps(bin, bamContent))
		# } else {
		# 	overlapCount = as.data.frame(bin)
		# 	overlapCount[, ] = 0
		# }
		
		# # Thresholding
		# overlapCount[overlapCount < THRESHOLD] <- 0
		bamName <- strsplit(tail(strsplit(bamFile, "/")[[1]], n=1), "[.]")[[1]][1]
		# print(overlapCount)
		colnames(overlapCount) <- bamName
		binChriDataframe <- cbind(binChriDataframe, overlapCount)
		# print("here")
		if (k %% 10 == 0) {
			print(object.size(binChriDataframe))
		}
	}

	# if (length(empty_target_pair_list) == 0) {
	# 	print(paste0("All target pairs at chr", chr_i, " are NOT empty!"))
	# } else {
	# 	for (empty_target_pair in empty_target_pair_list) {
	# 		overlapCount[, empty_target_pair] = 0
	# 	}
	# 	print(paste0("The following target pairs at chr", chr_i, "are empty!"))
	# 	print(empty_target_pair_list)
	# }
	print("overlap computation finished")
	
	# remove 0 ranges
	# binChriDataframe <- binChriDataframe[rowSums(binChriDataframe[, colnames(binChriDataframe)[4: length(colnames(binChriDataframe))]]) > 0, ]
	chr_df = data.frame(first_column=rep(c(paste0("chr", grpName)), times=dim(binChriDataframe)[1]))
	colnames(chr_df) = "CHR"
	binChriDataframe = cbind(chr_df, binChriDataframe)
	tmp_pos = binChriDataframe[, c("CHR", "start", "end")]
	tmp_pos$pos = paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
	pos_df = data.frame("pos"=tmp_pos$pos)
	tmp_wgc = binChriDataframe[ , -which(names(binChriDataframe) %in% c("CHR", "start", "end"))]
	# tmp_libnorm = sweep(tmp_wgc, 2, colSums(tmp_wgc), FUN = '/') * norm_scale_factor
	# tmp_libnorm_log2 = log2(tmp_libnorm + 1)

	binChriDataframe = cbind(pos_df, tmp_wgc)
	# binChriDataframe_libnorm = cbind(pos_df, tmp_libnorm)
	# binChriDataframe_libnorm_log2 = cbind(pos_df, tmp_libnorm_log2)

	binChriDataframe_list[[chr_idx]] = binChriDataframe
	# binChriDataframe_libnorm_list[[chr_idx]] = binChriDataframe_libnorm
	# binChriDataframe_libnorm_log2_list[[chr_idx]] = binChriDataframe_libnorm_log2
	
	# datasetName <- paste0((rangesBin[j] + 1), "-", rangesBin[j + 1])
}

binChriDataframe_full = as.data.frame(do.call(rbind, binChriDataframe_list))
# binChriDataframe_libnorm_full = as.data.frame(do.call(rbind, binChriDataframe_libnorm_list))
# binChriDataframe_libnorm_log2_full = as.data.frame(do.call(rbind, binChriDataframe_libnorm_log2_list))

# if (frag_type == "mixed") {
# 	datasetName_full <- paste0(scen, "_", frag_type, "_", BINSIZE, "_colQC-", target_qc_type)
# } else {
# 	datasetName_full <- paste0(scen, "_", BINSIZE, "_colQC-", target_qc_type)
# }
datasetName_full <- paste0(scen, "_", BINSIZE, "_colQC-", target_qc_type)

print(warnings())
print(datasetName_full)

preprocess_time <- Sys.time()
preprocess_time_taken <- round(preprocess_time - start_time, 2)
print(c("prepocess time taken: ", preprocess_time_taken))

datasetName_full_filename = paste0(datasetName_full, "_orig.feather")
datasetName_full_dir_filename = file.path(save_dir, datasetName_full_filename)
write_feather(binChriDataframe_full, datasetName_full_dir_filename)

saving_time_0 <- Sys.time()
saving_time_taken_0 <- round(saving_time_0 - preprocess_time, 2)
print(c("saving time original taken: ", saving_time_taken_0))

binChriDataframe_full = column_to_rownames(binChriDataframe_full, var=pos_colname)
col2idx_time <- Sys.time()
col2idx_time_taken <- round(col2idx_time - saving_time_0, 2)
print(c("col2idx time taken: ", col2idx_time_taken))

wgc_libnorm = sweep(binChriDataframe_full, 2, colSums(binChriDataframe_full), FUN = '/') * norm_scale_factor
libnorm_time <- Sys.time()
libnorm_time_taken <- round(libnorm_time - col2idx_time, 2)
print(c("libnorm time taken: ", libnorm_time_taken))

wgc_libnorm_log2 = log2(wgc_libnorm + 1)
libnorm_log2_time <- Sys.time()
libnorm_log2_time_taken <- round(libnorm_log2_time - libnorm_time, 2)
print(c("libnorm log2 time taken: ", libnorm_log2_time_taken))

wgc_libnorm_log2_minmaxnorm = sweep(sweep(wgc_libnorm_log2, 2, apply(wgc_libnorm_log2, 2, min), FUN = '-'), 2, apply(wgc_libnorm_log2, 2, function(x) max(x) - min(x)), FUN = '/')
libnorm_log2_minmaxnorm_time <- Sys.time()
libnorm_log2_minmaxnorm_time_taken <- round(libnorm_log2_minmaxnorm_time - libnorm_log2_time, 2)
print(c("libnorm log2 minmaxnorm time taken: ", libnorm_log2_minmaxnorm_time_taken))

wgc_libnorm_sqrt = sqrt(wgc_libnorm)
libnorm_sqrt_time <- Sys.time()
libnorm_sqrt_time_taken <- round(libnorm_sqrt_time - libnorm_log2_minmaxnorm_time, 2)
print(c("libnorm sqrt time taken: ", libnorm_sqrt_time_taken))

wgc_libnorm_sqrt_minmaxnorm = sweep(sweep(wgc_libnorm_sqrt, 2, apply(wgc_libnorm_sqrt, 2, min), FUN = '-'), 2, apply(wgc_libnorm_sqrt, 2, function(x) max(x) - min(x)), FUN = '/')
libnorm_sqrt_minmaxnorm_time <- Sys.time()
libnorm_sqrt_minmaxnorm_time_taken <- round(libnorm_sqrt_minmaxnorm_time - libnorm_sqrt_time, 2)
print(c("libnorm sqrt minmaxnorm time taken: ", libnorm_sqrt_minmaxnorm_time_taken))

# wgc_libnorm_sqrt_standard_col = scale(wgc_libnorm_sqrt, center=TRUE, scale=TRUE)
# libnorm_sqrt_standard_col_time <- Sys.time()
# libnorm_sqrt_standard_col_time_taken <- round(libnorm_sqrt_standard_col_time - libnorm_sqrt_minmaxnorm_time, 2)
# print(c("libnorm sqrt standard col time taken: ", libnorm_sqrt_standard_col_time_taken))

# wgc_libnorm_log2_standard_col = scale(wgc_libnorm_log2, center=TRUE, scale=TRUE)
# libnorm_log2_standard_col_time <- Sys.time()
# libnorm_log2_standard_col_time_taken <- round(libnorm_log2_standard_col_time - libnorm_sqrt_standard_col_time, 2)
# print(c("libnorm log2 standard col time taken: ", libnorm_log2_standard_col_time_taken))

# wgc_libnorm_sqrt_standard_row = t(scale(t(wgc_libnorm_sqrt), center=TRUE, scale=TRUE))
# libnorm_sqrt_standard_row_time <- Sys.time()
# libnorm_sqrt_standard_row_time_taken <- round(libnorm_sqrt_standard_row_time - libnorm_log2_standard_col_time, 2)
# print(c("libnorm sqrt standard row time taken: ", libnorm_sqrt_standard_row_time_taken))

# wgc_libnorm_log2_standard_row = t(scale(t(wgc_libnorm_log2), center=TRUE, scale=TRUE))
# libnorm_log2_standard_row_time <- Sys.time()
# libnorm_log2_standard_row_time_taken <- round(libnorm_log2_standard_row_time - libnorm_sqrt_standard_row_time, 2)
# print(c("libnorm log2 standard row time taken: ", libnorm_log2_standard_row_time_taken))

# # wgc_libnorm_log2_qnorm_matrix = normalize.quantiles(as.matrix(wgc_libnorm_log2), copy=TRUE)
# # wgc_libnorm_log2_qnorm = as.data.frame(wgc_libnorm_log2_qnorm_matrix)
# # rownames(wgc_libnorm_log2_qnorm) = rownames(wgc_libnorm_log2)
# # colnames(wgc_libnorm_log2_qnorm) = colnames(wgc_libnorm_log2)

wgc_nonzero = binChriDataframe_full[rowSums(binChriDataframe_full) != 0, ]
# nonzero_time <- Sys.time()
# nonzero_time_taken <- round(nonzero_time - libnorm_log2_standard_row_time, 2)
# print(c("nonzero time taken: ", nonzero_time_taken))

wgc_nonzero_libnorm = sweep(wgc_nonzero, 2, colSums(wgc_nonzero), FUN = '/') * norm_scale_factor
# nonzero_libnorm_time <- Sys.time()
# nonzero_libnorm_time_taken <- round(nonzero_libnorm_time - nonzero_time, 2)
# print(c("nonzero libnorm time taken: ", nonzero_libnorm_time_taken))

wgc_nonzero_libnorm_log2 = log2(wgc_nonzero_libnorm + 1)
# nonzero_libnorm_log2_time <- Sys.time()
# nonzero_libnorm_log2_time_taken <- round(nonzero_libnorm_log2_time - nonzero_libnorm_time, 2)
# print(c("nonzero libnorm log2 time taken: ", nonzero_libnorm_log2_time_taken))

wgc_nonzero_libnorm_log2_minmaxnorm = sweep(sweep(wgc_nonzero_libnorm_log2, 2, apply(wgc_nonzero_libnorm_log2, 2, min), FUN = '-'), 2, apply(wgc_nonzero_libnorm_log2, 2, function(x) max(x) - min(x)), FUN = '/')
# nonzero_libnorm_log2_minmaxnorm_time <- Sys.time()
# nonzero_libnorm_log2_minmaxnorm_time_taken <- round(nonzero_libnorm_log2_minmaxnorm_time - nonzero_libnorm_log2_time, 2)
# print(c("nonzero libnorm log2 minmaxnorm time taken: ", nonzero_libnorm_log2_minmaxnorm_time_taken))

wgc_nonzero_libnorm_sqrt = sqrt(wgc_nonzero_libnorm)
# nonzero_libnorm_sqrt_time <- Sys.time()
# nonzero_libnorm_sqrt_time_taken <- round(nonzero_libnorm_sqrt_time - nonzero_libnorm_log2_minmaxnorm_time, 2)
# print(c("nonzero libnorm sqrt time taken: ", nonzero_libnorm_sqrt_time_taken))

wgc_nonzero_libnorm_sqrt_minmaxnorm = sweep(sweep(wgc_nonzero_libnorm_sqrt, 2, apply(wgc_nonzero_libnorm_sqrt, 2, min), FUN = '-'), 2, apply(wgc_nonzero_libnorm_sqrt, 2, function(x) max(x) - min(x)), FUN = '/')
# nonzero_libnorm_sqrt_minmaxnorm_time <- Sys.time()
# nonzero_libnorm_sqrt_minmaxnorm_time_taken <- round(nonzero_libnorm_sqrt_minmaxnorm_time - nonzero_libnorm_sqrt_time, 2)
# print(c("nonzero libnorm sqrt minmaxnorm time taken: ", nonzero_libnorm_sqrt_minmaxnorm_time_taken))

# wgc_nonzero_libnorm_sqrt_standard_col = scale(wgc_nonzero_libnorm_sqrt, center=TRUE, scale=TRUE)
# nonzero_libnorm_sqrt_standard_col_time <- Sys.time()
# nonzero_libnorm_sqrt_standard_col_time_taken <- round(nonzero_libnorm_sqrt_standard_col_time - nonzero_libnorm_sqrt_minmaxnorm_time, 2)
# print(c("nonzero libnorm sqrt standard col time taken: ", nonzero_libnorm_sqrt_standard_col_time_taken))

# wgc_nonzero_libnorm_log2_standard_col = scale(wgc_nonzero_libnorm_log2, center=TRUE, scale=TRUE)
# nonzero_libnorm_log2_standard_col_time <- Sys.time()
# nonzero_libnorm_log2_standard_col_time_taken <- round(nonzero_libnorm_log2_standard_col_time - nonzero_libnorm_sqrt_standard_col_time, 2)
# print(c("nonzero libnorm log2 standard col time taken: ", nonzero_libnorm_log2_standard_col_time_taken))

# wgc_nonzero_libnorm_sqrt_standard_row = t(scale(t(wgc_nonzero_libnorm_sqrt), center=TRUE, scale=TRUE))
# nonzero_libnorm_sqrt_standard_row_time <- Sys.time()
# nonzero_libnorm_sqrt_standard_row_time_taken <- round(nonzero_libnorm_sqrt_standard_row_time - nonzero_libnorm_log2_standard_col_time, 2)
# print(c("nonzero libnorm sqrt standard row time taken: ", nonzero_libnorm_sqrt_standard_row_time_taken))

# wgc_nonzero_libnorm_log2_standard_row = t(scale(t(wgc_nonzero_libnorm_log2), center=TRUE, scale=TRUE))
# nonzero_libnorm_log2_standard_row_time <- Sys.time()
# nonzero_libnorm_log2_standard_row_time_taken <- round(nonzero_libnorm_log2_standard_row_time - nonzero_libnorm_sqrt_standard_row_time, 2)
# print(c("nonzero libnorm log2 standard row time taken: ", nonzero_libnorm_log2_standard_row_time_taken))
# # wgc_nonzero_libnorm_log2_qnorm_matrix = normalize.quantiles(as.matrix(wgc_nonzero_libnorm_log2), copy=TRUE)
# # wgc_nonzero_libnorm_log2_qnorm = as.data.frame(wgc_nonzero_libnorm_log2_qnorm_matrix)
# # rownames(wgc_nonzero_libnorm_log2_qnorm) = rownames(wgc_nonzero_libnorm_log2)
# # colnames(wgc_nonzero_libnorm_log2_qnorm) = colnames(wgc_nonzero_libnorm_log2)

wgc_libnorm = rownames_to_column(as.data.frame(wgc_libnorm), var=pos_colname)
wgc_libnorm_log2 = rownames_to_column(as.data.frame(wgc_libnorm_log2), var=pos_colname)
wgc_libnorm_log2_minmaxnorm = rownames_to_column(as.data.frame(wgc_libnorm_log2_minmaxnorm), var=pos_colname)
wgc_libnorm_sqrt = rownames_to_column(as.data.frame(wgc_libnorm_sqrt), var=pos_colname)
wgc_libnorm_sqrt_minmaxnorm = rownames_to_column(as.data.frame(wgc_libnorm_sqrt_minmaxnorm), var=pos_colname)
# idx2col_time_0 <- Sys.time()
# idx2col_time_0_taken <- round(idx2col_time_0 - nonzero_libnorm_log2_standard_row_time, 2)
# print(c("idx2col time taken: ", idx2col_time_0_taken))

# wgc_libnorm_sqrt_standard_col = rownames_to_column(as.data.frame(wgc_libnorm_sqrt_standard_col), var=pos_colname)
# wgc_libnorm_log2_standard_col = rownames_to_column(as.data.frame(wgc_libnorm_log2_standard_col), var=pos_colname)
# wgc_libnorm_sqrt_standard_row = rownames_to_column(as.data.frame(wgc_libnorm_sqrt_standard_row), var=pos_colname)
# wgc_libnorm_log2_standard_row = rownames_to_column(as.data.frame(wgc_libnorm_log2_standard_row), var=pos_colname)
# idx2col_time_1 <- Sys.time()
# idx2col_time_1_taken <- round(idx2col_time_1 - idx2col_time_0, 2)
# print(c("idx2col time taken: ", idx2col_time_1_taken))

wgc_nonzero = rownames_to_column(as.data.frame(wgc_nonzero), var=pos_colname)
wgc_nonzero_libnorm = rownames_to_column(as.data.frame(wgc_nonzero_libnorm), var=pos_colname)
wgc_nonzero_libnorm_log2 = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_log2), var=pos_colname)
wgc_nonzero_libnorm_log2_minmaxnorm = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_log2_minmaxnorm), var=pos_colname)
wgc_nonzero_libnorm_sqrt = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_sqrt), var=pos_colname)
# idx2col_time_2 <- Sys.time()
# idx2col_time_2_taken <- round(idx2col_time_2 - idx2col_time_1, 2)
# print(c("idx2col time taken: ", idx2col_time_2_taken))

wgc_nonzero_libnorm_sqrt_minmaxnorm = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_sqrt_minmaxnorm), var=pos_colname)
# wgc_nonzero_libnorm_sqrt_standard_col = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_sqrt_standard_col), var=pos_colname)
# wgc_nonzero_libnorm_log2_standard_col = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_log2_standard_col), var=pos_colname)
# wgc_nonzero_libnorm_sqrt_standard_row = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_sqrt_standard_row), var=pos_colname)
# wgc_nonzero_libnorm_log2_standard_row = rownames_to_column(as.data.frame(wgc_nonzero_libnorm_log2_standard_row), var=pos_colname)
# idx2col_time_3 <- Sys.time()
# idx2col_time_3_taken <- round(idx2col_time_3 - idx2col_time_2, 2)
# print(c("idx2col time taken: ", idx2col_time_3_taken))

# wgc_libnorm_log2 = rownames2first_col(wgc_libnorm_log2, pos_colname)
# wgc_libnorm_log2_minmaxnorm = rownames2first_col(wgc_libnorm_log2_minmaxnorm, pos_colname)
# wgc_libnorm_sqrt = rownames2first_col(wgc_libnorm_sqrt, pos_colname)
# wgc_libnorm_sqrt_minmaxnorm = rownames2first_col(wgc_libnorm_sqrt_minmaxnorm, pos_colname)
# wgc_libnorm_sqrt_standard_col = rownames2first_col(wgc_libnorm_sqrt_standard_col, pos_colname)
# wgc_libnorm_log2_standard_col = rownames2first_col(wgc_libnorm_log2_standard_col, pos_colname)
# wgc_libnorm_sqrt_standard_row = rownames2first_col(wgc_libnorm_sqrt_standard_row, pos_colname)
# wgc_libnorm_log2_standard_row = rownames2first_col(wgc_libnorm_log2_standard_row, pos_colname)
# # wgc_libnorm_log2_qnorm = rownames2first_col(wgc_libnorm_log2_qnorm, pos_colname)
# wgc_nonzero = rownames2first_col(wgc_nonzero, pos_colname)
# wgc_nonzero_libnorm = rownames2first_col(wgc_nonzero_libnorm, pos_colname)
# wgc_nonzero_libnorm_log2 = rownames2first_col(wgc_nonzero_libnorm_log2, pos_colname)
# wgc_nonzero_libnorm_log2_minmaxnorm = rownames2first_col(wgc_nonzero_libnorm_log2_minmaxnorm, pos_colname)
# wgc_nonzero_libnorm_sqrt = rownames2first_col(wgc_nonzero_libnorm_sqrt, pos_colname)
# wgc_nonzero_libnorm_sqrt_minmaxnorm = rownames2first_col(wgc_nonzero_libnorm_sqrt_minmaxnorm, pos_colname)
# wgc_nonzero_libnorm_sqrt_standard_col = rownames2first_col(wgc_nonzero_libnorm_sqrt_standard_col, pos_colname)
# wgc_nonzero_libnorm_log2_standard_col = rownames2first_col(wgc_nonzero_libnorm_log2_standard_col, pos_colname)
# wgc_nonzero_libnorm_sqrt_standard_row = rownames2first_col(wgc_nonzero_libnorm_sqrt_standard_row, pos_colname)
# wgc_nonzero_libnorm_log2_standard_row = rownames2first_col(wgc_nonzero_libnorm_log2_standard_row, pos_colname)
# # wgc_nonzero_libnorm_log2_qnorm = rownames2first_col(wgc_nonzero_libnorm_log2_qnorm, pos_colname)

wgc_libnorm_filename = paste0(datasetName_full, "_", libnorm_type, ".feather")
wgc_libnorm_dir_filename = file.path(save_dir, wgc_libnorm_filename)
write_feather(wgc_libnorm, wgc_libnorm_dir_filename)

wgc_libnorm_log2_filename = paste0(datasetName_full, "_", libnorm_type, "_log2.feather")
wgc_libnorm_log2_dir_filename = file.path(save_dir, wgc_libnorm_log2_filename)
write_feather(wgc_libnorm_log2, wgc_libnorm_log2_dir_filename)

wgc_libnorm_log2_minmaxnorm_filename = paste0(datasetName_full, "_", libnorm_type, "_log2_minmax.feather")
wgc_libnorm_log2_minmaxnorm_dir_filename = file.path(save_dir, wgc_libnorm_log2_minmaxnorm_filename)
write_feather(wgc_libnorm_log2_minmaxnorm, wgc_libnorm_log2_minmaxnorm_dir_filename)

wgc_libnorm_sqrt_filename = paste0(datasetName_full, "_", libnorm_type, "_sqrt.feather")
wgc_libnorm_sqrt_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_filename)
write_feather(wgc_libnorm_sqrt, wgc_libnorm_sqrt_dir_filename)

wgc_libnorm_sqrt_minmaxnorm_filename = paste0(datasetName_full, "_", libnorm_type, "_sqrt_minmax.feather")
wgc_libnorm_sqrt_minmaxnorm_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_minmaxnorm_filename)
write_feather(wgc_libnorm_sqrt_minmaxnorm, wgc_libnorm_sqrt_minmaxnorm_dir_filename)

# wgc_libnorm_sqrt_standard_col_filename = paste0(datasetName_full, "_", libnorm_type, "_sqrt_standard_col.feather")
# wgc_libnorm_sqrt_standard_col_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_standard_col_filename)
# write_feather(wgc_libnorm_sqrt_standard_col, wgc_libnorm_sqrt_standard_col_dir_filename)

# wgc_libnorm_log2_standard_col_filename = paste0(datasetName_full, "_", libnorm_type, "_log2_standard_col.feather")
# wgc_libnorm_log2_standard_col_dir_filename = file.path(save_dir, wgc_libnorm_log2_standard_col_filename)
# write_feather(wgc_libnorm_log2_standard_col, wgc_libnorm_log2_standard_col_dir_filename)

# wgc_libnorm_sqrt_standard_row_filename = paste0(datasetName_full, "_", libnorm_type, "_sqrt_standard_row.feather")
# wgc_libnorm_sqrt_standard_row_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_standard_row_filename)
# write_feather(wgc_libnorm_sqrt_standard_row, wgc_libnorm_sqrt_standard_row_dir_filename)

# wgc_libnorm_log2_standard_row_filename = paste0(datasetName_full, "_", libnorm_type, "_log2_standard_row.feather")
# wgc_libnorm_log2_standard_row_dir_filename = file.path(save_dir, wgc_libnorm_log2_standard_row_filename)
# write_feather(wgc_libnorm_log2_standard_row, wgc_libnorm_log2_standard_row_dir_filename)

# saving_time_1 <- Sys.time()
# saving_time_taken_1 <- round(saving_time_1 - idx2col_time_3, 2)
# print(c("saving time taken: ", saving_time_taken_1))

# wgc_libnorm_log2_qnorm_filename = paste0(datasetName_full, "_libnorm_log2_qnorm.feather")
# wgc_libnorm_log2_qnorm_dir_filename = file.path(save_dir, wgc_libnorm_log2_qnorm_filename)
# write_feather(wgc_libnorm_log2_qnorm, wgc_libnorm_log2_qnorm_dir_filename)

wgc_nonzero_filename = paste0(datasetName_full, "_noAllZero.feather")
wgc_nonzero_dir_filename = file.path(save_dir, wgc_nonzero_filename)
write_feather(wgc_nonzero, wgc_nonzero_dir_filename)

wgc_nonzero_libnorm_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, ".feather")
wgc_nonzero_libnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_filename)
write_feather(wgc_nonzero_libnorm, wgc_nonzero_libnorm_dir_filename)

wgc_nonzero_libnorm_log2_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_log2.feather")
wgc_nonzero_libnorm_log2_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_filename)
write_feather(wgc_nonzero_libnorm_log2, wgc_nonzero_libnorm_log2_dir_filename)

wgc_nonzero_libnorm_log2_minmaxnorm_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_log2_minmax.feather")
wgc_nonzero_libnorm_log2_minmaxnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_minmaxnorm_filename)
write_feather(wgc_nonzero_libnorm_log2_minmaxnorm, wgc_nonzero_libnorm_log2_minmaxnorm_dir_filename)

wgc_nonzero_libnorm_sqrt_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_sqrt.feather")
wgc_nonzero_libnorm_sqrt_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_filename)
write_feather(wgc_nonzero_libnorm_sqrt, wgc_nonzero_libnorm_sqrt_dir_filename)

wgc_nonzero_libnorm_sqrt_minmaxnorm_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_sqrt_minmax.feather")
wgc_nonzero_libnorm_sqrt_minmaxnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_minmaxnorm_filename)
write_feather(wgc_nonzero_libnorm_sqrt_minmaxnorm, wgc_nonzero_libnorm_sqrt_minmaxnorm_dir_filename)

# wgc_nonzero_libnorm_sqrt_standard_col_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_sqrt_standard_col.feather")
# wgc_nonzero_libnorm_sqrt_standard_col_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_standard_col_filename)
# write_feather(wgc_nonzero_libnorm_sqrt_standard_col, wgc_nonzero_libnorm_sqrt_standard_col_dir_filename)

# wgc_nonzero_libnorm_log2_standard_col_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_log2_standard_col.feather")
# wgc_nonzero_libnorm_log2_standard_col_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_standard_col_filename)
# write_feather(wgc_nonzero_libnorm_log2_standard_col, wgc_nonzero_libnorm_log2_standard_col_dir_filename)

# wgc_nonzero_libnorm_sqrt_standard_row_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_sqrt_standard_row.feather")
# wgc_nonzero_libnorm_sqrt_standard_row_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_standard_row_filename)
# write_feather(wgc_nonzero_libnorm_sqrt_standard_row, wgc_nonzero_libnorm_sqrt_standard_row_dir_filename)

# wgc_nonzero_libnorm_log2_standard_row_filename = paste0(datasetName_full, "_noAllZero_", libnorm_type, "_log2_standard_row.feather")
# wgc_nonzero_libnorm_log2_standard_row_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_standard_row_filename)
# write_feather(wgc_nonzero_libnorm_log2_standard_row, wgc_nonzero_libnorm_log2_standard_row_dir_filename)

# saving_time_2 <- Sys.time()
# saving_time_taken_2 <- round(saving_time_2 - saving_time_1, 2)
# print(c("saving time taken: ", saving_time_taken_2))

# wgc_nonzero_libnorm_log2_qnorm_filename = paste0(datasetName_full, "_noAllZero_libnorm_log2_qnorm.feather")
# wgc_nonzero_libnorm_log2_qnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_qnorm_filename)
# write_feather(wgc_nonzero_libnorm_log2_qnorm, wgc_nonzero_libnorm_log2_qnorm_dir_filename)

print(paste0(scen, " ", target_qc_type, " ", "saved!"))
