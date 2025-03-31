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
	# library(data.table)
})
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")

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
BINSIZE <- args[3]
target_qc_type <- args[4]

save_dir <- args[5]
align_dir <- args[6]

norm_scale_factor = 1E6
pos_colname = "pos"

# target_list = c('cJun', 'H3K9ac', 'H3K27me3', 'USF2', 'YY1', 'H3K4me3', 'cFos', 'POLR2AphosphoS2', 'NRF1',
# 				'H3K27ac', 'H3K14ac', 'Myc', 'H3K9me3', 'USF1', 'Max', 'CTCF')

# hm_list = c('H3K9ac', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K14ac', 'H3K9me3')
# hm_writer_list = c('H3K9ac', 'H3K27me3', 'H3K4me3', 'H3K27ac', 'H3K14ac', 'H3K9me3', 'CTCF', 'POLR2AphosphoS2')

target_list_dict = list("all"=c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
								"H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac", "SETD2", "MLL4_MLL2_KMT2B",
								"CBP_CREBBP", "EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2",
								"SuVar39_SUV39H1", "EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "POLR2AphosphoS2", "cJun",
								"cFos", "Max", "Myc", "USF1", "USF2", "NRF1", "YY1"),
						"all-qc"=c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2",
									"H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac", "SETD2", "MLL4_MLL2_KMT2B",
									"CBP_CREBBP", "EP300", "MSK1", "MSK2", "PIM1", "CDK8", "AURORA_Aurora_B", "EHMT2",
									"SuVar39_SUV39H1", "EHMT1", "EZH2", "MLL1_KMT2A", "CTCF", "POLR2AphosphoS2", "cJun",
									"cFos", "Max", "Myc", "USF1", "USF2", "NRF1", "YY1"),
						"writer-sptf"=list("writer"=c("SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP", "EP300", "MSK1", "MSK2", "PIM1", "CDK8",
														"AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1", "EHMT1", "EZH2", "MLL1_KMT2A"),
											"special_tf"=c("CTCF", "EZH2", "POLR2AphosphoS2")),
						"hm-only"=c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2", "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac"),
						"hm-writer"=list("hm"=c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2", "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac", "CTCF", "POLR2AphosphoS2"),
										"writer"=c("SETD2", "MLL4_MLL2_KMT2B", "CBP_CREBBP", "EP300", "MSK1", "MSK2", "PIM1", "CDK8",
													"AURORA_Aurora_B", "EHMT2", "SuVar39_SUV39H1", "EHMT1", "EZH2", "MLL1_KMT2A")),
						"hm-tf"=list("hm"=c("H3K36me3", "H3K4me1", "H3K27ac", "H3S10ph", "H2A_XS139ph", "H3K79me3", "H3K9me2", "H3K9me3", "H3K14ac", "H3K27me3", "H3K4me3", "H3K9ac", "CTCF", "POLR2AphosphoS2"),
									"tf"=c("cJun", "cFos", "Max", "Myc", "USF1", "USF2", "NRF1", "YY1"))
						)

target_list = target_list_dict[[target_qc_type]]
tags = target_pair_generation(target_list)
target_pairs_remained = filter_target_pairs(percentage_cutoff = 0.25)
if (grepl("-qc", target_qc_type)) {
	tags = intersect(tags, target_pairs_remained)
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

# promoter_dir = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/promoter_region"
promoter_dir = "/dcs05/hongkai/data/next_cutntag/bulk/RNA-seq/ccre_region"
promoter_filename = paste0(BINSIZE, ".tsv")
promoter_dir_filename = file.path(promoter_dir, promoter_filename)
promoter_df = read.table(promoter_dir_filename, header=TRUE, sep="\t", row.names=NULL)
print(c(BINSIZE, dim(promoter_df)))

# Extract promoter indicator and gene name and gene ID
promoter_info = promoter_df[, c("pos", "gene_id", "gene_name")]
bin_all_chr_df = promoter_df[, c("seqnames", "start", "end", "width", "strand")]
binChriDataframe_all_chr = promoter_df[, c("gene_id", "seqnames", "start", "end")]

chr_list = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X", "Y")
# chr_list = c("22", "Y")
chr_num = length(chr_list)
binChriDataframe_list = list()
# binChriDataframe_log2_list = list()
for (chr_idx in 1:chr_num) {
	chr_i = chr_list[chr_idx]

    bin_df = bin_all_chr_df[(bin_all_chr_df$seqnames) == chr_i, ]
    bin = makeGRangesFromDataFrame(bin_df)
	binChriDataframe = binChriDataframe_all_chr[(binChriDataframe_all_chr$seqnames) == chr_i, ]
	binChriDataframe <- binChriDataframe[c("gene_id", "start", "end")]
	seqlevels(bin) <- paste0("chr", seqlevels(bin))
	print(paste0("compute at chr ", chr_i, " with ", dim(binChriDataframe)[1], " bins"))
	
	grpName <- chr_i
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

	# chr_df = data.frame(first_column=rep(c(paste0("chr", grpName)), times=dim(binChriDataframe)[1]))
	# colnames(chr_df) = "CHR"
	# binChriDataframe = cbind(chr_df, binChriDataframe)
	# tmp_pos = binChriDataframe[, c("CHR", "start", "end")]
	# tmp_pos$pos = paste0(tmp_pos$CHR, "_", tmp_pos$start, "_", tmp_pos$end)
	# pos_df = data.frame("pos"=tmp_pos$pos)
	# tmp_wgc = binChriDataframe[ , -which(names(binChriDataframe) %in% c("CHR", "start", "end"))]
	# tmp_log2 = log2(tmp_wgc + 1)
	pos_df = data.frame("pos"=binChriDataframe$gene_id)
	tmp_wgc = binChriDataframe[ , -which(names(binChriDataframe) %in% c("gene_id", "start", "end"))]
	print(c(chr_i, dim(tmp_wgc)))
	# tmp_log2 = log2(tmp_wgc + 1)

	binChriDataframe = cbind(pos_df, tmp_wgc)
	# binChriDataframe_log2 = cbind(pos_df, tmp_log2)

	binChriDataframe_list[[chr_idx]] = binChriDataframe
	# binChriDataframe_log2_list[[chr_idx]] = binChriDataframe_log2
	
	# datasetName <- paste0((rangesBin[j] + 1), "-", rangesBin[j + 1])
}

binChriDataframe_full = as.data.frame(do.call(rbind, binChriDataframe_list))
print(dim(binChriDataframe_full))
# binChriDataframe_log2_full = as.data.frame(do.call(rbind, binChriDataframe_log2_list))

# if (frag_type == "mixed") {
# 	datasetName_full <- paste0(scen, "_", frag_type, "_", BINSIZE, "_colQC-", target_qc_type)
# } else {
# 	datasetName_full <- paste0(scen, "_", BINSIZE, "_colQC-", target_qc_type)
# }
datasetName_full <- paste0(scen, "_", BINSIZE, "_colQC-", target_qc_type)

print(warnings())
print(datasetName_full)
datasetName_full_filename = paste0(datasetName_full, "_orig.feather")
datasetName_full_dir_filename = file.path(save_dir, datasetName_full_filename)
write_feather(binChriDataframe_full, datasetName_full_dir_filename)

pos_list <- binChriDataframe_full[, "pos"]
binChriDataframe_full = binChriDataframe_full[, -1]
rownames(binChriDataframe_full) = pos_list

# for each column, divide each element by the maximum value of the column
wgc_libnorm = sweep(binChriDataframe_full, 2, colSums(binChriDataframe_full), FUN = '/') * norm_scale_factor
# wgc_libnorm_log2 = log2(wgc_libnorm + 1)
# wgc_libnorm_log2_minmaxnorm = sweep(sweep(wgc_libnorm_log2, 2, apply(wgc_libnorm_log2, 2, min), FUN = '-'), 2, apply(wgc_libnorm_log2, 2, function(x) max(x) - min(x)), FUN = '/')
# wgc_libnorm_sqrt = sqrt(wgc_libnorm)
# wgc_libnorm_sqrt_minmaxnorm = sweep(sweep(wgc_libnorm_sqrt, 2, apply(wgc_libnorm_sqrt, 2, min), FUN = '-'), 2, apply(wgc_libnorm_sqrt, 2, function(x) max(x) - min(x)), FUN = '/')

# wgc_nonzero = binChriDataframe_full[rowSums(binChriDataframe_full) != 0, ]
# wgc_nonzero_libnorm = sweep(wgc_nonzero, 2, colSums(wgc_nonzero), FUN = '/') * norm_scale_factor
# wgc_nonzero_libnorm_log2 = log2(wgc_nonzero_libnorm + 1)
# wgc_nonzero_libnorm_log2_minmaxnorm = sweep(sweep(wgc_nonzero_libnorm_log2, 2, apply(wgc_nonzero_libnorm_log2, 2, min), FUN = '-'), 2, apply(wgc_nonzero_libnorm_log2, 2, function(x) max(x) - min(x)), FUN = '/')
# wgc_nonzero_libnorm_sqrt = sqrt(wgc_nonzero_libnorm)
# wgc_nonzero_libnorm_sqrt_minmaxnorm = sweep(sweep(wgc_nonzero_libnorm_sqrt, 2, apply(wgc_nonzero_libnorm_sqrt, 2, min), FUN = '-'), 2, apply(wgc_nonzero_libnorm_sqrt, 2, function(x) max(x) - min(x)), FUN = '/')

wgc_libnorm = rownames2first_col(wgc_libnorm, pos_colname)
# wgc_libnorm_log2 = rownames2first_col(wgc_libnorm_log2, pos_colname)
# wgc_libnorm_log2_minmaxnorm = rownames2first_col(wgc_libnorm_log2_minmaxnorm, pos_colname)
# wgc_libnorm_sqrt = rownames2first_col(wgc_libnorm_sqrt, pos_colname)
# wgc_libnorm_sqrt_minmaxnorm = rownames2first_col(wgc_libnorm_sqrt_minmaxnorm, pos_colname)

# wgc_nonzero = rownames2first_col(wgc_nonzero, pos_colname)
# wgc_nonzero_libnorm = rownames2first_col(wgc_nonzero_libnorm, pos_colname)
# wgc_nonzero_libnorm_log2 = rownames2first_col(wgc_nonzero_libnorm_log2, pos_colname)
# wgc_nonzero_libnorm_log2_minmaxnorm = rownames2first_col(wgc_nonzero_libnorm_log2_minmaxnorm, pos_colname)
# wgc_nonzero_libnorm_sqrt = rownames2first_col(wgc_nonzero_libnorm_sqrt, pos_colname)
# wgc_nonzero_libnorm_sqrt_minmaxnorm = rownames2first_col(wgc_nonzero_libnorm_sqrt_minmaxnorm, pos_colname)

wgc_libnorm_filename = paste0(datasetName_full, "_libnorm.feather")
wgc_libnorm_dir_filename = file.path(save_dir, wgc_libnorm_filename)
write_feather(wgc_libnorm, wgc_libnorm_dir_filename)

# wgc_libnorm_log2_filename = paste0(datasetName_full, "_libnorm_log2.feather")
# wgc_libnorm_log2_dir_filename = file.path(save_dir, wgc_libnorm_log2_filename)
# write_feather(wgc_libnorm_log2, wgc_libnorm_log2_dir_filename)

# wgc_libnorm_log2_minmaxnorm_filename = paste0(datasetName_full, "_libnorm_log2_minmax.feather")
# wgc_libnorm_log2_minmaxnorm_dir_filename = file.path(save_dir, wgc_libnorm_log2_minmaxnorm_filename)
# write_feather(wgc_libnorm_log2_minmaxnorm, wgc_libnorm_log2_minmaxnorm_dir_filename)

# wgc_libnorm_sqrt_filename = paste0(datasetName_full, "_libnorm_sqrt.feather")
# wgc_libnorm_sqrt_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_filename)
# write_feather(wgc_libnorm_sqrt, wgc_libnorm_sqrt_dir_filename)

# wgc_libnorm_sqrt_minmaxnorm_filename = paste0(datasetName_full, "_libnorm_sqrt_minmax.feather")
# wgc_libnorm_sqrt_minmaxnorm_dir_filename = file.path(save_dir, wgc_libnorm_sqrt_minmaxnorm_filename)
# write_feather(wgc_libnorm_sqrt_minmaxnorm, wgc_libnorm_sqrt_minmaxnorm_dir_filename)

# wgc_nonzero_filename = paste0(datasetName_full, "_noAllZero.feather")
# wgc_nonzero_dir_filename = file.path(save_dir, wgc_nonzero_filename)
# write_feather(wgc_nonzero, wgc_nonzero_dir_filename)

# wgc_nonzero_libnorm_filename = paste0(datasetName_full, "_noAllZero_libnorm.feather")
# wgc_nonzero_libnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_filename)
# write_feather(wgc_nonzero_libnorm, wgc_nonzero_libnorm_dir_filename)

# wgc_nonzero_libnorm_log2_filename = paste0(datasetName_full, "_noAllZero_libnorm_log2.feather")
# wgc_nonzero_libnorm_log2_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_filename)
# write_feather(wgc_nonzero_libnorm_log2, wgc_nonzero_libnorm_log2_dir_filename)

# wgc_nonzero_libnorm_log2_minmaxnorm_filename = paste0(datasetName_full, "_noAllZero_libnorm_log2_minmax.feather")
# wgc_nonzero_libnorm_log2_minmaxnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_log2_minmaxnorm_filename)
# write_feather(wgc_nonzero_libnorm_log2_minmaxnorm, wgc_nonzero_libnorm_log2_minmaxnorm_dir_filename)

# wgc_nonzero_libnorm_sqrt_filename = paste0(datasetName_full, "_noAllZero_libnorm_sqrt.feather")
# wgc_nonzero_libnorm_sqrt_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_filename)
# write_feather(wgc_nonzero_libnorm_sqrt, wgc_nonzero_libnorm_sqrt_dir_filename)

# wgc_nonzero_libnorm_sqrt_minmaxnorm_filename = paste0(datasetName_full, "_noAllZero_libnorm_sqrt_minmax.feather")
# wgc_nonzero_libnorm_sqrt_minmaxnorm_dir_filename = file.path(save_dir, wgc_nonzero_libnorm_sqrt_minmaxnorm_filename)
# write_feather(wgc_nonzero_libnorm_sqrt_minmaxnorm, wgc_nonzero_libnorm_sqrt_minmaxnorm_dir_filename)

print(paste0(scen, " ", target_qc_type, " ", "saved!"))
