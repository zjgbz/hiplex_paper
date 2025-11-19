library(rtracklayer)
library(GenomicRanges)
library(glue)

options <- commandArgs(trailingOnly = TRUE)
if (length(options) != 0) {
  input_dir <- options[1]
  save_dir <- options[2]
}
control_bed_files <- Sys.glob(glue("{input_dir}/*/control/*_output_matchcontrol.bed"))
print(control_bed_files)
all_gr <- c()
for (bed_file in control_bed_files) {
  bed_content <- read.table(bed_file, sep = "\t", skip = 2)
  bed_grange <- makeGRangesFromDataFrame(
    bed_content,
    seqnames.field = "V1",
    start.field = "V2",
    end.field = "V3"
  )
  print(bed_file)
  print(bed_grange[start(bed_grange)==0,])
  all_gr <- append(all_gr, bed_grange)
}
ctrl_merged_bed <- reduce(all_gr)
export.bed(ctrl_merged_bed, glue("{save_dir}/ctrl.bed"))
# export.bed(ctrl_merged_bed, "/dcs05/hongkai/data/next_cutntag/script/motif_analysis/reorganized/control/ctrl.bed")
