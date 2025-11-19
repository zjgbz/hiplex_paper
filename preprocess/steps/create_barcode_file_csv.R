#usage Rscript create_barcode_file.R "/dcs05/hongkai/data/zwang5/NathansLab/pipeline/barcode.xlsx"   ${output_dir}/data/
args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_dir <- args[2]

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
})


barcode <- read.csv(input_file) %>%
	mutate(target = gsub(" ","_",target)) %>%
	mutate(target = gsub("\\/","_",target)) %>%
	mutate(target = gsub("\\.","_",target)) %>%
	mutate(target = gsub("-","",target))

print("The barcode being used is")
colnames(barcode)[2] <- "BC"
barcode

barcode_fasta <- lapply(c(1:nrow(barcode)), function(x){
	c(glue(">{barcode$target[x]}"), barcode$BC[x])
})

barcode_fasta <- Reduce(c,barcode_fasta)

##same barcode for both reads
write.table(barcode_fasta, file=paste0(output_dir,"Next_CUT_Tag_barcode_rev_next.fasta"), row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(barcode_fasta, file=paste0(output_dir,"Next_CUT_Tag_barcode_fwd_next.fasta"), row.names=FALSE, col.names=FALSE, quote=FALSE)
