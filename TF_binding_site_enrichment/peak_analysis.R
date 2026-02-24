install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    # Bioconductor packages
    bioc_packages <- c("plyranges", "ChIPseeker")
    if (pkg %in% bioc_packages) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

packages <- c("plyranges", "ChIPseeker")
lapply(packages, install_if_missing)

suppressPackageStartupMessages({
  library(plyranges)
  library(ChIPseeker)
})

source("../hiplex_paper/TF_binding_site_enrichment/peak_enrichment.R")
options <- commandArgs(trailingOnly = TRUE)
target_bed_file <- options[1]
control_bed_file <- options[2]
out_dir <- options[3]
num_skip <- as.numeric(options[4])
# target_bed_file <- "/dcs05/hongkai/data/next_cutntag/bulk/wgc/mixed/800/A.bed"
# control_bed_file <- "/dcs05/hongkai/data/next_cutntag/script/motif_analysis/test_run_result/archive/background/control/A_output_matchcontrol.bed"

if (grepl("stringent", target_bed_file, fixed = TRUE)) {
    target_region <- readPeakFile(target_bed_file)
} else {
    target_region <- read_bed(target_bed_file)
}

control_df <- as.data.frame(read.table(control_bed_file,
                                        header = FALSE, 
                                        sep="\t",
                                        stringsAsFactors=FALSE, 
                                        quote="", 
                                        skip= num_skip))
control_region <- makeGRangesFromDataFrame(control_df, seqnames.field="V1", start.field="V2", end.field="V3", strand.field="V6", keep.extra.columns=TRUE)
enrich_motif_res <- enrich_motif(target_region, control_region)
enrich_motif_res <- enrich_motif_res[order(enrich_motif_res$FDR),]
out_file <- out_dir

write.table(enrich_motif_res, file = out_file, quote = FALSE, sep = "\t", row.names = TRUE)
