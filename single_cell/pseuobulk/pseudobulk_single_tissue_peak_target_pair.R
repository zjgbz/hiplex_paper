rm(list=ls())
library(GenomicAlignments)
library(ChIPseeker)
library(Biostrings)
library(glue)

tissue_name = "K562"
target_pair_name_list = c("H3K27me3-H3K27me3", "H3K4me3-H3K4me3", "H3K27me3-H3K4me3")

sc_root_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
predefined_dir = file.path(sc_root_dir, "RNA-seq", "ccre_region")
tissue_root_dir = file.path(sc_root_dir, tissue_name)
pseudobulk_root_dir = file.path(tissue_root_dir, "homotone_heterotone_merged", "pseudobulk")
pseudobulk_peak_dir = file.path(pseudobulk_root_dir, "peak", "data_peak_auc_0.05")

for (target_pair_name in target_pair_name_list) {
    pseudobulk_peak_filename = glue("{target_pair_name}.stringent.bed")
    pseudobulk_peak_dir_filename = file.path(pseudobulk_peak_dir, pseudobulk_peak_filename)

    pseudobulk_peak <- ChIPseeker::readPeakFile(pseudobulk_peak_dir_filename, as="GRanges")
    pseudobulk_peak_df <- as.data.frame(pseudobulk_peak)

    pseudobulk_peak_df$pos <- with(pseudobulk_peak_df, paste0(seqnames, "_", start, "_", end))
    pseudobulk_peak_df$gene_id <- pseudobulk_peak_df$pos
    pseudobulk_peak_df$gene_name <- pseudobulk_peak_df$pos
    pseudobulk_peak_df$seqnames <-gsub("chr", "", pseudobulk_peak_df$seqnames)
    pseudobulk_peak_df <- pseudobulk_peak_df[, c("pos", "seqnames", "start", "end", "width", "strand", "gene_id", "gene_name")]

    pseudobulk_peak_filename = glue("{tissue_name}_{target_pair_name}_peak.tsv")
    pseudobulk_peak_dir_filename = file.path(predefined_dir, pseudobulk_peak_filename)
    write.table(pseudobulk_peak_df, file=pseudobulk_peak_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)
}
