rm(list=ls())
library(GenomicAlignments)
library(ChIPseeker)
library(Biostrings)
library(glue)

tissue_name_list = c("brain", "iPSC")
final_tissue_name = "brain_iPSC"
target_pair_name_list = c("H3K27me3-H3K27me3", "H3K4me3-H3K4me3", "H3K27me3-H3K4me3")
target_pair_name_final = "tpc-01"

sc_root_dir = "/projects/foundation_model_for_single_cell_multiomics_data/hiplex/single_cell"
predefined_dir = file.path(sc_root_dir, "RNA-seq", "ccre_region")

combine_pseudobulk_peak_list = c()
for (tissue_name in tissue_name_list) {
    tissue_root_dir = file.path(sc_root_dir, tissue_name)
    tissue_pseudobulk_root_dir = file.path(tissue_root_dir, "homotone_heterotone_merged", "pseudobulk")
    tissue_pseudobulk_peak_dir = file.path(tissue_pseudobulk_root_dir, "peak", "data_peak_auc_0.05")
    for (target_pair_name in target_pair_name_list) {
        pseudobulk_peak_filename = glue("{target_pair_name}.stringent.bed")
        pseudobulk_peak_dir_filename = file.path(tissue_pseudobulk_peak_dir, pseudobulk_peak_filename)
        pseudobulk_peak <- ChIPseeker::readPeakFile(pseudobulk_peak_dir_filename, as="GRanges")
        combine_pseudobulk_peak_list = append(combine_pseudobulk_peak_list, pseudobulk_peak)
    }
}

combine_pseudobulk_peak_grange <- GenomicRanges::reduce(combine_pseudobulk_peak_list)
combine_pseudobulk_peak_df <- as.data.frame(combine_pseudobulk_peak_grange)

combine_pseudobulk_peak_df$pos <- with(combine_pseudobulk_peak_df, paste0(seqnames, "_", start, "_", end))
combine_pseudobulk_peak_df$gene_id <- combine_pseudobulk_peak_df$pos
combine_pseudobulk_peak_df$gene_name <- combine_pseudobulk_peak_df$pos
combine_pseudobulk_peak_df$seqnames <-gsub("chr", "", combine_pseudobulk_peak_df$seqnames)
combine_pseudobulk_peak_df <- combine_pseudobulk_peak_df[, c("pos", "seqnames", "start", "end", "width", "strand", "gene_id", "gene_name")]

combine_pseudobulk_peak_filename = glue("{final_tissue_name}_{target_pair_name_final}_peak_union.tsv")
combine_pseudobulk_peak_dir_filename = file.path(predefined_dir, combine_pseudobulk_peak_filename)
write.table(combine_pseudobulk_peak_df, file=combine_pseudobulk_peak_dir_filename, sep="\t", quote=FALSE, row.names=FALSE)