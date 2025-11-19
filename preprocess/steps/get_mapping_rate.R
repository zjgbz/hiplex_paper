#usage:
#Rscript /dcs05/hongkai/data/next_cutntag/workflows/get_mapping_rate.R "/dcs05/hongkai/data/next_cutntag/workflows/test/ZL174/logs/step3_err.log" "/dcs05/hongkai/data/next_cutntag/workflows/test/ZL174/analysis" "ZL174"
#Rscript /dcs05/hongkai/data/next_cutntag/workflows/get_mapping_rate.R "/dcs05/hongkai/data/kyu/Andy_lab/Masa_project/log/alignment/A13D_Parental_S22_L006.stderr" "/dcs05/hongkai/data/next_cutntag/workflows/test/ZL175/analysis" "ZL175"
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
dir_analysis <- args[2]
sample_id <- args[3]

#
#filename <-  "/dcs05/hongkai/data/next_cutntag/workflows/test/ZL174/logs/step3_err.log" 
#dir_analysis<- "/dcs05/hongkai/data/next_cutntag/workflows/test/ZL174/analysis"
#sample_id <- "ZL174"
#
suppressPackageStartupMessages({
  library(glue)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
})

mapping_rate_all <- list()

data <- readLines(filename)

mapping_rate <- data[grep("overall alignment rate",data)]
mapping_rate <- as.vector(unlist(sapply(mapping_rate,function(x)sub("% overall alignment rate","",x))))
mapping_rate <- as.vector(unlist(sapply(mapping_rate,function(x)as.numeric(x))))/100

total_read <- data[grep("reads",data)]
total_read <- as.vector(unlist(sapply(total_read,function(x)sub(" reads; of these:","",x))))
total_read <- as.vector(unlist(sapply(total_read,function(x)sub(" reads","",x))))
total_read <- as.vector(unlist(sapply(total_read,function(x)as.numeric(x))))

pre_sample_name <- data[grep("_log_file.txt",data)]
pre_sample_name <- as.vector(unlist(sapply(pre_sample_name,function(x)sub(".*METRICS_FILE=","",x))))
pre_sample_name <- as.vector(unlist(sapply(pre_sample_name,function(x)sub("_log_file\\.txt.*","",x))))
# sample_name <- pre_sample_name[c(FALSE, TRUE)]
sample_name <- pre_sample_name


mapping_rate_all[[1]] <- data.frame(cell=sample_name,rate=mapping_rate,read=total_read)

mapping_rate_combine <- Reduce(rbind,mapping_rate_all)

##paired-end sequencing
dulication_rate <- sapply(mapping_rate_combine[,1], function(x){
    data <- data.frame(fread(glue("{x}_log_file.txt")))
    if ("V3" %in% colnames(data)) {
        return(cbind(0,0,0))
    }
    return(cbind(data[1,"READ_PAIRS_EXAMINED"],data[1,"PERCENT_DUPLICATION"],data[1,"READ_PAIRS_EXAMINED"] - data[1,"READ_PAIR_DUPLICATES"]))
})

align_info <- data.frame(mapping_rate_combine, t(dulication_rate))
colnames(align_info) <- c("sample_name","alignment_rate","sequence_read","map_read","duplication_rate","unique_read")
# align_info$sample_name <- basename(align_info$sample_name)
align_info <- align_info %>%
    mutate(sample_name = sub(".*/","", sample_name)) %>%
    separate(sample_name, c("read1_bc","read2_bc"),sep="-", remove = FALSE)

# User Input: Replace "{sample_id}" with the actual sample ID used by the user
write.csv(align_info, file=glue("{dir_analysis}/mapping_rate_pe_{sample_id}.csv"), row.names=FALSE)

align_info_2 <- align_info[align_info$read1_bc != align_info$read2_bc, ]
tmp_bc <- align_info_2$read1_bc
align_info_2$read1_bc <- align_info_2$read2_bc
align_info_2$read2_bc <- tmp_bc
if (!paste0(align_info_2[1, "read1_bc"], "-", align_info_2[1, "read2_bc"]) %in% align_info$sample_name) {
  align_info <- rbind(align_info, align_info_2)
}

# get total reads matrix
align_info_seq_read <- align_info %>%
    select(read1_bc, read2_bc, sequence_read) %>%
    pivot_wider(names_from = read2_bc, values_from = sequence_read)


align_info_seq_read_mat <- as.data.frame(align_info_seq_read[,-1])
row.names(align_info_seq_read_mat) <- align_info_seq_read$read1_bc
target_pair_order <- sort(rownames(align_info_seq_read_mat), method = "radix")
align_info_seq_read_mat <- align_info_seq_read_mat[target_pair_order, target_pair_order]

write.csv(align_info_seq_read_mat, file=glue("{dir_analysis}/sequence_read_mat_{sample_id}.csv"))

ht <- Heatmap(align_info_seq_read_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/sequence_read_mat_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

# get log10 total reads matrix
align_info_seq_read_mat_log10 <- log10(align_info_seq_read_mat+1)
# write.csv(align_info_seq_read_mat_log10, file=glue("{dir_analysis}/sequence_read_mat_log10_{sample_id}.csv"))

ht <- Heatmap(align_info_seq_read_mat_log10, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/sequence_read_mat_log10_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()



# get unique reads
align_info_unique_read <- align_info %>%
    select(read1_bc, read2_bc, unique_read) %>%
    pivot_wider(names_from = read2_bc, values_from = unique_read)
# align_info_part_1 <- align_info %>% select(read1_bc, read2_bc, unique_read)
# align_info_part_2 <- align_info_part_1[align_info_part_1$read1_bc != align_info_part_1$read2_bc, ]
# align_info_part_2$read1_bc <- align_info_part_1$read2_bc
# align_info_part_2$read2_bc <- align_info_part_1$read1_bc
# align_info_parts <- rbind(align_info_part_1, align_info_part_2)
# align_info_unique_read <- align_info_parts %>% 
#   pivot_wider(names_from = read2_bc, values_from = unique_read)

align_info_unique_read_mat <- as.data.frame(align_info_unique_read[,-1])
row.names(align_info_unique_read_mat) <- align_info_unique_read$read1_bc
align_info_unique_read_mat <- align_info_unique_read_mat[target_pair_order, target_pair_order]

write.csv(align_info_unique_read_mat, file=glue("{dir_analysis}/unique_read_mat_{sample_id}.csv"))

ht <- Heatmap(align_info_unique_read_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/unique_read_mat_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

# get log10 unique reads
align_info_unique_read_mat_log10 <- log10(align_info_unique_read_mat+1)
# write.csv(align_info_unique_read_mat_log10, file=glue("{dir_analysis}/unique_read_mat_log10_{sample_id}.csv"))

ht <- Heatmap(align_info_unique_read_mat_log10, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/unique_read_mat_log10_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

# mapping rate
align_info_map_rate <- align_info %>%
    select(read1_bc, read2_bc, alignment_rate) %>%
    pivot_wider(names_from = read2_bc, values_from = alignment_rate)

align_info_map_rate_mat <- as.data.frame(align_info_map_rate[,-1])
row.names(align_info_map_rate_mat) <- align_info_map_rate$read1_bc
align_info_map_rate_mat <- align_info_map_rate_mat[target_pair_order, target_pair_order]
write.csv(align_info_map_rate_mat, file=glue("{dir_analysis}/alignment_rate_mat_{sample_id}.csv"))

ht <- Heatmap(align_info_map_rate_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/alignment_rate_mat_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

# duplication rate
align_info_dup_rate <- align_info %>%
    select(read1_bc, read2_bc, duplication_rate) %>%
    pivot_wider(names_from = read2_bc, values_from = duplication_rate)

align_info_dup_rate_mat <- as.data.frame(align_info_dup_rate[,-1])
row.names(align_info_dup_rate_mat) <- align_info_dup_rate$read1_bc
align_info_dup_rate_mat <- align_info_dup_rate_mat[target_pair_order, target_pair_order]
write.csv(align_info_dup_rate_mat, file=glue("{dir_analysis}/duplication_rate_mat_{sample_id}.csv"))

ht <- Heatmap(align_info_dup_rate_mat, cluster_rows = F, cluster_columns = F, row_names_side = "left")

pdf(glue("{dir_analysis}/duplication_rate_mat_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()
