
args <- commandArgs(trailingOnly = TRUE)
filename <- args[1]
dir_analysis <- args[2]
sample_id <- args[3]

print(sample_id)
 


suppressPackageStartupMessages({
library(glue)
library(dplyr)
library(ComplexHeatmap)
library(tidyr)
})


temp <- read.csv(file=glue("{dir_analysis}/read_stat_{sample_id}.csv"), header=FALSE)
row_name <- sub(".*/([^/]+)\\..*", "\\1", temp[,1]) %>%
              sub("\\.fastq$", "", .) 
read_stat <- data.frame(barcode = row_name, count = temp[,2]/4)
row.names(read_stat) <- row_name

write.csv(read_stat, file=glue("{dir_analysis}/read_stat_counts_{sample_id}.csv"))

##combine
temp <- read.csv(file=glue("{dir_analysis}/read_stat_{sample_id}.csv"), header=FALSE)
row_name <- sub(".*/([^/]+)\\..*", "\\1", temp[,1]) %>%
              sub("\\.fastq$", "", .)
              
row_name_seq_R1 <- data.frame(row_name = row_name) %>%
    separate(row_name, c("barcode1", "barcode2", "read"), sep="-|\\.") %>%
    mutate(count = temp[,2]/4) %>%
    filter(read == "R1") %>%
    group_by(barcode1) %>%
    summarize(count = sum(count)) %>%
    rename(barcode = barcode1)

row_name_seq_R2 <- data.frame(row_name = row_name) %>%
    separate(row_name, c("barcode1", "barcode2", "read"), sep="-|\\.") %>%
    mutate(count = temp[,2]/4) %>%
    filter(read == "R2") %>%
    group_by(barcode2) %>%
    summarize(count = sum(count)) %>%
    rename(barcode = barcode2)

read_stat_total <- left_join(row_name_seq_R1, row_name_seq_R2, by = "barcode") %>%
    mutate(count = count.x + count.y) %>%
    select(-count.x, -count.y)

write.csv(read_stat_total, file=glue("{dir_analysis}/read_stat_total_{sample_id}.csv"))

##heatmap
read_stat_plot_R1 <- data.frame(row_name = row_name) %>%
    separate(row_name, c("barcode1", "barcode2", "read"), sep="-|\\.") %>%
    mutate(count = temp[,2]/4) %>%
    filter(read == "R1") %>%
    select(-read) %>%
    pivot_wider(names_from = barcode2, values_from = count)

read_stat_plot_R1_mat <- as.matrix(read_stat_plot_R1[,-1])
row.names(read_stat_plot_R1_mat) <- read_stat_plot_R1$barcode1

ht <- Heatmap(read_stat_plot_R1_mat, cluster_rows = T, cluster_columns = T, row_names_side = "left")

pdf(glue("{dir_analysis}/read_mat_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

## log10 heatmap
read_stat_plot_R1_mat_log10 <- log10(read_stat_plot_R1_mat+1)
ht <- Heatmap(read_stat_plot_R1_mat_log10, cluster_rows = T, cluster_columns = T, row_names_side = "left")

pdf(glue("{dir_analysis}/read_mat_log10_{sample_id}.pdf"), width = 10, height = 10)
print(ht)
dev.off()

write.csv(read_stat_plot_R1_mat, glue("{dir_analysis}/read_stat_mat_{sample_id}.csv"))
