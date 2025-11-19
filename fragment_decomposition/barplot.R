library(ggplot2)
library("dplyr")
library(glue)
library(latex2exp)
library(cowplot)

source("/dcs05/hongkai/data/next_cutntag/script/utils/map_target_pair_names.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/utils.R")
source("/dcs05/hongkai/data/next_cutntag/script/utils/filter_targets.R")

target_pair_list = filter_target_pairs(0.25)

list_to_dataframe <- function(dataList) {
  if (is.null(names(dataList)))
    return(do.call('rbind', dataList))
  
  cn <- lapply(dataList, colnames) %>% unlist %>% unique
  cn <- c('.id', cn)
  dataList2 <- lapply(seq_along(dataList), function(i) {
    data = dataList[[i]]
    data$.id = names(dataList)[i]
    idx <- ! cn %in% colnames(data)
    if (sum(idx) > 0) {
      for (i in cn[idx]) {
        data[, i] <- NA
      }
    }
    return(data[,cn])
  })
  res <- do.call('rbind', dataList2)
  res$.id <- factor(res$.id, levels=rev(names(dataList)))
  return(res)
}

frag_decomp_types <- c("valley-all-qc", "valley-V-qc")
plts <- list()
for (frag_decomp_type in frag_decomp_types) {
  frag_len_percent_file <- glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/frag_split_fastq-demux_per_V_{frag_decomp_type}.tsv")
  frag_len_percent <- read.table(frag_len_percent_file, sep = "\t", header = TRUE, row.names = 1)
  frag_len_percent_filtered <- frag_len_percent[target_pair_list,]
  colnames(frag_len_percent_filtered)[3] <- "dimer+"
  plts[[frag_decomp_type]] <- list()
  # Plot each top20 frags sorted on each fragment size_category
  for (size_category in colnames(frag_len_percent_filtered)) {
    frag_len_percent_filtered_ordered <- frag_len_percent_filtered[order(frag_len_percent_filtered[[size_category]], decreasing = TRUE), ]
    frag_len_percent_filtered_ordered_top50 <- frag_len_percent_filtered_ordered[1:20, ]
    frag_len_list <- list()
    for (tag in rownames(frag_len_percent_filtered_ordered_top50)) {
      percentage <- as.numeric(as.vector(frag_len_percent_filtered_ordered_top50[tag,]))
      frag_len_list[[tag]] <- data.frame(type=colnames(frag_len_percent_filtered_ordered_top50), percentage=percentage)
    }
    barplot_df <- list_to_dataframe(frag_len_list)
    barplot_df$.id <- map_target_names(barplot_df$.id , target_pair_mapping_df)
    barplot_df$.id <- factor(barplot_df$.id, levels=rev(map_target_names(rownames(frag_len_percent_filtered_ordered_top50),target_pair_mapping_df)))
    barplot_df$type <- factor(barplot_df$type, levels = rev(c("subnucleo", "monomer", "dimer+")))
    # data <- data.frame(x = age, y = hours, group = city)
    p <- ggplot(barplot_df, aes(fill=type, y=percentage, x=.id)) + 
      scale_x_discrete(labels = TeX) +
      geom_bar(position="fill", stat="identity")+
      theme_gray(base_size = 18) +
      scale_fill_manual(values = c("#00b0be","#ff8ca1","#9ac9db", "#FBE7C6")) + 
      theme_classic() + xlab("type") + ylab("Peak Number") +
      theme_bw() +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_rect(colour = "black", size=2),
            axis.title.y=element_blank()) +
      ggtitle(glue("Fragment length distribution with \nmost {size_category}")) + labs(color='') +
      guides(col = guide_legend(override.aes = list(linetype=2))) + 
      labs(fill= "",
           y = "Reads percentage",
           x = "condition") + coord_flip()
    plts[[frag_decomp_type]][[size_category]] <- p
    pdf(glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/barplot_frag_split_fastq-demux_per_V_{frag_decomp_type}_sort_on_{size_category}.pdf"), height=5, width=8)
    print(p)
    dev.off()
  }
}
legend <- get_legend(p)
for (frag_decomp_type in frag_decomp_types) {
  pdf(glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/barplot_frag_split_fastq-demux_per_V_{frag_decomp_type}_summary.pdf"), height=6, width=14)
  print(plot_grid(
    plts[[frag_decomp_type]][["subnucleo"]]+ theme(legend.position="none"),
    plts[[frag_decomp_type]][["monomer"]]+ theme(legend.position="none"),
    plts[[frag_decomp_type]][["dimer+"]]+ theme(legend.position="none"),
    legend,
    align = "h",
    axis = "bt",
    ncol = 4,
    rel_widths = c(1, 1, 1, 0.3)
  ))
  dev.off()
  
}




frag_decomp_types <- c("valley-all-qc", "valley-V-qc")
plts <- list()
for (frag_decomp_type in frag_decomp_types) {
  frag_len_percent_file <- glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/frag_split_fastq-demux_per_V_{frag_decomp_type}.tsv")
  frag_len_percent <- read.table(frag_len_percent_file, sep = "\t", header = TRUE, row.names = 1)
  frag_len_percent_filtered <- frag_len_percent[target_pair_list,]
  colnames(frag_len_percent_filtered)[3] <- "dimer+"
  plts[[frag_decomp_type]] <- list()
  # Plot each top20 frags sorted on each fragment size_category
  size_category=colnames(frag_len_percent_filtered)[1]
  
  frag_len_percent_filtered_ordered <- frag_len_percent_filtered[order(frag_len_percent_filtered[[size_category]], decreasing = TRUE), ]

  frag_len_list <- list()
  for (tag in rownames(frag_len_percent_filtered_ordered)) {
    percentage <- as.numeric(as.vector(frag_len_percent_filtered_ordered[tag,]))
    frag_len_list[[tag]] <- data.frame(type=colnames(frag_len_percent_filtered_ordered), percentage=percentage)
  }
  barplot_df <- list_to_dataframe(frag_len_list)
  barplot_df$.id <- map_target_names(barplot_df$.id , target_pair_mapping_df)
  barplot_df$.id <- factor(barplot_df$.id, levels=rev(map_target_names(rownames(frag_len_percent_filtered_ordered),target_pair_mapping_df)))
  barplot_df$type <- factor(barplot_df$type, levels = rev(c("subnucleo", "monomer", "dimer+")))
  write.csv(barplot_df, glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/barplot_frag_split_fastq-demux_per_V_{frag_decomp_type}_sort_on_all.csv"), quote = FALSE, row.names = FALSE)
  # data <- data.frame(x = age, y = hours, group = city)
  p <- ggplot(barplot_df, aes(fill=type, y=percentage, x=.id)) + 
    scale_x_discrete(labels = TeX) +
    geom_bar(position="fill", stat="identity")+
    theme_gray(base_size = 18) +
    scale_fill_manual(values = c("#00b0be","#ff8ca1","#9ac9db", "#FBE7C6")) + 
    theme_classic() + xlab("type") + ylab("Peak Number") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=2),
          axis.title.y=element_blank()) +
    ggtitle(glue("Fragment length distribution")) + labs(color='') +
    guides(col = guide_legend(override.aes = list(linetype=2))) + 
    labs(fill= "",
         y = "Reads percentage",
         x = "condition") + coord_flip()
  plts[[frag_decomp_type]][[size_category]] <- p
  pdf(glue("/dcs05/hongkai/data/next_cutntag/bulk/frag_len/barplot_frag_split_fastq-demux_per_V_{frag_decomp_type}_sort_on_all.pdf"), height=120, width=8)
  print(p)
  dev.off()
  
}
