library(httr)
library(jsonlite)
library(ChIPseeker)
library(rtracklayer)
library(GenomicRanges)
library(glue)

ENCODE_BASE_URL = 'https://www.encodeproject.org'
# Force return from the server in JSON format
headers <- add_headers('accept'='application/json')
ENCODE_BASE_URL <- "https://www.encodeproject.org/"


download_metadata <- read.csv("/dcs05/hongkai/data/next_cutntag/public_data/ChIPSeq/TFs/peak/download_file_info.csv")
peak_type <- "conservative IDR thresholded peaks"
download_metadata_selected <- download_metadata[download_metadata$output_type==peak_type,]

nopeak_type_target <- unique(download_metadata$target)[!unique(download_metadata$target) %in% unique(download_metadata_selected$target)]
download_metadata_selected_2 <- download_metadata[download_metadata$target %in% nopeak_type_target,]
download_metadata_selected_2 <- download_metadata_selected_2[download_metadata_selected_2$output_type=="IDR thresholded peaks",]
download_metadata_selected <- rbind(download_metadata_selected, download_metadata_selected_2)

encode4_vector <- c()
conditions <- c()
for (i in seq_along(download_metadata_selected$file_name)) {
  file_dir <- download_metadata_selected$file_name[i]
  split_array <- strsplit(file_dir, "[/]")[[1]]
  filename <- split_array[length(split_array)]
  filename <- gsub(".bed.gz", "", filename)
  file <- glue("files/{filename}")
  url <- paste0(ENCODE_BASE_URL, file, "?format=json")
  response <- GET(url, headers)
  if (status_code(response) == 200) {
    # Parse the JSON content
    search_file_json <- content(response, as = "parsed", type = "application/json")
    if (search_file_json$award$rfa != "ENCODE4") {
      print(search_file_json$award$rfa)
    } else {
      encode4_vector <- append(encode4_vector, i)
      conditions <- append(conditions, search_file_json$"simple_biosample_summary")
    }
  } else {
    # Handle errors
    stop("Failed to retrieve data: ", status_code(response))
  }
}
download_metadata_selected_new <- download_metadata_selected[encode4_vector,]
download_metadata_selected <- download_metadata_selected_new
rownames(download_metadata_selected) <- NULL
conditions_short <- gsub("genetically modified (insertion) using CRISPR targeting H. sapiens", "CRISPR targeting", conditions, fixed = TRUE)
conditions_short <- gsub("C-terminal eGFP-tagged ", "", conditions_short, fixed = TRUE)
conditions_short <- gsub("N-terminal eGFP-tagged ", "", conditions_short, fixed = TRUE)
download_metadata_selected$cond <- conditions_short
new_targets <- c()
for (i in 1: nrow(download_metadata_selected)) {
  target <- download_metadata_selected[i, "target"]
  cond_single <- download_metadata_selected[i, "cond"]
  if (cond_single != "") {
    new_target <- glue("({cond_single}){target}")
    new_targets <- append(new_targets, new_target)
  } else {
    new_targets <- append(new_targets, target)
  }
}
download_metadata_selected$new_targets <- new_targets

save_parent_dir <- "/dcs05/hongkai/data/next_cutntag/public_data/ChIPSeq/TFs/peak_renamed/"
for (tf in unique(download_metadata_selected$new_targets)) {
  tf_files <- download_metadata_selected[download_metadata_selected$new_targets==tf, "file_name"]
  if (length(tf_files) > 1) {
    peak <- c()
    for (tf_file in tf_files) {
      tf_peak <- readPeakFile(tf_file)
      peak <- append(peak, tf_peak)
    }
    peak <- reduce(peak)
  } else {
    peak <- readPeakFile(tf_files)
  }
  export.bed(peak, glue("{save_parent_dir}/{tf}.bed"))
}

allmotif <- list.files(glue("{save_parent_dir}"),pattern = ".bed")
gr <- GRangesList()

for (m in allmotif) {
  tf_gr <- import.bed(glue("{save_parent_dir}/{m}"))
  mcols(tf_gr) <- NULL
  gr[[gsub(".bed","",m)]] <- tf_gr
}

saveRDS(gr,file="/dcs05/hongkai/data/next_cutntag/public_data/ChIPSeq/TFs/peak.rds",compress="xz")
