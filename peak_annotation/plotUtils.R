# ==============================================================================
# Plot Utilities for Hi-Plex CUT&Tag Peak Annotation Visualization
# ==============================================================================
# This script provides utility functions for visualizing peak annotation
# results, creating bar plots, and extracting annotation statistics.
# Code adapted from ChIPSeeker
# Yu G, Wang L, He Q (2015). “ChIPseeker: an R/Bioconductor package for ChIP peak annotation, comparison and visualization.” Bioinformatics, 31(14), 2382-2383. doi:10.1093/bioinformatics/btv145.
# ==============================================================================

# Load Required Libraries ======================================================
library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(ComplexHeatmap)
library(grid)
library(stats)
library(latex2exp)

# Define Feature Categories ====================================================

# Standard features without enhancer annotation
allFeaturesNoEnhancer <- c(
  "Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron",
  "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)"
)

# Features with enhancer annotation
allFeaturesEnhancer <- c(
  "Enhancer", "Other Distal Intergenic", "Promoter", "5' UTR", "1st Exon",
  "1st Intron", "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)"
)

# Default feature set
allFeatures <- c(
  "Distal Intergenic", "Promoter", "5' UTR", "1st Exon", "1st Intron",
  "Other Exon", "Other Intron", "3' UTR", "Downstream (<=300)"
)

# Accessor Functions ===========================================================

#' Extract annotation statistics from csAnno object
#'
#' @param x csAnno object
#' @return Data frame with annotation statistics
getAnnoStat <- function(x) {
  if (!is(x, "csAnno")) {
    stop("Input must be a csAnno object")
  }
  return(x@annoStat)
}

# Data Transformation Functions ================================================

#' Convert list of data frames to single data frame with ID column
#'
#' @param dataList Named list of data frames
#' @return Combined data frame with .id column
list_to_dataframe <- function(dataList) {
  # Handle unnamed lists
  if (is.null(names(dataList))) {
    return(do.call('rbind', dataList))
  }
  
  # Get all unique column names
  cn <- lapply(dataList, colnames) %>% unlist() %>% unique()
  cn <- c('.id', cn)
  
  # Process each data frame to ensure consistent columns
  dataList2 <- lapply(seq_along(dataList), function(i) {
    data <- dataList[[i]]
    data$.id <- names(dataList)[i]
    
    # Add missing columns as NA
    idx <- !cn %in% colnames(data)
    if (sum(idx) > 0) {
      for (col in cn[idx]) {
        data[, col] <- NA
      }
    }
    return(data[, cn])
  })
  
  # Combine all data frames
  res <- do.call('rbind', dataList2)
  res$.id <- factor(res$.id, levels = rev(names(dataList)))
  
  return(res)
}

# Color Palette Functions ======================================================

#' Get color palette for plotting
#'
#' @param n Number of colors needed
#' @param option Integer (1-6) specifying color palette option
#' @param missingIdx Indices to exclude from palette
#' @return Vector of n colors
getCols <- function(n, option = 3, missingIdx = NULL) {
  # Define color palettes
  col1 <- c(
    "#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3", "#fdb462",
    "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f"
  )
  
  col2 <- c(
    "#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c", "#a6cee3",
    "#006d2c", "#d73027", "#4d4d4d", "#8c510a", "#78c679", "#7f0000",
    "#41b6c4", "#e7298a", "#54278f"
  )
  
  col3 <- c(
    "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c",
    "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"
  )
  
  # Composite palettes
  col4 <- c(col3, col2, col1)
  col5 <- rainbow(120)
  col6 <- c(col2, col1, col3)
  
  # Select palette
  cols <- list(col1, col2, col3, col4, col5, col6)
  col <- cols[[option]]
  
  # Remove missing indices if specified
  if (!is.null(missingIdx)) {
    col <- col[-missingIdx]
  }
  
  print(col)
  
  return(col[1:n])
}

# Plotting Functions ===========================================================

#' Create bar plot for annotation distribution
#'
#' @param anno.df Data frame with Feature and Frequency columns
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param title Plot title
#' @param categoryColumn Column name to use for categories
#' @param colorOption Integer (1-6) specifying color palette
#' @param width Bar width (optional)
#' @return ggplot object
plotAnnoBar.data.frame <- function(
    anno.df,
    xlab = "",
    ylab = "Percentage(%)",
    title = "Feature Distribution",
    categoryColumn,
    colorOption = 3,
    width = NULL
) {
  # Reverse factor levels for proper stacking
  anno.df$Feature <- factor(
    anno.df$Feature,
    levels = rev(levels(anno.df$Feature))
  )
  
  # Create base plot
  p <- ggplot(
    anno.df,
    aes_string(x = categoryColumn, fill = "Feature", y = "Frequency")
  ) +
    scale_x_discrete(labels = TeX)
  
  # Add bar geometry and formatting
  p <- p +
    geom_bar(stat = "identity", width = width) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)
  
  # Add color scale
  if (categoryColumn == 1) {
    p <- p +
      scale_x_continuous(breaks = NULL) +
      scale_fill_manual(
        values = rev(getCols(nrow(anno.df), option = colorOption)),
        guide = guide_legend(reverse = TRUE)
      )
  } else {
    p <- p +
      scale_fill_manual(
        values = rev(getCols(length(unique(anno.df$Feature)), option = colorOption)),
        guide = guide_legend(reverse = TRUE)
      )
  }
  
  return(p)
}

#' Create bar plot with custom color palette
#'
#' @param anno.df Data frame with Feature and Frequency columns
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param title Plot title
#' @param categoryColumn Column name to use for categories
#' @param colors Vector of colors to use
#' @param width Bar width (optional)
#' @return ggplot object
plotAnnoBar.data.frame.given.color <- function(
    anno.df,
    xlab = "",
    ylab = "Percentage(%)",
    title = "Feature Distribution",
    categoryColumn,
    colors,
    width = NULL
) {
  # Reverse factor levels for proper stacking
  anno.df$Feature <- factor(
    anno.df$Feature,
    levels = rev(levels(anno.df$Feature))
  )
  
  # Create base plot
  p <- ggplot(
    anno.df,
    aes_string(x = categoryColumn, fill = "Feature", y = "Frequency")
  )
  
  # Add bar geometry and formatting
  p <- p +
    geom_bar(stat = "identity", width = width) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)
  
  # Add custom color scale
  if (categoryColumn == 1) {
    p <- p +
      scale_x_continuous(breaks = NULL) +
      scale_fill_manual(
        values = rev(colors),
        guide = guide_legend(reverse = TRUE)
      )
  } else {
    p <- p +
      scale_fill_manual(
        values = rev(colors),
        guide = guide_legend(reverse = TRUE)
      )
  }
  
  return(p)
}

#' Create bar plot for single target with feature filtering
#'
#' @param anno.df Data frame with Feature and Frequency columns
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param title Plot title
#' @param categoryColumn Column name to use for categories
#' @param colorOption Integer (1-6) specifying color palette
#' @param features Vector of expected features
#' @return ggplot object
plotAnnoBar.data.frame.one.target <- function(
    anno.df,
    xlab = "",
    ylab = "Percentage(%)",
    title = "Feature Distribution",
    categoryColumn = ".id",
    colorOption = 3,
    features
) {
  # Reverse factor levels for proper stacking
  anno.df$Feature <- factor(
    anno.df$Feature,
    levels = rev(levels(anno.df$Feature))
  )
  
  # Identify missing features
  missingIndices <- c()
  print(anno.df)
  
  for (i in seq_along(features)) {
    feature <- features[i]
    if (!feature %in% anno.df$Feature) {
      missingIndices <- append(missingIndices, i)
    }
  }
  
  # Remove zero-frequency features
  anno.df <- anno.df[anno.df$Frequency != 0, ]
  print(missingIndices)
  
  # Create base plot
  p <- ggplot(
    anno.df,
    aes_string(x = categoryColumn, fill = "Feature", y = "Frequency")
  )
  
  # Add bar geometry and formatting
  p <- p +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title)
  
  # Add color scale
  if (categoryColumn == 1) {
    p <- p +
      scale_x_continuous(breaks = NULL) +
      scale_fill_manual(
        values = rev(getCols(nrow(anno.df), option = colorOption, missingIndices)),
        guide = guide_legend(reverse = TRUE)
      )
  } else {
    col <- getCols(
      length(unique(anno.df$Feature)),
      option = colorOption,
      missingIndices
    )
    print(col)
    col <- col[!is.na(col)]
    print(col)
    
    p <- p +
      scale_fill_manual(
        values = rev(col),
        guide = guide_legend(reverse = TRUE)
      )
  }
  
  return(p)
}

#' Create bar plot with red outline for specific rows
#'
#' @param anno.df Data frame with Feature and Frequency columns
#' @param xlab X-axis label
#' @param ylab Y-axis label
#' @param title Plot title
#' @param categoryColumn Column name to use for categories
#' @param colorOption Integer (1-6) specifying color palette
#' @return ggplot object
plotAnnoBar.data.frame.red.even.rows <- function(
    anno.df,
    xlab = "",
    ylab = "Percentage(%)",
    title = "Feature Distribution",
    categoryColumn,
    colorOption = 3
) {
  # Reverse factor levels for proper stacking
  anno.df$Feature <- factor(
    anno.df$Feature,
    levels = rev(levels(anno.df$Feature))
  )
  
  # Add outline color for treated samples
  anno.df$outline <- ifelse(
    grepl("treated", anno.df[[categoryColumn]]),
    "red",
    NA
  )
  
  # Create base plot
  p <- ggplot(anno.df) +
    geom_bar(
      aes(x = .id, fill = Feature, y = Frequency, col = outline),
      stat = "identity"
    ) +
    coord_flip() +
    theme_bw() +
    theme(axis.text = element_text(size = 16)) +
    ylab(ylab) +
    xlab(xlab) +
    ggtitle(title) +
    geom_bar(
      aes(x = .id, y = Frequency, group = .id),
      stat = "summary",
      fun = sum,
      fill = "transparent"
    ) +
    scale_color_identity()
  
  # Add color scale
  if (categoryColumn == 1) {
    p <- p +
      scale_x_continuous(breaks = NULL) +
      scale_fill_manual(
        values = rev(getCols(nrow(anno.df), option = colorOption)),
        guide = guide_legend(reverse = TRUE)
      )
  } else {
    p <- p +
      scale_fill_manual(
        values = rev(getCols(length(unique(anno.df$Feature)), option = colorOption)),
        guide = guide_legend(reverse = TRUE)
      )
  }
  
  return(p)
}

# Matrix Conversion Functions ==================================================

#' Convert annotation list to matrix format
#'
#' @param anno Named list of annotation data frames
#' @param allFeatures Vector of all possible features
#' @return Matrix with features as columns and samples as rows
getAnnotationMatrix <- function(anno, allFeatures) {
  # Initialize result list
  resultList <- list()
  for (f in allFeatures) {
    resultList[[f]] <- c()
  }
  
  # Process each combination
  for (comb in names(anno)) {
    features <- as.character(anno[[comb]]$Feature)
    frequencies <- as.numeric(anno[[comb]]$Frequency)
    names(frequencies) <- features
    frequencies <- as.list(frequencies)
    
    # Add zeros for missing features
    missingFeatures <- allFeatures[!allFeatures %in% features]
    for (f in missingFeatures) {
      frequencies[[f]] <- 0
    }
    
    # Append frequencies to result
    for (f in allFeatures) {
      resultList[[f]] <- append(resultList[[f]], frequencies[[f]])
    }
  }
  
  # Convert to matrix
  resultDF <- as.data.frame(resultList, check.names = FALSE)
  rownames(resultDF) <- names(anno)
  resultMat <- as.matrix(resultDF)
  
  return(resultMat)
}

# Enhancer Annotation Functions ================================================

#' Get annotation statistics with enhancer classification
#'
#' @param xT Named list of csAnno objects
#' @param ehancerDb GRanges object with enhancer regions
#' @return List of modified annotation statistics
getEnhancerAnnoStats <- function(xT, ehancerDb) {
  annoStats <- list()
  
  for (tag in names(xT)) {
    # Extract distal intergenic regions
    distalElem <- xT[[tag]]@anno[xT[[tag]]@anno$annotation == "Distal Intergenic"]
    
    # Calculate enhancer overlap percentage
    enhancerPercentage <- length(
      unique(queryHits(findOverlaps(distalElem, ehancerDb)))
    ) / length(distalElem)
    
    annoStat <- getAnnoStat(xT[[tag]])
    
    # Split distal intergenic into enhancer and other
    if (!is.nan(enhancerPercentage)) {
      distalFreq <- annoStat$Frequency[annoStat$Feature == "Distal Intergenic"]
      enhancerFreq <- distalFreq * enhancerPercentage
      otherDistalFreq <- distalFreq * (1 - enhancerPercentage)
      
      # Add new categories
      annoStat <- rbind(
        annoStat,
        data.frame(Feature = "Enhancer", Frequency = enhancerFreq)
      )
      annoStat <- rbind(
        annoStat,
        data.frame(Feature = "Other Distal Intergenic", Frequency = otherDistalFreq)
      )
      
      # Remove original distal intergenic category
      annoStat <- annoStat[annoStat$Feature != "Distal Intergenic", ]
      annoStats[[tag]] <- annoStat
    }
  }
  
  return(annoStats)
}

# Statistical Testing Functions ================================================

#' Filter tags by chi-square test significance
#'
#' @param resultMatT Matrix of frequencies for treatment group
#' @param resultMatV Matrix of frequencies for control group
#' @param xT Named list of treatment csAnno objects
#' @param xV Named list of control csAnno objects
#' @param cutoff P-value cutoff (not currently used)
#' @return Vector of significant tag names
filterByChisq <- function(resultMatT, resultMatV, xT, xV, cutoff) {
  tags <- c()
  
  for (tag in rownames(resultMatT)) {
    # Convert percentages to counts
    resultTagT <- resultMatT[tag, ] / 100 * length(xT[[tag]]@anno)
    resultTagV <- resultMatV[tag, ] / 100 * length(xV[[tag]]@anno)
    
    # Perform chi-square test
    a <- chisq.test(rbind(resultTagT, resultTagV))
    
    if (!is.nan(a$p.value)) {
      if (a$p.value < 0.05) {
        tags <- append(tags, tag)
        cat(tag, " ", a$p.value, "\n")
      }
    }
  }
  
  return(tags)
}

# Heterotone Analysis Functions ================================================

#' Get reversed heterotone tag name
#'
#' @param tag Character string in format "tag1-tag2"
#' @return Character string in format "tag2-tag1"
getRevHeterotone <- function(tag) {
  tag <- as.character(tag)
  tags <- strsplit(tag, "-")[[1]]
  tag1 <- tags[1]
  tag2 <- tags[2]
  return(paste0(tag2, "-", tag1))
}

#' Get unique heterotone combinations for a given tag
#'
#' @param anno.df Data frame with .id column containing tag names
#' @param tag Character string of tag to match
#' @return Vector of unique heterotone combinations
getUniqueHeterotones <- function(anno.df, tag) {
  uniqueHeterotones <- c()
  
  for (comb in unique(anno.df$.id)) {
    if (substring(comb, 1, nchar(tag)) == tag) {
      uniqueHeterotones <- append(uniqueHeterotones, comb)
    } else {
      uniqueHeterotones <- append(uniqueHeterotones, getRevHeterotone(comb))
    }
  }
  
  uniqueHeterotones <- unique(uniqueHeterotones)
  return(uniqueHeterotones)
}

#' Get combined annotation statistics for unique heterotones
#'
#' @param x Named list of csAnno objects
#' @param allFeatures Vector of all possible features
#' @param uniqueHeterotones Vector of unique heterotone combinations
#' @return List of combined annotation statistics
getAnnoStatForUniqueHeterotones <- function(x, allFeatures, uniqueHeterotones) {
  newAnnoStatV <- list()
  
  for (tag in uniqueHeterotones) {
    # Get forward annotation
    annoStat <- x[[tag]]@annoStat
    peakNum <- x[[tag]]@peakNum
    annoStat$Frequency <- annoStat$Frequency * peakNum
    
    # Get reverse annotation
    revAnnoStat <- x[[getRevHeterotone(tag)]]@annoStat
    revPeakNum <- x[[getRevHeterotone(tag)]]@peakNum
    revAnnoStat$Frequency <- revAnnoStat$Frequency * revPeakNum
    
    # Combine frequencies
    freq <- list()
    for (f in allFeatures) {
      freq[[f]] <- 0
    }
    
    for (f in annoStat$Feature) {
      freq[[f]] <- freq[[f]] + annoStat$Frequency[annoStat$Feature == f]
    }
    
    for (f in revAnnoStat$Feature) {
      freq[[f]] <- freq[[f]] + revAnnoStat$Frequency[revAnnoStat$Feature == f]
    }
    
    # Create new annotation statistics
    newannoStat <- data.frame(
      Feature = names(freq),
      Frequency = (as.numeric(freq) / sum(as.numeric(freq)) * 100)
    )
    
    newAnnoStatV[[tag]] <- newannoStat
  }
  
  return(newAnnoStatV)
}

# Legend Extraction Function ===================================================

#' Extract legend from ggplot object
#'
#' @param my_ggp ggplot object
#' @return grob containing the legend
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# ==============================================================================
# End of Script
# ==============================================================================