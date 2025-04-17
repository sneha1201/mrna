#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(ggvenn)
  library(tidyverse)
})

# Argument parser
option_list <- list(
  make_option(c("--normalized"), type = "character", help = "Path to normalized expression matrix CSV"),
  make_option(c("--metadata"), type = "character", help = "Path to metadata TSV file"),
  make_option(c("--comparisons"), type = "character", help = "Path to comparisons CSV file"),
  make_option(c("--outdir"), type = "character", help = "Output directory to save Venn diagrams")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load data
data <- read.csv(opt$normalized, row.names = 1, check.names = FALSE)
metadata <- read.delim(opt$metadata, header = TRUE, stringsAsFactors = FALSE)
comparisons_df <- read.csv(opt$comparisons)

# Construct sample_groups dynamically
sample_groups <- list()
for (condition in unique(metadata$Condition)) {
  sample_ids <- metadata$SampleID[metadata$Condition == condition]
  sample_ids <- sample_ids[sample_ids %in% colnames(data)]  # Ensure valid columns
  if (length(sample_ids) > 0) {
    sample_groups[[condition]] <- rowMeans(data[, sample_ids])
  }
}

# Create output directory if it doesn't exist
if (!dir.exists(opt$outdir)) {
  dir.create(opt$outdir, recursive = TRUE)
}

# Function to generate and save Venn diagram
plot_venn <- function(sample_name1, sample_name2, output_file) {
  if (!(sample_name1 %in% names(sample_groups)) | !(sample_name2 %in% names(sample_groups))) {
    warning(paste("Skipping", sample_name1, "vs", sample_name2, "- no valid samples found."))
    return()
  }

  genes1 <- names(sample_groups[[sample_name1]][sample_groups[[sample_name1]] > 0])
  genes2 <- names(sample_groups[[sample_name2]][sample_groups[[sample_name2]] > 0])

  if (length(genes1) == 0 | length(genes2) == 0) {
    warning(paste("Skipping", sample_name1, "vs", sample_name2, "- no expressed genes found."))
    return()
  }

  gene_sets <- setNames(list(genes1, genes2), c(sample_name1, sample_name2))

  venn_plot <- ggvenn(
    gene_sets,
    fill_color = c("blue", "red"),
    set_name_size = 6,
    text_size = 5
  ) +
    theme_void() +
    theme(plot.background = element_rect(fill = "white", color = NA))

  ggsave(filename = file.path(opt$outdir, output_file),
         plot = venn_plot, width = 8, height = 8, dpi = 300, bg = "white")
}

# Loop through comparisons
for (i in 1:nrow(comparisons_df)) {
  numerator <- comparisons_df$numerator[i]
  denominator <- comparisons_df$denominator[i]
  filename <- paste0(numerator, "_vs_", denominator, ".png")
  plot_venn(numerator, denominator, filename)
}

cat("✅ Venn diagrams saved successfully in", opt$outdir, "\n")

