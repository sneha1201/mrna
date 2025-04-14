#!/usr/bin/env Rscript

# Suppress startup messages and load required packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(optparse)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Path to counts matrix"),
  make_option(c("-m", "--metadata"), type="character", help="Path to metadata file"),
  make_option(c("-C", "--comparisons"), type="character", help="Path to comparison file (2 columns: numerator, denominator)"),
  make_option(c("-l", "--lfc"), type="numeric", default=1, help="Log2 fold change cutoff [default: %default]"),
  make_option(c("-p", "--padj"), type="numeric", default=0.05, help="Adjusted p-value cutoff [default: %default]"),
  make_option(c("-s", "--condition"), type="character", default="Condition", help="Condition column in metadata [default: %default]"),
  make_option(c("-o", "--outdir"), type="character", default="DESeq2_results", help="Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Load count matrix
cts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
des <- read.table(opt$metadata, header = TRUE, row.names = 1, check.names = FALSE)

# Clean up any whitespace
colnames(cts) <- trimws(colnames(cts))
rownames(des) <- trimws(rownames(des))

# Match and subset samples
cts <- cts[, rownames(des), drop = FALSE]

# Build and run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts, colData = des, design = as.formula(paste("~", opt$condition)))
dds <- DESeq(dds)

# Export normalized counts
ncts <- as.data.frame(counts(dds, normalized = TRUE))
ncts$IDS <- rownames(ncts)
write.csv(ncts, file = file.path(opt$outdir, "Normalized_expression_values.csv"), quote = FALSE)

# Load comparisons
comparisons <- read.table(opt$comparisons, header = TRUE, sep = ",", stringsAsFactors = FALSE)
if (!all(c("numerator", "denominator") %in% colnames(comparisons))) {
  stop("❌ Comparison file must have two columns: numerator, denominator")
}

# DE and plotting
DElist <- function(deseqObject, numerator, denominator, lfcCutoff, padjCutoff, condCol) {
  res <- results(deseqObject,
                 contrast = c(condCol, numerator, denominator),
                 alpha = padjCutoff,
                 lfcThreshold = lfcCutoff,
                 altHypothesis = "greaterAbs",
                 pAdjustMethod = "none")

  res <- as.data.frame(res)
  res$IDS <- rownames(res)

  summaryFile <- file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_summary.txt"))
  sink(summaryFile); print(summary(res)); sink()

  volcanoPlot <- EnhancedVolcano(res, lab = res$IDS,
                                  x = 'log2FoldChange', y = 'padj',
                                  title = paste0(numerator, " VS ", denominator),
                                  pCutoff = padjCutoff, FCcutoff = lfcCutoff,
                                  drawConnectors = TRUE, widthConnectors = 0.5)

  ggsave(file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_volcano_plot.png")),
         plot = volcanoPlot, width = 12, height = 10, units = "in", dpi = 300)

  res_filtered <- res[!is.na(res$padj) & res$padj < padjCutoff & abs(res$log2FoldChange) > lfcCutoff, ]

  write.csv(res, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Unfiltered_DE.csv")), quote = FALSE)
  write.csv(res_filtered, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Filtered_LFC_", lfcCutoff, "_FDR_", padjCutoff, ".csv")), quote = FALSE)

  # Heatmap of top DE genes (averaged by condition and Z-scored)
if (nrow(res_filtered) > 0) {
  top_genes <- head(order(res_filtered$padj), min(50, nrow(res_filtered)))
  valid_ids <- res_filtered$IDS[top_genes]

  # Subset normalized counts to selected genes and samples
  ncts_subset <- ncts[valid_ids, , drop = FALSE]
  ncts_subset$IDS <- rownames(ncts_subset)

  # Filter only samples from numerator and denominator
samples_to_use <- rownames(des)[des[[condCol]] %in% c(numerator, denominator)]
sample_conditions <- des[samples_to_use, condCol, drop = TRUE]
names(sample_conditions) <- samples_to_use

# Subset normalized counts to selected genes and selected samples
ncts_subset <- ncts[valid_ids, samples_to_use, drop = FALSE]

# Average expression by condition (numerator and denominator only)
averaged_expr <- sapply(c(numerator, denominator), function(cond) {
  samples_in_cond <- names(sample_conditions[sample_conditions == cond])
  rowMeans(ncts_subset[, samples_in_cond, drop = FALSE], na.rm = TRUE)
})


  # Compute Z-scores (row-wise)
  z_scores <- t(scale(t(averaged_expr)))
  rownames(z_scores) <- rownames(ncts_subset)
  colnames(z_scores) <- unique(sample_conditions)

  # Remove NA or infinite rows
  z_scores <- z_scores[apply(z_scores, 1, function(x) all(is.finite(x))), , drop = FALSE]

  if (nrow(z_scores) >= 2) {
    pheatmap(z_scores, cluster_cols = FALSE, cluster_rows = FALSE,
             filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_Genes_Heatmap_Zscore_Avg.png")),
             height = 10, width = 14)
  } else if (nrow(z_scores) == 1) {
    pheatmap(z_scores, cluster_cols = FALSE, cluster_rows = FALSE,
             filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_1_Gene_Heatmap_Zscore_Avg.png")),
             height = 4, width = 10)
  } else {
    warning(paste("⚠ Skipped heatmap for", numerator, "vs", denominator, ": No valid Z-score expression data."))
  }
}


  message(paste0("✓ Completed: ", numerator, " vs ", denominator))
}

# Run comparisons
for (i in seq_len(nrow(comparisons))) {
  DElist(dds, comparisons$numerator[i], comparisons$denominator[i], opt$lfc, opt$padj, opt$condition)
}

