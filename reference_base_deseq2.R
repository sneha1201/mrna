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
                                  subtitle = NULL,
                                  caption = "",  # Hides the EnhancedVolcano caption
                                  pCutoff = padjCutoff, FCcutoff = lfcCutoff,
                                  drawConnectors = TRUE, widthConnectors = 0.5)


  ggsave(file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_volcano_plot.png")),
         plot = volcanoPlot, width = 12, height = 10, units = "in", dpi = 300)

  res_filtered <- res[!is.na(res$padj) & res$padj < padjCutoff & abs(res$log2FoldChange) > lfcCutoff, ]

  write.csv(res, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Unfiltered_DE.csv")), quote = FALSE)
  write.csv(res_filtered, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Filtered_LFC_", lfcCutoff, "_FDR_", padjCutoff, ".csv")), quote = FALSE)

  # Heatmap of top DE genes
  if (nrow(res_filtered) > 0) {
    top_genes <- head(order(res_filtered$padj), min(50, nrow(res_filtered)))
    valid_ids <- res_filtered$IDS[top_genes]
    samples <- rownames(des[des[[condCol]] %in% c(numerator, denominator), , drop = FALSE])

    heat_mat <- ncts[valid_ids, samples, drop = FALSE]
    heat_mat <- as.matrix(heat_mat)
    heat_mat <- apply(heat_mat, 2, as.numeric)
    rownames(heat_mat) <- valid_ids
    colnames(heat_mat) <- samples

    # Ensure finite values and remove all-NA rows
    heat_mat[!is.finite(heat_mat)] <- NA
    heat_mat <- heat_mat[rowSums(!is.na(heat_mat)) > 0, , drop = FALSE]

    if (nrow(heat_mat) >= 2) {
      pheatmap(heat_mat, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_Genes_Heatmap.png")),
               height = 10, width = 14)
    } else if (nrow(heat_mat) == 1) {
      pheatmap(heat_mat, scale = "row", cluster_cols = TRUE, cluster_rows = FALSE,
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_1_Gene_Heatmap.png")),
               height = 4, width = 10)
    } else {
      warning(paste("⚠ Skipped heatmap for", numerator, "vs", denominator, ": No valid expression data."))
    }
  }

  message(paste0("✓ Completed: ", numerator, " vs ", denominator))
}

# Run comparisons
for (i in seq_len(nrow(comparisons))) {
  DElist(dds, comparisons$numerator[i], comparisons$denominator[i], opt$lfc, opt$padj, opt$condition)
}

