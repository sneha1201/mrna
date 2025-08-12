#!/usr/bin/env Rscript

# Suppress startup messages and load required packages
suppressPackageStartupMessages({
  library(DESeq2)
  library(EnhancedVolcano)
  library(pheatmap)
  library(optparse)
  library(cowplot)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Path to counts matrix"),
  make_option(c("-m", "--metadata"), type="character", help="Path to metadata file"),
  make_option(c("-C", "--comparisons"), type="character", help="Path to comparison file (2 columns: numerator, denominator)"),
  make_option(c("-l", "--lfc"), type="numeric", default=1, help="Log2 fold change cutoff [default: %default]"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05, help="P-value cutoff [default: %default]"),
  make_option(c("-s", "--condition"), type="character", default="Condition", help="Condition column in metadata [default: %default]"),
  make_option(c("-o", "--outdir"), type="character", default="DESeq2_results", help="Output directory [default: %default]"),
  make_option(c("-r", "--rdata"), type="character", default="data.Rdata", help="Output R session file [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Create output directory
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Load count matrix and metadata
cts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
des <- read.table(opt$metadata, header = TRUE, row.names = 1, check.names = FALSE)

# Clean whitespace
colnames(cts) <- trimws(colnames(cts))
rownames(des) <- trimws(rownames(des))

# Match and subset samples
cts <- cts[, rownames(des), drop = FALSE]

# Run DESeq2
dds <- DESeqDataSetFromMatrix(countData = cts, colData = des, design = as.formula(paste("~", opt$condition)))
dds <- DESeq(dds)

# Save dispersion plot
png(file.path(opt$outdir, "Dispersion_plot.png"), width = 1200, height = 1000, res = 150)
plotDispEsts(dds, main = "Dispersion Estimates")
dev.off()

# Export normalized counts
ncts <- as.data.frame(counts(dds, normalized = TRUE))
ncts$IDS <- rownames(ncts)
write.csv(ncts, file.path(opt$outdir, "Normalized_expression_values.csv"), quote = FALSE)

# Load comparisons
comparisons <- read.table(opt$comparisons, header = TRUE, sep = ",", stringsAsFactors = FALSE)
if (!all(c("numerator", "denominator") %in% colnames(comparisons))) {
  stop("❌ Comparison file must have two columns: numerator, denominator")
}

# Function for DE and plots
DElist <- function(deseqObject, numerator, denominator, lfcCutoff, pvalueCutoff, condCol) {
  res <- results(deseqObject,
                 contrast = c(condCol, numerator, denominator),
                 alpha = pvalueCutoff,
                 lfcThreshold = lfcCutoff,
                 altHypothesis = "greaterAbs",
                 pAdjustMethod = "none")

  res <- as.data.frame(res)
  res$IDS <- rownames(res)
  res$Regulation <- "Not Significant"
  res$Regulation[res$log2FoldChange > lfcCutoff & res$pvalue < pvalueCutoff] <- "Upregulated"
  res$Regulation[res$log2FoldChange < -lfcCutoff & res$pvalue < pvalueCutoff] <- "Downregulated"

  # Summary file
  summaryFile <- file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_summary.txt"))
  sink(summaryFile); print(summary(res)); sink()

  # Volcano plot
  lfc_max <- max(abs(res$log2FoldChange), na.rm = TRUE)
  lfc_limit <- ceiling(lfc_max)
  base_plot <- EnhancedVolcano(
    res,
    lab = res$IDS,
    x = 'log2FoldChange',
    y = 'pvalue',
    title = NULL,
    subtitle = NULL,
    caption = "",
    pCutoff = pvalueCutoff,
    FCcutoff = lfcCutoff,
    drawConnectors = FALSE,
    widthConnectors = 0.5,
    pointSize = 2.0,
    labSize = 3,
    colAlpha = 0.4,
    xlim = c(-lfc_limit, lfc_limit)
  )
  ggsave(file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_volcano_plot.png")),
         plot = base_plot, width = 12, height = 10, dpi = 300)

  # Filtered results
  res_filtered <- res[!is.na(res$pvalue) & res$pvalue < pvalueCutoff & abs(res$log2FoldChange) > lfcCutoff, ]
  write.csv(res, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Unfiltered_DE.csv")), quote = FALSE)
  write.csv(res_filtered, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Filtered_LFC_", lfcCutoff, "_pvalue_", pvalueCutoff, ".csv")), quote = FALSE)

  # Proceed only if filtered DE genes exist
  if (nrow(res_filtered) > 0) {
    top_genes <- head(order(res_filtered$pvalue), min(50, nrow(res_filtered)))
    valid_ids <- res_filtered$IDS[top_genes]

    # ==== Sample-wise heatmap ====
    samples <- rownames(des[des[[condCol]] %in% c(numerator, denominator), , drop = FALSE])
    heat_mat <- as.matrix(ncts[valid_ids, samples, drop = FALSE])
    heat_mat <- t(apply(heat_mat, 1, as.numeric))
    rownames(heat_mat) <- valid_ids
    colnames(heat_mat) <- samples
    heat_mat[!is.finite(heat_mat)] <- NA
    heat_mat <- heat_mat[rowSums(!is.na(heat_mat)) > 0, , drop = FALSE]

    if (nrow(heat_mat) >= 2) {
      pheatmap(heat_mat, scale = "row", cluster_cols = TRUE, cluster_rows = TRUE,
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_Genes_Sample_Heatmap.png")),
               height = 10, width = 14)
    } else if (nrow(heat_mat) == 1) {
      pheatmap(heat_mat, scale = "row", cluster_cols = TRUE, cluster_rows = FALSE,
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_1_Gene_Sample_Heatmap.png")),
               height = 4, width = 10)
    }

    # ==== Condition-wise heatmap ====
    samples_numerator <- rownames(des[des[[condCol]] == numerator, , drop = FALSE])
    samples_denominator <- rownames(des[des[[condCol]] == denominator, , drop = FALSE])
    heat_mat_raw <- ncts[valid_ids, c(samples_numerator, samples_denominator), drop = FALSE]
    heat_mat_raw <- t(apply(heat_mat_raw, 1, as.numeric))
    rownames(heat_mat_raw) <- valid_ids
    colnames(heat_mat_raw) <- c(samples_numerator, samples_denominator)

    # Compute average per condition
    condition_means <- data.frame(
      numerator_group = rowMeans(heat_mat_raw[, samples_numerator, drop = FALSE], na.rm = TRUE),
      denominator_group = rowMeans(heat_mat_raw[, samples_denominator, drop = FALSE], na.rm = TRUE)
    )
    colnames(condition_means) <- c(numerator, denominator)
    heat_mat_avg <- as.matrix(condition_means)
    heat_mat_avg <- heat_mat_avg[complete.cases(heat_mat_avg), , drop = FALSE]

    if (nrow(heat_mat_avg) >= 2) {
      pheatmap(heat_mat_avg, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
               angle_col = 45, fontsize_row = 10, fontsize_col = 12,
               main = paste(numerator, "vs", denominator),
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_Genes_Condition_Heatmap.png")),
               height = 10, width = 6)
    } else if (nrow(heat_mat_avg) == 1) {
      pheatmap(heat_mat_avg, scale = "row", cluster_cols = FALSE, cluster_rows = FALSE,
               angle_col = 45,
               filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Top_1_Gene_Condition_Heatmap.png")),
               height = 4, width = 5)
    }
  }

  message(paste0("✓ Completed: ", numerator, " vs ", denominator))
}

# Run all comparisons
for (i in seq_len(nrow(comparisons))) {
  DElist(dds, comparisons$numerator[i], comparisons$denominator[i], opt$lfc, opt$pvalue, opt$condition)
}

# Save R session
save.image(file = opt$rdata)

