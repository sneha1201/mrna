#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
})

option_list <- list(
  make_option(c("-c", "--counts"), type="character", help="Path to counts matrix"),
  make_option(c("-m", "--metadata"), type="character", help="Path to metadata file"),
  make_option(c("-C", "--comparisons"), type="character", help="CSV file with two columns: numerator,denominator"),
  make_option(c("-l", "--lfc"), type="numeric", default=1, help="Log2 fold change cutoff [default: %default]"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.05, help="P-value cutoff [default: %default]"),
  make_option(c("-s", "--condition"), type="character", default="Condition", help="Condition column in metadata [default: %default]"),
  make_option(c("-o", "--outdir"), type="character", default="edgeR_results", help="Output directory [default: %default]"),
  make_option(c("-r", "--rdata"), type="character", default="edgeR_session.Rdata", help="RData file to save session [default: %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Load data
cts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
meta <- read.table(opt$metadata, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cts <- cts[, rownames(meta), drop = FALSE]

# Create DGEList
group <- factor(meta[[opt$condition]])
dge <- DGEList(counts = cts, group = group)
dge <- dge[filterByExpr(dge), , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Save normalized CPM
norm_cpm <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
write.csv(norm_cpm, file = file.path(opt$outdir, "Normalized_expression_values.csv"), quote = FALSE)

# Load comparisons
comps <- read.csv(opt$comparisons, header = TRUE)
if (!all(c("numerator", "denominator") %in% colnames(comps))) {
  stop("❌ Comparison file must contain 'numerator' and 'denominator' columns")
}

# Function to perform edgeR DE analysis
run_edgeR <- function(numerator, denominator) {
  message("Running: ", numerator, " vs ", denominator)

  samples <- rownames(meta[meta[[opt$condition]] %in% c(numerator, denominator), , drop = FALSE])
  sub_dge <- dge[, samples]
  sub_dge$samples$group <- factor(meta[samples, opt$condition])

  res <- exactTest(sub_dge, pair = c(denominator, numerator), dispersion = 0.01)
  tab <- topTags(res, n = Inf, sort.by = "PValue")$table
  tab$GeneID <- rownames(tab)

  # Add Regulation label based on raw p-value
  tab$Regulation <- "Not Significant"
  tab$Regulation[tab$logFC > opt$lfc & tab$PValue < opt$pvalue] <- "Upregulated"
  tab$Regulation[tab$logFC < -opt$lfc & tab$PValue < opt$pvalue] <- "Downregulated"

  # Save unfiltered results
  write.csv(tab, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Unfiltered_DE.csv")), row.names = FALSE)

  # Filtered results
  tab_filt <- tab[abs(tab$logFC) > opt$lfc & tab$PValue < opt$pvalue, ]
  write.csv(tab_filt, file = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_Filtered_LFC_", opt$lfc, "_pvalue_", opt$pvalue, ".csv")), row.names = FALSE)

  # Volcano plot
  vol <- EnhancedVolcano(tab, lab = tab$GeneID,
                         x = 'logFC', y = 'PValue',
                         title = paste(numerator, "VS", denominator),
                         pCutoff = opt$pvalue, FCcutoff = opt$lfc,
                         pointSize = 3, labSize = 3, colAlpha = 0.4,
                         drawConnectors = TRUE, subtitle = NULL, caption = "")
  ggsave(file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_volcano.png")),
         plot = vol, width = 10, height = 8)

#Heatmap
 if (nrow(tab_filt) > 0) {
  top <- head(tab_filt$GeneID[order(tab_filt$PValue)], min(50, nrow(tab_filt)))

  logcpm_mat <- cpm(sub_dge, log = TRUE)   # log2 CPM normalization
  mat_top <- logcpm_mat[top, , drop = FALSE]

  # Remove rows with NA/Inf
  mat_top_clean <- mat_top[apply(mat_top, 1, function(x) all(is.finite(x))), , drop = FALSE]

  if (nrow(mat_top_clean) > 1) {  # Need at least 2 genes to cluster rows
    pheatmap(mat_top_clean,
             scale = "none",
             cluster_rows = TRUE,
             cluster_cols = FALSE,  # no clustering on samples
             show_rownames = TRUE,
             annotation_col = meta[samples, , drop = FALSE][, opt$condition, drop = FALSE],
             color = colorRampPalette(brewer.pal(11, "RdYlBu"))(100),
             filename = file.path(opt$outdir, paste0(numerator, "_VS_", denominator, "_heatmap_log2CPM.png")),
             main = paste0(numerator, " vs ", denominator, ": Top 50  DEGs"))
  } else {
    message("Not enough valid genes to plot heatmap for ", numerator, " vs ", denominator)
  }
}




  message("✓ Done: ", numerator, " vs ", denominator)
}

# Run all comparisons
for (i in seq_len(nrow(comps))) {
  run_edgeR(comps$numerator[i], comps$denominator[i])
}

save.image(file = opt$rdata)

