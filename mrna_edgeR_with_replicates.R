#!/usr/bin/env Rscript

# Load required libraries
suppressPackageStartupMessages({
  library(edgeR)
  library(optparse)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
})

# Define command-line options
option_list <- list(
  make_option(c("-c", "--counts"), type = "character", help = "Path to raw counts matrix"),
  make_option(c("-m", "--metadata"), type = "character", help = "Path to metadata file"),
  make_option(c("-C", "--comparisons"), type = "character", help = "Path to comparisons file (CSV with numerator,denominator)"),
  make_option(c("-l", "--lfc"), type = "numeric", default = 1, help = "Log2 fold change cutoff [default: %default]"),
  make_option(c("-p", "--pvalue"), type = "numeric", default = 0.05, help = "Raw p-value cutoff [default: %default]"),
  make_option(c("-s", "--condition"), type = "character", default = "Condition", help = "Condition column in metadata [default: %default]"),
  make_option(c("-o", "--outdir"), type = "character", default = "edgeR_results", help = "Output directory [default: %default]"),
  make_option(c("-r", "--rdata"), type = "character", default = "edgeR_workspace.RData", help = "File to save R environment")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Load data
cts <- read.table(opt$counts, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
des <- read.table(opt$metadata, header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
cts <- cts[, rownames(des), drop = FALSE]

# Match column names
colnames(cts) <- trimws(colnames(cts))
rownames(des) <- trimws(rownames(des))
group <- factor(des[[opt$condition]])

# Create DGEList
dge <- DGEList(counts = cts, group = group)
dge <- calcNormFactors(dge)

# Save normalized counts
norm_counts <- cpm(dge, normalized.lib.sizes = TRUE, log = FALSE)
write.table(norm_counts, file = file.path(opt$outdir, "Normalized_counts.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# Filter low-expressed genes
keep <- filterByExpr(dge)
dge <- dge[keep,, keep.lib.sizes = FALSE]

# Estimate dispersions
dge <- estimateDisp(dge)

# Load comparisons
comparisons <- read.csv(opt$comparisons, header = TRUE, stringsAsFactors = FALSE)
if (!all(c("numerator", "denominator") %in% colnames(comparisons))) {
  stop("❌ Comparison file must have two columns: numerator, denominator")
}

# DEG function
edgeR_DE <- function(dge, group_info, numerator, denominator, outdir, lfc_cutoff, pval_cutoff) {
  samples_idx <- group_info %in% c(numerator, denominator)
  dge_sub <- dge[, samples_idx]
  dge_sub$samples$group <- droplevels(group_info[samples_idx])
  dge_sub <- estimateDisp(dge_sub)

  res <- exactTest(dge_sub, pair = c(denominator, numerator))
  tab <- topTags(res, n = Inf, sort.by = "PValue")$table
  tab$GeneID <- rownames(tab)

  # Label regulation based on raw p-value
  tab$Regulation <- "Not Significant"
  tab$Regulation[tab$logFC > lfc_cutoff & tab$PValue < pval_cutoff] <- "Upregulated"
  tab$Regulation[tab$logFC < -lfc_cutoff & tab$PValue < pval_cutoff] <- "Downregulated"

  # Save unfiltered result
  write.csv(tab, file = file.path(outdir, paste0(numerator, "_vs_", denominator, "_Unfiltered.csv")), row.names = FALSE)

  # Save filtered result
  filtered <- tab[abs(tab$logFC) > lfc_cutoff & tab$PValue < pval_cutoff, ]
  write.csv(filtered, file = file.path(outdir, paste0(numerator, "_vs_", denominator, "_Filtered_LFC_", lfc_cutoff, "_pvalue_", pval_cutoff, ".csv")), row.names = FALSE)

  # Volcano plot
  volcano <- EnhancedVolcano(tab,
                              lab = tab$GeneID,
                              x = "logFC",
                              y = "PValue",
                              title = paste(numerator, "vs", denominator),
                              pCutoff = pval_cutoff,
                              FCcutoff = lfc_cutoff,
                              subtitle = NULL,
                              caption = "",
                              colAlpha = 0.4,
                              drawConnectors = TRUE)
  ggsave(file.path(outdir, paste0(numerator, "_vs_", denominator, "_Volcano.png")),
         plot = volcano, width = 10, height = 8)

  # Heatmap
  top_genes <- head(filtered$GeneID[order(filtered$PValue)], min(50, nrow(filtered)))
  if (length(top_genes) >= 1) {
    heat_mat <- norm_counts[top_genes, colnames(dge_sub$counts), drop = FALSE]
    heat_mat <- heat_mat[rowSums(is.finite(heat_mat)) > 0, , drop = FALSE]

    if (nrow(heat_mat) >= 1) {
      pheatmap(heat_mat,
               scale = "row",
               cluster_rows = TRUE,
               cluster_cols = TRUE,
               color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(255),
               fontsize_row = 8,
               main = paste(numerator, "vs", denominator, "- Top DE Genes"),
               filename = file.path(outdir, paste0(numerator, "_vs_", denominator, "_Heatmap.png")),
               width = 10, height = 10)
    }
  }

  message(paste("✓ Completed:", numerator, "vs", denominator))
}

# Loop through comparisons
for (i in seq_len(nrow(comparisons))) {
  edgeR_DE(dge, group, comparisons$numerator[i], comparisons$denominator[i], opt$outdir, opt$lfc, opt$pvalue)
}

# Save session
save.image(file = opt$rdata)

