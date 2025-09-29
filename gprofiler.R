#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(gprofiler2)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
  library(readr)
})

# -----------------------------
# Command-line options
# -----------------------------
option_list <- list(
  make_option(c("-f", "--gene_file"), type = "character", help = "CSV file with DESeq2 results (must have 'IDS' column)"),
  make_option(c("-o", "--outdir"), type = "character", default = "enrichment_results", help = "Output directory [default: %default]"),
  make_option(c("-s", "--species"), type = "character", help = "g:Profiler species code (required)")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Check required inputs
if (is.null(opt$gene_file)) stop("❌ --gene_file is required. Use -f <filename.csv>")
if (is.null(opt$species)) stop("❌ --species is required. Use -s <species_code>")
if (!file.exists(opt$gene_file)) stop("❌ Gene file not found: ", opt$gene_file)
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

message("🧬 Using species: ", opt$species)

# -----------------------------
# Read gene IDs from CSV
# -----------------------------
df <- read_csv(opt$gene_file, show_col_types = FALSE)

if (!"IDS" %in% colnames(df)) stop("❌ CSV must contain an 'IDS' column with gene IDs.")

gene_ids <- unique(na.omit(df$IDS))
gene_ids <- gene_ids[nzchar(gene_ids)]  # remove empty
if (length(gene_ids) < 3) stop("❌ Not enough valid gene IDs in input file.")
message("Number of (unique, non-blank) gene IDs read: ", length(gene_ids))

# Save gene IDs to file
prefix <- "genes"
writeLines(gene_ids, file.path(opt$outdir, paste0(prefix, "_gene_ids.txt")))

# -----------------------------
# g:Profiler enrichment function
# -----------------------------
gprof <- function(IDS, prefix, title_prefix, source = c("GO", "KEGG"), species, outdir) {
  
  gp <- tryCatch({
    gost(IDS, organism = species, sources = source, multi_query = FALSE, evcodes = TRUE,
         user_threshold = 1, correction_method = "gSCS")
  }, error = function(e) {
    message("❌ g:Profiler error: ", e$message)
    return(NULL)
  })
  
  if (is.null(gp) || is.null(gp$result) || nrow(gp$result) == 0) {
    message(paste0("❌ No ", paste(source, collapse = ","), " terms found for ", prefix))
    return(NULL)
  }
  
  # Save all results table
  outfile_all <- file.path(outdir, paste0(prefix, "_", paste(source, collapse = "_"), "_enrichment_all.csv"))
  gp_result <- as.data.frame(lapply(gp$result, function(x) {
    if (is.list(x)) sapply(x, function(elem) paste(elem, collapse = ";")) else x
  }), stringsAsFactors = FALSE)
  write.csv(gp_result, file = outfile_all, row.names = FALSE)
  
  # Filter significant terms (p_value <= 0.05)
  sig <- subset(gp_result, as.logical(significant) & p_value <= 0.05)
  if (nrow(sig) == 0) {
    message(paste0("⚠️ No significant ", paste(source, collapse = ","), " terms found for ", prefix))
    return(NULL)
  }
  
  # Prepare data for enrichplot
  gp_mod <- sig[, c("query", "source", "term_id", "term_name", "p_value",
                    "query_size", "intersection_size", "term_size",
                    "effective_domain_size", "intersection")]
  gp_mod$GeneRatio <- paste0(gp_mod$intersection_size, "/", gp_mod$query_size)
  gp_mod$BgRatio   <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
  colnames(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust",
                        "query_size", "Count", "term_size", "effective_domain_size",
                        "geneID", "GeneRatio", "BgRatio")
  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)
  rownames(gp_mod) <- gp_mod$ID
  
  # Save significant results
  outfile_csv <- file.path(outdir, paste0(prefix, "_", paste(source, collapse = "_"), "_enrichment.csv"))
  write.csv(gp_mod, file = outfile_csv, row.names = FALSE)
  
  # Dotplot visualization
  enrich_obj <- new("enrichResult", result = gp_mod)
  p <- dotplot(enrich_obj, showCategory = min(15, nrow(gp_mod)), font.size = 14)
  if (!is.null(title_prefix)) p <- p + ggtitle(title_prefix)
  
  ggsave(file = file.path(outdir, paste0(prefix, "_", paste(source, collapse = "_"), "_dotplot.png")),
         plot = p, width = 12, height = 8, units = "in", dpi = 300)
  
  message(paste0("✓ Enrichment complete for ", paste(source, collapse = ","), " - ", prefix))
}

# -----------------------------
# Run enrichment
# -----------------------------
title_prefix <- paste("Enrichment for", prefix)

# GO enrichment
gprof(gene_ids, prefix = prefix, title_prefix = title_prefix, source = "GO",
      species = opt$species, outdir = opt$outdir)

# KEGG enrichment
gprof(gene_ids, prefix = prefix, title_prefix = title_prefix, source = "KEGG",
      species = opt$species, outdir = opt$outdir)
