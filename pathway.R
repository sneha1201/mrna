#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(gprofiler2)
  library(enrichplot)
  library(ggplot2)
})

# Define command-line options
option_list <- list(
  make_option(c("-f", "--filtered"), type = "character", help = "Path to filtered DE gene list (CSV with 'IDS' column)"),
  make_option(c("-p", "--prefix"), type = "character", help = "Prefix for output filenames"),
  make_option(c("-o", "--outdir"), type = "character", default = "enrichment_results", help = "Directory to save enrichment results [default: %default]"),
  make_option(c("-s", "--species"), type = "character", default = "hsapiens", help = "Organism for gProfiler (e.g., hsapiens, mmusculus) [default: %default]"),
  make_option(c("--split_ids"), action = "store_true", default = FALSE, help = "Split 'IDS' on underscore (_) and use only the first part [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# Read filtered DE genes
filtered <- read.csv(opt$filtered)

if (!"IDS" %in% colnames(filtered)) {
  stop("❌ The filtered file must contain an 'IDS' column.")
}

# Process IDS column
if (opt$split_ids) {
  gene_ids <- sub("_.*", "", filtered$IDS)  # e.g. ENSG0000012345
  gene_names <- sub(".*_", "", filtered$IDS)  # e.g. GENE1
} else {
  gene_ids <- filtered$IDS
  gene_names <- filtered$IDS
}

# Create a mapping vector: gene_id -> gene_name
id_name_map <- setNames(gene_names, gene_ids)

# Function to run gProfiler enrichment
gprof <- function(IDS, prefix, source = c("GO", "KEGG"), species = "hsapiens") {
  gp <- gost(IDS, organism = species, sources = source,
             multi_query = FALSE, evcodes = TRUE, user_threshold = 0.05,
             correction_method = "gSCS")

  if (is.null(gp)) {
    message(paste0("❌ No enriched ", source, " terms found."))
    return(NULL)
  }

  gp_mod <- gp$result[, c("query", "source", "term_id",
                          "term_name", "p_value", "query_size",
                          "intersection_size", "term_size",
                          "effective_domain_size", "intersection")]

  gp_mod$GeneRatio <- paste0(gp_mod$intersection_size, "/", gp_mod$query_size)
  gp_mod$BgRatio   <- paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)

  colnames(gp_mod) <- c("Cluster", "Category", "ID", "Description", "p.adjust",
                        "query_size", "Count", "term_size", "effective_domain_size",
                        "geneID", "GeneRatio", "BgRatio")

  gp_mod$geneID <- gsub(",", "/", gp_mod$geneID)

  # Add gene names based on ID mapping
  gp_mod$geneName <- sapply(strsplit(gp_mod$geneID, "/"), function(ids) {
    ids_clean <- sub("\\..*", "", ids)  # Remove version numbers
    paste(unique(na.omit(id_name_map[ids_clean])), collapse = "/")
  })

  rownames(gp_mod) <- gp_mod$ID

  # Save table
  outfile_csv <- file.path(opt$outdir, paste0(prefix, "_", source, "_enrichment.csv"))
  write.csv(gp_mod, file = outfile_csv, row.names = FALSE)

  # Enrichplot object
  enrich_obj <- new("enrichResult", result = gp_mod)

  # Save dotplot
  p <- enrichplot::dotplot(enrich_obj, showCategory = 15, font.size = 16) +
    ggtitle(paste(source, "Enrichment -", prefix))

  png(file = file.path(opt$outdir, paste0(prefix, "_", source, "_dotplot.png")),
      width = 15, height = 10, units = "in", res = 300)
  print(p)
  dev.off()

  message(paste0("✓ Enrichment complete for ", source))
  return(list(gprofiler_object = gp,
              enrichplot_object = enrich_obj,
              enrichment_dotplot = p))
}

# Run enrichment
go <- gprof(gene_ids, opt$prefix, source = "GO", species = opt$species)
kegg <- gprof(gene_ids, opt$prefix, source = "KEGG", species = opt$species)

