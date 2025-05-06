# mRNA
mRNA Reference-Based RNA-seq Pipeline
Overview
The pipeline includes the following key steps:

1)Quality Control and Quantification using mrna_reference_base_quant.py

2)Conversion to DESeq2-Compatible Format using deseq2_compatable_featurecountfile.py

3)Differential Expression Analysis using reference_base_deseq2.R

4)Pathway Enrichment Analysis (optional) using pathway.R

5)For Venn Diagram venn_diagram.R (optional)

 Requirements
Python 3.8+
R (≥ 4.0)

Tools used internally: fastp,fastp,seqkit,multiqc, HISAT2,samtools, featureCounts, MultiQC, DESeq2, gprofiler2,EnhancedVolcano,pheatmap

Input File Formats

sample_info.tsv

Sample	File Name	Pair

Control_1	Control_1_R1_001.fastq.gz	Pair1

Control_1	Control_1_R2_001.fastq.gz	Pair2

metadata.tsv

SampleID	Condition

Control_1	control

PTX_12_1	PTX_12

DOXO_15_2	DOXO_15

comparisions.csv

numerator,denominator
PTX_12,control

PTX_15,control

DOXO_12,control

DOXO_15,control

Pipeline Steps

FASTQ to Feature Counts

Run the main quantification script:

python mrna_reference_base_quant.py \
  --sample_info sample_info.tsv \
  --fastq_folder input_folder_having_fastq \
  --output_folder path_to_output_folder \
  --genome_index /path_to_genome_index_name \
  --annotation_file /path_to_gtf_or_gff \
  --gene_attribute ID_or_gene_id_for_feature_count \
  --threads 20

Convert to DESeq2-Compatible Format

python deseq2_compatable_featurecountfile.py \
  -i path_to_output_folder/feature_counts_gene_names.txt \
  -o path_to_output_folder/feature_counts_gene_names_deseq2.txt

Differential Expression Analysis with DESeq2

Rscript reference_base_deseq2.R \
  --counts path_to_output_folder/feature_counts_gene_names_deseq2.txt \
  --metadata 5_deseq2/metadata.tsv \
  --comparisons 5_deseq2/comparisions.csv \
  --outdir 5_deseq2/

 Optional: Pathway Enrichment (via g:Profiler)
 
Rscript pathway.R \
  --filtered 5_deseq2/DOXO_15_VS_control_Filtered_LFC_1_FDR_0.05.csv \
  -p DOXO_15_VS_control \
  --outdir 5_deseq2 \
  --species hsapiens

  -p the name for the output files
Ensure filenames and sample names are consistent across all files.
The gene_attribute must match the attribute in your GTF/GFF used by featureCounts.
Set --species appropriately for g:Profiler (hsapiens, mmusculus, etc.)

  Rscript venn_diagram.R \
  --normalized 5_Deseq2/Normalized_expression_values.csv \
  --metadata /5_Deseq2/metadata.tsv \
  --comparisons /5_Deseq2/comparision.csv 
  --outdir /5_Deseq2

 for sample file
 
>     ls *.fastq.gz | sort | while read file; do
>     sample=$(echo "$file" | cut -d'_' -f1)  # Extract sample name like C1, C2, etc.
>     pair=$(echo "$file" | grep -o '_R[12]_001' | sed 's/_R1_001/Pair1/; s/_R2_001/Pair2/')
>     echo -e "${sample}\t${file}\t${pair}"
> done > sample_file.tsv
  


