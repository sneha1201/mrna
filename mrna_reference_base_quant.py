import os
import subprocess
import argparse
import pandas as pd

def run_command(command, log_file):
    """Run a shell command and log output."""
    print(f"Running: {command}")
    with open(log_file, "a") as log:
        log.write(f"\n[COMMAND] {command}\n")
        subprocess.run(command, shell=True, stdout=log, stderr=subprocess.STDOUT)

def create_dir(directory):
    """Create output directory if it doesn't exist."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def main():
    parser = argparse.ArgumentParser(description="RNA-Seq Analysis Pipeline")
    parser.add_argument("--sample_info", required=True, help="TSV file containing sample information (Sample, File Name, Pair)")
    parser.add_argument("--fastq_folder", required=True, help="Folder containing FASTQ files")
    parser.add_argument("--output_folder", required=True, help="Folder to store output results")
    parser.add_argument("--genome_index", required=True, help="HISAT2 genome index base name (path without file extension)")
    parser.add_argument("--annotation_file", required=True, help="GTF annotation file for featureCounts")
    parser.add_argument("--gene_attribute", required=True, help="Gene attribute to use in featureCounts (e.g., ID, gene_id, gene_name)")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads to use for processing")
    args = parser.parse_args()

    sample_info = pd.read_csv(args.sample_info, sep="\t")
    fastq_folder = os.path.abspath(args.fastq_folder)
    output_folder = os.path.abspath(args.output_folder)
    genome_index = os.path.abspath(args.genome_index)
    annotation_file = os.path.abspath(args.annotation_file)
    gene_attribute = args.gene_attribute
    threads = args.threads

    log_file = os.path.join(output_folder, "pipeline.log")
    with open(log_file, "w") as log:
        log.write("RNA-Seq Analysis Pipeline Log\n")

    # Define output directories
    qc_folder = os.path.join(output_folder, "1_QC")
    qc_before_folder = os.path.join(qc_folder, "qc_before")
    qc_after_folder = os.path.join(qc_folder, "qc_after")
    trimmed_folder = os.path.join(output_folder, "2_Trimmed")
    alignment_folder = os.path.join(output_folder, "3_Alignment")
    counts_folder = os.path.join(output_folder, "4_Quantification")

    for folder in [qc_folder, qc_before_folder, qc_after_folder, trimmed_folder, alignment_folder, counts_folder]:
        create_dir(folder)

    # FastQC (Pre-Processing QC)
    run_command(f"fastqc {fastq_folder}/*.fastq.gz -o {qc_before_folder} -t {threads}", log_file)
    
    # SeqKit Stats (Pre-Processing QC)
    run_command(f"seqkit stats -a {fastq_folder}/*.fastq.gz -o {qc_before_folder}/stats.tsv", log_file)

    # MultiQC for Pre-Processing QC
    run_command(f"multiqc {qc_before_folder} -o {qc_before_folder}", log_file)

    # Fastp (Trimming and QC)
    for sample in sample_info["Sample"].unique():
        sample_files = sample_info[sample_info["Sample"] == sample]
        
        pair1 = sample_files[sample_files["Pair"] == "Pair1"]["File Name"].values[0]
        pair2 = sample_files[sample_files["Pair"] == "Pair2"]["File Name"].values[0]
        
        input1 = os.path.join(fastq_folder, pair1)
        input2 = os.path.join(fastq_folder, pair2)
        output1 = os.path.join(trimmed_folder, f"{sample}_1_trimmed.fastq.gz")
        output2 = os.path.join(trimmed_folder, f"{sample}_2_trimmed.fastq.gz")
        html_report = os.path.join(trimmed_folder, f"{sample}_fastp.html")
        json_report = os.path.join(trimmed_folder, f"{sample}_fastp.json")
        
        run_command(
            f"fastp --in1 {input1} --in2 {input2} --out1 {output1} --out2 {output2} "
            f"--detect_adapter_for_pe -w {threads} -h {html_report} -j {json_report}",
            log_file
        )

    # FastQC (Post-Processing QC)
    run_command(f"fastqc {trimmed_folder}/*.fastq.gz -o {qc_after_folder} -t {threads}", log_file)
    
    # SeqKit Stats (Post-Processing QC)
    run_command(f"seqkit stats -a {trimmed_folder}/*.fastq.gz -o {qc_after_folder}/stats.tsv", log_file)
    
    # MultiQC for Post-Processing QC
    run_command(f"multiqc {qc_after_folder} -o {qc_after_folder}", log_file)

    # HISAT2 Alignment
    for sample in sample_info["Sample"].unique():
        left_reads = os.path.join(trimmed_folder, f"{sample}_1_trimmed.fastq.gz")
        right_reads = os.path.join(trimmed_folder, f"{sample}_2_trimmed.fastq.gz")
        sorted_bam = os.path.join(alignment_folder, f"{sample}_sorted.bam")
        
        run_command(
            f"hisat2 -p {threads} -x {genome_index} -1 {left_reads} -2 {right_reads} | "
            f"samtools view -bS | samtools sort -o {sorted_bam}",
            log_file
        )
        run_command(f"samtools index {sorted_bam}", log_file)

    # BAM Stats with Sample Header
    bam_stats_file = os.path.join(alignment_folder, "bam_stats.txt")
    with open(bam_stats_file, "w") as stats_out:
        for sample in sample_info["Sample"].unique():
            sorted_bam = os.path.join(alignment_folder, f"{sample}_sorted.bam")
            stats_out.write(f"\n=== Sample: {sample} ===\n")
            subprocess.run(
                f"samtools flagstat {sorted_bam}",
                shell=True,
                stdout=stats_out,
                stderr=subprocess.STDOUT
            )

    # FeatureCounts Quantification
    count_output = os.path.join(counts_folder, "feature_counts.txt")
    bam_files = " ".join([
        os.path.join(alignment_folder, f"{sample}_sorted.bam")
        for sample in sample_info["Sample"].unique()
    ])
    
    run_command(
        f"featureCounts -T {threads} -p -t gene -g {gene_attribute} "
        f"-a {annotation_file} -o {count_output} {bam_files}",
        log_file
    )

if __name__ == "__main__":
    main()

