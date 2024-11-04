# Snakefile
# modules
import os
import glob
import numpy
import pandas as pd
from pathlib import Path

# Load sample metadata
samples = pd.read_csv("./data/kirsten_samples.csv", sep="\t")

# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))

# Set the output directory path
output_dir_kirsten = "/home/workspace/jogrady/ML4TB/data/kirsten/individual/"


sample_ids = [
    "A6511_W-1", "A6514_W-1", "A6520_W-1", "A6526_W-1", "A6635_W-1", 
    "A6636_W-1", "A6637_W-1", "A6644_W-1", "A6698_W-1", "A6511_W1", 
    "A6514_W1", "A6520_W1", "A6526_W1", "A6635_W1", "A6636_W1", 
    "A6637_W1", "A6644_W1", "A6698_W1", "A6511_W2", "A6514_W2", 
    "A6635_W2", "A6636_W2", "A6637_W2", "A6644_W2", "A6698_W2", 
    "A6511_W6", "A6514_W6", "A6520_W6", "A6526_W6", "A6635_W6", 
    "A6636_W6", "A6637_W6", "A6644_W6", "A6698_W6", "A6511_W10", 
    "A6514_W10", "A6520_W10", "A6526_W10", "A6635_W10", "A6636_W10", 
    "A6637_W10", "A6644_W10", "A6698_W10", "A6511_W12", "A6514_W12", 
    "A6520_W12", "A6526_W12", "A6635_W12", "A6636_W12", "A6637_W12", 
    "A6644_W12", "A6698_W12"
]

lanes = ["001","002","003","004","005"]


# functions
def get_fq1(wildcards):
    return sorted(glob.glob("/home/workspace/jogrady/ML4TB/data/kirsten/individual/" + wildcards.sample_id + "_*_R1.fastq.gz"))
def get_fq2(wildcards):
    return sorted(glob.glob("/home/workspace/jogrady/ML4TB/data/kirsten/individual/" + wildcards.sample_id + "_*_R2.fastq.gz"))


# Rule to specify the final target files for all samples
rule all:
    input:
        expand(
            f"{output_dir_kirsten}{{sample_code}}_R{{run}}.fastq.gz",
            sample_code=samples["sample_code"],
            run=[1, 2]
        ),
        expand("/home/workspace/jogrady/ML4TB/data/kirsten/combined/{sample_id}_R1.fastq.gz", sample_id = sample_ids, lane = lanes),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids, N=[1, 2]),
        #"/home/workspace/jogrady/ML4TB/work/RNA_seq/multiqc_report.html"

# Rule to download paired-end files from ENA
# Note do not specify an input file as snakemake will give out cause it cant find a file that is not there
rule download_paired_end_rna_seq:
    output:
        r1=f"{output_dir_kirsten}{{sample_code}}_R1.fastq.gz",
        r2=f"{output_dir_kirsten}{{sample_code}}_R2.fastq.gz"
    params:
        ena_run=lambda wildcards: samples[samples["sample_code"] == wildcards.sample_code]["ena_run"].values[0]
    conda:
        "/home/workspace/jogrady/ML4TB/envs/fastqdl.yml"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_kirsten}

        # Rename files to match the sample code
        mv {output_dir_kirsten}/{params.ena_run}_1.fastq.gz {output.r1}
        mv {output_dir_kirsten}/{params.ena_run}_2.fastq.gz {output.r2}
        """


# Rule to concatenate lanes for each sample
rule combineFiles:
    input:
        reads_1 = get_fq1,
        reads_2 = get_fq2
    output:
        r1="/home/workspace/jogrady/ML4TB/data/kirsten/combined/{sample_id}_R1.fastq.gz",
        r2="/home/workspace/jogrady/ML4TB/data/kirsten/combined/{sample_id}_R2.fastq.gz"
    shell:
        """
        cat {input.reads_1} > {output.r1}
        cat {input.reads_2} > {output.r2}
        """


rule fastqc:
    input:
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/data/kirsten/combined/{wildcards.sample_id}_R1.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/data/kirsten/combined/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/'


rule mutltiqc:
    input:
        reads = expand('/home/workspace/jogrady/ML4TB/data/kirsten/combined/{sample_id}_R{N}.fastq.gz', sample_id = sample_ids, N = (1,2))
    output:
        report='/home/workspace/jogrady/ML4TB/work/RNA_seq/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/ML4TB/work/RNA_seq/
        """