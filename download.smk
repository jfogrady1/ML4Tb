# Snakefile

import pandas as pd

# Load sample metadata
samples = pd.read_csv("./data/kirsten_samples.csv", sep="\t")

# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))

# Set the output directory path
output_dir_kirsten = "/home/workspace/jogrady/ML4TB/data/kirsten/individual/"

# Rule to specify the final target files for all samples
rule all:
    input:
        expand(
            f"{output_dir_kirsten}{{sample_code}}_R{{run}}.fastq.gz",
            sample_code=samples["sample_code"],
            run=[1, 2]
        )

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

