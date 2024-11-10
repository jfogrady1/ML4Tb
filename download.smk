# Snakefile
# modules
import os
import glob
import numpy
import pandas as pd
from pathlib import Path

# Load sample metadata
samples = pd.read_csv("./data/kirsten/kirsten_samples.csv", sep="\t")
samples_wiarda = pd.read_csv("./data/wiarda/wiarda_samples.csv", sep = "\t")

# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))
sample_mapping_wiarda = dict(zip(samples_wiarda["Sample_SRR"], samples_wiarda["Animal_Code"]))
# Set the output directory path
output_dir_kirsten = "/home/workspace/jogrady/ML4TB/data/kirsten/individual/"
output_dir_wiarda = "/home/workspace/jogrady/ML4TB/data/wiarda/"


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


sample_ids_wiarda = [
    "864_Infected_W0", "851_Infected_W0", "865_Infected_W0", "867_Infected_W0", 
    "858_Infected_W0", "846_Control_W0", "849_Control_W0", "855_Control_W0", 
    "856_Control_W0", "863_Control_W0", "847_Infected_W0", "861_Infected_W0", 
    "869_Infected_W0", "846_Control_W4", "849_Control_W4", "847_Infected_W4", 
    "861_Infected_W4", "869_Infected_W4", "864_Infected_W4", "851_Infected_W4", 
    "865_Infected_W4", "855_Control_W4", "867_Infected_W4", "858_Infected_W4", 
    "856_Control_W4", "863_Control_W4", "847_Infected_W10", "861_Infected_W10", 
    "869_Infected_W10", "864_Infected_W10", "851_Infected_W10", "865_Infected_W10", 
    "867_Infected_W10", "858_Infected_W10", "846_Control_W10", "849_Control_W10", 
    "855_Control_W10", "856_Control_W10", "863_Control_W10"
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
         expand(
            f"{output_dir_wiarda}{{sample_code}}_R{{run}}.fastq.gz",
            sample_code=samples_wiarda["Animal_Code"],
            run=[1, 2]
        ),
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/multiqc_report.html',
        expand("/home/workspace/jogrady/ML4TB/data/kirsten/combined/{sample_id}_R1.fastq.gz", sample_id = sample_ids, lane = lanes),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids, N=[1, 2]),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids_wiarda, N=[1, 2]),
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids, N = [1,2]),
         expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids_wiarda, N = [1,2]),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids, N=[1, 2]),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids_wiarda, N=[1, 2]),
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/multiqc_report.html',
        directory("/home/workspace/jogrady/ML4TB/data/kirsten/star-genome/"),
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/multiqc_report.html',
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids),
        expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids_wiarda),
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/gene_counts.txt',
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt',
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt'
        
        

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
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids, N = (1,2))
    output:
        report='/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/fastqc/
        """

rule trimming:
    input:
        untrimmed_reads = lambda wildcards: [f"/home/workspace/jogrady/ML4TB/data/kirsten/combined/{wildcards.sample_id}_R1.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/data/kirsten/combined/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "/home/workspace/jogrady/ML4TB/data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """


rule fastqc_trimmed:
    input:
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/'


rule mutltiqc_trim:
    input:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids, N = (1,2))
    output:
        report='/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/fastqc/
        """


# Build this for 100 bps
rule build_reference:
    input:
        fa="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf" 
    output:
        STAR_dir = directory("/home/workspace/jogrady/ML4TB/data/kirsten/star-genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir /home/workspace/jogrady/ML4TB/data/kirsten/star-genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir /home/workspace/jogrady/ML4TB/data/kirsten/star-genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99
        ''' 



rule Alignment:
    input:
        genome = "/home/workspace/jogrady/ML4TB/data/kirsten/star-genome/",
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Log.final.out',
        interlog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''


rule featureCounts:
    input:
        bam = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC:
    input:
        count_matrix = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/gene_counts_temp.txt',
        cleaned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''

#################################
## Wiarda
#################################
rule download_paired_end_rna_seq_wiarda:
    output:
        r1=f"{output_dir_wiarda}{{sample_code}}_R1.fastq.gz",
        r2=f"{output_dir_wiarda}{{sample_code}}_R2.fastq.gz"
    params:
        ena_run=lambda wildcards: samples_wiarda[samples_wiarda["Animal_Code"] == wildcards.sample_code]["Sample_SRR"].values[0]
    conda:
        "/home/workspace/jogrady/ML4TB/envs/fastqdl.yml"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_wiarda}

        # Rename files to match the sample code
        mv {output_dir_wiarda}/{params.ena_run}_1.fastq.gz {output.r1}
        mv {output_dir_wiarda}/{params.ena_run}_2.fastq.gz {output.r2}
        """


rule fastqc_wiarda:
    input:
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/data/wiarda/{wildcards.sample_id}_R1.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/data/wiarda/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/'


rule mutltiqc_wiarda:
    input:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids_wiarda, N = (1,2))
    output:
        report='/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/fastqc/
        """

rule trimming_wiarda:
    input:
        untrimmed_reads = lambda wildcards: [f"/home/workspace/jogrady/ML4TB/data/wiarda/{wildcards.sample_id}_R1.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/data/wiarda/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "/home/workspace/jogrady/ML4TB/data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """
rule fastqc_trimmed_wiarda:
    input:
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o /home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/'

rule mutltiqc_trimmed_wiarda:
    input:
        reads = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids_wiarda, N = (1,2))
    output:
        report='/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o /home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/fastqc/
        """

rule Alignment_wiarda:
    input:
        genome = "/home/workspace/jogrady/ML4TB/data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Log.final.out',
        interlog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Log.progress.out',
        initiallog = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix /home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_wiarda:
    input:
        bam = expand('/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_wiarda),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC_wiarda:
    input:
        count_matrix = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/gene_counts_temp.txt',
        cleaned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''


