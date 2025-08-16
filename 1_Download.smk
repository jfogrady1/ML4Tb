# Snakefile
# modules
import os
import glob
import numpy
import pandas as pd
from pathlib import Path

# Load sample metadata
samples = pd.read_csv("./data/kirsten/kirsten_samples.csv", sep="\t")
samples_wiarda = pd.read_csv("./data/wiarda/wiarda_samples.csv", sep="\t")
samples_johnston = pd.read_csv("./data/johnston/johnston_samples.tsv", sep="\t")
samples_odonoghue = pd.read_csv("./data/odonoghue/odonoghue_samples.tsv", sep="\t")
samples_kirsten_pbl = pd.read_csv("./data/kirsten_pbl/kirsten_pbl_samples.csv", sep = "\t")
samples_abdelaal = pd.read_csv("./data/abdelaal/abdelaal_samples.csv", sep="\t")
samples_alonso = pd.read_csv("./data/alonso/alonso_samples.tsv", sep="\t")

# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))
sample_mapping_wiarda = dict(zip(samples_wiarda["Sample_SRR"], samples_wiarda["Animal_Code"]))
sample_mapping_johnston = dict(zip(samples_johnston["Sample_SRR"], samples_johnston["Animal_Code"]))
sample_mapping_odonoghue = dict(zip(samples_odonoghue["Sample_SRR"], samples_odonoghue["Animal_Code"]))
sample_mapping_abdelaal = dict(zip(samples_abdelaal["Sample_SRR"], samples_abdelaal["Run_Code"]))
sample_mapping_alonso = dict(zip(samples_alonso["Sample_SRR"], samples_alonso["Animal_Code"]))

# Set the output directory path
output_dir_kirsten = "data/kirsten/individual/"
output_dir_wiarda = "data/wiarda/"
output_dir_johnston = "data/johnston/"
output_dir_odonoghue = "data/odonoghue/"
output_dir_kirsten_pbl = "data/kirsten_pbl/"
output_dir_alonso = "data/alonso/"



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



sample_ids_abdelaal = [
    'Infected_1_20', 'Infected_1_8',
    'Infected_2_20', 'Infected_2_8',
    'Infected_3_20', 'Infected_3_8',
    'Infected_4_20', 'Infected_4_8',
    'Infected_5_20', 'Infected_5_8',
    'Infected_6_20', 'Infected_6_8',
    'Uninfected_1_20', 'Uninfected_1_8',
    'Uninfected_2_20', 'Uninfected_2_8',
    'Uninfected_3_20', 'Uninfected_3_8',
    'Uninfected_4_20', 'Uninfected_4_8',
    'Uninfected_5_20', 'Uninfected_5_8',
    'Uninfected_6_20', 'Uninfected_6_8'
]

alonso_run_codes = [
"CON_001", "CON_002", "CON_003", "INFEC_FOCAL_001", "INFEC_FOCAL_002", "INFEC_FOCAL_003",
"INFEC_FOCAL_004", "INFEC_FOCAL_005", "INFEC_DIFFUSE_001", "INFEC_DIFFUSE_002", "INFEC_DIFFUSE_003",
"INFEC_DIFFUSE_004", "INFEC_FOCAL_006", "INFEC_DIFFUSE_005"]

sample_ids_alonso = [
"CON_001", "CON_002", "CON_003", "INFEC_FOCAL_001", "INFEC_FOCAL_002", "INFEC_FOCAL_003",
"INFEC_FOCAL_004", "INFEC_FOCAL_005", "INFEC_DIFFUSE_001", "INFEC_DIFFUSE_002", "INFEC_DIFFUSE_003",
"INFEC_DIFFUSE_004", "INFEC_FOCAL_006", "INFEC_DIFFUSE_005"]


sample_ids_johnston = [
"001_BRSV_CON","002_BRSV_CON","003_BRSV_CON","004_BRSV_CON","005_BRSV_CON","006_BRSV_CON",
"007_BRSV_INFECTED","008_BRSV_INFECTED","009_BRSV_INFECTED","010_BRSV_INFECTED","011_BRSV_INFECTED","012_BRSV_INFECTED",
"013_BRSV_INFECTED","014_BRSV_INFECTED","015_BRSV_INFECTED","016_BRSV_INFECTED","017_BRSV_INFECTED","018_BRSV_INFECTED"]

sample_ids_odonoghue = ["001_BOHV1_INFECTED",
"002_BOHV1_INFECTED",
"003_BOHV1_INFECTED",
"004_BOHV1_INFECTED",
"005_BOHV1_INFECTED",
"006_BOHV1_INFECTED",
"007_BOHV1_INFECTED",
"008_BOHV1_INFECTED",
"009_BOHV1_INFECTED",
"010_BOHV1_INFECTED",
"011_BOHV1_INFECTED",
"012_BOHV1_INFECTED",
"013_BOHV1_CON",
"014_BOHV1_CON",
"015_BOHV1_CON",
"016_BOHV1_CON",
"017_BOHV1_CON",
"018_BOHV1_CON"]



abdelaal_run_codes = [
    "Infected_1_20_1", "Infected_1_20_2", "Infected_1_8_1", "Infected_1_8_2",
    "Infected_2_20_1", "Infected_2_20_2", "Infected_2_8_1", "Infected_2_8_2",
    "Infected_3_20_1", "Infected_3_20_2", "Infected_3_8_1", "Infected_3_8_2",
    "Infected_4_20_1", "Infected_4_20_2", "Infected_4_8_1", "Infected_4_8_2",
    "Infected_5_20_1", "Infected_5_20_2", "Infected_5_8_1", "Infected_5_8_2",
    "Infected_6_20_1", "Infected_6_20_2", "Infected_6_8_1", "Infected_6_8_2",
    "Uninfected_1_20_1", "Uninfected_1_20_2", "Uninfected_1_8_1", "Uninfected_1_8_2",
    "Uninfected_2_20_1", "Uninfected_2_20_2", "Uninfected_2_8_1", "Uninfected_2_8_2",
    "Uninfected_3_20_1", "Uninfected_3_20_2", "Uninfected_3_8_1", "Uninfected_3_8_2",
    "Uninfected_4_20_1", "Uninfected_4_20_2", "Uninfected_4_8_1", "Uninfected_4_8_2",
    "Uninfected_5_20_1", "Uninfected_5_20_2", "Uninfected_5_8_1", "Uninfected_5_8_2",
    "Uninfected_6_20_1", "Uninfected_6_20_2", "Uninfected_6_8_1", "Uninfected_6_8_2"
]



sample_ids_kirsten_pbl = ['A016_CON','A017_CON','A018_CON','A020_CON','A022_CON','A023_CON','A024_CON','A029_CON',
                          'A032_TB','A033_TB','A034_TB','A035_TB','A036_TB','A037_TB','A038_TB','A039_TB']


# functions
def get_fq1(wildcards):
    return sorted(glob.glob("data/kirsten/individual/" + wildcards.sample_id + "_*_R1.fastq.gz"))
def get_fq2(wildcards):
    return sorted(glob.glob("data/kirsten/individual/" + wildcards.sample_id + "_*_R2.fastq.gz"))


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
         expand(
            f"{output_dir_johnston}{{sample_code}}_R{{run}}.fastq.gz",
            sample_code=samples_johnston["Animal_Code"],
            run=[1, 2]
        ),
            expand(
            f"{output_dir_odonoghue}{{sample_code}}_R{{run}}.fastq.gz",
            sample_code=samples_odonoghue["Animal_Code"],
            run=[1, 2]
        ),
        'work/RNA_seq/wiarda/fastqc/multiqc_report.html',
        'work/RNA_seq/johnston/fastqc/multiqc_report.html',
        'work/RNA_seq/odonoghue/fastqc/multiqc_report.html',
        expand("data/alonso/{sample_code_alonso}.fastq.gz", sample_code_alonso = alonso_run_codes), # not paired end
        expand("data/kirsten/combined/{sample_id}_R1.fastq.gz", sample_id = sample_ids, lane = lanes),
        expand('work/RNA_seq/kirsten/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids, N=[1, 2]),
        expand('work/RNA_seq/wiarda/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids_wiarda, N=[1, 2]),
        expand('work/RNA_seq/johnston/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids_johnston, N=[1, 2]),
        expand('work/RNA_seq/odonoghue/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id=sample_ids_odonoghue, N=[1, 2]),
        'work/RNA_seq/kirsten/fastqc/multiqc_report.html',
        expand('work/RNA_seq/kirsten/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids, N = [1,2]),
         expand('work/RNA_seq/wiarda/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids_wiarda, N = [1,2]),
         expand('work/RNA_seq/johnston/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids_johnston, N = [1,2]),
         expand('work/RNA_seq/odonoghue/trimmed/{sample_id}_R{N}_trimmed.fastq.gz', sample_id = sample_ids_odonoghue, N = [1,2]),
        expand('work/RNA_seq/kirsten/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids, N=[1, 2]),
        expand('work/RNA_seq/wiarda/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids_wiarda, N=[1, 2]),
        expand('work/RNA_seq/johnston/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids_johnston, N=[1, 2]),
        expand('work/RNA_seq/odonoghue/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id=sample_ids_odonoghue, N=[1, 2]),
        'work/RNA_seq/kirsten/trimmed/fastqc/multiqc_report.html',
        directory("data/kirsten/star-genome/"),
        'work/RNA_seq/wiarda/trimmed/fastqc/multiqc_report.html',
        'work/RNA_seq/johnston/trimmed/fastqc/multiqc_report.html',
        'work/RNA_seq/odonoghue/trimmed/fastqc/multiqc_report.html',
        expand('work/RNA_seq/kirsten/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids),
        expand('work/RNA_seq/wiarda/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids_wiarda),
        expand('work/RNA_seq/johnston/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids_johnston),
        expand('work/RNA_seq/odonoghue/Alignment/{sample_id}_Log.final.out', sample_id = sample_ids_odonoghue),
        'work/RNA_seq/kirsten/Quantification/gene_counts.txt',
        'work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt',
        'work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt',
        'work/RNA_seq/johnston/Quantification/johnston_count_matrix_clean.txt',
        'work/RNA_seq/odonoghue/Quantification/odonoghue_count_matrix_clean.txt',
        expand("data/kirsten_pbl/{sample_code_kirsten_pbl}.fastq.gz", sample_code_kirsten_pbl = sample_ids_kirsten_pbl), # not paired end
        expand('work/RNA_seq/kirsten_pbl/fastqc/{sample_id}_fastqc.zip', sample_id=sample_ids_kirsten_pbl),
        'work/RNA_seq/kirsten_pbl/fastqc/multiqc_report.html',
        expand('work/RNA_seq/kirsten_pbl/trimmed/{sample_id}_trimmed.fastq.gz', sample_id = sample_ids_kirsten_pbl),
        expand('work/RNA_seq/kirsten_pbl/trimmed/fastqc/{sample_id}_trimmed_fastqc.zip', sample_id = sample_ids_kirsten_pbl),
        'work/RNA_seq/kirsten_pbl/trimmed/fastqc/multiqc_report.html',
        'work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt',
        expand('work/RNA_seq/alonso/trimmed/{sample_id}_trimmed.fastq.gz', sample_id = sample_ids_alonso),
        'work/RNA_seq/alonso/fastqc/multiqc_report.html',
        'work/RNA_seq/alonso/trimmed/fastqc/multiqc_report.html',
        'work/RNA_seq/alonso/Quantification/alonso_count_matrix_clean.txt',
        
        

# Rule to download paired-end files from ENA
# Note do not specify an input file as snakemake will give out cause it cant find a file that is not there
rule download_paired_end_rna_seq:
    output:
        r1=f"{output_dir_kirsten}{{sample_code}}_R1.fastq.gz",
        r2=f"{output_dir_kirsten}{{sample_code}}_R2.fastq.gz"
    params:
        ena_run=lambda wildcards: samples[samples["sample_code"] == wildcards.sample_code]["ena_run"].values[0]
    conda:
        "envs/fastqdl.yml"
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
        r1="data/kirsten/combined/{sample_id}_R1.fastq.gz",
        r2="data/kirsten/combined/{sample_id}_R2.fastq.gz"
    shell:
        """
        cat {input.reads_1} > {output.r1}
        cat {input.reads_2} > {output.r2}
        """


rule fastqc:
    input:
        reads = lambda wildcards:[f"data/kirsten/combined/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/kirsten/combined/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/kirsten/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/kirsten/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/kirsten/fastqc/'


rule mutltiqc:
    input:
        reads = expand('work/RNA_seq/kirsten/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids, N = (1,2))
    output:
        report='work/RNA_seq/kirsten/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/kirsten/fastqc/
        """

rule trimming:
    input:
        untrimmed_reads = lambda wildcards: [f"data/kirsten/combined/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/kirsten/combined/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('work/RNA_seq/kirsten/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/kirsten/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
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
        reads = lambda wildcards:[f"work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/kirsten/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/kirsten/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/kirsten/trimmed/fastqc/'


rule mutltiqc_trim:
    input:
        reads = expand('work/RNA_seq/kirsten/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids, N = (1,2))
    output:
        report='work/RNA_seq/kirsten/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/kirsten/trimmed/fastqc/
        """


# Build this for 100 bps
rule build_reference:
    input:
        fa="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf" 
    output:
        STAR_dir = directory("data/kirsten/star-genome/")  # Please provide the output files names from this step.
    threads: 20
    shell:
        '''
        mkdir data/kirsten/star-genome/ # STAR cannot make directory
        STAR-2.7.1a --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir data/kirsten/star-genome/ \
        --genomeFastaFiles {input.fa} \
        --sjdbGTFfile {input.gtf} \
        --sjdbOverhang 99
        ''' 



rule Alignment:
    input:
        genome = "data/kirsten/star-genome/",
        reads = lambda wildcards:[f"work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/kirsten/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/kirsten/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/kirsten/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/kirsten/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/kirsten/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/kirsten/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''


rule featureCounts:
    input:
        bam = expand('work/RNA_seq/kirsten/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/kirsten/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC:
    input:
        count_matrix = 'work/RNA_seq/kirsten/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/kirsten/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/kirsten/Alignment/##g' {output.count_matrix_temp}
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
        "envs/fastqdl.yml"
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
        reads = lambda wildcards:[f"data/wiarda/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/wiarda/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/wiarda/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/wiarda/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/wiarda/fastqc/'


rule mutltiqc_wiarda:
    input:
        reads = expand('work/RNA_seq/wiarda/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids_wiarda, N = (1,2))
    output:
        report='work/RNA_seq/wiarda/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/wiarda/fastqc/
        """

rule trimming_wiarda:
    input:
        untrimmed_reads = lambda wildcards: [f"data/wiarda/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/wiarda/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('work/RNA_seq/wiarda/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/wiarda/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
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
        reads = lambda wildcards:[f"work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/wiarda/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/wiarda/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/wiarda/trimmed/fastqc/'

rule mutltiqc_trimmed_wiarda:
    input:
        reads = expand('work/RNA_seq/wiarda/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids_wiarda, N = (1,2))
    output:
        report='work/RNA_seq/wiarda/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/wiarda/trimmed/fastqc/
        """

rule Alignment_wiarda:
    input:
        genome = "data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/wiarda/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/wiarda/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/wiarda/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/wiarda/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/wiarda/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/wiarda/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_wiarda:
    input:
        bam = expand('work/RNA_seq/wiarda/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_wiarda),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/wiarda/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC_wiarda:
    input:
        count_matrix = 'work/RNA_seq/wiarda/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/wiarda/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/wiarda/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''

########################################
### Abdelaal samples
###
########################################
rule download_paired_end_rna_seq_johnston:
    output:
        r1=f"{output_dir_johnston}{{sample_code}}_R1.fastq.gz",
        r2=f"{output_dir_johnston}{{sample_code}}_R2.fastq.gz"
    params:
        ena_run=lambda wildcards: samples_johnston[samples_johnston["Animal_Code"] == wildcards.sample_code]["Sample_SRR"].values[0]
    conda:
        "envs/fastqdl.yml"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_johnston}

        # Rename files to match the sample code
        mv {output_dir_johnston}/{params.ena_run}_1.fastq.gz {output.r1}
        mv {output_dir_johnston}/{params.ena_run}_2.fastq.gz {output.r2}
        """



rule fastqc_johnston:
    input:
        reads = lambda wildcards:[f"data/johnston/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/johnston/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/johnston/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/johnston/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/johnston/fastqc/'


rule mutltiqc_johnston:
    input:
        reads = expand('work/RNA_seq/johnston/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids_johnston, N = (1,2))
    output:
        report='work/RNA_seq/johnston/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/johnston/fastqc/
        """

rule trimming_johnston:
    input:
        untrimmed_reads = lambda wildcards: [f"data/johnston/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/johnston/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('work/RNA_seq/johnston/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/johnston/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """
rule fastqc_trimmed_johnston:
    input:
        reads = lambda wildcards:[f"work/RNA_seq/johnston/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/johnston/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/johnston/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/johnston/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/johnston/trimmed/fastqc/'

rule mutltiqc_trimmed_johnston:
    input:
        reads = expand('work/RNA_seq/johnston/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids_johnston, N = (1,2))
    output:
        report='work/RNA_seq/johnston/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/johnston/trimmed/fastqc/
        """

rule Alignment_johnston:
    input:
        genome = "data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"work/RNA_seq/johnston/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/johnston/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/johnston/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/johnston/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/johnston/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/johnston/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/johnston/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_johnston:
    input:
        bam = expand('work/RNA_seq/johnston/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_johnston),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/johnston/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC_johnston:
    input:
        count_matrix = 'work/RNA_seq/johnston/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/johnston/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/johnston/Quantification/johnston_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/johnston/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''


########################################
### Abdelaal samples
###
########################################
rule download_paired_end_rna_seq_odonoghue:
    output:
        r1=f"{output_dir_odonoghue}{{sample_code}}_R1.fastq.gz",
        r2=f"{output_dir_odonoghue}{{sample_code}}_R2.fastq.gz"
    params:
        ena_run=lambda wildcards: samples_odonoghue[samples_odonoghue["Animal_Code"] == wildcards.sample_code]["Sample_SRR"].values[0]
    conda:
        "envs/fastqdl.yml"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_odonoghue}

        # Rename files to match the sample code
        mv {output_dir_odonoghue}/{params.ena_run}_1.fastq.gz {output.r1}
        mv {output_dir_odonoghue}/{params.ena_run}_2.fastq.gz {output.r2}
        """



rule fastqc_odonoghue:
    input:
        reads = lambda wildcards:[f"data/odonoghue/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/odonoghue/{wildcards.sample_id}_R2.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/odonoghue/fastqc/{{sample_id}}_R{N}_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/odonoghue/fastqc/{{sample_id}}_R{N}_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/odonoghue/fastqc/'


rule mutltiqc_odonoghue:
    input:
        reads = expand('work/RNA_seq/odonoghue/fastqc/{sample_id}_R{N}_fastqc.zip', sample_id = sample_ids_odonoghue, N = (1,2))
    output:
        report='work/RNA_seq/odonoghue/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/odonoghue/fastqc/
        """

rule trimming_odonoghue:
    input:
        untrimmed_reads = lambda wildcards: [f"data/odonoghue/{wildcards.sample_id}_R1.fastq.gz",
                                  f"data/odonoghue/{wildcards.sample_id}_R2.fastq.gz"],
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads=expand('work/RNA_seq/odonoghue/trimmed/{{sample_id}}_R{N}_trimmed.fastq.gz', N = (1,2)),
        unpaired_reads=expand('work/RNA_seq/odonoghue/trimmed/{{sample_id}}_R{N}_unpaired.fastq.gz', N = (1,2))
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic PE -threads {threads} -phred33 {input.untrimmed_reads[0]} {input.untrimmed_reads[1]} {output.trimmed_reads[0]} {output.unpaired_reads[0]} {output.trimmed_reads[1]} {output.unpaired_reads[1]} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """
rule fastqc_trimmed_odonoghue:
    input:
        reads = lambda wildcards:[f"work/RNA_seq/odonoghue/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/odonoghue/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    output:
        reads = expand('work/RNA_seq/odonoghue/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.zip', N = (1,2)),
        html = expand('work/RNA_seq/odonoghue/trimmed/fastqc/{{sample_id}}_R{N}_trimmed_fastqc.html', N = (1,2))
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads[0]} {input.reads[1]} -t {threads} -o work/RNA_seq/odonoghue/trimmed/fastqc/'

rule mutltiqc_trimmed_odonoghue:
    input:
        reads = expand('work/RNA_seq/odonoghue/trimmed/fastqc/{sample_id}_R{N}_trimmed_fastqc.zip', sample_id = sample_ids_odonoghue, N = (1,2))
    output:
        report='work/RNA_seq/odonoghue/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/odonoghue/trimmed/fastqc/
        """

rule Alignment_odonoghue:
    input:
        genome = "data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"work/RNA_seq/odonoghue/trimmed/{wildcards.sample_id}_R1_trimmed.fastq.gz",
                                  f"work/RNA_seq/odonoghue/trimmed/{wildcards.sample_id}_R2_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/odonoghue/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/odonoghue/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/odonoghue/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/odonoghue/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/odonoghue/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_odonoghue:
    input:
        bam = expand('work/RNA_seq/odonoghue/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_odonoghue),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/odonoghue/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -B -p -C -R BAM -T {threads} -s 0 -t gene -g gene_id
        '''

rule cleanup_FC_odonoghue:
    input:
        count_matrix = 'work/RNA_seq/odonoghue/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/odonoghue/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/odonoghue/Quantification/odonoghue_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/odonoghue/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''





########################################
### Kirsten_PBL samples
###
########################################

rule download_single_end_kirsten_pbl:
    output:
        r1=f"{output_dir_kirsten_pbl}{{sample_code_kirsten_pbl}}.fastq.gz"
    params:
        ena_run= lambda wildcards: samples_kirsten_pbl[samples_kirsten_pbl["Run_Code"] == wildcards.sample_code_kirsten_pbl]["Sample_SRR"].values[0]
    conda:
        "envs/fastqdl.yml"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_kirsten_pbl}

        # Rename files to match the sample code
        mv data/kirsten_pbl/{params.ena_run}.fastq.gz {output.r1}
        """

rule fastqc_kirsten_pbl:
    input:
        reads = lambda wildcards: f"data/kirsten_pbl/{wildcards.sample_id}.fastq.gz"
    output:
        reads = 'work/RNA_seq/kirsten_pbl/fastqc/{sample_id}_fastqc.zip',
        html = 'work/RNA_seq/kirsten_pbl/fastqc/{sample_id}_fastqc.html'
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads} -t {threads} -o work/RNA_seq/kirsten_pbl/fastqc/'



rule mutltiqc_kirsten_pbl:
    input:
        reads = expand('work/RNA_seq/kirsten_pbl/fastqc/{sample_id}_fastqc.zip', sample_id = sample_ids_kirsten_pbl)
    output:
        report='work/RNA_seq/kirsten_pbl/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/kirsten_pbl/fastqc/
        """

rule trimming_kirsten_pbl:
    input:
        untrimmed_reads = "data/kirsten_pbl/{sample_id}.fastq.gz",
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads = 'work/RNA_seq/kirsten_pbl/trimmed/{sample_id}_trimmed.fastq.gz'
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic SE -threads {threads} -phred33 {input.untrimmed_reads} {output.trimmed_reads} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """

rule fastqc_kirsten_pbl_trimmed:
    input:
        reads = lambda wildcards: f"work/RNA_seq/kirsten_pbl/trimmed/{wildcards.sample_id}_trimmed.fastq.gz"
    output:
        reads = 'work/RNA_seq/kirsten_pbl/trimmed/fastqc/{sample_id}_trimmed_fastqc.zip',
        html = 'work/RNA_seq/kirsten_pbl/trimmed/fastqc/{sample_id}_trimmed_fastqc.html'
    threads:
        40
    resources:
        mem_mb = 16000
    shell:
        'fastqc {input.reads} -t {threads} -o work/RNA_seq/kirsten_pbl/trimmed/fastqc/'


rule mutltiqc_kirsten_pbl_trimmed:
    input:
        reads = expand('work/RNA_seq/kirsten_pbl/trimmed/fastqc/{sample_id}_trimmed_fastqc.zip', sample_id = sample_ids_kirsten_pbl)
    output:
        report='work/RNA_seq/kirsten_pbl/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/kirsten_pbl/trimmed/fastqc/
        """

rule Alignment_kirsten_pbl:
    input:
        genome = "data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"work/RNA_seq/kirsten_pbl/trimmed/{wildcards.sample_id}_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/kirsten_pbl/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/kirsten_pbl/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/kirsten_pbl/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/kirsten_pbl/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/kirsten_pbl/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_kirsten_pbl:
    input:
        bam = expand('work/RNA_seq/kirsten_pbl/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_kirsten_pbl),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/kirsten_pbl/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -T {threads} -s 0 -t gene -g gene_id
        '''


rule cleanup_FC_kirsten_pbl:
    input:
        count_matrix = 'work/RNA_seq/kirsten_pbl/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/kirsten_pbl/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/kirsten_pbl/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''



rule download_single_end_alonso:
    output:
        r1=f"{output_dir_alonso}{{sample_code_alonso}}.fastq.gz"
    params:
        ena_run= lambda wildcards: samples_alonso[samples_alonso["Animal_Code"] == wildcards.sample_code_alonso]["Sample_SRR"].values[0]
    conda:
        "envs/fastqdl.yml"
    #wildcard_constraints: 
     #   sample_code_alonso = "^[CI].*"
    shell:
        """
        # Download paired-end reads using fastq-dl
        fastq-dl -a {params.ena_run} --cpus 20 -o {output_dir_alonso}

        # Rename files to match the sample code
        mv data/alonso/{params.ena_run}.fastq.gz {output.r1}
        """

rule fastqc_alonso:
    input:
        reads = lambda wildcards: f"data/alonso/{wildcards.sample_id}.fastq.gz"
    output:
        reads = 'work/RNA_seq/alonso/fastqc/{sample_id}_fastqc.zip',
        html = 'work/RNA_seq/alonso/fastqc/{sample_id}_fastqc.html'
    threads:
        40
    resources:
        mem_mb = 4000
    shell:
        'fastqc {input.reads} -t {threads} -o work/RNA_seq/alonso/fastqc/'


rule mutltiqc_alonso:
    input:
        reads = expand('work/RNA_seq/alonso/fastqc/{sample_id}_fastqc.zip', sample_id = sample_ids_alonso)
    output:
        report='work/RNA_seq/alonso/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/alonso/fastqc/
        """

rule trimming_alonso:
    input:
        untrimmed_reads = "data/alonso/{sample_id}.fastq.gz",
        adapters = "data/adapters/Illumina_adpters.fa"
    output:
        trimmed_reads = 'work/RNA_seq/alonso/trimmed/{sample_id}_trimmed.fastq.gz'
    threads:
        50
    singularity:
        "docker://staphb/trimmomatic:latest"
    shell:
        """
        /Trimmomatic-0.39/trimmomatic SE -threads {threads} -phred33 {input.untrimmed_reads} {output.trimmed_reads} ILLUMINACLIP:{input.adapters}:2:30:10 TRAILING:30 MINLEN:36 
        """

rule fastqc_alonso_trimmed:
    input:
        reads = lambda wildcards: f"work/RNA_seq/alonso/trimmed/{wildcards.sample_id}_trimmed.fastq.gz"
    output:
        reads = 'work/RNA_seq/alonso/trimmed/fastqc/{sample_id}_trimmed_fastqc.zip',
        html = 'work/RNA_seq/alonso/trimmed/fastqc/{sample_id}_trimmed_fastqc.html'
    threads:
        40
    resources:
        mem_mb = 16000
    shell:
        'fastqc {input.reads} -t {threads} -o work/RNA_seq/alonso/trimmed/fastqc/'

rule mutltiqc_alonso_trimmed:
    input:
        reads = expand('work/RNA_seq/alonso/trimmed/fastqc/{sample_id}_trimmed_fastqc.zip', sample_id = sample_ids_alonso)
    output:
        report='work/RNA_seq/alonso/trimmed/fastqc/multiqc_report.html'
    shell:
        """
        multiqc {input} -f -o work/RNA_seq/alonso/trimmed/fastqc/
        """
    

rule Alignment_alonso:
    input:
        genome = "data/kirsten/star-genome/", # use the same
        reads = lambda wildcards:[f"work/RNA_seq/alonso/trimmed/{wildcards.sample_id}_trimmed.fastq.gz"]
    params:
        prefix = lambda wildcards: f'{wildcards.sample_id}'
    output:
        aligned = 'work/RNA_seq/alonso/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam',
        finallog = 'work/RNA_seq/alonso/Alignment/{sample_id}_Log.final.out',
        interlog = 'work/RNA_seq/alonso/Alignment/{sample_id}_Log.progress.out',
        initiallog = 'work/RNA_seq/alonso/Alignment/{sample_id}_Log.out'
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad LoadAndKeep --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads} --readFilesCommand gunzip -c \
        --outFileNamePrefix work/RNA_seq/alonso/Alignment/{params.prefix}_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''

rule featureCounts_alonso:
    input:
        bam = expand('work/RNA_seq/alonso/Alignment/{sample_id}_Aligned.sortedByCoord.out.bam', sample_id  = sample_ids_alonso),
        annotation="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf"
    output:
        count_matrix = 'work/RNA_seq/alonso/Quantification/gene_counts.txt'
    threads: 40
    shell:
        '''
        # use new version of feature counts
        featureCounts -a {input.annotation} -o {output.count_matrix} {input.bam} -T {threads} -s 0 -t gene -g gene_id
        '''


rule cleanup_FC_alonso:
    input:
        count_matrix = 'work/RNA_seq/alonso/Quantification/gene_counts.txt'
    output:
        count_matrix_temp =  'work/RNA_seq/alonso/Quantification/gene_counts_temp.txt',
        cleaned = 'work/RNA_seq/alonso/Quantification/alonso_count_matrix_clean.txt',
    shell:
        ''' 
        tail -n+2 {input.count_matrix} | cut -f 1,7-58  > {output.count_matrix_temp}
        sed -i 's#work/RNA_seq/alonso/Alignment/##g' {output.count_matrix_temp}
        sed -i 's/'"_Aligned\.sortedByCoord\.out\.bam"'//g' {output.count_matrix_temp} 
        cat {output.count_matrix_temp} > {output.cleaned} 
        '''