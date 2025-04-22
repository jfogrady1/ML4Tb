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
samples_abdelaal = pd.read_csv("./data/abdelaal/abdelaal_samples.csv", sep="\t")
samples_kirsten_pbl = pd.read_csv("./data/kirsten_pbl/kirsten_pbl_samples.csv", sep = "\t")
samples_ogrady = pd.read_csv("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt", sep="\t")["Sample"].to_list()  
# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))
sample_mapping_wiarda = dict(zip(samples_wiarda["Sample_SRR"], samples_wiarda["Animal_Code"]))
sample_mapping_abdelaal = dict(zip(samples_abdelaal["Sample_SRR"], samples_abdelaal["Run_Code"]))

# Set the output directory path
output_dir_kirsten = "/home/workspace/jogrady/ML4TB/data/kirsten/individual/"
output_dir_wiarda = "/home/workspace/jogrady/ML4TB/data/wiarda/"
output_dir_abdelaal = "/home/workspace/jogrady/ML4TB/data/abdelaal/individual/"
output_dir_kirsten_pbl = "/home/workspace/jogrady/ML4TB/data/kirsten_pbl/"




Kirsten_sample_ids = [
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



# Extract unique IDs to merge
unique_kirsten_ids = set([s.split("_")[0] for s in Kirsten_sample_ids])
unique_wiarda_ids = set([s.split("_")[0] for s in sample_ids_wiarda])

print(unique_kirsten_ids)
rule all:
    input:
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/{sample_ogrady}.filtered.renamed.vcf.gz", sample_ogrady = samples_ogrady),
        "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.filtered.vcf.gz",
        "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/DNA_ref_alt.txt",
        '/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_filtered_Aligned.sortedByCoord.out.bam',
        "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/call/selected.T064.filtered.vcf.gz"


rule rename_vcf_sample_ogrady:
    input:
        vcf="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/reviewer/RNA_seq/variant_call/call/selected.{sample_ogrady}.filtered.vcf.gz"
    output:
        renamed_vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/{sample_ogrady}.filtered.renamed.vcf.gz"
    params:
        sample="{sample_ogrady}"
    shell:
        "echo '{params.sample}' > {params.sample}.txt && "
        "bcftools reheader -s {params.sample}.txt {input.vcf} -o {output.renamed_vcf} && tabix -p vcf {output.renamed_vcf} && "
        "rm {params.sample}.txt"

rule merge_vcfs_ogrady:
    input:
        vcfs = expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/{sample_ogrady}.filtered.renamed.vcf.gz", sample_ogrady = samples_ogrady)

    output:
        merged_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.vcf.gz"
    shell:
        """
        bcftools merge {input.vcfs} -O z -o {output.merged_vcf} && tabix -p vcf {output.merged_vcf}
        """

rule filter_merged_vcfs_ogrady:
    input:
        merged_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.vcf.gz"
    output:
        filtered_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.merged_vcf} -Oz | bcftools view -q 0.01:minor -Oz | bcftools view -e 'HWE < 0.00001 || F_MISSING > 0.2 || N_MISSING > 0.05' -Oz -o {output.filtered_vcf} && tabix -p vcf {output.filtered_vcf}
        """

rule extract_pos:
    input:
        merged_vcf ="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.vcf.gz",
        imputed_vcf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/imputed/SNP_Data_ALL_CHR.Renamed.IMPUTED.FILTERED.dose.vcf.gz"
    output:
        pos_ref_alt_RNA = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/RNA_ref_alt.txt",
        pos_ref_alt_impute = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/DNA_ref_alt.txt"
    
    shell:
        """
        bcftools query -f '%CHROM %POS %REF %ALT\n' {input.merged_vcf} > {output.pos_ref_alt_RNA}
        bcftools query -f '%CHROM %POS %REF %ALT\n' {input.imputed_vcf} > {output.pos_ref_alt_impute}
        """
###################

# WASP filtering
##################


rule Alignment:
    input:
        genome = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/star-genome/",
        reads = lambda wildcards: expand(f'/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/T064_{{N}}.fq.gz', N=(1,2)),
        vcf =  "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/Dutta_WGS_panel_filtered_WASP.vcf"
    params:
        prefix = "T064_"
    output:
        aligned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_Aligned.sortedByCoord.out.bam',
    threads: 40
    shell:
        '''
        STAR-2.7.1a  --genomeLoad NoSharedMemory --genomeDir {input.genome} --runThreadN {threads} \
        --readFilesIn {input.reads[0]} {input.reads[1]} --readFilesCommand gunzip -c \
        --waspOutputMode SAMtag \
        --twopassMode Basic \
        --varVCFfile {input.vcf} \
        --outFileNamePrefix /home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_ --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000
        '''
        
rule samtools_filter:
  input:
    aligned = '/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_Aligned.sortedByCoord.out.bam',
    samtools = "/home/workspace/jogrady/heQTL/software/samtools-1.15.1/samtools"
  output:
    filtered ='/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_filtered_Aligned.sortedByCoord.out.bam'
  
  shell:
    """
     {input.samtools} view -h -b -e '[vW]==1' -q 255 --threads 25 -o {output.filtered} {input.aligned}
    """
        
# 1. Add read groups, sort, mark duplicates, and create index
rule task1:
    input:
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/Alignment/T064_filtered_Aligned.sortedByCoord.out.bam"
    output:
        read_group = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/rg_added_T064-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/dedupped_T064-STARAligned.sortedByCoord.out.bam",
        duplicated_metrics = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/T064_marked_dup_metrix.txt"
    threads: 30
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/"
    shell:
        '''        
        java -jar {input.picard} AddOrReplaceReadGroups I={input.in_bam} O={output.read_group} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        {input.gatk} MarkDuplicates -I {output.read_group} -O {output.dedupped} \
        -CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M {output.duplicated_metrics}
        '''

# 2. Split N trim and reassign mapping qualities
rule task2:
    input:
        reference="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        dedupped ="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/dedupped_T064-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/SplitNCigarReads_T064-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/"
    threads: 30
    shell:
        '''
        {input.gatk} SplitNCigarReads -R {input.reference} --tmp-dir {params.TMPDIR} -I {input.dedupped} -O {output.split}  
        '''


# 3. Base recalibration (BQSR)
rule task3:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        cigar = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/SplitNCigarReads_T064-STARAligned.sortedByCoord.out.bam"
    output:
        recalibration_table="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/BQSR/T064-recal.table",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/BQSR/BQSR_T064-STARAligned.sortedByCoord.out.bam"

    shell:
        '''

        {input.gatk}  BaseRecalibrator -I {input.cigar} -R {input.reference} --known-sites {input.dbSNP} -O {output.recalibration_table}
        {input.gatk}  ApplyBQSR -R {input.reference} -I {input.cigar} --bqsr-recal-file {output.recalibration_table} -O {output.BQSR}
        '''


# 4. Run the haplotypecaller
rule task4:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/BQSR/BQSR_T064-STARAligned.sortedByCoord.out.bam",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
    
    output:
        vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/call/T064.vcf.gz"
    wildcard_constraints:
        sample="[^.]+"
    shell:
        '''
        {input.gatk} HaplotypeCaller -R {input.reference} -I {input.BQSR} -O {output.vcf} -dbsnp {input.dbSNP} --dont-use-soft-clipped-bases --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 0
        '''

rule filter_snps:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        vcf = rules.task4.output.vcf,
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",

    output:
        out_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/call/T064.filtered.vcf.gz"
    shell:
        " {input.gatk} " +
        " VariantFiltration " +
        " -R {input.reference} " +
        " -V {input.vcf} " +
        " -O {output.out_filter} " +
        " --verbosity ERROR " +
        " --filter \"QD < 2.0\" --filter-name \"snp_QD2\" "
        " --filter \"FS > 30.0\" --filter-name \"snp_FS30\" " 
        " --filter \"DP < 4.0\" --filter-name \"snp_DP4\" "




# 6. Selected variants
rule select_variants:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        in_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/call/T064.filtered.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk"

    output:
        out_filter2="/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/call/selected.T064.filtered.vcf.gz",
        
    shell:
        '''
        {input.gatk} SelectVariants -R {input.reference} -V {input.in_filter} --select-type-to-include SNP --exclude-filtered -O {output.out_filter2}
        '''
        
