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
samples_kirsten_pbl = pd.read_csv("./data/kirsten_pbl/kirsten_pbl_samples.csv", sep = "\t")

# Define a dictionary for mapping run IDs to sample codes
sample_mapping = dict(zip(samples["ena_run"], samples["sample_code"]))
sample_mapping_wiarda = dict(zip(samples_wiarda["Sample_SRR"], samples_wiarda["Animal_Code"]))

# Set the output directory path
output_dir_kirsten = "/home/workspace/jogrady/ML4TB/data/kirsten/individual/"
output_dir_wiarda = "/home/workspace/jogrady/ML4TB/data/wiarda/"
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







sample_ids_kirsten_pbl = ['A016_CON','A017_CON','A018_CON','A020_CON','A022_CON','A023_CON','A024_CON','A029_CON',
                          'A032_TB','A033_TB','A034_TB','A035_TB','A036_TB','A037_TB','A038_TB','A039_TB']



# Extract unique IDs to merge
unique_kirsten_ids = set([s.split("_")[0] for s in Kirsten_sample_ids])
unique_wiarda_ids = set([s.split("_")[0] for s in sample_ids_wiarda])
print(unique_kirsten_ids)
rule all:
    input:
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/merged/{s}.bam", s=unique_kirsten_ids),
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/merged/{s}.bam", s=unique_wiarda_ids),
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/selected.{sample_wiarda}.filtered.vcf.gz", sample_wiarda = unique_wiarda_ids),
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/selected.{sample}.filtered.vcf.gz", sample = unique_kirsten_ids),
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/selected.{sample_kpbl}.filtered.vcf.gz", sample_kpbl = sample_ids_kirsten_pbl),
        "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_ALL_Pruned.vcf.gz", # un comment when running
        expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/{study}_filtered_PRUNED.vcf.gz", study = ["kirsten", "kirsten_pbl", "wiarda", "ogrady"])
        

rule merge_bam:
    input:
        reads = lambda wildcards: [f"/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Alignment/{s}_Aligned.sortedByCoord.out.bam" for s in Kirsten_sample_ids if s.startswith(wildcards.s)]
    output:
        bam = multiext("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/merged/{s}", ".bam", ".bam.bai")
    threads: 30
    shell:
        """
        samtools merge {output.bam[0]} {input.reads} --threads {threads}
        samtools index {output.bam[0]}
        """

rule merge_bam_wiarda:
    input:
        reads = lambda wildcards: [f"/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Alignment/{s}_Aligned.sortedByCoord.out.bam" for s in sample_ids_wiarda if s.startswith(wildcards.s)]
    output:
        bam = multiext("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/merged/{s}", ".bam", ".bam.bai")
    threads: 30
    shell:
        """
        samtools merge {output.bam[0]} {input.reads} --threads {threads}
        samtools index {output.bam[0]}
        """
        

# Kirsten

# 1. Add read groups, sort, mark duplicates, and create index
rule task1:
    input:
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/merged/{sample}.bam"
    output:
        read_group = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/rg_added_{sample}-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/dedupped_{sample}-STARAligned.sortedByCoord.out.bam",
        duplicated_metrics = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/{sample}_marked_dup_metrix.txt"
    threads: 30
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/"
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
        dedupped ="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/dedupped_{sample}-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/"
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
        cigar = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/SplitNCigarReads_{sample}-STARAligned.sortedByCoord.out.bam"
    output:
        recalibration_table="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/BQSR/{sample}-recal.table",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/BQSR/BQSR_{sample}-STARAligned.sortedByCoord.out.bam"

    shell:
        '''

        {input.gatk}  BaseRecalibrator -I {input.cigar} -R {input.reference} --known-sites {input.dbSNP} -O {output.recalibration_table}
        {input.gatk}  ApplyBQSR -R {input.reference} -I {input.cigar} --bqsr-recal-file {output.recalibration_table} -O {output.BQSR}
        '''


# 4. Run the haplotypecaller
rule task4:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/BQSR/BQSR_{sample}-STARAligned.sortedByCoord.out.bam",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
    
    output:
        vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/{sample}.vcf.gz"
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
        out_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/{sample}.filtered.vcf.gz"
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
        in_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/{sample}.filtered.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk"

    output:
        out_filter2="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/selected.{sample}.filtered.vcf.gz",
        
    shell:
        '''
        {input.gatk} SelectVariants -R {input.reference} -V {input.in_filter} --select-type-to-include SNP --exclude-filtered -O {output.out_filter2}
        '''
        
        
        
        

# Wiarda

# 1. Add read groups, sort, mark duplicates, and create index
rule task1_wiarda:
    input:
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/merged/{sample_wiarda}.bam"
    output:
        read_group = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/rg_added_{sample_wiarda}-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/dedupped_{sample_wiarda}-STARAligned.sortedByCoord.out.bam",
        duplicated_metrics = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/{sample_wiarda}_marked_dup_metrix.txt"
    threads: 30
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/"
    shell:
        '''        
        java -jar {input.picard} AddOrReplaceReadGroups I={input.in_bam} O={output.read_group} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        {input.gatk} MarkDuplicates -I {output.read_group} -O {output.dedupped} \
        -CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M {output.duplicated_metrics}
        '''

# 2. Split N trim and reassign mapping qualities
rule task2_wiarda:
    input:
        reference="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        dedupped ="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/dedupped_{sample_wiarda}-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/SplitNCigarReads_{sample_wiarda}-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/"
    threads: 30
    shell:
        '''
        {input.gatk} SplitNCigarReads -R {input.reference} --tmp-dir {params.TMPDIR} -I {input.dedupped} -O {output.split}  
        '''


# 3. Base recalibration (BQSR)
rule task3_wiarda:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        cigar = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/SplitNCigarReads_{sample_wiarda}-STARAligned.sortedByCoord.out.bam"
    output:
        recalibration_table="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/BQSR/{sample_wiarda}-recal.table",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/BQSR/BQSR_{sample_wiarda}-STARAligned.sortedByCoord.out.bam"

    shell:
        '''

        {input.gatk}  BaseRecalibrator -I {input.cigar} -R {input.reference} --known-sites {input.dbSNP} -O {output.recalibration_table}
        {input.gatk}  ApplyBQSR -R {input.reference} -I {input.cigar} --bqsr-recal-file {output.recalibration_table} -O {output.BQSR}
        '''


# 4. Run the haplotypecaller
rule task4_wiarda:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/BQSR/BQSR_{sample_wiarda}-STARAligned.sortedByCoord.out.bam",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
    
    output:
        vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/{sample_wiarda}.vcf.gz"
    wildcard_constraints:
        sample_wiarda="[^.]+"
    shell:
        '''
        {input.gatk} HaplotypeCaller -R {input.reference} -I {input.BQSR} -O {output.vcf} -dbsnp {input.dbSNP} --dont-use-soft-clipped-bases --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 0
        '''

rule filter_snps_wiarda:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        vcf = rules.task4_wiarda.output.vcf,
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",

    output:
        out_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/{sample_wiarda}.filtered.vcf.gz"
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
rule select_variants_wiarda:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        in_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/{sample_wiarda}.filtered.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk"

    output:
        out_filter2="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/selected.{sample_wiarda}.filtered.vcf.gz",
        
    shell:
        '''
        {input.gatk} SelectVariants -R {input.reference} -V {input.in_filter} --select-type-to-include SNP --exclude-filtered -O {output.out_filter2}
        '''



# Kirsten_pbl

# 1. Add read groups, sort, mark duplicates, and create index
rule task1_kpbl:
    input:
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        picard = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/picard.jar",
        in_bam="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Alignment/{sample_kpbl}_Aligned.sortedByCoord.out.bam"
    output:
        read_group = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/rg_added_{sample_kpbl}-STARAligned.sortedByCoord.out.bam",
        dedupped = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/dedupped_{sample_kpbl}-STARAligned.sortedByCoord.out.bam",
        duplicated_metrics = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/{sample_kpbl}_marked_dup_metrix.txt"
    threads: 30
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/"
    shell:
        '''        
        java -jar {input.picard} AddOrReplaceReadGroups I={input.in_bam} O={output.read_group} \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=run RGSM=20 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT SORT_ORDER=coordinate
        
        {input.gatk} MarkDuplicates -I {output.read_group} -O {output.dedupped} \
        -CREATE_INDEX true --VALIDATION_STRINGENCY SILENT -M {output.duplicated_metrics}
        '''

# 2. Split N trim and reassign mapping qualities
rule task2_kpbl:
    input:
        reference="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        dedupped ="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/dedupped_{sample_kpbl}-STARAligned.sortedByCoord.out.bam"
    output:
        split = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/SplitNCigarReads_{sample_kpbl}-STARAligned.sortedByCoord.out.bam"
    params: TMPDIR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/"
    threads: 30
    shell:
        '''
        {input.gatk} SplitNCigarReads -R {input.reference} --tmp-dir {params.TMPDIR} -I {input.dedupped} -O {output.split}  
        '''


# 3. Base recalibration (BQSR)
rule task3_kpbl:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
        cigar = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/SplitNCigarReads_{sample_kpbl}-STARAligned.sortedByCoord.out.bam"
    output:
        recalibration_table="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/BQSR/{sample_kpbl}-recal.table",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/BQSR/BQSR_{sample_kpbl}-STARAligned.sortedByCoord.out.bam"

    shell:
        '''

        {input.gatk}  BaseRecalibrator -I {input.cigar} -R {input.reference} --known-sites {input.dbSNP} -O {output.recalibration_table}
        {input.gatk}  ApplyBQSR -R {input.reference} -I {input.cigar} --bqsr-recal-file {output.recalibration_table} -O {output.BQSR}
        '''


# 4. Run the haplotypecaller
rule task4_kpbl:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        BQSR = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/BQSR/BQSR_{sample_kpbl}-STARAligned.sortedByCoord.out.bam",
        dbSNP="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/TWAS/ARS1.2PlusY_BQSR.vcf",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",
    
    output:
        vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/{sample_kpbl}.vcf.gz"
    wildcard_constraints:
        sample_kpbl="[^.]+"
    shell:
        '''
        {input.gatk} HaplotypeCaller -R {input.reference} -I {input.BQSR} -O {output.vcf} -dbsnp {input.dbSNP} --dont-use-soft-clipped-bases --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 0
        '''

rule filter_snps_kpbl:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        vcf = rules.task4_kpbl.output.vcf,
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk",

    output:
        out_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/{sample_kpbl}.filtered.vcf.gz"
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
rule select_variants_kpbl:
    input:
        reference ="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa",
        in_filter="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/{sample_kpbl}.filtered.vcf.gz",
        gatk = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/bin/Response/gatk-4.3.0.0/gatk"

    output:
        out_filter2="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/selected.{sample_kpbl}.filtered.vcf.gz",
        
    shell:
        '''
        {input.gatk} SelectVariants -R {input.reference} -V {input.in_filter} --select-type-to-include SNP --exclude-filtered -O {output.out_filter2}
        '''
        
        
        


        

###########################

# Rename and merge

###########################




# Kirsten
rule rename_vcf_sample_kirsten:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/call/selected.{sample}.filtered.vcf.gz",
    output:
        renamed_vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/{sample}.filtered.renamed.vcf.gz"
    params:
        sample="{sample}"
    shell:
        "echo '{params.sample}' > {params.sample}.txt && "
        "bcftools reheader -s {params.sample}.txt {input.vcf} -o {output.renamed_vcf} && tabix -p vcf {output.renamed_vcf} && "
        "rm {params.sample}.txt"



rule merge_vcfs_kirsten:
    input:
        vcfs = expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/{sample}.filtered.renamed.vcf.gz", sample = unique_kirsten_ids)

    output:
        merged_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/kirsten_RNA_merged.vcf.gz"
    shell:
        """
        bcftools merge {input.vcfs} -O z -o {output.merged_vcf} && tabix -p vcf {output.merged_vcf}
        """



# Kirsten_ PBL

rule rename_vcf_sample_kirsten_pbl:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/call/selected.{sample_kpbl}.filtered.vcf.gz",
    output:
        renamed_vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/{sample_kpbl}.filtered.renamed.vcf.gz"
    params:
        sample="{sample_kpbl}"
    shell:
        "echo '{params.sample}' > {params.sample}.txt && "
        "bcftools reheader -s {params.sample}.txt {input.vcf} -o {output.renamed_vcf} && tabix -p vcf {output.renamed_vcf} && "
        "rm {params.sample}.txt"



rule merge_vcfs_kirsten_pbl:
    input:
        vcfs = expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/{sample_kpbl}.filtered.renamed.vcf.gz", sample_kpbl = sample_ids_kirsten_pbl)

    output:
        merged_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/kirsten_pbl_RNA_merged.vcf.gz"
    shell:
        """
        bcftools merge {input.vcfs} -O z -o {output.merged_vcf} && tabix -p vcf {output.merged_vcf}
        """



# Wiarda
rule rename_vcf_sample_wiarda:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/call/selected.{sample_wiarda}.filtered.vcf.gz",
    output:
        renamed_vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/{sample_wiarda}.filtered.renamed.vcf.gz"
    params:
        sample="{sample_wiarda}"
    shell:
        "echo '{params.sample}' > {params.sample}.txt && "
        "bcftools reheader -s {params.sample}.txt {input.vcf} -o {output.renamed_vcf} && tabix -p vcf {output.renamed_vcf} && "
        "rm {params.sample}.txt"


rule merge_vcfs_wiarda:
    input:
        vcfs = expand("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/{sample_wiarda}.filtered.renamed.vcf.gz", sample_wiarda = unique_wiarda_ids)

    output:
        merged_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/wiarda_RNA_merged.vcf.gz"
    shell:
        """
        bcftools merge {input.vcfs} -O z -o {output.merged_vcf} && tabix -p vcf {output.merged_vcf}
        """



rule intersect_vcfs:
    input:
        kirsten="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/variant_call/kirsten_RNA_merged.vcf.gz",
        kirsten_pbl="/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/variant_call/kirsten_pbl_RNA_merged.vcf.gz",
        wiarda="/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/variant_call/wiarda_RNA_merged.vcf.gz",
        ogrady_WGS="/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/EQTL/formatting/ALL_IMPUTED_UPDATED.vcf.gz"

    output:
        intersect="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0000.vcf.gz",
        intersect1="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0001.vcf.gz",
        intersect2="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0002.vcf.gz",
        intersect3="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0003.vcf.gz"
    params:
        prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2"

    shell:
        """
        bcftools isec -n=4 -c none -Oz -p {params.prefix} {input.kirsten} {input.kirsten_pbl} {input.wiarda} {input.ogrady_WGS}
        """

rule merge_all_vcfs:
    input:
        kirsten="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0000.vcf.gz",
        kirsten_pbl="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0001.vcf.gz",
        wiarda="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0002.vcf.gz",
        ogrady_WGS="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0003.vcf.gz"

    output:
        intersect="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_merged.vcf.gz"

    shell:
        """
        bcftools merge {input.kirsten} {input.kirsten_pbl} {input.wiarda} {input.ogrady_WGS} -Oz -o {output.intersect}
        tabix -p vcf {output.intersect}
        """


rule filter_all:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_merged.vcf.gz",
    output:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz | bcftools view -q 0.05:minor -Oz | bcftools view -e 'HWE < 0.000001 || F_MISSING > 0.2' -Oz -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
        """
rule make_plink_all:
    input:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_filtered.vcf.gz"    
    output:
        plink = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_filtered.bed",
        wiarda = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_ALL_Pruned.vcf.gz",
        kirsten = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_ALL_Pruned.vcf.gz",
        kirsten_pbl = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_ALL_Pruned.vcf.gz",
        ogrady = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_ALL_Pruned.vcf.gz",

    params:
        plink_prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_filtered",
        wiarda_samples = ",".join(f"{id}_{id}" for id in unique_wiarda_ids),
        kirsten_samples = ",".join(f"{id}_{id}" for id in unique_kirsten_ids),
        kirsten_pbL_samples = ",".join(sample_ids_kirsten_pbl),
        # PCA
        wiarda_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_ALL_Pruned",
        kirsten_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_ALL_Pruned",
        kirsten_pbl_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_ALL_Pruned",
        ogrady_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_ALL_Pruned"

    shell:
        """
        plink --cow --keep-allele-order --vcf {input.vcf_filtered} --autosome --allow-extra-chr --make-bed --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --id-delim _ --indep-pairwise 1000 5 0.1 --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --id-delim _ --extract {params.plink_prefix}.prune.in --make-bed --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --id-delim _ --pca --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --id-delim _ --recode vcf bgz --out {params.plink_prefix}_PRUNED_ALL

        tabix -f -p vcf {params.plink_prefix}_PRUNED_ALL.vcf.gz

        bcftools view -s {params.wiarda_samples} {params.plink_prefix}_PRUNED_ALL.vcf.gz -Oz -o  {output.wiarda} && tabix -p vcf {output.wiarda}
        bcftools view -s {params.kirsten_samples} {params.plink_prefix}_PRUNED_ALL.vcf.gz -Oz -o  {output.kirsten} && tabix -p vcf {output.kirsten}
        bcftools view -s {params.kirsten_pbL_samples} {params.plink_prefix}_PRUNED_ALL.vcf.gz -Oz -o  {output.kirsten_pbl} && tabix -p vcf {output.kirsten_pbl}
        # ogrady
        bcftools view -s ^{params.kirsten_pbL_samples},{params.kirsten_samples},{params.wiarda_samples} {params.plink_prefix}_PRUNED_ALL.vcf.gz -Oz -o  {output.ogrady} && tabix -p vcf {output.ogrady}

        # Now do the pca on each data set
        plink --cow --keep-allele-order --vcf {output.wiarda} --autosome --allow-extra-chr --pca --out {params.wiarda_pca}
        plink --cow --keep-allele-order --vcf {output.kirsten} --autosome --allow-extra-chr --pca --out {params.kirsten_pca}
        plink --cow --keep-allele-order --vcf {output.kirsten_pbl} --autosome --allow-extra-chr --pca --out {params.kirsten_pbl_pca}
        plink --cow --keep-allele-order --vcf {output.ogrady} --autosome --allow-extra-chr --pca --out {params.ogrady_pca}
        """



rule filter_kirsten:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0000.vcf.gz",
    output:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz | bcftools view -q 0.1:minor -Oz | bcftools view -e 'HWE < 0.000001' -Oz -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
        """


rule filter_kirsten_pbl:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0001.vcf.gz",
    output:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz | bcftools view -q 0.1:minor -Oz | bcftools view -e 'HWE < 0.000001' -Oz -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
        """


rule filter_wiarda:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0002.vcf.gz",
    output:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz | bcftools view -q 0.1:minor -Oz | bcftools view -e 'HWE < 0.000001' -Oz -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
        """

rule filter_ogrady:
    input:
        vcf="/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/0003.vcf.gz",
    output:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered.vcf.gz"
    shell:
        """
        bcftools +fill-tags {input.vcf} -Oz | bcftools view -q 0.05:minor -Oz | bcftools view -e 'HWE < 0.000001' -Oz -o {output.vcf_filtered}
        tabix -p vcf {output.vcf_filtered}
        """
        
rule make_plink_kirsten:
    input:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered.vcf.gz"    
    output:
        plink = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_PRUNED.vcf.gz"

    params:
        plink_prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered"
    shell:
        """
        plink --cow --keep-allele-order --vcf {input.vcf_filtered} --autosome --allow-extra-chr --make-bed --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --indep-pairwise 1000 5 0.1 --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --extract {params.plink_prefix}.prune.in --make-bed --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --recode vcf bgz --out {params.plink_prefix}_PRUNED && tabix -p vcf {params.plink_prefix}_PRUNED.vcf.gz
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --pca --out {params.plink_prefix}_PRUNED 
        """
        
rule make_plink_kirsten_pbl:
    input:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered.vcf.gz"    
    output:
        plink = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_PRUNED.vcf.gz"

    params:
        plink_prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered"
    shell:
        """
        plink --cow --keep-allele-order --vcf {input.vcf_filtered} --autosome --allow-extra-chr --make-bed --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --indep-pairwise 1000 5 0.1 --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --extract {params.plink_prefix}.prune.in --make-bed --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --recode vcf bgz --out {params.plink_prefix}_PRUNED && tabix -p vcf {params.plink_prefix}_PRUNED.vcf.gz
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --pca --out {params.plink_prefix}_PRUNED 
        """
        

rule make_plink_ogrady:
    input:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered.vcf.gz"    
    output:
        plink = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_PRUNED.vcf.gz"

    params:
        plink_prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered"
    shell:
        """
        plink --cow --keep-allele-order --vcf {input.vcf_filtered} --autosome --allow-extra-chr --make-bed --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --indep-pairwise 1000 5 0.1 --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --extract {params.plink_prefix}.prune.in --make-bed --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --recode vcf bgz --out {params.plink_prefix}_PRUNED && tabix -p vcf {params.plink_prefix}_PRUNED.vcf.gz
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --pca --out {params.plink_prefix}_PRUNED 
        """


rule make_plink_wiarda:
    input:
        vcf_filtered = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered.vcf.gz"    
    output:
        plink = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_PRUNED.vcf.gz"
    params:
        plink_prefix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered"
    shell:
        """
        plink --cow --keep-allele-order --vcf {input.vcf_filtered} --autosome --allow-extra-chr --make-bed --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --indep-pairwise 1000 5 0.1 --out {params.plink_prefix}
        plink --cow --bfile {params.plink_prefix} --keep-allele-order --extract {params.plink_prefix}.prune.in --make-bed --out {params.plink_prefix}_PRUNED
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --recode vcf bgz --out {params.plink_prefix}_PRUNED && tabix -p vcf {params.plink_prefix}_PRUNED.vcf.gz
        plink --cow --bfile {params.plink_prefix}_PRUNED --keep-allele-order --pca --out {params.plink_prefix}_PRUNED 
        """
