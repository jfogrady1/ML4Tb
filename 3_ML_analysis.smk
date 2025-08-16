# Snakefile
# modules
import os
import glob
import numpy
import pandas as pd
from pathlib import Path


rule all:
  input:
      "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/PC2_RNA_SNP_wilcox.pdf",
      "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_SNPRELATE_eigenvec.txt"
        
rule PCA_plot:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/PCA_plot.R",
    wiarda_samples = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    wiarda_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_PRUNED.vcf.gz",
    kpbl_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_PRUNED.vcf.gz",
    kirsten_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_PRUNED.vcf.gz",
    ogrady_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_PRUNED.vcf.gz",
    ogrady_snp_eigenvec = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec",
    ogrady_snp_eigenval = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenval"
    
  output:
    wiarda_pdf_pve = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_PVE.pdf", #output
    wiarda_pdf_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_PCA.pdf", # output
    wiarda_gds = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda2.gds", # output
    kpbl_pve = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_PVE.pdf", # output
    kpbl_pve_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_PCA.pdf", #output
    kpbl_gds = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl2.gds",
    kirsten_pve = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_PVE.pdf", # output
    kirsten_pve_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_PCA.pdf", #output
    kirsten_gds = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten.gds",
    ogrady_pve =  "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_PVE.pdf", #output
    ogrady_pve_pca = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_PCA.pdf", #output
    ogrady_gds = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady.gds",
    PC1_correlation = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/PC1_RNA_SNP_correlation.pdf", #output - FigS02
    PC2_correlation = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/PC2_RNA_SNP_correlation.pdf", #output -FigS02
    
    PC1_wilcox = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/PC1_RNA_SNP_wilcox.pdf", #output -FigS02
    PC2_wilcox = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/PC2_RNA_SNP_wilcox.pdf", #output -FigS02
    wiarda_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_SNPRELATE_eigenvec.txt",
    kirsten_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_SNPRELATE_eigenvec.txt",
    kirsten_pbl_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_SNPRELATE_eigenvec.txt",
    ogrady_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_SNPRELATE_eigenvec.txt"

  shell:
    """
    Rscript {input.script} {input.wiarda_samples} {output.wiarda_pdf_pve} {output.wiarda_pdf_pca} {input.wiarda_vcf} {output.wiarda_gds} \
    {output.kpbl_pve} {output.kpbl_pve_pca} {input.kpbl_vcf} {output.kpbl_gds} {output.kirsten_pve} {output.kirsten_pve_pca} {input.kirsten_vcf} \
    {output.kirsten_gds} {output.ogrady_pve} {output.ogrady_pve_pca} {input.ogrady_vcf} {output.ogrady_gds} {input.ogrady_snp_eigenvec} {input.ogrady_snp_eigenval} \
    {output.PC1_correlation} {output.PC2_correlation} {output.PC1_wilcox} {output.PC2_wilcox} {output.wiarda_eigenvec} {output.kirsten_eigenvec} {output.kirsten_pbl_eigenvec} {output.ogrady_eigenvec}
    """
  