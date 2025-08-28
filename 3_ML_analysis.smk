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
      "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_SNPRELATE_eigenvec.txt",
      "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_8_models.pdf",
      "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/Upset_plot.pdf",
      "/home/workspace/jogrady/ML4TB/work/normalisation/Figures/kirsten_PCA.pdf",
      "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/ROC_iteration_curve.pdf",
      "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_OD_AUROC.pdf",
      "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_GFS_models_FINAL.pdf",
      "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_8_models_FINAL.pdf"



rule exp_PCA_plot:
  input:
    script =  "/home/workspace/jogrady/ML4TB/bin/Preprocessing/1_data_preprocessing_visualisation.R",
    ogrady_data_raw = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt",
    wiarda_data_raw = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt",
    kirsten_data_raw = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt",
    kirsten_pbl_raw = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt",
    ogrady_covariate = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt",
    wiarda_covariate = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    kirsten_covariate = "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_samples.csv",
    kirsten_pbl_covariate = "/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv",
    

  output:
    dds_ogrady_full_vst = "/home/workspace/jogrady/ML4TB/work/normalisation/ogrady_dds_vst.rds",
    ogrady_pca = "/home/workspace/jogrady/ML4TB/work/normalisation/Figures/ogrady_PCA.pdf",
    dds_wiarda_full_vst = "/home/workspace/jogrady/ML4TB/work/normalisation/wiarda_dds_vst.rds",
    wiarda_pca = "/home/workspace/jogrady/ML4TB/work/normalisation/Figures/wiarda_PCA.pdf",
    dds_kirsten_full_vst = "/home/workspace/jogrady/ML4TB/work/normalisation/kirsten_dds_vst.rds",
    kirsten_pca = "/home/workspace/jogrady/ML4TB/work/normalisation/Figures/kirsten_PCA.pdf",
    dds_kirsten_pbl_full_vst = "/home/workspace/jogrady/ML4TB/work/normalisation/kirsten_pbl_dds_vst.rds",
    kirsten_pbl_pca = "/home/workspace/jogrady/ML4TB/work/normalisation/Figures/kirsten_pbl_PCA.pdf"
  shell:
    """
    Rscript {input.script} {input.ogrady_data_raw} {input.wiarda_data_raw} {input.kirsten_data_raw} {input.kirsten_pbl_raw} \
    {input.ogrady_covariate} {input.wiarda_covariate} {input.kirsten_covariate} {input.kirsten_pbl_covariate} \
    {output.dds_ogrady_full_vst} {output.ogrady_pca} {output.dds_wiarda_full_vst} {output.wiarda_pca} {output.dds_kirsten_full_vst} {output.kirsten_pca} {output.dds_kirsten_pbl_full_vst} {output.kirsten_pbl_pca}
    """


rule geno_PCA_plot:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/PCA_plot.R",
    wiarda_samples = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    wiarda_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_ALL.vcf.gz",
    kpbl_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_ALL.vcf.gz",
    kirsten_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_ALL.vcf.gz",
    ogrady_vcf = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_ALL.vcf.gz",
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

rule DE_analysis_all_sets:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/DE_analysis_all_sets.R",
    functions = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Funcitons.R",
    ogr25_data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt", # From previous study (O'grady et al., 2025)
    ensemble_gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
    ogr_labels = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt",
    ogr_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_SNPRELATE_eigenvec.txt",
    wia20_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt",
    wia_labels = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    wiarda_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_SNPRELATE_eigenvec.txt",
    mcl21_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt",
    mcl_labels = "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt",
    mcl_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_SNPRELATE_eigenvec.txt",
    mcl14_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt",
    mcl14_labels = "/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv",
    mcl14_eigenvec = "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_SNPRELATE_eigenvec.txt"

  output:
    upset = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/Upset_plot.pdf",
    heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/Heatmap_DE_genes.pdf",
    ogrady_gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/ogrady_gprofiler.pdf",
    mcl14_gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/mcloughlin_pbl_gprofiler.pdf",
    mcl21_gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/mcloughlin_gprofiler.pdf",
    wiarda_gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/wiarda_gprofiler.pdf",
    text_gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/gprofiler_results.txt",
    DE_analysis_all_sets = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/all_de_results.txt",
    heatmap_text = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/DE_analysis/Heatmap_lfc.txt"

  shell:
    """
    Rscript {input.script} {input.functions} {input.ogr25_data} {input.ensemble_gtf} {input.ogr_labels} {input.ogr_eigenvec} {input.wia20_data} \
    {input.wia_labels} {input.wiarda_eigenvec} {input.mcl21_data} {input.mcl_labels} {input.mcl_eigenvec} {input.mcl14_data} {input.mcl14_labels} {input.mcl14_eigenvec} \
    {output.upset} {output.heatmap} {output.ogrady_gprofiler} {output.mcl14_gprofiler} {output.mcl21_gprofiler} {output.wiarda_gprofiler} {output.text_gprofiler} {output.DE_analysis_all_sets} {output.heatmap_text}
    """
  
rule Model_training:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Merged_analysis_smk.R",
    functions = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Funcitons.R",
    ogr25_data = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt", # From previous study (O'grady et al., 2025)
    ensemble_gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
    ogr_labels = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt",
    wia20_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt",
    wia_labels = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    mcl21_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt",
    mcl_labels = "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt",
    mcl14_data = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt",
    mcl14_labels = "/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv"


  output:
    dds_test = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_RAW_DESEQ2_matrix.rds",
    dds_train = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Train_RAW_DESEQ2_matrix.rds",
    train_PCA = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/PCA_batch/PCA_train_unadjusted.pdf",
    test_PCA = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/PCA_batch/PCA_test_unadjusted.pdf",
    train_geno_pve = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/vcf_isec2/Train_filtered_PRUNED_PVE.pdf",
    train_geno_pca = "/home/workspace/jogrady/ML4TB/snktest/work/RNA_seq/vcf_isec2/Train_filtered_PRUNED_PCA.pdf",
    de_results = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/DE_results_integrated.txt",
    volcano = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/DESeq2_results_top15_annotated.pdf",
    gprofiler = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Merged_gprofiler.pdf",
    Train_adjusted_matrix = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Train_adjusted_DESEQ2_matrix.rds",
    test_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_counts.txt",
    final_de_training_results = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Training_DE_results.txt",
    ora_full_results = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ORA_training_set.txt",
    train_metadata = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_data_manuscript.txt",
    train_vst_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_normalised_filtered_counts.txt",
    train_folds = expand("/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_normalised_filtered_counts_Fold{n}.txt", n = range(1, 6)),
    test_vst_counts =  "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/test_normalised_filtered_counts.txt",
    folds_integers = expand("/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Fold{n}_integers.txt", n = range(1, 6)),
    test_metadata = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_set_manuscript.txt",

    glm = "/home/workspace/jogrady/ML4TB/models/GLM_models.rds",
    glmridge = "/home/workspace/jogrady/ML4TB/models/GLMRIDGE_model.rds",
    glmlasso = "/home/workspace/jogrady/ML4TB/models/GLMLASSO_model.rds",
    glmenet = "/home/workspace/jogrady/ML4TB/models/GLMENET_model.rds",
    rf = "/home/workspace/jogrady/ML4TB/models/Random_forest_list_models.rds",
    rf_et = "/home/workspace/jogrady/ML4TB/models/Random_forest_extratrees_list_models.rds",
    nb = "/home/workspace/jogrady/ML4TB/models/NB_model.rds",
    mlp = "/home/workspace/jogrady/ML4TB/models/MLP_model.rds",


    # ROC training
    glm_train_roc = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/GLM_CV_ROC.pdf",
    glmridge_train_roc = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/GLMRIDGE_CV_ROC.pdf",
    glmlasso_train_roc = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/GLMLASSO_CV_ROC.pdf",
    glmenet_train_roc = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/GLMENET_CV_ROC.pdf",

    # RF
    rf_paramaters = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/RF_Gini_node_mtry_variables.pdf",
    rf_best_performer = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/RF_CV_ROC.pdf",
    rf_et_paramaters = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/RF_EXTRATREES_Gini_node_mtry_variables.pdf",
    rf_et_best_performer = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/RF_ET_CV_ROC.pdf",

    # nb
    nb_best_performer = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/NB_train_ROC.pdf",

    # mlp
    mlp_best_performer = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/MLP_train_ROC.pdf",

    dotplot = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/Train_CV_ROC_dotplot.pdf",
    table_dotplot = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Train_dot_plot_results.txt",

    test_model_results = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/TEST1_ROC.pdf",

    natural_V_experimental = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Train_CV_results/NATURAL_V_Infected_ROCs.pdf",

    natural_v_experimental_p = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/NATURAL_V_Infected_results_wilcox.txt",

    natural_v_experimental_txt = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/NATURAL_V_Infected_raw_results.txt",

    heatmap_txt = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Heatmap_8_models.txt",

    heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_8_models.pdf",

  shell:
    """
    Rscript {input.script} {input.functions} {input.ogr25_data} {input.ensemble_gtf} {input.ogr_labels} {input.wia20_data} \
    {input.wia_labels} {input.mcl21_data} {input.mcl_labels} {input.mcl14_data} {input.mcl14_labels} \
    {output.dds_test} {output.dds_train} {output.train_PCA} {output.test_PCA} {output.train_geno_pve} \
    {output.train_geno_pca} {output.de_results} {output.volcano} {output.gprofiler} {output.Train_adjusted_matrix} \
    {output.test_counts} {output.final_de_training_results} {output.ora_full_results} {output.train_metadata} {output.train_vst_counts} \
    {output.train_folds} {output.test_vst_counts} {output.folds_integers} {output.test_metadata} \
    {output.glm} {output.glmridge} {output.glmlasso} {output.glmenet} {output.glm_train_roc} {output.glmridge_train_roc} {output.glmlasso_train_roc} {output.glmenet_train_roc} \
    {output.rf} {output.rf_paramaters} {output.rf_best_performer} {output.rf_et} {output.rf_et_paramaters} {output.rf_et_best_performer} \
    {output.nb} {output.nb_best_performer} {output.mlp} {output.mlp_best_performer} {output.dotplot} {output.table_dotplot} \
    {output.test_model_results} {output.natural_v_experimental_p} {output.natural_v_experimental_txt} {output.natural_V_experimental} {output.heatmap_txt} {output.heatmap}
    """

rule GFS:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/GFS_search.R",
    train_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_normalised_filtered_counts.txt",
    ensemble_gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
    ogrady_labs = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt",
    wiarda_counts = "/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt",
    wiarda_labs = "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv",
    kirsten_counts = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt",
    kirsten_labs = "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt",
    kirsten_pbl_counts = "/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt",
    kirsten_pbl_labs = "/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv",
    DE_results_integrated = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/DE_results_integrated.txt",
    functions_script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Funcitons.R",
    train_metadata = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_data_manuscript.txt",
    IGRA_results = "/home/workspace/jogrady/ML4TB/data/ogrady/Ogrady_IGRA_results.txt",
    test_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/test_normalised_filtered_counts.txt",
    test_metadata = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_set_manuscript.txt"

    
  output:
    ROC_set_1 = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC1_genes.txt",
    ROC_set_2 = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC2_genes.txt",
    TB_score_training = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_score_training.pdf",
    scores_wilcox_train = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_score_wilcox_training.txt",
    GFS_training_performance ="/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_AUC_performance_training.txt",
    GFS_dot_plot = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/AUC_0.95_3Sets_GFS.pdf",
    GFS_Train_ROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_Train_ROC_curves.pdf",
    IGRA_correlation = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/IGRA_TB_Score_correlation.pdf",
    GFS_testing_score = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_score_testing.pdf",
    test_score_wilcox = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_score_wilcox_testing.txt",
    heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_3_models.pdf",
    GFS_ROC_test = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_Test_ROC_curves.pdf",
    GFS_ROC_iteration_curve = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/ROC_iteration_curve.pdf"

  shell:
    """
    Rscript {input.script} {input.train_counts} {input.ensemble_gtf} {input.ogrady_labs} {input.wiarda_counts} {input.wiarda_labs} \
    {input.kirsten_counts} {input.kirsten_labs} {input.kirsten_pbl_counts} {input.kirsten_pbl_labs} {input.DE_results_integrated} {input.functions_script} \
    {input.train_metadata} {output.ROC_set_1} {output.ROC_set_2} {output.TB_score_training} {output.scores_wilcox_train} {output.GFS_training_performance} {output.GFS_dot_plot} {output.GFS_Train_ROC} \
    {input.IGRA_results} {output.IGRA_correlation} {input.test_counts} {output.GFS_testing_score} {output.test_score_wilcox} {input.test_metadata} {output.heatmap} {output.GFS_ROC_test} {output.GFS_ROC_iteration_curve}
    """

rule External_evaluation:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/External_evaluation_1.R",
    functions = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Funcitons.R",
    ensemble_gtf = "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf",
    bohv1_matrix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/odonoghue/Quantification/odonoghue_count_matrix_clean.txt",
    brsv_matrix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/johnston/Quantification/johnston_count_matrix_clean.txt",
    map_matrix = "/home/workspace/jogrady/ML4TB/work/RNA_seq/alonso/Quantification/alonso_count_matrix_clean.txt",
    bohv1_metadata = "/home/workspace/jogrady/ML4TB/data/odonoghue/odonoghue_samples.tsv",
    brsv_metadata = "/home/workspace/jogrady/ML4TB/data/johnston/johnston_samples.tsv",
    map_metadata = "/home/workspace/jogrady/ML4TB/data/alonso/alonso_samples.tsv",
    dds_train = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Train_adjusted_DESEQ2_matrix.rds",
    DE_results_integrated = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/DE_results_integrated.txt",
    glm_model = "/home/workspace/jogrady/ML4TB/models/GLM_models.rds",
    glmridge_model = "/home/workspace/jogrady/ML4TB/models/GLMRIDGE_model.rds",
    glmlasso_model = "/home/workspace/jogrady/ML4TB/models/GLMLASSO_model.rds",
    glmenet_model = "/home/workspace/jogrady/ML4TB/models/GLMENET_model.rds",
    rf_model = "/home/workspace/jogrady/ML4TB/models/Random_forest_list_models.rds",
    rf_et_model = "/home/workspace/jogrady/ML4TB/models/Random_forest_extratrees_list_models.rds",
    nb_model = "/home/workspace/jogrady/ML4TB/models/NB_model.rds",
    mlp_model = "/home/workspace/jogrady/ML4TB/models/MLP_model.rds",
    ROC1_genes = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC1_genes.txt",
    ROC2_genes = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC2_genes.txt",
    train_normalised_filtered_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_normalised_filtered_counts.txt",
    train_labels = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_data_manuscript.txt",
    test_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_counts.txt", # Need Raw
    test_labels = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_set_manuscript.txt"


    
  output:
    map_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/MAP_AUROC.pdf",
    brsv_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/BRSV_AUROC.pdf",
    bohv1_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/BOHV1_AUROC.pdf",
    GFS_OD_Score = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_OD_Score.pdf",
    wilcox_OD_Score = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Wilcox_External_GFS_comparison.txt",
    GFS_MAP_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_MAP_AUROC.pdf",
    GFS_BRSV_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_BRSV_AUROC.pdf",
    GFS_BOHV1_AUROC = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/GFS_BOHV1_AUROC.pdf",
    ROC_OD_curve = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/TB_OD_AUROC.pdf"
    
  shell:
    """
    Rscript {input.script} {input.functions} {input.ensemble_gtf} {input.bohv1_matrix} {input.brsv_matrix} {input.map_matrix} {input.bohv1_metadata} {input.brsv_metadata} {input.map_metadata} \
    {input.dds_train} {input.DE_results_integrated} {input.glm_model} {input.glmridge_model} {input.glmlasso_model} {input.glmenet_model} {input.rf_model} {input.rf_et_model} {input.nb_model} {input.mlp_model} \
    {output.map_AUROC} {output.brsv_AUROC} {output.bohv1_AUROC} {input.ROC1_genes} {input.ROC2_genes} {input.train_normalised_filtered_counts} {input.train_labels} {output.GFS_OD_Score} {output.wilcox_OD_Score} \
    {output.GFS_MAP_AUROC} {output.GFS_BRSV_AUROC} {output.GFS_BOHV1_AUROC} {input.test_counts} {input.test_labels} {output.ROC_OD_curve} 
    
    """

rule GFS_Heatmap:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Heatmap_GFS_metadata.R",
    functions = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Funcitons.R",
    DE_results_integrated = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/DE_results_integrated.txt",
    test_normalised_filtered_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/test_normalised_filtered_counts.txt",
    test_labels = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_set_manuscript.txt",
    ROC1_genes = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC1_genes.txt",
    ROC2_genes = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/ROC2_genes.txt",
    train_normalised_filtered_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_normalised_filtered_counts.txt",
    train_labels = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_data_manuscript.txt",
    
  output:
    GFS_final_threshold = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/GFS_thresholds_FINAL.txt",
    GFS_text_heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/GFS_HEATMAP_FINAL.txt",
    heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_GFS_models_FINAL.pdf"

  shell:
    """
    Rscript {input.script} {input.functions} {input.DE_results_integrated} {input.test_normalised_filtered_counts} {input.test_labels} {input.ROC1_genes} {input.ROC2_genes} \
    {input.train_normalised_filtered_counts} {input.train_labels} {output.GFS_final_threshold} {output.GFS_text_heatmap} {output.heatmap}
    """

rule Heatmap_models:
  input:
    script = "/home/workspace/jogrady/ML4TB/bin/Merged_analysis/Remake_heatmap.R",
    GLM_model = "/home/workspace/jogrady/ML4TB/models/GLM_models.rds",
    ENET_model = "/home/workspace/jogrady/ML4TB/models/GLMENET_model.rds",
    LASSO_model = "/home/workspace/jogrady/ML4TB/models/GLMLASSO_model.rds",
    RIDGE_model = "/home/workspace/jogrady/ML4TB/models/GLMRIDGE_model.rds",
    MLP_model = "/home/workspace/jogrady/ML4TB/models/MLP_model.rds",
    NB_model = "/home/workspace/jogrady/ML4TB/models/NB_model.rds",
    RF_model = "/home/workspace/jogrady/ML4TB/models/Random_forest_list_models.rds",
    RF_ET_model = "/home/workspace/jogrady/ML4TB/models/Random_forest_extratrees_list_models.rds",
    test_normalised_filtered_counts = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/test_normalised_filtered_counts.txt",
    test_labels = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Test_set_manuscript.txt",

  output:
    heatmap_text = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Heatmap_8_models_FINAL.txt",
    heatmap = "/home/workspace/jogrady/ML4TB/snktest/work/merged/figures/Heatmap_8_models_FINAL.pdf",
    thresholds = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Threshold.txt"


  shell:
    """
    Rscript {input.script} {input.GLM_model} {input.ENET_model} {input.LASSO_model} {input.RIDGE_model} {input.MLP_model} {input.NB_model} {input.RF_model} {input.RF_ET_model} \
    {input.test_normalised_filtered_counts} {input.test_labels} {output.heatmap_text} {output.heatmap} {output.thresholds}
    """

