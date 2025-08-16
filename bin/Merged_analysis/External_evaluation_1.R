
# load libraries
library(tidyverse)
library(ggplot2)
library(caret)
library(sva)
library(DESeq2)
library(pROC)
library("ochRe")
library("ranger")
library(ggsci)
library(plotROC)
library(nnet)
library(pROC)
library(RColorBrewer)
library(grid)
library("naivebayes")
library(ggrepel)
library(data.table)
library(pROC)
library(doParallel)
library(cowplot)
library(kernelshap)
library(matrixStats)
library(ggpubr)
library(gprofiler2)
library(viridis)
library(magrittr)
library(dplyr)
source("~/ML4TB/bin/Merged_analysis/Funcitons.R")
set.seed(42) # For reproducibility


# Ensemble
ensemble <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/Bos_taurus.ARS-UCD1.2.110.gtf")
ensemble <- ensemble %>% filter(V3 == "gene")
head(ensemble)
ensemble <- ensemble %>% separate(., V9, into = c("gene_id", "gene_version", "gene_name"), sep = ";")
ensemble$gene_id <- gsub("^gene_id ", "", ensemble$gene_id)
ensemble$gene_id <- gsub('"', '', ensemble$gene_id)
ensemble$gene_name <- gsub("gene_name ", "", ensemble$gene_name)
ensemble$gene_name <- gsub("gene_source ", "", ensemble$gene_name)
ensemble$gene_name <- gsub('"', '', ensemble$gene_name)

ensemble$gene_name <- if_else(ensemble$gene_name == " ensembl", ensemble$gene_id, ensemble$gene_name)
ensemble$gene_name <- if_else(ensemble$gene_name == " 5S_rRNA", ensemble$gene_id, ensemble$gene_name)
colnames(ensemble)[1] <- "chr"
ensemble <- ensemble %>% dplyr::select(gene_id, gene_name, chr, V4)
colnames(ensemble)[4] <- "pos"
ensemble <- ensemble %>% select(1:2)



ensemble$gene_name <- if_else(duplicated(ensemble$gene_name), ensemble$gene_id, ensemble$gene_name)
ensemble$gene_name <- gsub(" ", "", ensemble$gene_name)


# Read in the external evaluation data

bohv1_matrix = read.csv("/home/workspace/jogrady/ML4TB/work/RNA_seq/odonoghue/Quantification/odonoghue_count_matrix_clean.txt", sep = "\t", row.names = 1, check.names = FALSE)
brsv_matrix = read.csv("/home/workspace/jogrady/ML4TB/work/RNA_seq/johnston/Quantification/johnston_count_matrix_clean.txt", sep = "\t", row.names = 1, check.names = FALSE)
map_matrix = read.csv("/home/workspace/jogrady/ML4TB/work/RNA_seq/alonso/Quantification/alonso_count_matrix_clean.txt", sep = "\t", row.names = 1, check.names = FALSE)

all(rownames(bohv1_matrix) == ensemble$gene_id)
rownames(bohv1_matrix) <- ensemble$gene_name 
rownames(brsv_matrix) <- ensemble$gene_name
rownames(map_matrix) <- ensemble$gene_name
rownames(map_matrix)

bohv1_metadata = read.csv("/home/workspace/jogrady/ML4TB/data/odonoghue/odonoghue_samples.tsv", sep = "\t", row.names = "Animal_Code", check.names = FALSE)
brsv_metadata = read.csv("/home/workspace/jogrady/ML4TB/data/johnston/johnston_samples.tsv", sep = "\t", row.names = "Animal_Code", check.names = FALSE)
map_metadata = read.csv("/home/workspace/jogrady/ML4TB/data/alonso/alonso_samples.tsv", sep = "\t", row.names = "Animal_Code", check.names = FALSE)

head(brsv_metadata)

bohv1_metadata$Status = c(rep("Infected", 12), rep("Control", 6))
brsv_metadata$Status = c(rep("Control", 6), rep("Infected", 12))
map_metadata$Status = c(rep("Control", 3), rep("Infected", 11))
bohv1_metadata$TB_status = "Control"
brsv_metadata$TB_status = "Control"
map_metadata$TB_status = "Control"

head(bohv1_matrix)
# Generate DESEQ2 object of external evaluation data
bohv1_dds <- DESeqDataSetFromMatrix(countData = bohv1_matrix,
                                     colData = bohv1_metadata,
                                     design = ~ 1)
brsv_dds <- DESeqDataSetFromMatrix(countData = brsv_matrix,
                                     colData = brsv_metadata,
                                     design = ~ 1)
map_dds <- DESeqDataSetFromMatrix(countData = map_matrix,
                                     colData = map_metadata,
                                     design = ~ 1)

ddsTrain_adjusted <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ddsTrain_adjusted.rds")
DE_genes <- fread("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/DE_results_integrated.txt") %>%
 filter(Retained != "Excluded") %>% select(Symbol)
DE_genes$Symbol

bohv1_dds <- DESeq(bohv1_dds)
brsv_dds <- DESeq(brsv_dds)
map_dds <- DESeq(map_dds)

dispersionFunction(bohv1_dds) <- dispersionFunction(ddsTrain_adjusted)
dispersionFunction(brsv_dds) <- dispersionFunction(ddsTrain_adjusted)
dispersionFunction(map_dds) <- dispersionFunction(ddsTrain_adjusted)

bohv1_vst <- varianceStabilizingTransformation(bohv1_dds, blind = FALSE) # Now perform the normalisation
brsv_vst <- varianceStabilizingTransformation(brsv_dds, blind = FALSE) # Now perform the normalisation
map_vst <- varianceStabilizingTransformation(map_dds, blind = FALSE) # Now perform the normalisation
ensemble$gene_name
rownames(bohv1_vst) <- gsub(" ", "", rownames(bohv1_vst))
rownames(brsv_vst) <- gsub(" ", "", rownames(brsv_vst))
rownames(map_vst) <- gsub(" ", "", rownames(map_vst)) 


rownames(DE_genes) <- DE_genes$Symbol
rownames(DE_genes)
bohv1_filtered_counts <- as.data.frame(assay(bohv1_vst))
bohv1_filtered_counts <- bohv1_filtered_counts[DE_genes$Symbol, ]

brsv_filtered_counts <- as.data.frame(assay(brsv_vst))
brsv_filtered_counts <- brsv_filtered_counts[DE_genes$Symbol, ]
map_filtered_counts <- as.data.frame(assay(map_vst))
map_filtered_counts <- map_filtered_counts[DE_genes$Symbol, ]

bohv1_filtered_counts
map_filtered_counts




dim(brsv_filtered_counts)

train_normalised_filtered_counts["IFNAR1",]
brsv_filtered_counts["IFNAR1",]

# read in the models
GLM_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/GLM_models.rds")
GLMRIDGE_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/GLMRIDGE_model.rds")
GLMLASSO_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/GLMLASSO_model.rds")
GLMENET_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/GLMENET_model.rds")
RF_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/Random_forest_list_models.rds")

resamples_rf <- resamples(RF_model[1:6])

summary(resamples_rf)
# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF_model <- RF_model[[best_rf_model_name]]


RF_ET = readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/Random_forest_extratrees_list_models.rds")
resamples_rf <- resamples(RF_ET[1:6])

summary(resamples_rf)
# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF_ET_model <- RF_ET[[best_rf_model_name]]


NB_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/NB_model.rds")
MLP_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/MLP_model.rds")


NB_model
MLP_model

bohv1_filtered_counts["PALM2",]
bohv1_glm_results <- predict(GLM_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "GLM")
bohv1_glmridge_results <- predict(GLMRIDGE_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "RIDGE")
bohv1_glmlasso_results <- predict(GLMLASSO_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "LASSO")
bohv1_glmenet_results <- predict(GLMENET_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "ENET")
bohv1_RF_results <- predict(RF_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "RF")
bohv1_RF_ET_results <- predict(RF_ET_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "RF_ET")
bohv1_NB_results <- predict(NB_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "NB")
bohv1_MLP_results <- predict(MLP_model, newdata = t(bohv1_filtered_counts), type = "prob") %>% mutate(Model = "MLP")
bohv1_glmlasso_results

map_glm_results <- predict(GLM_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "GLM") 
map_glmridge_results <- predict(GLMRIDGE_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "RIDGE")
map_glmlasso_results <- predict(GLMLASSO_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "LASSO")
map_glmenet_results <- predict(GLMENET_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "ENET")
map_RF_results <- predict(RF_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "RF")
map_RF_ET_results <- predict(RF_ET_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "RF_ET")
map_NB_results <- predict(NB_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "NB")
map_MLP_results <- predict(MLP_model, newdata = t(map_filtered_counts), type = "prob")%>% mutate(Model = "MLP")


brsv_glm_results <- predict(GLM_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "GLM")
brsv_glmridge_results <- predict(GLMRIDGE_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "RIDGE")
brsv_glmlasso_results <- predict(GLMLASSO_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "LASSO")
brsv_glmenet_results <- predict(GLMENET_model, newdata = t(brsv_filtered_counts), type = "prob")  %>% mutate(Model = "ENET")
brsv_RF_results <- predict(RF_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "RF")
brsv_RF_ET_results <- predict(RF_ET_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "RF_ET")
brsv_NB_results <- predict(NB_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "NB")
brsv_MLP_results <- predict(MLP_model, newdata = t(brsv_filtered_counts), type = "prob") %>% mutate(Model = "MLP")



df_brsv_results = data.frame(rbind(brsv_glm_results,
                        brsv_glmridge_results,
                        brsv_glmlasso_results,
                        brsv_glmenet_results,
                        brsv_RF_results,
                        brsv_RF_ET_results,
                        brsv_MLP_results,
                        brsv_NB_results))

#colnames(df_brsv_results)[4] <- "TB_status"
#colnames(df_brsv_results)[5] <- "brsv_status"

ROC_test_combined(df_brsv_results, as.character(brsv_metadata$Status))

df_map_results = data.frame(rbind(map_glm_results,
                        map_glmridge_results,
                        map_glmlasso_results,
                        map_glmenet_results,
                        map_RF_results,
                        map_RF_ET_results,
                        map_MLP_results,
                        map_NB_results))

df_bohv1_results = data.frame(rbind(bohv1_glm_results,
                       bohv1_glmridge_results,
                       bohv1_glmlasso_results,
                       bohv1_glmenet_results,
                       bohv1_RF_results,
                       bohv1_RF_ET_results,
                       bohv1_MLP_results,
                       bohv1_NB_results))




ROC_test_combined(df_map_results, as.character(map_metadata$Status))

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/MAP_AUROC.pdf", width = 12, height = 12, dpi = 600)

ROC_test_combined(df_brsv_results, as.character(brsv_metadata$Status))

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/BRSV_AUROC.pdf", width = 12, height = 12, dpi = 600)

ROC_test_combined(df_bohv1_results, as.character(bohv1_metadata$Status))

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/BOHV1_AUROC.pdf", width = 12, height = 12, dpi = 600)



# Now onto score

DE_results = fread("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/DE_results_integrated.txt") %>% filter(padj < 0.05) %>% filter(baseMean > 100) %>% select(-V1) %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))

# Genes


lines <- readLines("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC1_genes.txt")
lines_fixed <- gsub("/t", "\t", lines)
writeLines(lines_fixed, "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC1_genes.txt")
ROC1 <- fread("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC1_genes.txt", sep = "\t")


lines <- readLines("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC2_genes.txt")
lines_fixed <- gsub("/t", "\t", lines)
writeLines(lines_fixed, "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC2_genes.txt")
ROC2 <- fread("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC2_genes.txt", sep = "\t")

GFS_1_genes <- ROC1$combination[13]
GFS_1_genes <- str_split_1(GFS_1_genes, pattern = "_")


GFS_2_genes <- ROC2$combination[17]
GFS_2_genes <- str_split_1(GFS_2_genes, pattern = "_")


Pos_genes1 = DE_results %>% filter(Symbol %in% GFS_1_genes & log2FoldChange > 0) %>% select(Symbol)
Neg_genes1 = DE_results %>% filter(Symbol %in% GFS_1_genes & log2FoldChange < 0) %>% select(Symbol)

Pos_genes2 = DE_results %>% filter(Symbol %in% GFS_2_genes & log2FoldChange > 0) %>% select(Symbol)
Neg_genes2 = DE_results %>% filter(Symbol %in% GFS_2_genes & log2FoldChange < 0) %>% select(Symbol)

Pos_genes3 = DE_results %>% filter(Symbol %in% c(GFS_1_genes,GFS_2_genes) & log2FoldChange > 0) %>% select(Symbol)
Neg_genes3 = DE_results %>% filter(Symbol %in% c(GFS_1_genes,GFS_2_genes) & log2FoldChange < 0) %>% select(Symbol)



map_vst <- assay(map_vst)
bohv1_vst <- assay(bohv1_vst)
brsv_vst <- assay(brsv_vst)

if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = map_vst[rownames(map_vst) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =map_vst[rownames(map_vst) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =map_vst[rownames(map_vst) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =map_vst[rownames(map_vst) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = map_vst[rownames(map_vst) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =map_vst[rownames(map_vst) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =map_vst[rownames(map_vst) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =map_vst[rownames(map_vst) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = map_vst[rownames(map_vst) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =map_vst[rownames(map_vst) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}



if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =map_vst[rownames(map_vst) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =map_vst[rownames(map_vst) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}

map_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
map_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
map_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result




if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = brsv_vst[rownames(brsv_vst) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =brsv_vst[rownames(brsv_vst) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =brsv_vst[rownames(brsv_vst) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =brsv_vst[rownames(brsv_vst) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}



# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = brsv_vst[rownames(brsv_vst) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =brsv_vst[rownames(brsv_vst) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =brsv_vst[rownames(brsv_vst) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =brsv_vst[rownames(brsv_vst) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = brsv_vst[rownames(brsv_vst) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =brsv_vst[rownames(brsv_vst) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =brsv_vst[rownames(brsv_vst) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =brsv_vst[rownames(brsv_vst) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}

brsv_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
brsv_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
brsv_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result




if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = bohv1_vst[rownames(bohv1_vst) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =bohv1_vst[rownames(bohv1_vst) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = bohv1_vst[rownames(bohv1_vst) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =bohv1_vst[rownames(bohv1_vst) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = bohv1_vst[rownames(bohv1_vst) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =bohv1_vst[rownames(bohv1_vst) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =bohv1_vst[rownames(bohv1_vst) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}


bohv1_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
bohv1_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
bohv1_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result



map_Scores_for_plotting1$Score <- map_Scores_for_plotting1$pos_Score
map_Scores_for_plotting2$Score <- map_Scores_for_plotting2$pos_Score - map_Scores_for_plotting2$neg_Score
map_Scores_for_plotting3$Score <- map_Scores_for_plotting3$pos_Score - map_Scores_for_plotting3$neg_Score

brsv_Scores_for_plotting1$Score <- brsv_Scores_for_plotting1$pos_Score
brsv_Scores_for_plotting2$Score <- brsv_Scores_for_plotting2$pos_Score - brsv_Scores_for_plotting2$neg_Score
brsv_Scores_for_plotting3$Score <- brsv_Scores_for_plotting3$pos_Score - brsv_Scores_for_plotting3$neg_Score


bohv1_Scores_for_plotting1$Score <- bohv1_Scores_for_plotting1$pos_Score
bohv1_Scores_for_plotting2$Score <- bohv1_Scores_for_plotting2$pos_Score - bohv1_Scores_for_plotting2$neg_Score
bohv1_Scores_for_plotting3$Score <- bohv1_Scores_for_plotting3$pos_Score - bohv1_Scores_for_plotting3$neg_Score

map_Scores_for_plotting1$Sample <- rownames(map_Scores_for_plotting1)
map_Scores_for_plotting2$Sample <- rownames(map_Scores_for_plotting2)
map_Scores_for_plotting3$Sample <- rownames(map_Scores_for_plotting3)


brsv_Scores_for_plotting1$Sample <- rownames(brsv_Scores_for_plotting1)
brsv_Scores_for_plotting2$Sample <- rownames(brsv_Scores_for_plotting2)
brsv_Scores_for_plotting3$Sample <- rownames(brsv_Scores_for_plotting3)


bohv1_Scores_for_plotting1$Sample <- rownames(bohv1_Scores_for_plotting1)
bohv1_Scores_for_plotting2$Sample <- rownames(bohv1_Scores_for_plotting2)
bohv1_Scores_for_plotting3$Sample <- rownames(bohv1_Scores_for_plotting3)


map_metadata$Sample <- rownames(map_metadata)


map_Scores_for_plotting1 <- left_join(map_Scores_for_plotting1, map_metadata)
map_Scores_for_plotting2 <- left_join(map_Scores_for_plotting2, map_metadata)
map_Scores_for_plotting3 <- left_join(map_Scores_for_plotting3, map_metadata)

brsv_Scores_for_plotting1 <- left_join(brsv_Scores_for_plotting1, brsv_metadata)
brsv_Scores_for_plotting2 <- left_join(brsv_Scores_for_plotting2, brsv_metadata)
brsv_Scores_for_plotting3 <- left_join(brsv_Scores_for_plotting3, brsv_metadata)

bohv1_Scores_for_plotting1 <- left_join(bohv1_Scores_for_plotting1, bohv1_metadata)
bohv1_Scores_for_plotting2 <- left_join(bohv1_Scores_for_plotting2, bohv1_metadata)
bohv1_Scores_for_plotting3 <- left_join(bohv1_Scores_for_plotting3, bohv1_metadata)





View(df_bohv1_results)
View(df_map_results)
View(df_brsv_results)


bohv1_glm_results <- predict(GLM_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "GLM", Sample = colnames(bohv1_filtered_counts))
bohv1_glmridge_results <- predict(GLMRIDGE_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Ridge", Sample = colnames(bohv1_filtered_counts))
bohv1_glmlasso_results <- predict(GLMLASSO_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Lasso", Sample = colnames(bohv1_filtered_counts))
bohv1_glmenet_results <- predict(GLMENET_model, newdata = t(bohv1_filtered_counts), type = "raw")  %>% data.frame() %>% mutate(Model = "Enet", Sample = colnames(bohv1_filtered_counts))
bohv1_RF_results <- predict(RF_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF", Sample = colnames(bohv1_filtered_counts))
bohv1_RF_ET_results <- predict(RF_ET_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF-ET", Sample = colnames(bohv1_filtered_counts))
bohv1_NB_results <- predict(NB_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "NB", Sample = colnames(bohv1_filtered_counts))
bohv1_MLP_results <- predict(MLP_model, newdata = t(bohv1_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "MLP", Sample = colnames(bohv1_filtered_counts))

bohv1_results = rbind(bohv1_glm_results, bohv1_glmridge_results, bohv1_glmlasso_results, bohv1_glmenet_results, bohv1_RF_results, bohv1_RF_ET_results, bohv1_NB_results, bohv1_MLP_results)
colnames(bohv1_results)[1] <- "Bohv1_prediction"
bohv1_results <- bohv1_results %>% left_join(., rownames_to_column(bohv1_metadata[,c("Status","TB_status")]), by=c("Sample" = "rowname"))
bohv1_results$Model <- factor(bohv1_results$Model,levels = c("GLM","Ridge","Lasso","Enet", "RF", "RF-ET","MLP", "NB"))
bohv1_results$Status <- factor(bohv1_results$Status)
bohv1_results$Bohv1_prediction <- factor(bohv1_results$Bohv1_prediction)


bohv1_results
map_glm_results <- predict(GLM_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "GLM", Sample = colnames(map_filtered_counts))
map_glmridge_results <- predict(GLMRIDGE_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Ridge", Sample = colnames(map_filtered_counts))
map_glmlasso_results <- predict(GLMLASSO_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Lasso", Sample = colnames(map_filtered_counts))
map_glmenet_results <- predict(GLMENET_model, newdata = t(map_filtered_counts), type = "raw")  %>% data.frame() %>% mutate(Model = "Enet", Sample = colnames(map_filtered_counts))
map_RF_results <- predict(RF_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF", Sample = colnames(map_filtered_counts))
map_RF_ET_results <- predict(RF_ET_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF-ET", Sample = colnames(map_filtered_counts))
map_NB_results <- predict(NB_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "NB", Sample = colnames(map_filtered_counts))
map_MLP_results <- predict(MLP_model, newdata = t(map_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "MLP", Sample = colnames(map_filtered_counts))

map_results = rbind(map_glm_results, map_glmridge_results, map_glmlasso_results, map_glmenet_results, map_RF_results, map_RF_ET_results, map_NB_results, map_MLP_results)
colnames(map_results)[1] <- "Map_prediction"
map_metadata
map_results <- map_results %>% left_join(., rownames_to_column(map_metadata[,c("Status","TB_status")]), by=c("Sample" = "rowname"))
map_results$Model <- factor(map_results$Model,levels = c("GLM","Ridge","Lasso","Enet", "RF", "RF-ET","MLP", "NB"))
map_results$Status <- factor(map_results$Status)
map_results$Map_prediction <- factor(map_results$Map_prediction)


brsv_glm_results <- predict(GLM_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "GLM", Sample = colnames(brsv_filtered_counts))
brsv_glmridge_results <- predict(GLMRIDGE_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Ridge", Sample = colnames(brsv_filtered_counts))
brsv_glmlasso_results <- predict(GLMLASSO_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "Lasso", Sample = colnames(brsv_filtered_counts))
brsv_glmenet_results <- predict(GLMENET_model, newdata = t(brsv_filtered_counts), type = "raw")  %>% data.frame() %>% mutate(Model = "Enet", Sample = colnames(brsv_filtered_counts))
brsv_RF_results <- predict(RF_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF", Sample = colnames(brsv_filtered_counts))
brsv_RF_ET_results <- predict(RF_ET_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "RF-ET", Sample = colnames(brsv_filtered_counts))
brsv_NB_results <- predict(NB_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "NB", Sample = colnames(brsv_filtered_counts))
brsv_MLP_results <- predict(MLP_model, newdata = t(brsv_filtered_counts), type = "raw") %>% data.frame() %>% mutate(Model = "MLP", Sample = colnames(brsv_filtered_counts))
brsv_MLP_results
brsv_results = rbind(brsv_glm_results, brsv_glmridge_results, brsv_glmlasso_results, brsv_glmenet_results, brsv_RF_results, brsv_RF_ET_results, brsv_NB_results, brsv_MLP_results)
colnames(brsv_results)[1] <- "Brsv_prediction"
brsv_results <- brsv_results %>% left_join(., rownames_to_column(brsv_metadata[,c("Status","TB_status")]), by=c("Sample" = "rowname"))
brsv_results$Model <- factor(brsv_results$Model,levels = c("GLM","Ridge","Lasso","Enet", "RF", "RF-ET","MLP", "NB"))
brsv_results$Status <- factor(brsv_results$Status)
brsv_results$Brsv_prediction <- factor(brsv_results$Brsv_prediction)
brsv_results


df_brsv_results = data.frame(brsv_glm_results,
                             brsv_glmridge_results,
                             brsv_glmlasso_results,
                             brsv_glmenet_results,
                             brsv_RF_results,
                             brsv_RF_ET_results,
                             brsv_MLP_results,
                             brsv_NB_results,
                             brsv_metadata$TB_status,
                             brsv_metadata$Status)
colnames(df_brsv_results)[9] <- "TB_status"
colnames(df_brsv_results)[10] <- "brsv_status"

df_map_results = data.frame(map_glm_results,
                            map_glmridge_results,
                            map_glmlasso_results,
                            map_glmenet_results,
                            map_RF_results,
                            map_RF_ET_results,
                            map_MLP_results,
                            map_NB_results,
                            map_metadata$TB_status,
                            map_metadata$Status)
colnames(df_map_results)[9] <- "TB_status"
colnames(df_map_results)[10] <- "TB_status"

df_bohv1_results = data.frame(bohv1_glm_results,
                              bohv1_glmridge_results,
                              bohv1_glmlasso_results,
                              bohv1_glmenet_results,
                              bohv1_RF_results,
                              bohv1_RF_ET_results,
                              bohv1_MLP_results,
                              bohv1_NB_results,
                              bohv1_metadata$TB_status,
                              bohv1_metadata$Status)
colnames(df_bohv1_results)[9] <- "TB_status"
colnames(df_bohv1_results)[10] <- "bohv1_status"




brsv_results_summary = brsv_results %>% group_by(Model) %>% summarize(
  Dataset = "BRSV",
  Accuracy = as.vector(confusionMatrix(data = Brsv_prediction, reference = Status, positive = "Infected")$overall[1]),
  Sens = as.vector(confusionMatrix(data = Brsv_prediction, reference = Status, positive = "Infected")$byClass[1]),
  Spec = as.vector(confusionMatrix(data = Brsv_prediction, reference = Status, positive = "Infected")$byClass[2]),
  Accuracy_TB = as.vector(sum(Brsv_prediction == TB_status, na.rm = T) / n() *100))




map_results_summary = map_results %>% group_by(Model) %>% summarize(
  Dataset = "MAP",
  Accuracy = as.vector(confusionMatrix(data = Map_prediction, reference = Status, positive = "Infected")$overall[1]),
  Sens = as.vector(confusionMatrix(data = Map_prediction, reference = Status, positive = "Infected")$byClass[1]),
  Spec = as.vector(confusionMatrix(data = Map_prediction, reference = Status, positive = "Infected")$byClass[2]),
  Accuracy_TB = as.vector(sum(Map_prediction == TB_status, na.rm = T) / n() *100))


bohv1_results_summary = bohv1_results %>% group_by(Model) %>% summarize(
  Dataset = "BoHV-1",
  Accuracy = as.vector(confusionMatrix(data = Bohv1_prediction, reference = Status, positive = "Infected")$overall[1]),
  Sens = as.vector(confusionMatrix(data = Bohv1_prediction, reference = Status, positive = "Infected")$byClass[1]),
  Spec = as.vector(confusionMatrix(data = Bohv1_prediction, reference = Status, positive = "Infected")$byClass[2]),
  Accuracy_TB = as.vector(sum(Bohv1_prediction == TB_status, na.rm = T) / n() *100))




External_summary <- rbind(bohv1_results_summary, brsv_results_summary, map_results_summary)

External_summary <- External_summary %>% pivot_wider(values_from = c("Accuracy", "Sens", "Spec", "Accuracy_TB"), names_from = Model)

View(External_summary)

View(brsv_results)
brsv_results_summary
Accuracy = confusionMatrix(brsv_results$Brsv_prediction, reference = brsv_results$Status, positive = "Infected")
Accuracy$byClass

BOHV1 = ROC_test_combined(bohv1_results, bohv1_dds$Status)
BOHV1



AUC = as.numeric(pROC::auc(pROC::roc(bohv1_dds$Status, as.numeric(bohv1_glmenet_results$Infected, plot = TRUE))))
AUC

cbind(brsv_dds$Status, bohv1_glmenet_results$Infected)
