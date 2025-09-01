
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
library(rstatix)
args = commandArgs(trailingOnly = TRUE)
source(args[1])
set.seed(42) # For reproducibility


# Ensemble
ensemble <- fread(args[2])
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

bohv1_matrix = read.csv(args[3], sep = "\t", row.names = 1, check.names = FALSE)
brsv_matrix = read.csv(args[4], sep = "\t", row.names = 1, check.names = FALSE)
map_matrix = read.csv(args[5], sep = "\t", row.names = 1, check.names = FALSE)

all(rownames(bohv1_matrix) == ensemble$gene_id)
rownames(bohv1_matrix) <- ensemble$gene_name 
rownames(brsv_matrix) <- ensemble$gene_name
rownames(map_matrix) <- ensemble$gene_name
rownames(map_matrix)

bohv1_metadata = read.csv(args[6], sep = "\t", row.names = "Animal_Code", check.names = FALSE)
brsv_metadata = read.csv(args[7], sep = "\t", row.names = "Animal_Code", check.names = FALSE)
map_metadata = read.csv(args[8], sep = "\t", row.names = "Animal_Code", check.names = FALSE)

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

ddsTrain_adjusted <- readRDS(args[9])
DE_genes <- fread(args[10]) %>%
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






# read in the models
GLM_model <- readRDS(args[11])
GLMRIDGE_model <- readRDS(args[12])
GLMLASSO_model <- readRDS(args[13])
GLMENET_model <- readRDS(args[14])
RF_model <- readRDS(args[15])

resamples_rf <- resamples(RF_model[1:6])

summary(resamples_rf)
# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF_model <- RF_model[[best_rf_model_name]]


RF_ET = readRDS(args[16])
resamples_rf <- resamples(RF_ET[1:6])

summary(resamples_rf)
# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF_ET_model <- RF_ET[[best_rf_model_name]]


NB_model <- readRDS(args[17])
MLP_model <- readRDS(args[18])


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

ggsave(args[19], width = 12, height = 12, dpi = 600)

ROC_test_combined(df_brsv_results, as.character(brsv_metadata$Status))

ggsave(args[20], width = 12, height = 12, dpi = 600)

ROC_test_combined(df_bohv1_results, as.character(bohv1_metadata$Status))

ggsave(args[21], width = 12, height = 12, dpi = 600)



# Now onto score

DE_results = fread(args[10]) %>% filter(padj < 0.05) %>% filter(baseMean > 100) %>% select(-V1) %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))

# Genes



ROC1 <- fread(args[22], sep = "\t")


ROC2 <- fread(args[23], sep = "\t")

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


train_data = as.data.frame(t(read.table(args[24], sep = "\t", check.names=F)))
train_labels = fread(args[25])

if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = train_data[rownames(train_data) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =train_data[rownames(train_data) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =train_data[rownames(train_data) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes1$Symbol)>1) {
  Neg_counts_ROC1 =train_data[rownames(train_data) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = train_data[rownames(train_data) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =train_data[rownames(train_data) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =train_data[rownames(train_data) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =train_data[rownames(train_data) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = train_data[rownames(train_data) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =train_data[rownames(train_data) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =train_data[rownames(train_data) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =train_data[rownames(train_data) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}



Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result

Scores_for_plotting1$Score <- Scores_for_plotting1$pos_Score

# Calculate Score
Scores_for_plotting2$Score <- Scores_for_plotting2$pos_Score - Scores_for_plotting2$neg_Score
#Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))

Scores_for_plotting3$Score <- Scores_for_plotting3$pos_Score - Scores_for_plotting3$neg_Score
#Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))

Scores_for_plotting1
Scores_for_plotting2
Scores_for_plotting3

#fold_assignments

Scores_for_plotting1$Sample <- rownames(Scores_for_plotting1)
Scores_for_plotting2$Sample <- rownames(Scores_for_plotting2)
Scores_for_plotting3$Sample <- rownames(Scores_for_plotting3)


Scores_for_plotting1 <- left_join(Scores_for_plotting1, train_labels)
Scores_for_plotting2 <- left_join(Scores_for_plotting2, train_labels)
Scores_for_plotting3 <- left_join(Scores_for_plotting3, train_labels)


Scores_for_plotting1$Study <- gsub("1_", "", Scores_for_plotting1$Study)
Scores_for_plotting1$Study <- gsub("2_", "", Scores_for_plotting1$Study)


Scores_for_plotting2$Study <- gsub("1_", "", Scores_for_plotting2$Study)
Scores_for_plotting2$Study <- gsub("2_", "", Scores_for_plotting2$Study)


Scores_for_plotting3$Study <- gsub("1_", "", Scores_for_plotting3$Study)
Scores_for_plotting3$Study <- gsub("2_", "", Scores_for_plotting3$Study)



Scores_for_plotting1$Study <- factor(Scores_for_plotting1$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Scores_for_plotting2$Study <- factor(Scores_for_plotting2$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Scores_for_plotting3$Study <- factor(Scores_for_plotting3$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))



Scores_for_plotting1$Condition <- factor(Scores_for_plotting1$Condition, levels = c(0,1), labels = c("Control", "Infected"))
Scores_for_plotting2$Condition <- factor(Scores_for_plotting2$Condition, levels = c(0,1), labels = c("Control", "Infected"))
Scores_for_plotting3$Condition <- factor(Scores_for_plotting3$Condition, levels = c(0,1), labels = c("Control", "Infected"))

Scores_for_plotting1$Set = "Set1"
Scores_for_plotting2$Set = "Set2"
Scores_for_plotting3$Set = "Set3"

Scores_for_plotting1$Unscaled_Score <- as.numeric(Scores_for_plotting1$Score)
Scores_for_plotting2$Unscaled_Score <- as.numeric(Scores_for_plotting2$Score)
Scores_for_plotting3$Unscaled_Score <- as.numeric(Scores_for_plotting3$Score)



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
} else if (length(Neg_genes1$Symbol)>1) {
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
} else if (length(Neg_genes1$Symbol)>1) {
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
} else if (length(Neg_genes1$Symbol)>1) {
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

# Scale
map_Scores_for_plotting1$Score <- as.numeric((map_Scores_for_plotting1$Score - mean(Scores_for_plotting1$Unscaled_Score)) / sd(Scores_for_plotting1$Unscaled_Score))
map_Scores_for_plotting2$Score <- as.numeric((map_Scores_for_plotting2$Score - mean(Scores_for_plotting2$Unscaled_Score)) / sd(Scores_for_plotting2$Unscaled_Score))
map_Scores_for_plotting3$Score <- as.numeric((map_Scores_for_plotting3$Score - mean(Scores_for_plotting3$Unscaled_Score)) / sd(Scores_for_plotting3$Unscaled_Score))



brsv_Scores_for_plotting1$Score <- brsv_Scores_for_plotting1$pos_Score
brsv_Scores_for_plotting2$Score <- brsv_Scores_for_plotting2$pos_Score - brsv_Scores_for_plotting2$neg_Score
brsv_Scores_for_plotting3$Score <- brsv_Scores_for_plotting3$pos_Score - brsv_Scores_for_plotting3$neg_Score

brsv_Scores_for_plotting1$Score <- as.numeric((brsv_Scores_for_plotting1$Score - mean(Scores_for_plotting1$Unscaled_Score)) / sd(Scores_for_plotting1$Unscaled_Score))
brsv_Scores_for_plotting2$Score <- as.numeric((brsv_Scores_for_plotting2$Score - mean(Scores_for_plotting2$Unscaled_Score)) / sd(Scores_for_plotting2$Unscaled_Score))
brsv_Scores_for_plotting3$Score <- as.numeric((brsv_Scores_for_plotting3$Score - mean(Scores_for_plotting3$Unscaled_Score)) / sd(Scores_for_plotting3$Unscaled_Score))





bohv1_Scores_for_plotting1$Score <- bohv1_Scores_for_plotting1$pos_Score 
bohv1_Scores_for_plotting2$Score <- bohv1_Scores_for_plotting2$pos_Score - bohv1_Scores_for_plotting2$neg_Score
bohv1_Scores_for_plotting3$Score <- bohv1_Scores_for_plotting3$pos_Score - bohv1_Scores_for_plotting3$neg_Score

bohv1_Scores_for_plotting1$Score <- as.numeric((bohv1_Scores_for_plotting1$Score - mean(Scores_for_plotting1$Unscaled_Score)) / sd(Scores_for_plotting1$Unscaled_Score))
bohv1_Scores_for_plotting2$Score <- as.numeric((bohv1_Scores_for_plotting2$Score - mean(Scores_for_plotting2$Unscaled_Score)) / sd(Scores_for_plotting2$Unscaled_Score))
bohv1_Scores_for_plotting3$Score <- as.numeric((bohv1_Scores_for_plotting3$Score - mean(Scores_for_plotting3$Unscaled_Score)) / sd(Scores_for_plotting3$Unscaled_Score))





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
brsv_metadata$Sample <- rownames(brsv_metadata)
bohv1_metadata$Sample <- rownames(bohv1_metadata)



map_Scores_for_plotting1 <- left_join(map_Scores_for_plotting1, map_metadata)
map_Scores_for_plotting2 <- left_join(map_Scores_for_plotting2, map_metadata)
map_Scores_for_plotting3 <- left_join(map_Scores_for_plotting3, map_metadata)

brsv_Scores_for_plotting1 <- left_join(brsv_Scores_for_plotting1, brsv_metadata)
brsv_Scores_for_plotting2 <- left_join(brsv_Scores_for_plotting2, brsv_metadata)
brsv_Scores_for_plotting3 <- left_join(brsv_Scores_for_plotting3, brsv_metadata)

bohv1_Scores_for_plotting1 <- left_join(bohv1_Scores_for_plotting1, bohv1_metadata)
bohv1_Scores_for_plotting2 <- left_join(bohv1_Scores_for_plotting2, bohv1_metadata)
bohv1_Scores_for_plotting3 <- left_join(bohv1_Scores_for_plotting3, bohv1_metadata)


common_cols = intersect(colnames(bohv1_Scores_for_plotting1), intersect(colnames(brsv_Scores_for_plotting1), colnames(bohv1_Scores_for_plotting1)))

map_Scores_for_plotting1 <- map_Scores_for_plotting1 %>% select(common_cols)
map_Scores_for_plotting2 <- map_Scores_for_plotting2 %>% select(common_cols)
map_Scores_for_plotting3 <- map_Scores_for_plotting3 %>% select(common_cols)

brsv_Scores_for_plotting1 <- brsv_Scores_for_plotting1 %>% select(common_cols)
brsv_Scores_for_plotting2 <- brsv_Scores_for_plotting2 %>% select(common_cols)
brsv_Scores_for_plotting3 <- brsv_Scores_for_plotting3 %>% select(common_cols)


bohv1_Scores_for_plotting1 <- bohv1_Scores_for_plotting1 %>% select(common_cols)
bohv1_Scores_for_plotting2 <- bohv1_Scores_for_plotting2 %>% select(common_cols)
bohv1_Scores_for_plotting3 <- bohv1_Scores_for_plotting3 %>% select(common_cols)


map_Scores_for_plotting1$Disease = "MAP"
map_Scores_for_plotting2$Disease = "MAP"
map_Scores_for_plotting3$Disease = "MAP"


brsv_Scores_for_plotting1$Disease = "brsv"
brsv_Scores_for_plotting2$Disease = "brsv"
brsv_Scores_for_plotting3$Disease = "brsv"


bohv1_Scores_for_plotting1$Disease = "bohv1"
bohv1_Scores_for_plotting2$Disease = "bohv1"
bohv1_Scores_for_plotting3$Disease = "bohv1"


map_Scores_for_plotting1$Set = "Pass 1"
map_Scores_for_plotting2$Set = "Pass 2"
map_Scores_for_plotting3$Set = "Combined"


brsv_Scores_for_plotting1$Set = "Pass 1"
brsv_Scores_for_plotting2$Set = "Pass 2"
brsv_Scores_for_plotting3$Set = "Combined"


bohv1_Scores_for_plotting1$Set = "Pass 1"
bohv1_Scores_for_plotting2$Set = "Pass 2"
bohv1_Scores_for_plotting3$Set = "Combined"





bohv1_Scores_for_plotting1



External_scores_plotting <- rbind(map_Scores_for_plotting1,
                                  map_Scores_for_plotting2,
                                  map_Scores_for_plotting3,
                                  brsv_Scores_for_plotting1,
                                  brsv_Scores_for_plotting2,
                                  brsv_Scores_for_plotting3,
                                  bohv1_Scores_for_plotting1,
                                  bohv1_Scores_for_plotting2,
                                  bohv1_Scores_for_plotting3)

External_scores_plotting
External_scores_plotting$Set <- factor(External_scores_plotting$Set, levels = c("Pass 1", "Pass 2", "Combined"))
External_scores_plotting$Disease <- factor(External_scores_plotting$Disease, levels = c("MAP", "bohv1", "brsv"))


ggplot(External_scores_plotting, aes(y = Score, x = Disease, fill = Set, shape = Status, col = Status)) + 
  geom_boxplot(position = position_dodge2(width = 2, padding = 1.9), outlier.colour = NA, alpha = 1) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5, size = 3) + 
  scale_colour_manual(values = c("#2166ac", "#b2182b")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  xlab("Study") +
  ylab("Score") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        #legend.position = c(.9, .25),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        #legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14)) #change legend title font size

ggsave(args[26], width = 12, height = 12, dpi = 600)

External_scores_plotting %>%
  group_by(Disease, Set) %>%
  wilcox_test(Score ~ Status) %>%
  adjust_pvalue(method = "BH") %>%  # Optional: adjust p-values
  add_significance() %>% arrange(desc(Disease)) %>% write.table(args[27])  


#Disease Set      .y.   group1  group2      n1    n2 statistic        p    p.adj p.adj.signif
#<fct>   <fct>    <chr> <chr>   <chr>    <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#  1 brsv    Pass 1   Score Control Infected     6    12        33 0.82     0.82     ns          
#2 brsv    Pass 2   Score Control Infected     6    12         1 0.000215 0.000967 ***         
#  3 brsv    Combined Score Control Infected     6    12        28 0.494    0.635    ns          
#4 bohv1   Pass 1   Score Control Infected     6    12         1 0.000215 0.000967 ***         
#  5 bohv1   Pass 2   Score Control Infected     6    12         2 0.000431 0.00129  **          
#  6 bohv1   Combined Score Control Infected     6    12        10 0.0135   0.0304   *           
#  7 MAP     Pass 1   Score Control Infected     3    11        20 0.659    0.741    ns          
#8 MAP     Pass 2   Score Control Infected     3    11        10 0.368    0.635    ns          
#9 MAP     Combined Score Control Infected     3    11        11 0.456    0.635    ns   







GFS_map_roc <- External_scores_plotting %>% filter(Disease == "MAP") %>% group_by(Set) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Status, as.numeric(Score), direction = "<"))),
                                                                             AUC_CI_lower = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[1],
                                                                             AUC_CI_upper = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[3]) %>% mutate(set = "Set1")

GFS_map_roc
GFS_brsv_roc <- External_scores_plotting %>% filter(Disease == "brsv") %>% group_by(Set) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Status, as.numeric(Score), direction = "<"))),
                                                                              AUC_CI_lower = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[1],
                                                                              AUC_CI_upper = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[3]) %>% mutate(set = "Set2")




GFS_bohv1_roc <- External_scores_plotting %>% filter(Disease == "bohv1") %>% group_by(Set) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Status, as.numeric(Score), direction = "<"))),
                                                                                                       AUC_CI_lower = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[1],
                                                                                                       AUC_CI_upper = (pROC::ci(pROC::roc(Status, as.numeric(Score), direction = "<")))[3]) %>% mutate(set = "Set2")
GFS_map_roc

GFS_bohv1_roc

GFS_map_roc <- GFS_map_roc %>% arrange(desc(AUC))

roc_text <- paste(paste0("ROC ", seq_along(GFS_map_roc$Set), " = ", GFS_map_roc$AUC, "(95% CI:", round(GFS_map_roc$AUC_CI_lower, 2), " - ", round(GFS_map_roc$AUC_CI_upper,2), ")"), collapse = "\n")
resample_colors <- c("#fc8d62","#8da0cb","#66c2a5")

GFS_map_score = ggplot(External_scores_plotting %>% filter(Disease == "MAP"), aes(m=Score, d=factor(Status, levels = c("Control", "Infected")), colour = factor(Set, levels = GFS_map_roc$Set))) + # Note do not need to filter unless saving all predictions
  style_roc() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
  geom_roc(n.cuts=0) + 
  coord_equal() +
  scale_x_continuous("1 - Specificity (FPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous("Sensitivity (TPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + 
  scale_colour_manual(values = resample_colors)

for (i in seq_along(GFS_map_roc$AUC)) {
  GFS_map_score <- GFS_map_score + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("AUROC ", GFS_map_roc$Set[i], " = ", round(GFS_map_roc$AUC[i],3), " (", round(GFS_map_roc$AUC_CI_lower[i],2), " - ", round(GFS_map_roc$AUC_CI_upper[i],2), ")"), 
             hjust = .25, 
             size = 5,
             col = resample_colors[i])
}

GFS_map_score
ggsave(args[28], width = 12, height = 12, dpi = 600)


ggsave 


# brsv
GFS_brsv_roc <- GFS_brsv_roc %>% arrange(desc(AUC))
GFS_map_roc
roc_text <- paste(paste0("ROC ", seq_along(GFS_brsv_roc$Set), " = ", GFS_brsv_roc$AUC, "(95% CI:", round(GFS_brsv_roc$AUC_CI_lower, 2), " - ", round(GFS_brsv_roc$AUC_CI_upper,2), ")"), collapse = "\n")
resample_colors <- c("#fc8d62","#8da0cb","#66c2a5")

GFS_brsv_score = ggplot(External_scores_plotting %>% filter(Disease == "brsv"), aes(m=Score, d=factor(Status, levels = c("Control", "Infected")), colour = factor(Set, levels = GFS_brsv_roc$Set))) + # Note do not need to filter unless saving all predictions
  style_roc() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
  geom_roc(n.cuts=0) + 
  coord_equal() +
  scale_x_continuous("1 - Specificity (FPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous("Sensitivity (TPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + 
  scale_colour_manual(values = resample_colors)

for (i in seq_along(GFS_brsv_roc$AUC)) {
  GFS_brsv_score <- GFS_brsv_score + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("AUROC ", GFS_brsv_roc$Set[i], " = ", round(GFS_brsv_roc$AUC[i],3), " (", round(GFS_brsv_roc$AUC_CI_lower[i],2), " - ", round(GFS_brsv_roc$AUC_CI_upper[i],2), ")"), 
             hjust = .25, 
             size = 5,
             col = resample_colors[i])
}

GFS_brsv_score
ggsave(args[29], width = 12, height = 12, dpi = 600)



# bohv1

GFS_bohv1_roc <- GFS_bohv1_roc %>% arrange(desc(AUC))

roc_text <- paste(paste0("ROC ", seq_along(GFS_bohv1_roc$Set), " = ", GFS_bohv1_roc$AUC, "(95% CI:", round(GFS_bohv1_roc$AUC_CI_lower, 2), " - ", round(GFS_bohv1_roc$AUC_CI_upper,2), ")"), collapse = "\n")
resample_colors <- c("#66c2a5","#fc8d62","#8da0cb")

GFS_bohv1_score = ggplot(External_scores_plotting %>% filter(Disease == "bohv1"), aes(m=Score, d=factor(Status, levels = c("Control", "Infected")), colour = factor(Set, levels = GFS_bohv1_roc$Set))) + # Note do not need to filter unless saving all predictions
  style_roc() +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
  geom_roc(n.cuts=0) + 
  coord_equal() +
  scale_x_continuous("1 - Specificity (FPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  scale_y_continuous("Sensitivity (TPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + 
  scale_colour_manual(values = resample_colors)

for (i in seq_along(GFS_bohv1_roc$AUC)) {
  GFS_bohv1_score <- GFS_bohv1_score + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("AUROC ", GFS_bohv1_roc$Set[i], " = ", round(GFS_bohv1_roc$AUC[i],3), " (", round(GFS_bohv1_roc$AUC_CI_lower[i],2), " - ", round(GFS_bohv1_roc$AUC_CI_upper[i],2), ")"), 
             hjust = .25, 
             size = 5,
             col = resample_colors[i])
}

GFS_bohv1_score

ggsave(args[30], width = 12, height = 12, dpi = 600)








####################################
## M.bovis versus other diseases
###################################

test_data = fread(file = args[31])
test_labels = fread(args[32])

test_infec =  test_labels %>% filter(Sample %in% colnames(test_data)) %>% filter(Condition == "Infected") %>% mutate(TB_status = "Infected")
rownames(test_infec)
rownames(test_infec) <- test_infec$Sample
colnames(test_data)
test_data <- test_data %>% select(rownames(test_infec))

map_infec = map_metadata %>% filter(Status == "Infected") %>% select(Status, TB_status)
bohv1_infec = as.data.frame(bohv1_metadata) %>% select(Status, TB_status)
bohv1_infec = bohv1_infec %>% filter(Status == "Infected")
brsv_infec = brsv_metadata %>% filter(Status == "Infected") %>% select(Status, TB_status)

test_infec <- test_infec %>% select(Condition, TB_status)

tb_od = cbind(test_data, map_matrix[,rownames(map_infec)], bohv1_matrix[,rownames(bohv1_infec)], brsv_matrix[,rownames(brsv_infec)])

colnames(test_infec)[1] <- "Status"

tb_od_metadata = rbind(test_infec, map_infec, bohv1_infec, brsv_infec)

tb_od_metadata

ddsTest_adjusted  <- DESeqDataSetFromMatrix(countData = tb_od, 
                                            colData = tb_od_metadata, 
                                            design = ~ 1)



# training normalised dataset
ddsTrain_adjusted = readRDS(file = args[9])

#write.table(test_counts,file = "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/Test_counts.txt", row.names = F, quote = F)



ddsTest_adjusted <- DESeq(ddsTest_adjusted)




# Apply the dispersion function on dds test
dispersionFunction(ddsTest_adjusted) <- dispersionFunction(ddsTrain_adjusted)
vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest_adjusted, blind = FALSE) # Now perform the normalisation

tb_od_vst <- assay(vstNormalizedExpressionDataForTest)
head(tb_od_vst)

rownames(tb_od_vst) <- rownames(map_vst)

tb_od_vst_filtered_counts <- tb_od_vst[DE_genes$Symbol,]






tb_od_vst_glm_results <- predict(GLM_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "GLM")
tb_od_vst_glmridge_results <- predict(GLMRIDGE_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "RIDGE")
tb_od_vst_glmlasso_results <- predict(GLMLASSO_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "LASSO")
tb_od_vst_glmenet_results <- predict(GLMENET_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "ENET")
tb_od_vst_RF_results <- predict(RF_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "RF")
tb_od_vst_RF_ET_results <- predict(RF_ET_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "RF_ET")
tb_od_vst_NB_results <- predict(NB_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "NB")
tb_od_vst_MLP_results <- predict(MLP_model, newdata = t(tb_od_vst_filtered_counts), type = "prob") %>% mutate(Model = "MLP")





# Search


if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes1$Symbol)>1) {
  Neg_counts_ROC1 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}



# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =tb_od_vst_filtered_counts[rownames(tb_od_vst_filtered_counts) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}

tb_od_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
tb_od_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
tb_od_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result




tb_od_Scores_for_plotting1$Set = "Pass 1"
tb_od_Scores_for_plotting2$Set = "Pass 2"
tb_od_Scores_for_plotting3$Set = "Combined"

tb_od_Scores_for_plotting1$Score <- tb_od_Scores_for_plotting1$pos_Score
tb_od_Scores_for_plotting2$Score <- tb_od_Scores_for_plotting2$pos_Score - tb_od_Scores_for_plotting2$neg_Score
tb_od_Scores_for_plotting3$Score <- tb_od_Scores_for_plotting3$pos_Score - tb_od_Scores_for_plotting3$neg_Score



tb_od_Scores_for_plotting1$Score <- as.numeric((tb_od_Scores_for_plotting1$Score - mean(Scores_for_plotting1$Unscaled_Score)) / sd(Scores_for_plotting1$Unscaled_Score))
tb_od_Scores_for_plotting2$Score <- as.numeric((tb_od_Scores_for_plotting2$Score - mean(Scores_for_plotting2$Unscaled_Score)) / sd(Scores_for_plotting2$Unscaled_Score))
tb_od_Scores_for_plotting3$Score <- as.numeric((tb_od_Scores_for_plotting3$Score - mean(Scores_for_plotting3$Unscaled_Score)) / sd(Scores_for_plotting3$Unscaled_Score))

colnames(tb_od_Scores_for_plotting1)[1] <- "Control" # Dummy
colnames(tb_od_Scores_for_plotting1)[3] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(tb_od_Scores_for_plotting1)[2] <- "Model"

colnames(tb_od_Scores_for_plotting2)[1] <- "Control" # Dummy
colnames(tb_od_Scores_for_plotting2)[4] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(tb_od_Scores_for_plotting2)[3] <- "Model"

colnames(tb_od_Scores_for_plotting3)[1] <- "Control" # Dummy
colnames(tb_od_Scores_for_plotting3)[4] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(tb_od_Scores_for_plotting3)[3] <- "Model"

tb_od_Scores_for_plotting1 <- tb_od_Scores_for_plotting1 %>% select(Control, Infected, Model)
tb_od_Scores_for_plotting2 <- tb_od_Scores_for_plotting2 %>% select(Control, Infected, Model)
tb_od_Scores_for_plotting3 <- tb_od_Scores_for_plotting3 %>% select(Control, Infected, Model)





df_tb_od_vst_results = data.frame(rbind(tb_od_vst_glm_results,
                                   tb_od_vst_glmridge_results,
                                   tb_od_vst_glmlasso_results,
                                   tb_od_vst_glmenet_results,
                                   tb_od_vst_RF_results,
                                   tb_od_vst_RF_ET_results,
                                   tb_od_vst_MLP_results,
                                   tb_od_vst_NB_results,
                                   tb_od_Scores_for_plotting1,
                                   tb_od_Scores_for_plotting2,
                                   tb_od_Scores_for_plotting3))

df_tb_od_vst_results$Model = factor(df_tb_od_vst_results$Model)

ROC_test_combined(df_tb_od_vst_results, tb_od_metadata$TB_status)







ggsave(args[33], width = 12, height = 12, dpi = 600)


