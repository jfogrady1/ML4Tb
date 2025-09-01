# Remake heatmap 2 (GFS genes)


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



DE_results = fread(args[2]) %>% filter(padj < 0.05) %>% filter(baseMean > 100) %>% select(-V1) %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))

# Genes

Test_data = as.data.frame(t(read.table(args[3], sep = "\t", check.names=F)))
test_labels = fread(args[4])


ROC1 <- fread(args[5], sep = "\t")


ROC2 <- fread(args[6], sep = "\t")

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



##########

# Training

#########

train_data = as.data.frame(t(read.table(args[7], sep = "\t", check.names=F)))
train_labels = fread(args[8])

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
} else if (length(Neg_genes1)>1) {
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


Scores_for_plotting1$Score <- as.numeric(scale(Scores_for_plotting1$Score))
Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))
Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))


Scores_for_plotting2 <- Scores_for_plotting2[,colnames(Scores_for_plotting1)]
Scores_for_plotting3 <- Scores_for_plotting3[,colnames(Scores_for_plotting1)]

#######
#Test
#######




if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = Test_data[rownames(Test_data) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =Test_data[rownames(Test_data) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =Test_data[rownames(Test_data) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes1)>1) {
  Neg_counts_ROC1 =Test_data[rownames(Test_data) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = Test_data[rownames(Test_data) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =Test_data[rownames(Test_data) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =Test_data[rownames(Test_data) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =Test_data[rownames(Test_data) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = Test_data[rownames(Test_data) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =Test_data[rownames(Test_data) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =Test_data[rownames(Test_data) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =Test_data[rownames(Test_data) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}


Test_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
Test_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
Test_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result


Test_Scores_for_plotting1$Score <- Test_Scores_for_plotting1$pos_Score

# Calculate Score
Test_Scores_for_plotting2$Score <- Test_Scores_for_plotting2$pos_Score - Test_Scores_for_plotting2$neg_Score
#Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))

Test_Scores_for_plotting3$Score <- Test_Scores_for_plotting3$pos_Score - Test_Scores_for_plotting3$neg_Score
#Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))

Test_Scores_for_plotting1
Test_Scores_for_plotting2
Test_Scores_for_plotting3



Test_Scores_for_plotting1$Sample <- rownames(Test_Scores_for_plotting1)
Test_Scores_for_plotting2$Sample <- rownames(Test_Scores_for_plotting2)
Test_Scores_for_plotting3$Sample <- rownames(Test_Scores_for_plotting3)


Test_Scores_for_plotting1 <- left_join(Test_Scores_for_plotting1, test_labels)
Test_Scores_for_plotting2 <- left_join(Test_Scores_for_plotting2, test_labels)
Test_Scores_for_plotting3 <- left_join(Test_Scores_for_plotting3, test_labels)


Test_Scores_for_plotting1$Study <- gsub("1_", "", Test_Scores_for_plotting1$Study)
Test_Scores_for_plotting1$Study <- gsub("2_", "", Test_Scores_for_plotting1$Study)


Test_Scores_for_plotting2$Study <- gsub("1_", "", Test_Scores_for_plotting2$Study)
Test_Scores_for_plotting2$Study <- gsub("2_", "", Test_Scores_for_plotting2$Study)


Test_Scores_for_plotting3$Study <- gsub("1_", "", Test_Scores_for_plotting3$Study)
Test_Scores_for_plotting3$Study <- gsub("2_", "", Test_Scores_for_plotting3$Study)


Test_Scores_for_plotting1$Study <- factor(Test_Scores_for_plotting1$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Test_Scores_for_plotting2$Study <- factor(Test_Scores_for_plotting2$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Test_Scores_for_plotting3$Study <- factor(Test_Scores_for_plotting3$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))



Test_Scores_for_plotting1$Condition <- factor(Test_Scores_for_plotting1$Condition, levels = c("Control", "Infected"))
Test_Scores_for_plotting2$Condition <- factor(Test_Scores_for_plotting2$Condition, levels = c("Control", "Infected"))
Test_Scores_for_plotting3$Condition <- factor(Test_Scores_for_plotting3$Condition, levels = c("Control", "Infected"))

Test_Scores_for_plotting1$Set = "Set1"
Test_Scores_for_plotting2$Set = "Set2"
Test_Scores_for_plotting3$Set = "Set3"


# This is very important to scale based on training set
Test_Scores_for_plotting1$Scaled_Score <- as.numeric((Test_Scores_for_plotting1$Score - mean(Scores_for_plotting1$Unscaled_Score)) / sd(Scores_for_plotting1$Unscaled_Score))
Test_Scores_for_plotting2$Scaled_Score <- as.numeric((Test_Scores_for_plotting2$Score - mean(Scores_for_plotting2$Unscaled_Score)) / sd(Scores_for_plotting2$Unscaled_Score))
Test_Scores_for_plotting3$Scaled_Score <- as.numeric((Test_Scores_for_plotting3$Score - mean(Scores_for_plotting3$Unscaled_Score)) / sd(Scores_for_plotting3$Unscaled_Score))



Test_Scores_for_plotting2 <- Test_Scores_for_plotting2[,colnames(Test_Scores_for_plotting1)]
Test_Scores_for_plotting3 <- Test_Scores_for_plotting3[,colnames(Test_Scores_for_plotting1)]

Test_Scores_for_plotting1$Score <- Test_Scores_for_plotting1$Scaled_Score
Test_Scores_for_plotting2$Score <- Test_Scores_for_plotting2$Scaled_Score
Test_Scores_for_plotting3$Score <- Test_Scores_for_plotting3$Scaled_Score


Test_GFS_thresholds_Set1 = roc(Condition ~ Score, data = Test_Scores_for_plotting1) 
Test_GFS_thresholds_Set2 = roc(Condition ~ Score, data = Test_Scores_for_plotting2)
Test_GFS_thresholds_Set3 = roc(Condition ~ Score, data = Test_Scores_for_plotting3) 

Test_GFS_thresholds_Set2

Test_GFS_thresholds_Set1 <- cbind(Test_GFS_thresholds_Set1$thresholds, Test_GFS_thresholds_Set1$sensitivities, Test_GFS_thresholds_Set1$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set1")
Test_GFS_thresholds_Set2 <- cbind(Test_GFS_thresholds_Set2$thresholds, Test_GFS_thresholds_Set2$sensitivities, Test_GFS_thresholds_Set2$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set2")
Test_GFS_thresholds_Set3 <- cbind(Test_GFS_thresholds_Set3$thresholds, Test_GFS_thresholds_Set3$sensitivities, Test_GFS_thresholds_Set3$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set3")

Test_GFS_thresholds_Set1$Combined = rowSums(Test_GFS_thresholds_Set1[,c("Sensitivity", "Specificity")])
Test_GFS_thresholds_Set2$Combined = rowSums(Test_GFS_thresholds_Set2[,c("Sensitivity", "Specificity")])
Test_GFS_thresholds_Set3$Combined = rowSums(Test_GFS_thresholds_Set3[,c("Sensitivity", "Specificity")])


Test_GFS_combined = rbind(Test_GFS_thresholds_Set1,Test_GFS_thresholds_Set2,Test_GFS_thresholds_Set3)


write.table(Test_GFS_combined, args[9], sep = "\t", quote = F, row.names = F)

set1_thresh_vec = Test_GFS_thresholds_Set1 %>% filter(Sensitivity > 0.85) %>% filter(Combined == max(Combined)) %>% slice_min(order_by = Sensitivity) #  0.8837282   0.9117647   0.5294118 Set1 1.441176
set2_thresh_vec = Test_GFS_thresholds_Set2 %>% filter(Sensitivity > 0.85) %>% filter(Combined == max(Combined)) # -0.2379047   0.9411765   0.7647059 Set2 1.705882
set3_thresh_vec = Test_GFS_thresholds_Set3 %>% filter(Sensitivity > 0.85) %>% filter(Combined == max(Combined))#  -0.3740885   0.9411765   0.7647059 Set3 1.705882

set2_thresh_vec
# Dealing with scaled score now rememeber
Test_Scores_for_plotting1 = Test_Scores_for_plotting1 %>% mutate(Prediction = if_else(Score > set1_thresh_vec$Threshold, "Infected", "Control"))
Test_Scores_for_plotting2 = Test_Scores_for_plotting2 %>% mutate(Prediction = if_else(Score > set2_thresh_vec$Threshold, "Infected", "Control"))
Test_Scores_for_plotting3 = Test_Scores_for_plotting3 %>% mutate(Prediction = if_else(Score > set3_thresh_vec$Threshold, "Infected", "Control"))


thresholds_df <- as.matrix(cbind(Test_Scores_for_plotting1$Prediction, Test_Scores_for_plotting2$Prediction, Test_Scores_for_plotting3$Prediction))
thresholds_df
thresholds_df <- ifelse(thresholds_df == "Infected", 1, 0)
mode(thresholds_df) <- "numeric"
rownames(thresholds_df) <- Test_Scores_for_plotting1$Sample

colnames(thresholds_df) <- c("Set1", "Set2", "Set3")

test_set = test_labels

test_set <- test_set %>% mutate(Infection_administration = if_else(Study == "Mcloughlin_pbl", "Natural", Infection_administration))
test_set$Study <- if_else(test_set$Study == "1_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "2_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Wiarda", "Wiarda et al., (2020)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin", "McLoughlin et al., (2021)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin_pbl", "McLoughlin et al., (2014)", test_set$Study)

test_set$Study <- factor(test_set$Study, levels = c("O'Grady et al., (2025)", "McLoughlin et al., (2014)", "Wiarda et al., (2020)", "McLoughlin et al., (2021)"))
test_set <- test_set %>% arrange(Study)
test_set
test_set
rownames(test_set) <- test_set$Sample
thresholds_df <- thresholds_df[rownames(test_set),]

ann_colors <- list(
  Study = c("O'Grady et al., (2025)" = "#351338",
            "Wiarda et al., (2020)" = "#dad2ff",
            "McLoughlin et al., (2021)" = "#734500",
            "McLoughlin et al., (2014)" = "#69a3a5"),
  Condition = c("Control" = "#2166ac",
                "Infected" = "#b2182b"),
  Infection_administration = c(Natural = "#708238",
                               Experimental = "beige"))

annotation_columns = test_set %>% select(c(Condition, Study, Infection_administration))
annotation_columns
rownames(annotation_columns) <- test_set$Sample 
annotation_columns <- annotation_columns %>% arrange(desc(Infection_administration))
thresholds_df <- thresholds_df[rownames(annotation_columns),]
write.table(thresholds_df, file = args[10], sep = "\t")

library(pheatmap)
pdf(args[11], width = 15, height = 20)
pheatmap(t(thresholds_df),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_row = c(seq(2:3)),
         annotation_col = annotation_columns,
         annotation_colors = ann_colors,
         border_color = "black",
         show_colnames = F,
         color = c("lightblue", "orange"))
dev.off()
