
# This script is for manual reproduction of the threshold heatmap as the pheatmap package is difficult to output the pdf and you
# need to play around with the dev.off() setting etc.

library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(edgeR)
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
library(SNPRelate)
library(dplyr)
library(matrixStats)
library(rstatix)


args = commandArgs(trailingOnly = T)
# Remake the heatmap again
GLM = readRDS(args[1])
ENET = readRDS(args[2])
LASSO = readRDS(args[3])
RIDGE = readRDS(args[4])
MLP = readRDS(args[5])
NB = readRDS(args[6])
RF = readRDS(args[7])
RF_ET = readRDS(args[8])

test_data = read.table(args[9]) 

names(test_data)[names(test_data) == 'bta.mir.223'] <- 'bta-mir-223'
names(test_data)[names(test_data) == 'WC1.12'] <- 'WC1-12'
test_data[,"bta-mir-223"]


resamples_rf_et <- resamples(RF[1:6])
best_rf_et_model_name <- names(sort(summary(resamples_rf_et)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
RF <- RF[[best_rf_et_model_name]]
RF


resamples_rf_et <- resamples(RF_ET[1:6])
best_rf_et_model_name <- names(sort(summary(resamples_rf_et)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
RF_ET <- RF_ET[[best_rf_et_model_name]]
RF_ET

GLM_predict = predict(object = GLM, test_data, type = "prob") %>% mutate(Model = "GLM")
RIDGE_predict = predict(object = RIDGE, test_data, type = "prob") %>% mutate(Model = "RIDGE")
LASSO_predict = predict(object = LASSO, test_data, type = "prob") %>% mutate(Model = "LASSO")
ENET_predict = predict(object = ENET, test_data, type = "prob") %>% mutate(Model = "ENET")
MLP_predict = predict(object = MLP, test_data, type = "prob") %>% mutate(Model = "MLP")
NB_predict = predict(object = NB, test_data, type = "prob") %>% mutate(Model = "NB")
RF_predict = predict(object = RF, test_data, type = "prob") %>% mutate(Model = "RF")
RF_ET_predict = predict(object = RF_ET, test_data, type = "prob") %>% mutate(Model = "RF-ET")

test_set = read.table(args[10], header = T)

test_set$Sample == rownames(test_data)

GLM_predict <- cbind(GLM_predict, test_set) 
GLMRIDGE_predict <- cbind(RIDGE_predict, test_set) 
GLMLASSO_predict<- cbind(LASSO_predict, test_set) 
GLMENET_predict<- cbind(ENET_predict, test_set) 
RF_predict <- cbind(RF_predict, test_set) 
RF_ET_predict <- cbind(RF_ET_predict, test_set) 
nnet_predict<- cbind(MLP_predict, test_set) 
NB_predict<- cbind(NB_predict, test_set) 


GLM_ROC <- pROC::roc(GLM_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMRIDGE_ROC <- pROC::roc(GLMRIDGE_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMLASSO_ROC <- pROC::roc(GLMLASSO_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMENET_ROC <- pROC::roc(GLMENET_predict, response = Condition, predictor = Infected, quiet = TRUE)
RF_ROC <- pROC::roc(RF_predict, response = Condition, predictor = Infected, quiet = TRUE)
RF_ET_ROC <- pROC::roc(RF_ET_predict, response = Condition, predictor = Infected, quiet = TRUE)
nnet_ROC <- pROC::roc(nnet_predict, response = Condition, predictor = Infected, quiet = TRUE)
NB_ROC <- pROC::roc(NB_predict, response = Condition, predictor = Infected, quiet = TRUE)


GLM_thresholds = cbind(GLM_ROC$thresholds, GLM_ROC$specificities, GLM_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLM_thresholds$Combined = rowSums(GLM_thresholds[,c("Spec", "Sens")])

GLMRIDGE_thresholds = cbind(GLMRIDGE_ROC$thresholds, GLMRIDGE_ROC$specificities, GLMRIDGE_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMRIDGE_thresholds$Combined = rowSums(GLMRIDGE_thresholds[,c("Spec", "Sens")])

GLMLASSO_thresholds = cbind(GLMLASSO_ROC$thresholds, GLMLASSO_ROC$specificities, GLMLASSO_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMLASSO_thresholds$Combined = rowSums(GLMLASSO_thresholds[,c("Spec", "Sens")])

GLMENET_thresholds = cbind(GLMENET_ROC$thresholds, GLMENET_ROC$specificities, GLMENET_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMENET_thresholds$Combined = rowSums(GLMENET_thresholds[,c("Spec", "Sens")])


RF_thresholds = cbind(RF_ROC$thresholds, RF_ROC$specificities, RF_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
RF_thresholds$Combined = rowSums(RF_thresholds[,c("Spec", "Sens")])


RF_ET_thresholds = cbind(RF_ET_ROC$thresholds, RF_ET_ROC$specificities, RF_ET_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
RF_ET_thresholds$Combined = rowSums(RF_ET_thresholds[,c("Spec", "Sens")])


nnet_thresholds = cbind(nnet_ROC$thresholds, nnet_ROC$specificities, nnet_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
nnet_thresholds$Combined = rowSums(nnet_thresholds[,c("Spec", "Sens")])


NB_thresholds = cbind(NB_ROC$thresholds, NB_ROC$specificities, NB_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
NB_thresholds$Combined = rowSums(NB_thresholds[,c("Spec", "Sens")])

GLM_thresholds
NB_thresholds
GLM_final_threshold = GLM_thresholds %>% filter(Sens > 0.7) %>% filter(Spec > 0.01) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMRIDGE_final_threshold = GLMRIDGE_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMLASSO_final_threshold = GLMLASSO_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMENET_final_threshold = GLMENET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
RF_final_threshold = RF_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
RF_ET_final_threshold = RF_ET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
nnet_final_thresholds = nnet_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
NB_final_thresholds = NB_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()





GLM_thresholds_predict = predict(GLM, test_data, type="prob")[,"Infected"]
GLM_thresholds_predict <- ifelse(GLM_thresholds_predict > GLM_final_threshold, "Infected", "Control")

GLMENET_thresholds_predict = predict(ENET, test_data, type="prob")[,"Infected"]
GLMENET_thresholds_predict <- ifelse(GLMENET_thresholds_predict > GLMENET_final_threshold, "Infected", "Control")

GLMRIDGE_thresholds_predict = predict(RIDGE, test_data, type="prob")[,"Infected"]
GLMRIDGE_thresholds_predict <- ifelse(GLMRIDGE_thresholds_predict > GLMRIDGE_final_threshold, "Infected", "Control")

GLMLASSO_thresholds_predict = predict(LASSO, test_data, type="prob")[,"Infected"]
GLMLASSO_thresholds_predict <- ifelse(GLMLASSO_thresholds_predict > GLMLASSO_final_threshold, "Infected", "Control")

RF_thresholds_predict = predict(RF, test_data, type="prob")[,"Infected"]
RF_thresholds_predict <- ifelse(RF_thresholds_predict > RF_final_threshold, "Infected", "Control")

RF_ET_thresholds_predict = predict(RF_ET, test_data, type="prob")[,"Infected"]
RF_ET_thresholds_predict <- ifelse(RF_ET_thresholds_predict > RF_ET_final_threshold, "Infected", "Control")

nnet_thresholds_predict = predict(MLP, test_data, type="prob")[,"Infected"]
nnet_thresholds_predict <- ifelse(nnet_thresholds_predict > nnet_final_thresholds, "Infected", "Control")

NB_thresholds_predict = predict(NB, test_data, type="prob")[,"Infected"]
NB_thresholds_predict <- ifelse(NB_thresholds_predict > NB_final_thresholds, "Infected", "Control")





thresholds_df <- as.matrix(cbind(GLMENET_thresholds_predict, GLMLASSO_thresholds_predict, nnet_thresholds_predict, RF_ET_thresholds_predict, RF_thresholds_predict, GLMRIDGE_thresholds_predict, NB_thresholds_predict, GLM_thresholds_predict))
thresholds_df <- ifelse(thresholds_df == "Infected", 1, 0)
mode(thresholds_df) <- "numeric"
rownames(thresholds_df) <- test_set$Sample

test_set <- test_set %>% mutate(Infection_administration = if_else(Study == "Mcloughlin_pbl", "Natural", Infection_administration))
test_set$Study <- if_else(test_set$Study == "1_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "2_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Wiarda", "Wiarda et al., (2020)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin", "McLoughlin et al., (2021)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin_pbl", "McLoughlin et al., (2014)", test_set$Study)

test_set$Study <- factor(test_set$Study, levels = c("O'Grady et al., (2025)", "McLoughlin et al., (2014)", "Wiarda et al., (2020)", "McLoughlin et al., (2021)"))
test_set <- test_set %>% arrange(Study)
test_set

#thresholds_df <- thresholds_df[rownames(test_set),]
thresholds_df

thresholds_df
colnames(thresholds_df) <- c("ENET", "LASSO", "MLP", "RF(ET)", "RF", "RIDGE", "NB", "GLM")

ann_colors <- list(
  Study = c("O'Grady et al., (2025)" = "#351338",
            "Wiarda et al., (2020)" = "#dad2ff",
            "McLoughlin et al., (2021)" = "#734500",
            "McLoughlin et al., (2014)" = "#69a3a5"),
  Condition = c("Control" = "#2166ac",
                "Infected" = "#b2182b"),
  Infection_administration = c("Natural" = "#708238",
                               "Experimental" = "beige"))

thresholds_df

annotation_columns = test_set %>% select(c(Condition, Study, Infection_administration))
rownames(annotation_columns) <- test_set$Sample 
annotation_columns <- annotation_columns %>% arrange(desc(Infection_administration))
annotation_columns
thresholds_df <- thresholds_df[rownames(annotation_columns),]
write.table(thresholds_df, file = args[11], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
library(pheatmap)
pdf(args[12], width = 15, height = 20)
pheatmap(t(thresholds_df),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_row = c(seq(2:8)),
         annotation_col = annotation_columns,
         annotation_colors = ann_colors,
         border_color = "black",
         show_colnames = F,
         color = c("lightblue", "orange"))
dev.off()

threshold_summary = matrix(ncol = 5, nrow = 0)
colnames(threshold_summary) <- c("Threshold", "Spec", "Sens", "Combined", "Model")

threshold_summary <- rbind(threshold_summary, GLM_thresholds %>% filter(Sens > 0.7) %>% filter(Spec > 0.01) %>% filter(Combined == max(Combined)) %>% mutate(Model = "GLM"))
threshold_summary <- rbind(threshold_summary,GLMRIDGE_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "RIDGE"))
threshold_summary <- rbind(threshold_summary,GLMLASSO_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "LASSO"))
threshold_summary <- rbind(threshold_summary,GLMENET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "ENET"))
threshold_summary <- rbind(threshold_summary,RF_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "RF"))
threshold_summary <- rbind(threshold_summary,RF_ET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "RF-ET"))
threshold_summary <- rbind(threshold_summary,nnet_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens)  %>% mutate(Model = "MLP"))
threshold_summary <- rbind(threshold_summary,NB_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% mutate(Model = "NB"))
write.table(threshold_summary, args[13], sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

