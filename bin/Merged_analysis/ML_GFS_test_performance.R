
# test

test_counts <- read.table("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/test_normalised_filtered_counts.txt")
tb_test_vst_glm_results <- predict(GLM_model, newdata = test_counts, type = "prob") %>% mutate(Model = "GLM")
tb_test_vst_glmridge_results <- predict(GLMRIDGE_model, newdata = test_counts, type = "prob") %>% mutate(Model = "RIDGE")
tb_test_vst_glmlasso_results <- predict(GLMLASSO_model, newdata = test_counts, type = "prob") %>% mutate(Model = "LASSO")
tb_test_vst_glmenet_results <- predict(GLMENET_model, newdata = test_counts, type = "prob") %>% mutate(Model = "ENET")
tb_test_vst_RF_results <- predict(RF_model, newdata = test_counts, type = "prob") %>% mutate(Model = "RF")
tb_test_vst_RF_ET_results <- predict(RF_ET_model, newdata = test_counts, type = "prob") %>% mutate(Model = "RF_ET")
tb_test_vst_NB_results <- predict(NB_model, newdata = test_counts, type = "prob") %>% mutate(Model = "NB")
tb_test_vst_MLP_results <- predict(MLP_model, newdata = test_counts, type = "prob") %>% mutate(Model = "MLP")


test_counts

t(test_counts)
head(test_counts)

tb_test_vst_glmenet_results
if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = t(test_counts)[rownames(t(test_counts)) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =t(test_counts)[rownames(t(test_counts)) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}



# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = t(test_counts)[rownames(t(test_counts)) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =t(test_counts)[rownames(t(test_counts)) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = t(test_counts)[rownames(t(test_counts)) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =t(test_counts)[rownames(t(test_counts)) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =t(test_counts)[rownames(t(test_counts)) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}


test_counts_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
test_counts_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
test_counts_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result




test_counts_Scores_for_plotting1$Set = "Pass 1"
test_counts_Scores_for_plotting2$Set = "Pass 2"
test_counts_Scores_for_plotting3$Set = "Combined"

test_counts_Scores_for_plotting1$Score <- test_counts_Scores_for_plotting1$pos_Score
test_counts_Scores_for_plotting2$Score <- test_counts_Scores_for_plotting2$pos_Score - test_counts_Scores_for_plotting2$neg_Score
test_counts_Scores_for_plotting3$Score <- test_counts_Scores_for_plotting3$pos_Score - test_counts_Scores_for_plotting3$neg_Score



test_counts_Scores_for_plotting1$Score <- as.numeric(scale(test_counts_Scores_for_plotting1$Score))
test_counts_Scores_for_plotting2$Score <- as.numeric(scale(test_counts_Scores_for_plotting2$Score))
test_counts_Scores_for_plotting3$Score <- as.numeric(scale(test_counts_Scores_for_plotting3$Score))

colnames(test_counts_Scores_for_plotting1)[1] <- "Control" # Dummy
colnames(test_counts_Scores_for_plotting1)[3] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(test_counts_Scores_for_plotting1)[2] <- "Model"

colnames(test_counts_Scores_for_plotting2)[1] <- "Control" # Dummy
colnames(test_counts_Scores_for_plotting2)[4] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(test_counts_Scores_for_plotting2)[3] <- "Model"

colnames(test_counts_Scores_for_plotting3)[1] <- "Control" # Dummy
colnames(test_counts_Scores_for_plotting3)[4] <- "Infected" # This is column that is used for differentiating between ML models so need to rename for input to ROC function.
colnames(test_counts_Scores_for_plotting3)[3] <- "Model"

test_counts_Scores_for_plotting1 <- test_counts_Scores_for_plotting1 %>% select(Control, Infected, Model)
test_counts_Scores_for_plotting2 <- test_counts_Scores_for_plotting2 %>% select(Control, Infected, Model)
test_counts_Scores_for_plotting3 <- test_counts_Scores_for_plotting3 %>% select(Control, Infected, Model)



df_test_vst_results = data.frame(rbind(tb_test_vst_glm_results,
                                        tb_test_vst_glmridge_results,
                                        tb_test_vst_glmlasso_results,
                                        tb_test_vst_glmenet_results,
                                        tb_test_vst_RF_results,
                                        tb_test_vst_RF_ET_results,
                                        tb_test_vst_MLP_results,
                                        tb_test_vst_NB_results,
                                        test_counts_Scores_for_plotting1,
                                        test_counts_Scores_for_plotting2,
                                        test_counts_Scores_for_plotting3))

df_test_vst_results$Model = factor(df_test_vst_results$Model)
test
ROC_test_combined(df_test_vst_results, Test_data_score$Condition)
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/TB_test_all_ML_sets.pdf", width = 12, height = 12, dpi = 600)
