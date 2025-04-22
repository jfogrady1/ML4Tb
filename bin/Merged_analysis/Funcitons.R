####################################################
####################################################
#####################################################
####
#     Helper functions
####
####################################################
####################################################
####################################################

# Load necessary library
library(matrixStats)  # For row geometric mean calculation

# Function to calculate geometric mean
geomMean <- function (x, na.rm = FALSE) 
{
  # Remove NA values if requested
  if (na.rm) x <- x[!is.na(x)]
  
  # Filter out zero and negative values
  x <- x[x > 0] 
  # If all values were zero/negative, return 0 to avoid NaN
  if (length(x) == 0) return(0) 
  
  # Compute the geometric mean
  return(exp(mean(log(x))))
}

getGenesData <- function(genes.mtx, genes){
  tmp <- genes.mtx[genes, , drop = FALSE]
  tmp[is.na(tmp)] <- 1
  return(tmp)
}

ScoreGenesMtx <- function(GeneMtx, pos.genes, neg.genes){
  posScore <- 0 
  if (sum(pos.genes %in% rownames(GeneMtx)) > 0){
    posMatch <- getGenesData(GeneMtx, pos.genes)
    posScore <- apply(posMatch, 2, geomMean)
  }
  negScore <- 0
  if (sum(neg.genes %in% rownames(GeneMtx)) > 0){
    negMatch <- getGenesData(GeneMtx, neg.genes)
    negScore <- apply(negMatch, 2, geomMean)
  }
  
  ## Weight pos and neg by how many genes are being called
  if(length(posScore)>1){
    ratio <- length(neg.genes)/length(pos.genes); 
    negScore <- ratio*negScore
    totalScore <- scale(posScore - negScore)
  } else {
    totalScore <- scale(-negScore)
  }
  
  return(totalScore)
}


# CV ROC plot
CV_ROC_plot <- function (MODEL) {
  MODEL$resample <- MODEL$resample[order(MODEL$resample$Resample),]
  roc_values <- round(as.numeric(MODEL$resample$ROC), 3)
  roc_text <- paste(paste0("ROC Resample ", seq_along(roc_values), " = ", roc_values), collapse = "\n")
  resample_colors <- ggsci::pal_npg("nrc")(length(roc_values))
  
  MODEL_train_ROC <- MODEL$pred %>% ggplot(aes(m=Infected, d=factor(obs, levels = c("Control", "Infected")), colour = Resample)) + # Note do not need to filter unless saving all predictions
    style_roc() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
    geom_roc(n.cuts=0) + 
    coord_equal() +
    scale_colour_npg() +
    scale_x_continuous("1 - Specificity (FPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous("Sensitivity (TPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
          axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
          legend.text = element_text(size = 15, colour = "black"),
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.title.y = element_text(size = 18, color = "black"),
          legend.position="none") +
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (6-1) * 0.025, # Adjust y position for each line
             label = paste0("ROC Average CV = ", round(mean(MODEL$resample$ROC),3)), 
             hjust = 0,
             col = "black",
             size = 5,)
  # Add annotations for each resample
  for (i in seq_along(roc_values)) {
    MODEL_train_ROC <- MODEL_train_ROC + 
      annotate("text", 
               x = 0.75, 
               y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
               label = paste0("ROC Fold ", i, " = ", roc_values[i]), 
               hjust = 0, 
               size = 5, 
               color = resample_colors[i])
  }
  
  return(MODEL_train_ROC)
}


ROC_predict_individual <- function(Predictions, labels) {
  
  Predictions <- cbind(Predictions, labels)
  colnames(Predictions)[length(colnames(Predictions))] <- "obs"
  Predictions_plot <- Predictions %>% ggplot(aes(m=Infected, d=factor(obs, levels = c("Control", "Infected")))) + # Note do not need to filter unless saving all predictions
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
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.title.y = element_text(size = 18, color = "black"),
          legend.position="none")
  Predictions_plot <- Predictions_plot + annotate("text", x = .75, y = .25, 
                                                  label = paste("AUC =", round(calc_auc(Predictions_plot)$AUC, 3)))
}


# Calculate ROC values for all models
ROC_test_combined <- function(combined_predictions, labels) {
  combined_predictions$Model <- factor(combined_predictions$Model)
  combined_predictions <- cbind(combined_predictions, rep(labels, length(unique(combined_predictions$Model))))
  colnames(combined_predictions)[length(colnames(combined_predictions))] <- "Condition"
  roc_values <- combined_predictions %>%
    group_by(Model) %>%
    summarise(
      roc_values = as.numeric(pROC::roc(
        response = ifelse(Condition == "Infected", 1, 0),
        predictor = Infected,
        quiet = TRUE
      )$auc)
    )
  roc_values <- roc_values %>% arrange(desc(roc_values))
  
  roc_text <- paste(paste0("ROC Resample ", seq_along(roc_values$Model), " = ", roc_values$roc_values), collapse = "\n")
  resample_colors <- ggsci::pal_npg("nrc")(length(roc_values$roc_values))
  
  
  
  combined_predictions_ROC <- combined_predictions %>% ggplot(aes(m=Infected, d=factor(Condition, levels = c("Control", "Infected")), colour = factor(Model, levels = roc_values$Model))) + # Note do not need to filter unless saving all predictions
    style_roc() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
    geom_roc(n.cuts=0) + 
    coord_equal() +
    scale_colour_npg() +
    scale_x_continuous("1 - Specificity (FPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    scale_y_continuous("Sensitivity (TPR)", breaks = seq(0, 1, by = 0.1), limits = c(0,1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
          axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
          legend.text = element_text(size = 15, colour = "black"),
          axis.title.x = element_text(size = 18, colour = "black"),
          axis.title.y = element_text(size = 18, color = "black"),
          legend.position="none")
  
  # Add annotations for each resample
  for (i in seq_along(roc_values$Model)) {
    combined_predictions_ROC <- combined_predictions_ROC + 
      annotate("text", 
               x = 0.75, 
               y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
               label = paste0("ROC ", roc_values$Model[i], " = ", round(roc_values$roc_values[i],3)), 
               hjust = 0, 
               size = 5, 
               color = resample_colors[i])
  }
  
  
  return(combined_predictions_ROC)
}



pred_fun <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, 2]  # Return probabilities for the positive class
}



greedy_forward_search = function(gene.list, cpm_matrix, de_results, metadata) {
  train_set$Study
  metadata$study = metadata$Study
  metadata$study = if_else(metadata$study == "1_OGrady", "OGrady", metadata$study)
  metadata$study = if_else(metadata$study == "2_OGrady", "OGrady", metadata$study)
  metadata$study = factor(metadata$study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))

  
  Ogrady_n = as.numeric(sum(metadata$study == "OGrady"))
  Wiarda_n = as.numeric(sum(metadata$study == "Wiarda"))
  Mcloughlin_n = as.numeric(sum(metadata$study == "Mcloughlin"))
  Mcloughlin_pbl_n = as.numeric(sum(metadata$study == "Mcloughlin_pbl"))

  de_results$variable = rownames(de_results)
  starts = gene.list
  starts_df = data.frame(gene=character(), OGrady=numeric(), Wiarda=numeric(), Mcloughlin=numeric(), Mcloughlin_pbl=numeric(), Average=numeric(), W_Average=numeric(), stringsAsFactors=FALSE)
  for (s in starts){
    df = merge(train_counts_cpms[rownames(train_counts_cpms) == s,]  %>%
                 as.data.frame() %>%
                 mutate(Score = log2(. + 1)),  # Rename column properly
               metadata, 
               by.x = "row.names", by.y = "Sample")
    
    auc = df %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
    colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
    auc$Average <- rowMeans(auc)
    auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
    colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Average", "W_Average")
    print(auc)
    starts_df = rbind(starts_df, cbind("gene"=s, "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
    colnames(starts_df)[7] <- "W_Average"
  }
    
  starts_df = starts_df %>% as.data.frame() %>% dplyr::arrange(desc(W_Average)) 
  colnames(starts_df)[7] <- "W_Average"
  head(starts_df)
  plot_df = matrix(ncol=7, NA)
  colnames(plot_df)=c("combination","OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Average", "W_Average")

  plot_df = rbind(plot_df, cbind("combination"=starts_df[1,1], "OGrady"=starts_df[1,2], "Wiarda"=starts_df[1,3], "Mcloughlin"=starts_df[1,4], "Mcloughlin_pbl"=starts_df[1,5], "Average" = starts_df[1,6], "W_Average" = starts_df[1,7]))
  
  best_auc = plot_df %>% as.data.frame() %>% filter(row_number()==n())
  best_auc = best_auc$W_Average
  best_auc
  df2 = matrix(ncol=7,NA)
  colnames(df2)=c("combination","OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Average", "W_Average")
  start = plot_df[2,1]
  print(start)
  print(plot_df)
  print(starts_df)
  rest = starts[!(starts %in% start)]
  for (r in 1:length(rest)){
    print(r)
    comb = c(start,rest[r])
    pos = de_results[de_results$Direction == "Positive",] %>% filter(variable %in% comb) %>% pull(variable)
    neg =de_results[de_results$Direction == "Negative",] %>% filter(variable %in% comb)%>% pull(variable)
    
    if (length(pos)==1){
      pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]  %>% as.data.frame() %>%
        mutate(across(everything(), ~ log2(. + 1))) %>% t() %>% t() %>% set_colnames("pos_Score")
    } else if (length(pos)>1) {
      pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]   %>% 
        as.data.frame() %>% mutate_all(as.numeric) %>% 
        mutate(across(everything(), ~ log2(. + 1))) %>% 
        summarise_all(mean) %>% t() %>% set_colnames("pos_Score")
    } else{
      pos_counts=NA
    }
    
    if (length(neg)==1){
      neg_counts = cpm_matrix[rownames(cpm_matrix) %in% neg,]  %>% as.data.frame() %>% mutate_all(as.numeric) %>% 
        mutate(across(everything(), ~ log2(. + 1))) %>%
        t() %>% t() %>% set_colnames("neg_Score")
    } else if (length(neg)>1) {
      neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,] %>%
        as.data.frame() %>% mutate_all(as.numeric) %>%
        mutate(across(everything(), ~ log2(. + 1))) %>%
        summarise_all(mean) %>% t() %>%  set_colnames("neg_Score")
    } else {
      neg_counts = NA
    }

    
    if (all(!is.na(neg_counts)) & all(!is.na(pos_counts))){
      df_test = merge(pos_counts,neg_counts, by="row.names")
      df_test = merge(df_test, metadata, by.x="Row.names", by.y="Sample")
      df_test$Score = df_test$pos_Score - df_test$neg_Score
      
      
      
      auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
      colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
      auc$Average <- rowMeans(auc)
      auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
      print("HERE")
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
    } 
      else if (all(is.na(neg_counts)) & all(!is.na(pos_counts))){
      df_test = merge(pos_counts, metadata, by.x="row.names", by.y="Sample")
      df_test$Score = df_test$pos_Score
      auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
      colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
      auc$Average <- rowMeans(auc)
      auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
      print("HERE")
      print(auc)
      print(df2)
      df_temp = cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6])
      print(colnames(df_temp))
      print(colnames(df2))
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
    } 
      else {
      df_test = merge(neg_counts, metadata, by.x="row.names", by.y="Sample")
      df_test$Score = df_test$neg_Score
      auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
      colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
      auc$Average <- rowMeans(auc)
      auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
    }
  }
  
  while (!is.na(best_auc)){
    sorted_df = df2 %>% as.data.frame()%>% filter(as.numeric(W_Average) > best_auc) %>% arrange(desc(W_Average))
    new_start = unlist(str_split(sorted_df[1,1], "_"))
    new_rest = starts[!(starts %in% new_start)]
    
    df2 = matrix(ncol=7,NA)
    colnames(df2)=c("combination","OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Average", "W_Average")
    best_auc = as.numeric(sorted_df[1,6])
    plot_df = rbind(plot_df, cbind("combination"=sorted_df[1,1], "OGrady"=sorted_df[1,2], "Wiarda"=sorted_df[1,3], "Mcloughlin"=sorted_df[1,4], "Mcloughlin_pbl"=sorted_df[1,5], "Average" = sorted_df[1,6], "W_Average" = sorted_df[1,7]))
    print(best_auc)
    for (r in 1:length(new_rest)){
      
      comb = c(new_start,new_rest[r])
      print(r)
      pos = de_results[de_results$Direction == "Positive",] %>% filter(variable %in% comb) %>% pull(variable)
      neg =de_results[de_results$Direction == "Negative",] %>% filter(variable %in% comb)%>% pull(variable)
      if (length(pos)==1){
        pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]   %>% as.data.frame() %>% mutate_all(as.numeric) %>% 
          mutate(across(everything(), ~ log2(. + 1))) %>% t() %>% t() %>% set_colnames("pos_Score")
      } else if (length(pos)>1) {
        pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]   %>% 
          as.data.frame() %>% mutate_all(as.numeric) %>% 
          mutate(across(everything(), ~ log2(. + 1))) %>% 
          summarise_all(mean) %>% t()%>% set_colnames("pos_Score")    
      } else{
        pos_counts=NA
      }
      
      if (length(neg)==1){
        neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,]   %>% as.data.frame() %>% mutate_all(as.numeric) %>% 
          mutate(across(everything(), ~ log2(. + 1)))%>% 
          t() %>% t() %>% set_colnames("neg_Score")
      } else if (length(neg)>1) {
        neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,]   %>% 
          as.data.frame() %>% mutate_all(as.numeric) %>% 
          mutate(as.data.frame( . + 1))%>% 
          summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
      } else {
        neg_counts = NA
      }
      
      if (all(!is.na(neg_counts)) & all(!is.na(pos_counts))){
        df_test = merge(pos_counts,neg_counts, by="row.names")
        df_test = merge(df_test, metadata, by.x="Row.names", by.y="Sample")
        
        df_test$Score = df_test$pos_Score - df_test$neg_Score
        auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
        colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
        auc$Average <- rowMeans(auc)
        auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
      } 
      else if (all(is.na(neg_counts)) & all(!is.na(pos_counts))){
        df_test = merge(pos_counts, metadata, by.x="row.names", by.y="Sample")
        df_test$Score = df_test$pos_Score
        auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
        colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
        auc$Average <- rowMeans(auc)
        auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
      } 
      else {
        df_test = merge(neg_counts, metadata, by.x="row.names", by.y="Sample")
        df_test$Score = df_test$neg_Score
        auc = df_test %>% group_by(study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))), .groups = "drop") %>% pivot_wider(names_from = "study", values_from = "AUC")
        colnames(auc) <- c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl")
        auc$Average <- rowMeans(auc)
        auc$W_Average = as.numeric(((auc[1,1] * Ogrady_n) + (auc[1,2] * Wiarda_n) + (auc[1,3] * Mcloughlin_n) + (auc[1,4] * Mcloughlin_pbl_n)) / sum(Ogrady_n, Wiarda_n, Mcloughlin_n, Mcloughlin_pbl_n))
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "OGrady"=auc[1,1], "Wiarda"=auc[1,2], "Mcloughlin"=auc[1,3], "Mcloughlin_pbl"=auc[1,4], "Average" = auc[1,5], "W_Average" = auc[1,6]))
      }
    }
    plot.roc(df_test$Condition, df_test$Score, print.auc=TRUE)
  }
  return(as.data.frame(plot_df))
}





