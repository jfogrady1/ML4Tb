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
  return(2^(mean(x))) # already on log2 scale (vst counts)
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
      roc_values = as.numeric(pROC::auc(
        response = ifelse(Condition == "Infected", 1, 0),
        predictor = Infected,
        direction = "<",
        quiet = TRUE
      )[1]),
                AUC_CI_lower = (pROC::ci(pROC::roc(Condition, Infected, direction = "<")))[1],
                AUC_CI_upper = (pROC::ci(pROC::roc(Condition, Infected, direction = "<")))[3])
  roc_values <- roc_values %>% arrange(desc(roc_values))
  
  roc_text <- paste(paste0("ROC Resample ", seq_along(roc_values$Model), " = ", roc_values$roc_values, "(95% CI:", round(roc_values$AUC_CI_lower, 2), " - ", round(roc_values$AUC_CI_upper,2), ")"), collapse = "\n")
  resample_colors <- ggsci::pal_npg("nrc")(length(roc_values$roc_values))
  
  if (length(resample_colors) >= 10) {
    resample_colors = c(resample_colors[!is.na(resample_colors)], "darkblue")
  }
  print(resample_colors)
  combined_predictions_ROC <- combined_predictions %>% ggplot(aes(m=Infected, d=factor(Condition, levels = c("Control", "Infected")), colour = factor(Model, levels = roc_values$Model))) + # Note do not need to filter unless saving all predictions
    style_roc() +
    geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 2, colour = "black") +
    geom_roc(n.cuts=0) + 
    coord_equal() +
    scale_colour_manual(values = resample_colors) +
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
               label = paste0("ROC ", roc_values$Model[i], " = ", round(roc_values$roc_values[i],3)," (95% CI: ", round(roc_values$AUC_CI_lower[i], 2), " - ", round(roc_values$AUC_CI_upper[i],2), ")"), 
               hjust = 0, 
               size = 5, 
               color = resample_colors[i])
  }
  
  
  return(combined_predictions_ROC)
}



pred_fun <- function(model, newdata) {
  predict(model, newdata = newdata, type = "prob")[, 2]  # Return probabilities for the positive class
}




# Greedy forward search algorithm
greedy_forward_search = function(gene.list, cpm_matrix, de_results, metadata) {
  
  # Format input metadata to obtain a study specific batch
  metadata$study = metadata$Study
  metadata$Condition = factor(metadata$Condition, levels = c(0, 1), labels = c("Control", "Infected"))
  metadata$study = if_else(metadata$study == "1_OGrady", "OGrady", metadata$study)
  metadata$study = if_else(metadata$study == "2_OGrady", "OGrady", metadata$study)
  metadata$Fold = factor(metadata$Folds, labels = c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5"), levels = c(1,2,3,4,5))

  
  # Count all samples for weighted AUC
  Ogrady_n = as.numeric(sum(metadata$study == "OGrady"))
  Wiarda_n = as.numeric(sum(metadata$study == "Wiarda"))
  Mcloughlin_n = as.numeric(sum(metadata$study == "Mcloughlin"))
  Mcloughlin_pbl_n = as.numeric(sum(metadata$study == "Mcloughlin_pbl"))
  
  View(metadata)
  de_results$variable = rownames(de_results)
  starts = gene.list
  # Obtain a dataframe to collect results
  starts_df = data.frame(gene=character(), Fold1=numeric(), Fold2=numeric(), Fold2=numeric(), Fold3=numeric(), Fold4=numeric(), Fold5 = numeric(), Average=numeric(), stringsAsFactors=FALSE)
  
  for (s in starts){
    
    # For all genes, select the gene and derive the score
    df <- cpm_matrix[rownames(cpm_matrix) == s, ] %>%
      as.data.frame()  %>%
      rownames_to_column(var = "Sample") %>%
      left_join(metadata, by = c("Sample" = "Sample")) %>%
      as.data.frame()
  
    colnames(df)[2] <- "Score"
    
    # Calculate AUC for each of the first gene
    # Note, geometric mean is not calculated as only one gene considered
    auc = df %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(as.vector(unlist(Score))),direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
    colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
    auc$Average <- rowMeans(auc)
    colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5", "Average")
    starts_df = rbind(starts_df, cbind("gene"=s, "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5"=auc[1,5], "Average" = auc[1,6]))
    colnames(starts_df)[7] <- "Average"
  }
  
  # Order based on Weighted average
  starts_df = starts_df %>% as.data.frame() %>% dplyr::arrange(desc(Average)) 
  colnames(starts_df)[7] <- "Average"
  
  # Set up a new dataframe for plotting
  plot_df = matrix(ncol=7, NA)
  colnames(plot_df)=c("combination","Fold1", "Fold2", "Fold3", "Fold4", "Fold5", "Average")
  
  # Select best gene from ordered starts_df
  plot_df = rbind(plot_df, cbind("combination"=starts_df[1,1], "Fold1"=starts_df[1,2], "Fold2"=starts_df[1,3], "Fold3"=starts_df[1,4], "Fold4"=starts_df[1,5], "Fold5"=starts_df[1,6], "Average" = starts_df[1,7]))
  
  # Get the top AUC
  best_auc = plot_df %>% as.data.frame() %>% filter(row_number()==n())
  best_auc = best_auc$Average
  
  # Set up second dataframe for other genes
  df2 = matrix(ncol=7,NA)
  colnames(df2)=c("combination","Fold1", "Fold2", "Fold3", "Fold4", "Fold5", "Average")
  start = plot_df[2,1] # first gene - note NA is on first column
  # Isolate other genes
  rest = starts[!(starts %in% start)]
  for (r in 1:length(rest)){
    
    # Generate combination
    comb = c(start,rest[r])
    print(comb)
    pos = de_results[de_results$Direction == "Positive",] %>% filter(variable %in% comb) %>% pull(variable)
    neg = de_results[de_results$Direction == "Negative",] %>% filter(variable %in% comb)%>% pull(variable)
    

    if (length(pos)==1){
      pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]  %>% as.data.frame() %>% set_colnames("pos_Score")
      
    } else if (length(pos)>1) {
      pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]   %>% 
        as.data.frame() %>% mutate_all(as.numeric) %>%  summarise_all(mean) %>% t() %>% set_colnames("pos_Score")
    } else{
      pos_counts=NA
    }
    
    if (length(neg)==1){
      neg_counts = cpm_matrix[rownames(cpm_matrix) %in% neg,]  %>% as.data.frame() %>% mutate_all(as.numeric)%>% set_colnames("neg_Score")
    } else if (length(neg)>1) {
      neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,] %>%
        as.data.frame() %>% mutate_all(as.numeric) %>% summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
    } else {
      neg_counts = NA
    }
    

    if (any(!is.na(neg_counts)) && any(!is.na(pos_counts))) {
      df_test = merge(pos_counts,neg_counts, by="row.names")
      df_test = merge(df_test, metadata, by.x="Row.names", by.y="Sample")
      df_test$Score = df_test$pos_Score - df_test$neg_Score
      
      
      
      auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score),direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
      colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
      auc$Average <- rowMeans(auc)
      print("HERE1")
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
    } 
      else if (any(is.na(neg_counts)) & any(!is.na(pos_counts))) {
      df_test = merge(pos_counts, metadata, by.x="row.names", by.y="Sample")
      df_test$Score = df_test$pos_Score
      auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score),direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
      colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
      auc$Average <- rowMeans(auc)
      print("HERE2")
      print(auc)
      print(df2)
      df_temp = cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6])
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
      print(df2)
    } 
      else {
      df_test = merge(neg_counts, metadata, by.x="row.names", by.y="Sample")
      df_test$Score = df_test$neg_Score
      auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score), direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
      colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
      auc$Average <- rowMeans(auc)
      df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
    }
  }
  
  while (!is.na(best_auc)){
    sorted_df = df2 %>% as.data.frame()%>% filter(as.numeric(Average) > best_auc) %>% arrange(desc(Average))
    new_start = unlist(str_split(sorted_df[1,1], "_"))
    new_rest = starts[!(starts %in% new_start)]
    print(new_rest)
    df2 = matrix(ncol=7,NA)
    colnames(df2)=c("combination","Fold1", "Fold2", "Fold3", "Fold4", "Fold5", "Average")
    best_auc = as.numeric(sorted_df[1,7])
    plot_df = rbind(plot_df, cbind("combination"=sorted_df[1,1], "Fold1"=sorted_df[1,2], "Fold2"=sorted_df[1,3], "Fold3"=sorted_df[1,4], "Fold4"=sorted_df[1,5], "Fold5" = sorted_df[1,6], "Average" = sorted_df[1,7]))
    print(best_auc)
    for (r in 1:length(new_rest)){
      
      comb = c(new_start,new_rest[r])
      print(r)
      pos = de_results[de_results$Direction == "Positive",] %>% filter(variable %in% comb) %>% pull(variable)
      neg =de_results[de_results$Direction == "Negative",] %>% filter(variable %in% comb)%>% pull(variable)
      if (length(pos)==1){
        pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
      } else if (length(pos)>1) {
        pos_counts =cpm_matrix[rownames(cpm_matrix) %in% pos,]   %>% 
          as.data.frame() %>% mutate_all(as.numeric) %>%  summarise_all(mean) %>% t()%>% set_colnames("pos_Score")    
      } else{
        pos_counts=NA
      }
      
      if (length(neg)==1){
        neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,]   %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("neg_Score")
      } else if (length(neg)>1) {
        neg_counts =cpm_matrix[rownames(cpm_matrix) %in% neg,]   %>% 
          as.data.frame() %>% mutate_all(as.numeric) %>%
          summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
      } else {
        neg_counts = NA
      }
      
      if (any(!is.na(neg_counts)) && any(!is.na(pos_counts))){
        df_test = merge(pos_counts,neg_counts, by="row.names")
        df_test = merge(df_test, metadata, by.x="Row.names", by.y="Sample")
        df_test$Score = df_test$pos_Score - df_test$neg_Score
        auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score),direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
        colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
        auc$Average <- rowMeans(auc)
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
      } 
      else if (any(is.na(neg_counts)) & any(!is.na(pos_counts))){
        df_test = merge(pos_counts, metadata, by.x="row.names", by.y="Sample")
        df_test$Score = df_test$pos_Score
        auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score),direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
        colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
        auc$Average <- rowMeans(auc)
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
      } 
      else {
        df_test = merge(neg_counts, metadata, by.x="row.names", by.y="Sample")
        df_test$Score = df_test$neg_Score
        auc = df_test %>% group_by(Fold) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score), direction = "<"))), .groups = "drop") %>% pivot_wider(names_from = "Fold", values_from = "AUC")
        colnames(auc) <- c("Fold1", "Fold2", "Fold3", "Fold4", "Fold5")
        auc$Average <- rowMeans(auc)
        df2 = rbind(df2, cbind("combination" = paste(comb, collapse="_"), "Fold1"=auc[1,1], "Fold2"=auc[1,2], "Fold3"=auc[1,3], "Fold4"=auc[1,4], "Fold5" = auc[1,5], "Average" = auc[1,6]))
      }
    }
    plot.roc(df_test$Condition, df_test$Score, print.auc=TRUE)
  }
  return(as.data.frame(plot_df))
}

  




