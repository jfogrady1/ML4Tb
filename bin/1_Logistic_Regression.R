## Script to use logistic regression to build a model to classify animals as control and bTB+ animals.

# Load in our packages
library(data.table)
library(tidyverse)
library(caret)
library(glmnet)
library(DESeq2)
library(ggplot2)
library(magrittr)
library(vsn)
library(nmf)
library("boot")
set.seed(12789)


##########################################################################

# Quick function to estimate number of PCs to include from PCA4QTL #######
# Need to include as do not want to scale the vst transformed RNA-seq data
# The genes with high variance are typically those in which the samples are different for biological reasons, eg differential expression. If you squash all the genes to have equal variance you demote the biological signal and promote the noise. I’ve never understood why people suggest to do this. M.I.L
##########################################################################
JFOG_runBE<-function(X,B=20,alpha=0.05,
                mc.cores=min(B,parallel::detectCores()-1),
                verbose=FALSE){

  if(alpha<0 || alpha>1){
    stop("alpha must be between 0 and 1.")
  }

  n<-nrow(X) #Number of observations.
  p<-ncol(X) #Number of features.
  d<-min(n,p) #This is the total number of PCs.

  if(verbose) cat("Running PCA on permuted data...\n")
  results<-parallel::mclapply(1:B,FUN=function(b){
    # b<-3
    if(verbose) cat("b=",b," out of ",B," permutations...\n",sep="")

    #Permute each column of X. That is, permute the observations in each feature.
    XPermuted<-matrix(data=NA,nrow=n,ncol=p)
    for(j in 1:p){
      # j<-7
      XPermuted[,j]<-sample(x=X[,j],size=n,replace=FALSE)
    }

    prcompResultPerm<-prcomp(x=XPermuted,center=TRUE,scale. = FALSE) #Key step.
    importanceTablePerm<-summary(prcompResultPerm)$importance
    PVEsPerm<-importanceTablePerm[2,]
    return(PVEsPerm)
  },mc.cores=mc.cores) #results is a list of vectors.
  temp<-unlist(results)
  testStatsPerm<-matrix(data=temp,nrow=d,byrow=FALSE) #PC by permutation.

  if(verbose) cat("Running PCA on the unpermuted data...\n")
  prcompResult<-prcomp(x=X,center=TRUE,scale. = FALSE) #Key step.
  importanceTable<-summary(prcompResult)$importance
  PVEs<-importanceTable[2,]
  # Compare PVEs to testStatsPerm.
  # temp<-(testStatsPerm>=PVEs) #temp is calculated as desired.
  pValues<-(rowSums(testStatsPerm>=PVEs)+1)/(B+1) #The p-value for the jth PC is calculated as, roughly speaking, the proportion of permutations where the PVE of the jth PC is greater than or equal to PVE_j.

  for(j in 2:d){ #Enforce monotone increase of the p-values.
    if(pValues[j]<pValues[j-1]){
      pValues[j]<-pValues[j-1]
    }
  }

  numOfPCsChosen<-sum(pValues<=alpha)
  toReturn<-list(pValues=pValues,alpha=alpha,numOfPCsChosen=numOfPCsChosen)
  return(toReturn)
}




##########################################################################

# Read in data and load
##########################################################################

# raw counts from featurecounts
data_raw <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt")

# Labels and wrangling
labels <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt")
colnames(labels) <- labels[1,]
labels <- labels[-1,]
labels <- labels %>% select(Sample, Condition)

# Check to make sure everything is in order - labels need to be in same order
all(colnames(data_raw[,2:length(colnames(data_raw))]) == labels$Sample)

# Extract gene names as will be important later
genes <- data_raw %>% select(Geneid) %>% as.vector()


# Convert to matrix and name rows and columns
# Need to convert count matrix to a format which can be used by Caret (transpose)
data_caret_parition <- data_raw[,-1]
data_caret_parition <- as.matrix(data_caret_parition)
rownames(data_caret_parition) <- genes$Geneid
data_caret_parition <- as.data.frame(t(data_caret_parition))


# Bind the labels together and rename
data_caret_parition <- cbind(data_caret_parition, labels$Condition)
data_caret_parition[, "labels$Condition"]
colnames(data_caret_parition)[length(colnames(data_caret_parition))] <- "Condition"


# Some wrangling
data_caret_parition <- as.data.frame(data_caret_parition)

# Factorise the condition variable
data_caret_parition$Condition <- factor(data_caret_parition$Condition, levels = c("Control", "Infected"))


# Create the parition between the training and testing
trainIndex <- createDataPartition(data_caret_parition[["Condition"]], p = .75, 
                                  list = FALSE, 
                                  times = 1)


# Subset the data into training and testing sets
train_data_raw <- data_caret_parition[trainIndex, ]
test_data_raw <- data_caret_parition[-trainIndex, ]


train_counts <- t(train_data_raw[, 1:(ncol(train_data_raw) - 1)])  # Only count columns
train_labels <- train_data_raw[, "Condition", drop = FALSE]         # Label column as a dataframe


# Create folds
# This is required for cross validation
folds_train <- createFolds(train_labels$Condition, k = 5, list = TRUE, returnTrain = FALSE)

foldid_train <- rep(NA, length(train_labels))
for (i in seq_along(folds_train)) {
  foldid_train[folds_train[[i]]] <- i
}





##########################################################################

# DESEQ2 VST normalisation from DESEQ2

# Note: "The variance stabilizing transformation from a previous dataset can be 'frozen' and reapplied to new samples. The frozen VST is accomplished by saving the dispersion function accessible with dispersionFunction, assigning this to the DESeqDataSet with the new samples, and running varianceStabilizingTransformation with 'blind' set to FALSE." - M.I.L
##########################################################################


# Now make the DESEQ2 object
# Prepare DESeq2 dataset
ddsTrain  <- DESeqDataSetFromMatrix(countData = train_counts, 
                              colData = train_labels, 
                              design = ~ Condition)  # Adjust design based on the outcome variable

meanSdPlot(assay(ddsTrain))
# Normalise
ddsTrain <- DESeq(ddsTrain)
# variance stabilised transformation on the training data
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(ddsTrain, blind = FALSE)

meanSdPlot(assay(vstNormalizedExpressionDataForTrain))

"The variance stabilizing transformation from a previous dataset can be 'frozen' and reapplied to new samples. The frozen VST is accomplished by saving the dispersion function accessible with dispersionFunction, assigning this to the DESeqDataSet with the new samples, and running
varianceStabilizingTransformation with 'blind' set to FALSE."



# Get the testing set organised
test_counts <- t(test_data_raw[, 1:(ncol(test_data_raw) - 1)])  
test_labels <- test_data_raw[, "Condition", drop = FALSE]         




# Get into DESEQ2 object
ddsTest <- DESeqDataSetFromMatrix(countData = test_counts, colData = test_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddsTest <- DESeq(ddsTest)


# Apply the dispersion function on dds test
dispersionFunction(ddsTest) <- dispersionFunction(ddsTrain) # This is related to point above from Mike Love


vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation


# Visualise - should be homoskedasic
meanSdPlot(assay(vstNormalizedExpressionDataForTrain))
meanSdPlot(assay(vstNormalizedExpressionDataForTest))


# extract the counts form transformed object
train_counts_normalised  <- assay(vstNormalizedExpressionDataForTrain)
test_counts_normalised <- assay(vstNormalizedExpressionDataForTest)




# Now have both normalised train and test - get labels back in
train_counts_normalised <- t(train_counts_normalised)
test_counts_normalised <- t(test_counts_normalised)






results_logistic_genes_l.min <- data.frame()
results_logistic_genes_l.se <- data.frame()
# Generate list of all values
list_of_fits <- list()

# Perform logistic regression for various regularisation paramaters.
for(i in 0:10) {
  print(i)

  # Save name
  fit.name <- paste0("alpha", i/10)
  list_of_fits[[fit.name]] <- cv.glmnet(train_counts_normalised, train_labels$Condition, type.measure="class", alpha = i/10, family = "binomial", foldid = foldid_train)

  # Save model
  temp_model <- list_of_fits[[fit.name]]

  # Save model complexity with different lambda paramaters
  non_zero_count <- sum(coef(temp_model, s = "lambda.1se") != 0) - 1
  non_zero_count_min <- sum(coef(temp_model, s = "lambda.min") != 0) - 1

  # Predict back on the train
  predicted_1.se <- predict(list_of_fits[[fit.name]], s=list_of_fits[[fit.name]]$lambda.1se, newx=train_counts_normalised, type = "class")


  # Extract various metrics
  conf_matrix <- caret::confusionMatrix(table(predicted_1.se, train_labels$Condition), mode = "everything", positive = "Infected")
  Accuracy <- conf_matrix$overall["Accuracy"]
  Sensitivity <- conf_matrix$byClass["Sensitivity"]
  Specificity <- conf_matrix$byClass["Specificity"]
  Precision <- conf_matrix$byClass["Precision"]
  Recall <- conf_matrix$byClass["Recall"]
  Pos_pred_value <- conf_matrix$byClass["Pos Pred Value"]
  Neg_pred_value <- conf_matrix$byClass["Neg Pred Value"]
  F1 <- conf_matrix$byClass["F1"]
  temp <- data.frame(Comparison = "Training_Set_lambda_1se", alpha=i/10, lambda =list_of_fits[[fit.name]]$lambda.1se, fit.name=fit.name, non_zero_variables=non_zero_count, Accuracy = Accuracy, Sensitivity = Sensitivity, Specificity = Specificity, Precision = Precision, Recall = Recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  results_logistic_genes_l.se <- rbind(results_logistic_genes_l.se, temp)

  # Predict back on the train
  predicted_l.min <- predict(list_of_fits[[fit.name]], s=list_of_fits[[fit.name]]$lambda.min, newx=train_counts_normalised, type = "class")

  #Extract metrics for lambda .min
  conf_matrix <- caret::confusionMatrix(table(predicted_l.min, train_labels$Condition), mode = "everything", positive = "Infected")
  Accuracy <- conf_matrix$overall["Accuracy"]
  Sensitivity <- conf_matrix$byClass["Sensitivity"]
  Specificity <- conf_matrix$byClass["Specificity"]
  Precision <- conf_matrix$byClass["Precision"]
  Recall <- conf_matrix$byClass["Recall"]
  Pos_pred_value <- conf_matrix$byClass["Pos Pred Value"]
  Neg_pred_value <- conf_matrix$byClass["Neg Pred Value"]
  F1 <- conf_matrix$byClass["F1"]
  temp <- data.frame(Comparison = "Training_Set_lambda_min", alpha=i/10, lambda =list_of_fits[[fit.name]]$lambda.min, fit.name=fit.name, non_zero_variables=non_zero_count_min, Accuracy = Accuracy, Sensitivity = Sensitivity, Specificity = Specificity, Precision = Precision, Recall = Recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  results_logistic_genes_l.min <- rbind(results_logistic_genes_l.min, temp)
}


# Now predict on the testing set

for (i in 0:10) {
  fit.name <- paste0("alpha", i/10)

  temp_model <- list_of_fits[[fit.name]]
  non_zero_count <- sum(coef(temp_model, s = "lambda.1se") != 0) - 1
  non_zero_count_min <- sum(coef(temp_model, s = "lambda.min") != 0) - 1
  ## Use each model to predict 'y' given the Testing dataset
  predicted_1.se <- predict(list_of_fits[[fit.name]], s=list_of_fits[[fit.name]]$lambda.1se, newx=test_counts_normalised, type = "class")

  conf_matrix <- caret::confusionMatrix(table(predicted_1.se, test_labels$Condition), mode = "everything", positive = "Infected")
  Accuracy <- conf_matrix$overall["Accuracy"]
  Sensitivity <- conf_matrix$byClass["Sensitivity"]
  Specificity <- conf_matrix$byClass["Specificity"]
  Precision <- conf_matrix$byClass["Precision"]
  Recall <- conf_matrix$byClass["Recall"]
  Pos_pred_value <- conf_matrix$byClass["Pos Pred Value"]
  Neg_pred_value <- conf_matrix$byClass["Neg Pred Value"]
  F1 <- conf_matrix$byClass["F1"]


  ## Store the results
  temp <- data.frame(Comparison = "Testing_Set_lambda_1se", alpha=i/10, lambda =list_of_fits[[fit.name]]$lambda.1se, fit.name=fit.name, non_zero_variables=non_zero_count, Accuracy = Accuracy, Sensitivity = Sensitivity, Specificity = Specificity, Precision = Precision, Recall = Recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  results_logistic_genes_l.se <- rbind(results_logistic_genes_l.se, temp)

  predicted_l.min <- predict(list_of_fits[[fit.name]], s=list_of_fits[[fit.name]]$lambda.min, newx=test_counts_normalised, type = "class")

  #Extract metrics for lambda .min
  conf_matrix <- caret::confusionMatrix(table(predicted_l.min, test_labels$Condition), mode = "everything", positive = "Infected")
  Accuracy <- conf_matrix$overall["Accuracy"]
  Sensitivity <- conf_matrix$byClass["Sensitivity"]
  Specificity <- conf_matrix$byClass["Specificity"]
  Precision <- conf_matrix$byClass["Precision"]
  Recall <- conf_matrix$byClass["Recall"]
  Pos_pred_value <- conf_matrix$byClass["Pos Pred Value"]
  Neg_pred_value <- conf_matrix$byClass["Neg Pred Value"]
  F1 <- conf_matrix$byClass["F1"]
  temp <- data.frame(Comparison = "Training_Set_lambda_min", alpha=i/10, lambda =list_of_fits[[fit.name]]$lambda.min, fit.name=fit.name, non_zero_variables=non_zero_count_min, Accuracy = Accuracy, Sensitivity = Sensitivity, Specificity = Specificity, Precision = Precision, Recall = Recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  results_logistic_genes_l.min <- rbind(results_logistic_genes_l.min, temp)
}

View(results_logistic_genes_l.min)
View(results_logistic_genes_l.se)


#################################################################################
#### Latent variable analysis
#################################################################################


################################################################################
################################################################################

#### 1. PCA

################################################################################
################################################################################

# Apply the MAD function row-wise (across genes) - use 1 for this
mad_values <- apply(as.matrix(train_counts), 1, mad)

# Step 3: Create a data frame that combines gene names (row names) and MAD values
mad_df <- data.frame(gene = rownames(as.matrix(train_counts)), mad = mad_values)


# Step 4: Sort the data frame by MAD values in descending order
mad_df_sorted <- mad_df %>% arrange(desc(mad)) %>% as.data.frame()

top_genes <- mad_df_sorted[1:8000,]






# Perform PCA - center but do not scale
prcompResult <- prcomp(train_counts_normalised, center = TRUE, scale. = FALSE)
PCs<-prcompResult$x #The columns are PCs.


# Extract elbow and BE method
resultRunElbow<-PCAForQTL::runElbow(prcompResult=prcompResult)
resultRunElbow
RNGkind("L'Ecuyer-CMRG")
resultRunBE <- JFOG_runBE(train_counts_normalised,B=20,alpha=0.05, verbose = TRUE)
print(resultRunBE$numOfPCsChosen)

# Extract elbow and BE
K_elbow<-resultRunElbow
K_BE<-resultRunBE$numOfPCsChosen

# Scree plot
PCAForQTL::makeScreePlot(prcompResult,labels=c("Elbow","BE"),values=c(K_elbow,K_BE),
                         titleText="Training data scree plot")




# Just do it again
prcompResult <- prcomp(train_counts_normalised, center = TRUE, scale. = FALSE)

train_labels$Condition
# Predict the normalised
train_PCA_predict <- predict(prcompResult, train_counts_normalised)
train_labels$Condition
train_PCA_predict <- data.frame(train_PCA_predict, train_labels$Condition)
colnames(train_PCA_predict)[length(colnames(train_PCA_predict))] <- "Condition"

head(train_PCA_predict)

# Cannot do lasso/ridge regression on PCA as the variables are not correlated to each other, hence glmnet isn't appropriate
# Will therefore do CV.glm but we need to include our own folds etc


K = 5
cv_accuracy <- numeric(5)
cv_auc <- numeric(5)
cv_error <- numeric(5)



head(train_data)


# Make a list of all the results
PCA_cv_results <- list()

# Save the values for train data set

head(test_data)
dim(train_PCA_predict)
PCA_results <- data.frame()
train_PCA_predict$Condition
for (n in 1:dim(train_PCA_predict)[1]) {

  # Initialise the variables
  PCA_train_accuracy_values <- c()
  PCA_train_sensitivity_values <- c()
  PCA_train_specificity_values <- c()
  PCA_train_precision_values <- c()
  PCA_train_recall_values <- c()
  PCA_train_Pos_pred_value_values <- c()
  PCA_train_Neg_pred_value_values <- c()
  PCA_train_F1 <- c()
  for (i in 1:K) {
    # Train on all data except the ith fold, and test on the ith fold
    train_PCA_predict_temp <- train_PCA_predict # Select the max number of PCs
    test_indices <- which(foldid_train == i)
    train_data <- train_PCA_predict_temp[-test_indices, ]
    test_data <- train_PCA_predict_temp[test_indices, ]

    
    # Fit the logistic model on the training data
    glm_fit <- glm(Condition ~ ., data = train_data[,c(1:n, length(colnames(train_PCA_predict_temp)))], family = binomial (link='logit'))
    # Perform cross-validation (calculate error on test data)
    if (n == 1) {
      print("n == 1")
      test_data <- as.data.frame(test_data[,1])
      colnames(test_data) <- "PC1"
      PCA_cv_results_prob <- predict(glm_fit, newdata = test_data[,c(1:n)], type = "response")
    }
    else { PCA_cv_results_prob <- predict(glm_fit, newdata = test_data[,c(1:n)], type = "response")
    }
    
    PCA_cv_results[[paste0("PCA_Fold_", i)]] <- PCA_cv_results_prob

    PCA_cv_results_class <- if_else(PCA_cv_results_prob > 0.5, "Infected", "Control") 
    
    conf_matrix <- caret::confusionMatrix(table(PCA_cv_results_class, test_data$Condition), mode = "everything", positive = "Infected")
        
    PCA_cv_results[[paste0("PCA_Fold_Class_", i)]] <- PCA_cv_results_class
    PCA_train_sensitivity_values[i] <- conf_matrix$byClass['Sensitivity']
    PCA_train_specificity_values[i] <- conf_matrix$byClass['Specificity']
    PCA_train_accuracy_values[i] <- conf_matrix$overall['Accuracy']
    PCA_train_precision_values[i] <- conf_matrix$byClass['Precision']
    PCA_train_recall_values[i] <- conf_matrix$byClass['Recall']
    PCA_train_Pos_pred_value_values[i] <- conf_matrix$byClass["Pos Pred Value"]
    PCA_train_Neg_pred_value_values[i] <- conf_matrix$byClass["Neg Pred Value"]
    PCA_train_F1[i] <-  conf_matrix$byClass["F1"]
  }

  accuracy <- mean(PCA_train_accuracy_values)
  sensitivity <- mean(PCA_train_sensitivity_values)
  specificity <- mean(PCA_train_specificity_values)
  precision <- mean(PCA_train_precision_values)
  recall <- mean(PCA_train_recall_values)
  Pos_pred_value <- mean(PCA_train_Pos_pred_value_values)
  Neg_pred_value <- mean(PCA_train_Neg_pred_value_values)
  F1 <- mean(PCA_train_F1)

  temp <- data.frame(Comparison = "Training_PCA", number_of_PCs = n, Accuracy = accuracy, Sensitivity = sensitivity, Specificity = specificity, Precision = precision, Recall = recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  PCA_results <- rbind(PCA_results, temp)
}


train_PCA_predict <- predict(prcompResult, train_counts_normalised)

# Now predict on the actual dataset






View(PCA_results)

temp

c(1:n, length(colnames(train_PCA_predict_temp)))

avg_sensitivity <- mean(sensitivity_values)
avg_specificity <- mean(PCA_train_specificity_values)
avg_specificity









PCA_cv_results$PCA_Fold_Class_5
library(pROC)
roc_obj <- roc(test_labels$Condition, PCA_cv_results_prob, levels=c("Control", "Infected"))

PCA_cv_results
summary(PCA_cv_results)
PCA_cv_error

  list_of_fits[[fit.name]] <- cv.glmnet(train_PCA_predict, train_labels$Condition, type.measure="class", alpha = 1, family = "binomial", foldid = foldid_train)

  temp_model <- list_of_fits[[fit.name]]
  predicted <- predict(list_of_fits[[fit.name]], s=list_of_fits[[fit.name]]$lambda.1se, newx=train_PCA_predict, type = "class")
  
  conf_matrix <- caret::confusionMatrix(table(predicted, train_labels$Condition), mode = "everything", positive = "Infected")
  Accuracy <- conf_matrix$overall["Accuracy"]
  print(Accuracy)
  Sensitivity <- conf_matrix$byClass["Sensitivity"]
  Specificity <- conf_matrix$byClass["Specificity"]
  Precision <- conf_matrix$byClass["Precision"]
  Recall <- conf_matrix$byClass["Recall"]
  Pos_pred_value <- conf_matrix$byClass["Pos Pred Value"]
  Neg_pred_value <- conf_matrix$byClass["Neg Pred Value"]
  F1 <- conf_matrix$byClass["F1"]
  temp <- data.frame(Comparison = "Training Set", alpha=i/10, fit.name=fit.name, Accuracy = Accuracy, Sensitivity = Sensitivity, Specificity = Specificity, Precision = Precision, Recall = Recall, Pos_pred_value = Pos_pred_value, Neg_pred_value = Neg_pred_value, F1 = F1)
  results_PCA <- rbind(results_PCA, temp)
  if (Accuracy == 1) {

  }


View(cor(train_PCA_predict))

??cv.glm()

list_of_fits
prcompResult_test <- prcomp(test_counts_normalised, center = TRUE, scale. = FALSE)
PCs_test<-prcompResult_test$x #The columns are PCs.

  temp_model <- list_of_fits[[fit.name]]
  non_zero_count <- sum(coef(temp_model, s = "lambda.1se") != 0) - 1
  ## Use each model to predict 'y' given the Testing dataset
predicted <- predict(list_of_fits[["PC18"]], newx=test_counts_normalised, type = "prob")
predicted_classes <- predict(list_of_fits[["PC18"]], newdata =test_counts_normalised)
test_labels
predicted
predicted_classes
print(fit.name)
