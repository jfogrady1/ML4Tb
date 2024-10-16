## Script to use logistic regression to build a model to classify animals as control and bTB+ animals.

# Load in our packages
library(data.table)
library(tidyverse)
library(caret)
library(glmnet)
library(DESeq2)
library(ggplot2)
library(magrittr)
set.seed(12789)


# Load in data and labels
data_raw <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt")

# Labels and wrangling
labels <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt")
colnames(labels) <- labels[1,]
labels <- labels[-1,]
labels <- labels %>% select(Sample, Condition)

# Check to make sure everything is in order - labels need to be in same order
all(colnames(data[,2:length(colnames(data))]) == labels$Sample)
genes <- data %>% select(Geneid) %>% as.vector()

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
data_caret_parition$Condition

data_caret_parition <- as.data.frame(data_caret_parition)
data_caret_parition$Condition <- factor(data_caret_parition$Condition, levels = c("Control", "Infected"))
trainIndex <- createDataPartition(data_caret_parition[["Condition"]], p = .7, 
                                  list = FALSE, 
                                  times = 1)
head(trainIndex)

# Subset the data into training and testing sets
train_data <- data_caret_parition[trainIndex, ]
test_data <- data_caret_parition[-trainIndex, ]


train_counts <- t(train_data[, 1:(ncol(train_data) - 1)])  # Only count columns
head(train_data)
train_labels <- train_data[, "Condition", drop = FALSE]         # Label column as a dataframe
head(train_labels)


# Now make the DESEQ2 object
# Prepare DESeq2 dataset
ddsTrain  <- DESeqDataSetFromMatrix(countData = train_counts, 
                              colData = train_labels, 
                              design = ~ Condition)  # Adjust design based on the outcome variable

# Normalise
ddsTrain <- DESeq(ddsTrain)
# variance stabilised transformation on the training data
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(ddsTrain, blind = FALSE)



"The variance stabilizing transformation from a previous dataset can be 'frozen' and reapplied to new samples. The frozen VST is accomplished by saving the dispersion function accessible with dispersionFunction, assigning this to the DESeqDataSet with the new samples, and running
varianceStabilizingTransformation with 'blind' set to FALSE."

# Get the testing set organised
test_counts <- t(test_data[, 1:(ncol(test_data) - 1)])  
test_labels <- test_data[, "Condition", drop = FALSE]         
head(test_labels)

# Get into DESEQ2 object
ddsTest <- DESeqDataSetFromMatrix(countData = test_counts, colData = test_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddsTest <- DESeq(ddsTest)
dispersionFunction(ddsTest) <- dispersionFunction(ddsTrain) # This is related to point above from Mike Love
vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation

train_counts_normalised  <- assay(vstNormalizedExpressionDataForTrain)
test_counts_normalised <- assay(vstNormalizedExpressionDataForTest)

# Now have both normalised train and test - get labels back in
train_counts_normalised <- as.data.frame(t(train_counts_normalised))
train_counts_normalised <- cbind(train_counts_normalised, train_labels)

test_counts_normalised <- as.data.frame(t(test_counts_normalised))
test_counts_normalised <- cbind(test_counts_normalised, test_labels)

# Need to adjust for model.matrix function
options(expressions = 5e5)
train_test <- train_counts_normalised[,c(2000:2500, length(colnames(train_counts_normalised)))]

x <- model.matrix(Condition~., train_counts_normalised)[,-1]


default_glm_mod = train(
  form = Condition ~ .,
  data = train_test,
  trControl = trainControl(method = "cv", number = 5),
  method = "glm",
  family = "binomial" # using a factor
)
default_glm_mod
train()