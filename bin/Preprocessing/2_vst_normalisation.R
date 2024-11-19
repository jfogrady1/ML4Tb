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
library(fastICA)
library("boot")
set.seed(12789)



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
trainIndex <- createDataPartition(data_caret_parition[["Condition"]], p = .70, 
                                  list = FALSE, 
                                  times = 1)


# Subset the data into training and testing sets
train_data_raw <- data_caret_parition[trainIndex, ]
test_data_raw <- data_caret_parition[-trainIndex, ]


train_counts <- t(train_data_raw[, 1:(ncol(train_data_raw) - 1)])  # Only count columns
train_labels <- train_data_raw[, "Condition", drop = FALSE]         # Label column as a dataframe






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
saveRDS(ddsTrain, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/ogrady_train.rds")
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
saveRDS(ddsTest, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/ogrady_test.rds")


# Apply the dispersion function on dds test
dispersionFunction(ddsTest) <- dispersionFunction(ddsTrain) # This is related to point above from Mike Love


vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation


# Visualise - should be homoskedasic
meanSdPlot(assay(vstNormalizedExpressionDataForTrain))
meanSdPlot(assay(vstNormalizedExpressionDataForTest))


# extract the counts form transformed object
train_counts_normalised  <- assay(vstNormalizedExpressionDataForTrain)
test_counts_normalised <- assay(vstNormalizedExpressionDataForTest)


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


head(duplicated(ensemble$gene_name))
ensemble$gene_name <- if_else(duplicated(ensemble$gene_name), ensemble$gene_id, ensemble$gene_name)
head(ensemble,20)

table(duplicated(ensemble$gene_name))

tested_genes <- data.frame(rownames(train_counts_normalised))
colnames(tested_genes) <- "gene_id"
head(tested_genes$gene_id,20)
head(ensemble$gene_name,20)

all(tested_genes$gene_id == ensemble$gene_id)
ensemble$gene_name <- gsub(' ', '', ensemble$gene_name)
ensemble$gene_id <- gsub(' ', '', ensemble$gene_id)
head(tested_genes$gene_id,20)
head(ensemble$gene_id,20)

tested_genes <- left_join(tested_genes, ensemble, by = c("gene_id" = "gene_id"))

head(tested_genes,20)
all(tested_genes$gene_id == rownames(train_counts_normalised))

rownames(train_counts_normalised) <- tested_genes$gene_name


rownames(test_counts_normalised) <- tested_genes$gene_name
table(duplicated(rownames(test_counts_normalised)))


write.table(train_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/OGrady_train_vst_normalised_data.txt", quote = FALSE, sep = "\t")
write.table(test_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/OGrady_test_vst_normalised_data.txt", quote = FALSE, sep = "\t")
write.table(train_labels, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/Train_labels.txt", quote = FALSE, row.names = FALSE, sep = "\t")
write.table(test_labels, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/Test_labels.txt", quote = FALSE, sep = "\t")





# Read in other dataset - Here we will normalise the data with respect to the O'Grady train DESEQ2 transformation
# WIARDA
wiarda = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt", sep = "\t") %>% select(-1)
wiarda_labels = fread("/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv", sep = "\t")
rownames(wiarda_labels) <- wiarda_labels$Animal_Code
rownames(wiarda) <- tested_genes$gene_name
ddswiarda <- DESeqDataSetFromMatrix(countData = wiarda, colData = wiarda_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddswiarda <- DESeq(ddswiarda)
saveRDS(ddswiarda, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/wiarda.rds")
dispersionFunction(ddswiarda) <- dispersionFunction(ddsTrain)
vstNormalizedExpressionDataForwiarda <- varianceStabilizingTransformation(ddswiarda, blind = FALSE) # Now perform the normalisation
wiarda_counts_normalised <- assay(vstNormalizedExpressionDataForwiarda)
write.table(wiarda_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/wiarda_ogrady_vst_normalised_data.txt", quote = FALSE, sep = "\t")

# Kirsten data
Kirsten <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt", sep = "\t") %>% select(-1)
Kirsten <- as.matrix(Kirsten)
rownames(Kirsten) <- tested_genes$gene_name
kirsten_labels <- fread("/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt") %>% as.data.frame()
ddsKirsten <- DESeqDataSetFromMatrix(countData = Kirsten, colData = kirsten_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddsKirsten <- DESeq(ddsKirsten)
saveRDS(ddsKirsten, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/kirsten.rds")
dispersionFunction(ddsKirsten) <- dispersionFunction(ddsTrain) # This is related to point above from Mike Love
vstNormalizedExpressionDataForKirsten <- varianceStabilizingTransformation(ddsKirsten, blind = FALSE) # Now perform the normalisation
kirsten_counts_normalised <- assay(vstNormalizedExpressionDataForKirsten)
write.table(kirsten_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_ogrady_vst_normalised_data.txt", quote = FALSE, sep = "\t")

# Kirsten PBL data
Kirsten_pbl <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt", sep = "\t") %>% select(-1)
Kirsten_pbl <- as.matrix(Kirsten_pbl)
rownames(Kirsten_pbl) <- tested_genes$gene_name
kirsten_pbl_labels <- fread("/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv") %>% as.data.frame()
rownames(kirsten_pbl_labels) <- kirsten_pbl_labels$Run_Code
ddsKirsten_pbl <- DESeqDataSetFromMatrix(countData = Kirsten_pbl, colData = kirsten_pbl_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddsKirsten_pbl <- DESeq(ddsKirsten_pbl)
saveRDS(ddsKirsten_pbl, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/kirsten_pbl.rds")
dispersionFunction(ddsKirsten_pbl) <- dispersionFunction(ddsTrain) # This is related to point above from Mike Love
vstNormalizedExpressionDataForKirsten_pbl <- varianceStabilizingTransformation(ddsKirsten_pbl, blind = FALSE) # Now perform the normalisation
kirsten_pbl_counts_normalised <- assay(vstNormalizedExpressionDataForKirsten_pbl)
write.table(kirsten_pbl_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_pbl_ogrady_vst_normalised_data.txt", quote = FALSE, sep = "\t")

# Abdelaal
abdelaal = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/abdelaal/Quantification/abdelaal_count_matrix_clean.txt", sep = "\t") %>% select(-1)
abelaal <- as.matrix(abdelaal)
abdelaal_labels = fread("/home/workspace/jogrady/ML4TB/data/abdelaal/abdelaal_samples.csv", sep = "\t")
abdelaal_labels <- abdelaal_labels[seq(1,48,2),]
rownames(abdelaal_labels) <- unique(abdelaal_labels$Animal_Code)
rownames(abdelaal) <- tested_genes$gene_name
ddsabdelaal <- DESeqDataSetFromMatrix(countData = abdelaal, colData = abdelaal_labels, design = ~ 1) # keep to 1 so that it doesnt know the labels
ddsabdelaal <- DESeq(ddsabdelaal)
saveRDS(ddsabdelaal, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/abdelaal.rds")
dispersionFunction(ddsabdelaal) <- dispersionFunction(ddsTrain)
vstNormalizedExpressionDataForabdelaal <- varianceStabilizingTransformation(ddsabdelaal, blind = FALSE) # Now perform the normalisation
abdelaal_counts_normalised <- assay(vstNormalizedExpressionDataForabdelaal)
write.table(abdelaal_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/abdelaal_ogrady_vst_normalised_data.txt", quote = FALSE, sep = "\t")






# Now will normalise pgrady (test) for all the other data sets
# First get the dds objects vst

# Regenerate wiarda

ddswiarda <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/wiarda.rds")
ddsKirsten <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/kirsten.rds") 
ddsKirsten_pbl <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/kirsten_pbl.rds")
ddsabdelaal <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/abdelaal.rds") 


# Start with O'Grady
# regenerate
ddsTest <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/ogrady_test.rds")
ddsTrain <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/ogrady_train.rds")


# Now apply vst normalisation using frozen vst paramaters from each of the "training sets"
dispersionFunction(ddsTest) <- dispersionFunction(ddswiarda)
vstNormalizedExpressionDataForogrady <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation
ogrady_counts_normalised <- assay(vstNormalizedExpressionDataForogrady)
write.table(ogrady_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/ogradytest_wiarda_vst_normalised_data.txt", quote = FALSE, sep = "\t")




dispersionFunction(ddsTest) <- dispersionFunction(ddsKirsten)
vstNormalizedExpressionDataForogrady <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation
ogrady_counts_normalised <- assay(vstNormalizedExpressionDataForogrady)
write.table(ogrady_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/ogradytest_kirsten_vst_normalised_data.txt", quote = FALSE, sep = "\t")


dispersionFunction(ddsTest) <- dispersionFunction(ddsKirsten_pbl)
vstNormalizedExpressionDataForogrady <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation
ogrady_counts_normalised <- assay(vstNormalizedExpressionDataForogrady)
write.table(ogrady_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/ogradytest_kirsten_pbl_vst_normalised_data.txt", quote = FALSE, sep = "\t")




dispersionFunction(ddsTest) <- dispersionFunction(ddsabelaal)
vstNormalizedExpressionDataForogrady <- varianceStabilizingTransformation(ddsTest, blind = FALSE) # Now perform the normalisation
ogrady_counts_normalised <- assay(vstNormalizedExpressionDataForogrady)
write.table(ogrady_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/ogradytest_abdelaal_vst_normalised_data.txt", quote = FALSE, sep = "\t")




# Lets write a function

Test_normalisation <- function(ddsnormalised, ddsReference, output_file_path) {
  dispersionFunction(ddsnormalised) <- dispersionFunction(ddsReference)
  vstNormalizedExpressionData <- varianceStabilizingTransformation(ddsnormalised, blind = FALSE)
  counts_normalised <- assay(vstNormalizedExpressionData)
  write.table(counts_normalised, file = output_file_path, quote = FALSE, sep = "\t")
  print(head(counts_normalised))
}



# Now do wiarda toe everything else

#WIARDA
Test_normalisation(ddsnormalised = ddswiarda, ddsReference = ddsabdelaal , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/wiarda_abdelaal_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddswiarda, ddsReference = ddsKirsten , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/wiarda_kirsten_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddswiarda, ddsReference = ddsKirsten_pbl, output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/wiarda_kirsten_pbl_vst_normalised_data.txt")


#ABDELAAL
Test_normalisation(ddsnormalised = ddsabdelaal, ddsReference = ddswiarda , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/abdelaal_wiarda_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsabdelaal, ddsReference = ddsKirsten , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/abdelaal_kirsten_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsabdelaal, ddsReference = ddsKirsten_pbl, output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/abdelaal_kirsten_pbl_vst_normalised_data.txt")



# Kirsten
Test_normalisation(ddsnormalised = ddsKirsten, ddsReference = ddswiarda , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_wiarda_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsKirsten, ddsReference = ddsabdelaal , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_abdelaal_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsKirsten, ddsReference = ddsKirsten_pbl, output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_kirsten_pbl_vst_normalised_data.txt")




# Kirsten_pbl
Test_normalisation(ddsnormalised = ddsKirsten_pbl, ddsReference = ddswiarda , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_pbl_wiarda_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsKirsten_pbl, ddsReference = ddsabdelaal , output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_pbl_abdelaal_vst_normalised_data.txt")
Test_normalisation(ddsnormalised = ddsKirsten_pbl, ddsReference = ddsKirsten_pbl, output_file_path = "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/kirsten_pbl_kirsten_vst_normalised_data.txt")










#################################

# Output explanation
###############################
### Note files are written with first name being data set that IS normalised and second name that the dispersion information was derived from.
# Essentially, for comparison of models on external datasets, the idea is to take a model made by a data set (second name) and use the first name to identify the dataset in use that is normalised with respect to the dataset used to generate the model.

