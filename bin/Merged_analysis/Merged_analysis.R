# Merging bovine samples

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
library(pROC)
library(RColorBrewer)
library(grid)
set.seed(42) # For reproducibility



####################################################
####################################################
#####################################################
####
#     Helper functions
####
####################################################
####################################################
####################################################


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



####################################################
####################################################
#####################################################
####
#     Read in raw counts and metadata
####
####################################################
####################################################
####################################################

# raw counts from featurecounts
data_raw <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt")

tested_genes <- data_raw %>% select(Geneid) %>% as.vector()

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

all(tested_genes$Geneid == ensemble$gene_id)

rownames(data_raw) <- ensemble$gene_name





# Labels and wrangling
labels <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt")
colnames(labels) <- labels[1,]
labels
labels <- labels[-1,]
labels <- labels %>% select(Sample, Condition, Age, Batch)
labels$Study <- paste0(labels$Batch, "_OGrady")
labels <- labels %>% select(-c("Batch"))
rownames(labels) <- labels$Sample
labels$Infection_administration <- "Natural"
labels$Sex <- "M"
labels$Location <- "EU"
labels$Tissue <- "PB"
labels$Animal_Code <- labels$Sample




# WIARDA
wiarda = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt", sep = "\t") %>% select(-1)
wiarda_labels = fread("/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv", sep = "\t")
rownames(wiarda_labels) <- wiarda_labels$Animal_Code
rownames(wiarda) <- ensemble$gene_name
wiarda_labels <- wiarda_labels %>% select(Animal_Code, Status)
colnames(wiarda_labels)[2] <- "Condition"
wiarda_labels[1:5,2] <- "Control" # Sampled prior to infection with M. bovis in the time series
colnames(wiarda_labels)[1] <- "Sample"
wiarda_labels$Age <- c(rep(5, 13), rep(6, 13), rep(7,13))
wiarda_labels$Sex <- "M"
wiarda_labels$Infection_administration <- "Experimental"
wiarda_labels$Tissue <- "PBL"
wiarda_labels$Location <- "US"
wiarda_labels$Study <- "Wiarda"
wiarda_labels$Animal_Code <- gsub("", "", wiarda_labels$Sample)
wiarda_labels$Animal_Code <- gsub("_.*", "", wiarda_labels$Animal_Code)
rownames(wiarda_labels) <- wiarda_labels$Sample

# Kirsten data
Kirsten <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt", sep = "\t") %>% select(-1)
Kirsten <- as.matrix(Kirsten)
rownames(Kirsten) <- ensemble$gene_name
kirsten_labels <- fread("/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt") %>% as.data.frame() #%>% select(1,2)
kirsten_labels <- kirsten_labels %>% select(-c(Week))
kirsten_labels$Animal_Code <- gsub("_W.*", "", kirsten_labels$Sample)
rownames(labels) <- labels$Sample
kirsten_labels$Infection_administration <- "Experimental"
kirsten_labels$Sex <- "M"
kirsten_labels$Location <- "EU"
kirsten_labels$Tissue <- "PB"
kirsten_labels$Age <- c(rep(5,18), rep(6, 7), rep(7,9), rep(8,18)) # Animals get older as they are sampled etc
kirsten_labels$Study <- "Mcloughlin"
rownames(kirsten_labels) <- kirsten_labels$Sample


# Kirsten PBL data
Kirsten_pbl <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt", sep = "\t") %>% select(-1)
Kirsten_pbl <- as.matrix(Kirsten_pbl)
rownames(Kirsten_pbl) <- ensemble$gene_name
kirsten_pbl_labels <- fread("/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv") %>% as.data.frame()
rownames(kirsten_pbl_labels) <- kirsten_pbl_labels$Run_Code
kirsten_pbl_labels <- kirsten_pbl_labels %>% select(2,3)
colnames(kirsten_pbl_labels) <- c("Sample", "Condition")
kirsten_pbl_labels$Infection_administration <- "Experimental"
kirsten_pbl_labels$Sex <- "F"
kirsten_pbl_labels$Location <- "EU"
kirsten_pbl_labels$Tissue <- "PBL"
kirsten_pbl_labels$Age <- 12 # M.bovis at elast 12 months old, that is all we know.
kirsten_pbl_labels$Study <- "Mcloughlin_pbl"
kirsten_pbl_labels$Animal_Code <- kirsten_pbl_labels$Sample
rownames(kirsten_pbl_labels)


# Abdelaal
abdelaal = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/abdelaal/Quantification/abdelaal_count_matrix_clean.txt", sep = "\t") %>% select(-1)
abelaal <- as.matrix(abdelaal)
abdelaal_labels = fread("/home/workspace/jogrady/ML4TB/data/abdelaal/abdelaal_samples.csv", sep = "\t")
abdelaal_labels <- abdelaal_labels[seq(1,48,2),]
rownames(abdelaal_labels) <- unique(abdelaal_labels$Animal_Code)
rownames(abdelaal) <- ensemble$gene_name
abdelaal_labels <- abdelaal_labels %>% select(3,4)
colnames(abdelaal_labels) <- c("Sample", "Condition")
abdelaal_labels$Infection_administration <- "Experimental"
abdelaal_labels$Sex <- "M"
abdelaal_labels$Location <- "US"
abdelaal_labels$Tissue <- "PB"
abdelaal_labels$Animal_Code <- abdelaal_labels$Sample
abdelaal_labels$Animal_Code <- gsub("_8", "", abdelaal_labels$Animal_Code)
abdelaal_labels$Animal_Code <- gsub("_20", "", abdelaal_labels$Animal_Code)
abdelaal_labels$Age <- rep(c(14,11), 12)
abdelaal_labels$Study <- "Abdelaal"
rownames(abdelaal_labels)


labels



wiarda_labels <- wiarda_labels %>% select(colnames(labels))
wiarda_labels <- data.frame(wiarda_labels)
rownames(wiarda_labels) <- wiarda_labels$Sample

abdelaal_labels <- abdelaal_labels %>% select(colnames(labels))
abdelaal_labels <- data.frame(abdelaal_labels)
rownames(abdelaal_labels) <- abdelaal_labels$Sample


kirsten_labels <- kirsten_labels[,colnames(labels)]
kirsten_pbl_labels <- kirsten_pbl_labels[,colnames(labels)]
abdelaal_labels <- abdelaal_labels[,colnames(labels)]

data_raw <- data_raw %>% select(-1)
data_raw <- as.matrix(data_raw)
rownames(data_raw) <- ensemble$gene_name

abdelaal_labels
wiarda_labels

head(abdelaal_labels)

data_all <- cbind(data_raw, wiarda, Kirsten, Kirsten_pbl, abdelaal)
metadata_all <- rbind(labels, wiarda_labels, kirsten_labels, kirsten_pbl_labels, abdelaal_labels)

dim(metadata_all)




#####################################################
#####################################################
####
#     Test train split - make sure to ensure time series replicates are not in test and train
####
####################################################
####################################################
####################################################


# Step 1: Identify unique Animal_Codes
unique_animals <- unique(metadata_all$Animal_Code)

unique_animals
floor
# Step 2: Randomly partition Animal_Codes into training and testing
train_animals <- sample(unique_animals, size = floor(0.7 * length(unique_animals)))
length(train_animals)
test_animals <- setdiff(unique_animals, train_animals)

# Step 3: Assign samples based on Animal_Code
train_set <- metadata_all %>% filter(Animal_Code %in% train_animals)
test_set <- metadata_all %>% filter(Animal_Code %in% test_animals)
dim(metadata_all)
dim(train_set)
dim(test_set)
179 + 75
173 + 81
table(train_set$Animal_Code %in% test_set$Animal_Code)



train_counts <- data_all %>% select(train_set$Sample)
test_counts <- data_all %>% select(test_set$Sample)

head(train_counts)

####################################################
####################################################
####################################################
####
#     Batch Correction for test and train
####
####################################################
####################################################
####################################################


# First generate PCA of batches
ddsTrain <- DESeqDataSetFromMatrix(countData = train_counts, 
                                   colData = train_set, 
                                   design = ~ 1)

# Test
ddsTest <- DESeqDataSetFromMatrix(countData = test_counts, 
                                   colData = test_set, 
                                   design = ~ 1)




vsdTrain <- vst(ddsTrain)
vsdTest <- vst(ddsTest)

library(RColorBrewer)



pcaData <- plotPCA(vsdTrain, ntop=1500, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Train_PCA <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  xlim(c(-25,50)) + ylim(c(-45, 45))+ theme_bw() + scale_colour_brewer(palette = "Dark2")


Train_PCA



pcaData <- plotPCA(vsdTest, ntop=1500, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Test_PCA <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_brewer(palette = "Dark2") + scale_y_reverse(limits = c(45, -45))
Test_PCA



plot_grid(Train_PCA, Test_PCA)

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_batch/PCA_train_test_batch_unadjusted.pdf", width = 12, height = 12, dpi = 600)

# Perform the batch correction

batches= sapply(as.character(train_set$Study), switch, "1_OGrady" = 1, "2_OGrady" = 2, "Wiarda" = 3, "Abdelaal" = 4,  "Mcloughlin" = 5, "Mcloughlin_pbl" = 6, USE.NAMES = FALSE)
groups= sapply(as.character(train_set$Condition), switch, "Control" = 1, "Infected" = 2, USE.NAMES = FALSE)

train_adjusted_counts <- sva::ComBat_seq(as.matrix(train_counts), batch = train_set$Study, group = train_set$Condition)


batches= sapply(as.character(test_set$Study), switch, "1_OGrady" = 1, "2_OGrady" = 2, "Wiarda" = 3, "Abdelaal" = 4,  "Mcloughlin" = 5, "Mcloughlin_pbl" = 6, USE.NAMES = FALSE)
groups= sapply(as.character(test_set$Condition), switch, "Control" = 1, "Infected" = 2, USE.NAMES = FALSE)
test_adjusted_counts <- sva::ComBat_seq(as.matrix(test_counts), batch = batches, group = groups)



rownames(train_adjusted_counts) <- ensemble$gene_name
rownames(test_adjusted_counts) <- ensemble$gene_name
rownames(train_adjusted_counts) <- gsub(" ", "", rownames(train_adjusted_counts))
rownames(test_adjusted_counts) <- gsub(" ", "", rownames(test_adjusted_counts))



train_set$Age <- cut(as.numeric(train_set$Age), breaks=c(0,6,12,22, 50), labels=c("1-6","6-12","12-22", "20+"))



head(train_adjusted_counts[,1])

# Regenerate the DESeq2 object
ddsTrain_adjusted  <- DESeqDataSetFromMatrix(countData = train_adjusted_counts, 
                                  colData = train_set, 
                                  design = ~ Age + Infection_administration + Sex + Location + Tissue + Condition)


# Regenerate the DESeq2 object

ddsTest_adjusted  <- DESeqDataSetFromMatrix(countData = test_adjusted_counts, 
                                             colData = test_set, 
                                             design = ~ 1) # Note variables arent really important as not doing DE on them, just easy to copy and paste



vsdTrain_adjusted <- vst(ddsTrain_adjusted)
vsdTest_adjusted <- vst(ddsTest_adjusted)





pcaData <- plotPCA(vsdTrain_adjusted, ntop=1500, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Train_PCA_adjusted <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_brewer(palette = "Dark2") + ylim(-25,25)


Train_PCA_adjusted



pcaData <- plotPCA(vsdTest_adjusted, ntop=1500, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Test_PCA_adjusted <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_brewer(palette = "Dark2") + scale_y_reverse(limits = c(25, -25))

plot_grid(Train_PCA_adjusted, Test_PCA_adjusted)




## GGSAVE to file
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_batch/PCA_train_test_batch_adjusted.pdf", width = 12, height = 12, dpi = 600)



# Now colour by condition
pcaData <- plotPCA(vsdTrain_adjusted, ntop=1500, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Train_PCA_adjusted_condition <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_manual(values = c("#2166ac", "#b2182b"))+ ylim(-25,25)


Train_PCA_adjusted_condition



pcaData <- plotPCA(vsdTest_adjusted, ntop=1500, intgroup=c("Condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Test_PCA_adjusted_condition <- ggplot(pcaData, aes(PC1, PC2, color=Condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_manual(values = c("#2166ac", "#b2182b")) + scale_y_reverse(limits = c(25, -25))

plot_grid(Train_PCA_adjusted_condition, Test_PCA_adjusted_condition)



## GGSAVE to file
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_batch/PCA_train_test_Condition_adjusted.pdf", width = 12, height = 12, dpi = 600)





####################################################
####################################################
####################################################
####
#     Differential Expression analysis on training set
####
####################################################
####################################################
####################################################



ddsTrain_adjusted <- DESeq(ddsTrain_adjusted)
ddsTest_adjusted <- DESeq(ddsTest_adjusted)
results(ddsTrain_adjusted, alpha = 0.05, name="Condition_Infected_vs_Control")
res <- lfcShrink(ddsTrain_adjusted, coef="Condition_Infected_vs_Control", type="apeglm")


#out of 22540 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 786, 3.5%
#LFC < 0 (down)     : 768, 3.4%
#outliers [1]       : 4, 0.018%
#low counts [2]     : 4759, 21%
#(mean count < 1)

DE_genes_train <- as.data.frame(res)


DE_genes_train_expressed <- DE_genes_train %>% filter(baseMean >= 100 & padj < 0.05)

DE_genes_train_expressed_to_subset <- rownames(DE_genes_train_expressed)


####################################################
####################################################
####################################################
####
#     Now perform variance stabilising transformation (VST) on training and extract the genes of interest
#     Then use functions to VST the testing data and then extract genes of interest.
####
####################################################
####################################################
####################################################


# variance stabilised transformation on the training data
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(ddsTrain_adjusted, blind = FALSE)


# Apply the dispersion function on dds test
dispersionFunction(ddsTest_adjusted) <- dispersionFunction(ddsTrain_adjusted)
vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(ddsTest_adjusted, blind = FALSE) # Now perform the normalisation

train_normalised_filtered_counts <- assay(vstNormalizedExpressionDataForTrain)[DE_genes_train_expressed_to_subset,]

test_normalised_filtered_counts <- assay(vstNormalizedExpressionDataForTest)[DE_genes_train_expressed_to_subset,]


# At this stage could add in information about genotype / admixture / breed etc into these models
train_normalised_filtered_counts <- t(train_normalised_filtered_counts)

train_set$Condition <- factor(train_set$Condition, levels = c("Control", "Infected"), labels = c(0,1))
test_set$Condition <- factor(test_set$Condition, levels = c("Control", "Infected"))
train_normalised_filtered_counts <- cbind(train_normalised_filtered_counts, train_set$Condition)

colnames(train_normalised_filtered_counts)[length(colnames(train_normalised_filtered_counts))] <- "Condition"
train_normalised_filtered_counts[,length(colnames(train_normalised_filtered_counts))]



# Also transform the the sets
test_normalised_filtered_counts <- t(test_normalised_filtered_counts)



####################################################
####################################################
####################################################
####
#     For K fold cross validation, need to ensure no leaking
####
####################################################
####################################################
####################################################


# Need to generate folds by sample ID
num_folds = 5

unique_train_animals <- unique(train_set$Animal_Code)
shuffled_ids <- sample(unique_train_animals)
fold_assignments <- cut(seq_along(shuffled_ids), breaks = num_folds, labels = FALSE)
# Map fold assignments back to the original vector
train_set$Folds <- sapply(train_set$Animal_Code, function(id) fold_assignments[shuffled_ids == id])

# Compatability with caret
id_to_fold <- setNames(fold_assignments, shuffled_ids)

# Create caret-compatible folds
custom_folds <- lapply(1:num_folds, function(fold) {
  which(!(train_set$Animal_Code %in% names(id_to_fold[id_to_fold == fold])))
})
custom_folds

table(train_set$Folds, train_set$Condition)
# Fit into object
train_control <- trainControl(
  method = "cv",
  index = custom_folds,
  allowParallel = TRUE,
  classProbs = TRUE,
  number = 5,
  savePredictions = "final",
  verboseIter = TRUE,
  summaryFunction = twoClassSummary
)

train_normalised_filtered_counts <- as.data.frame(train_normalised_filtered_counts)
train_normalised_filtered_counts$Condition <- factor(train_normalised_filtered_counts$Condition, labels = c("Control","Infected"))

train_control

####################################################
####################################################
####################################################
####
#     Model development --> Simple - logistic regression, penalised (ridge, lasso, elastic net)
####
####################################################
####################################################
####################################################


GLM_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                   method = "glm",
                   trControl = train_control,
                   metric = "ROC")


GLMRIDGE_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
      method = "glmnet",
      tuneGrid = expand.grid(alpha = 0, lambda = 10 ^ seq(5, -3, length = 1000)),
      trControl = train_control,
      metric = "ROC",
      verbose = TRUE)


GLMRIDGE_model$results[425,]


GLMLASSO_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                        method = "glmnet",
                        tuneGrid = expand.grid(alpha = 1, lambda = 10 ^ seq(5, -3, length = 1000)),
                        trControl = train_control,
                        metric = "ROC",
                        verbose = TRUE)

GLMLASSO_model


GLMENET_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                        method = "glmnet",
                        tuneGrid = expand.grid(alpha = seq(0.05,0.95,0.05), lambda = 10 ^ seq(5, -3, length = 1000)),
                        trControl = train_control,
                        metric = "ROC",
                        verbose = TRUE)


####################################################
####################################################
####################################################


# ROC curves for logistic regression models
# And prediction on test data sets

####################################################
####################################################
####################################################

  
GLM_train_ROC <- CV_ROC_plot(GLM_model)

GLM_train_ROC
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLM_CV_ROC.pdf", width = 12, height = 12, dpi = 600)

GLMRIDGE_train_ROC <- CV_ROC_plot(GLMRIDGE_model)
GLMRIDGE_train_ROC
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLMRIDGE_CV_ROC.pdf", width = 12, height = 12, dpi = 600)

GLMLASSO_train_ROC <- CV_ROC_plot(GLMLASSO_model)
GLMLASSO_train_ROC
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLMLASSO_CV_ROC.pdf", width = 12, height = 12, dpi = 600)

GLMENET_train_ROC <- CV_ROC_plot(GLMENET_model)
GLMENET_train_ROC
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLMENET_CV_ROC.pdf", width = 12, height = 12, dpi = 600)





Test_AUC <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(Test_AUC) <- c("Model", "AUC")



GLM_predict <- predict(GLM_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLM")
GLM_predict_ROC <- ROC_predict_individual(GLM_predict, test_set$Condition)



GLMRIDGE_predict <- predict(GLMRIDGE_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMRIDGE")
Test_AUC <- rbind(Test_AUC, c("GLMRIDGE", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMRIDGE_predict[[2]])$auc))
GLMRIDGE_ROC <- ROC_predict_individual(GLMRIDGE_predict, test_set$Condition)


GLMLASSO_predict <- predict(GLMLASSO_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMLASSO")
Test_AUC <- rbind(Test_AUC, c("GLMRIDGE", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMLASSO_predict[[2]])$auc))
GLMLASSO_ROC <- ROC_predict_individual(GLMLASSO_predict, test_set$Condition)


GLMENET_predict <- predict(GLMENET_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMENET")
Test_AUC <- rbind(Test_AUC, c("GLMRIDGE", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMENET_predict[[2]])$auc))
GLMENET_ROC <- ROC_predict_individual(GLMENET_predict, test_set$Condition)



plot_grid(GLM_predict_ROC, GLMRIDGE_ROC, GLMLASSO_ROC, GLMENET_ROC)


#test <- rbind(GLM_predict, GLMRIDGE_predict, GLMLASSO_predict, GLMENET_predict)
#real_test <- ROC_test_combined(test, test_set$Condition)






twoClassSummary(test_set, lev = levels(GLM_ridge_predict$obs))



# save
#ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLM_CV_ROC.pdf", width = 12, height = 12, dpi = 600)


# How to get coefficients for model and filter for those that are non-zero
#test = exp(coef(GLMENET_model$finalModel, GLMENET_model$bestTune$lambda)[,1])
#test







library(doParallel) 
registerDoParallel(25) 
FOREST <- list()


####################################################
####################################################
####################################################


# RF model tuning with loops and parallelisation

####################################################
####################################################
####################################################

n_features <- length(setdiff(names(train_normalised_filtered_counts), "Condition"))
n_features

sqrt(n_features)
rf_grid = expand.grid(
  mtry=c(10, 25, 50, 75, 100, 150, 200, 250, 300, 350, 400, 450, 500),
  min.node.size = c(1:10),
  splitrule = c("gini", "extratrees", "hellinger")
)
colnames(rf_grid)

FOREST <- foreach(num.trees = c(51,101, 251, 501, 751, 1001, 1501), .packages = c("caret", "ranger"), .combine = 'c') %dopar% {
  set.seed(42)
  print(num.trees)
  # Fit into object
  fit <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                method = "ranger",
                trControl = train_control,
                metric = "ROC",
                tuneGrid = rf_grid,
                importance = "impurity",
                num.trees = num.trees,
                verbose = TRUE)
  list(fit)
}

names(FOREST) <- c(51,101, 251, 501, 751, 1001, 1501)



plot_grid(plot(FOREST[["51"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["101"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["251"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["501"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["751"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["1001"]], metric = "ROC", y = c(0.82,0.95)),plot(FOREST[["1501"]], metric = "ROC", y = c(0.82,0.95)))


# Compare models
resamples_rf <- resamples(FOREST)

summary(resamples_rf)


# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF <- FOREST[[best_rf_model_name]]



####################################################
####################################################
####################################################


# MLP model tuning

####################################################
####################################################
####################################################

library(pROC)

MLP_model <- train(Condition ~ ., data = train_normalised_filtered_counts[,500:length(colnames(train_normalised_filtered_counts))], 
                        method = "nnet",
                        trControl = train_control,
                        metric = "ROC",
                        maxit=1000,
                        na.action = "na.omit")

MLP_model
dim(train_normalised_filtered_counts)

# Step 3: Model Training
library(nnet)

# Train the neural network for regression
model <- nnet(Condition ~ colnames(train_normalised_filtered_counts[,-1]), data = train_data, size = 5, maxit = 1000, linout = TRUE)

table(train_normalised_filtered_counts$Condition)
dim(test_normalised_filtered_counts)




