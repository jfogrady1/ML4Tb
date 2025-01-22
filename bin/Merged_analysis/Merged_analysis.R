# Merging bovine samples

library(tidyverse)
library(ggplot2)
library(caret)
library(sva)
library(DESeq2)
library(pROC)
library("ochRe")


set.seed(42) # For reproducibility

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

# Step 2: Randomly partition Animal_Codes into training and testing
train_animals <- sample(unique_animals, size = floor(0.7 * length(unique_animals)))
test_animals <- setdiff(unique_animals, train_animals)

# Step 3: Assign samples based on Animal_Code
train_set <- metadata_all %>% filter(Animal_Code %in% train_animals)
test_set <- metadata_all %>% filter(Animal_Code %in% test_animals)


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

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_train_test_batch_unadjusted.pdf", width = 12, height = 12, dpi = 600)

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
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_train_test_batch_adjusted.pdf", width = 12, height = 12, dpi = 600)



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
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/PCA_train_test_Condition_adjusted.pdf", width = 12, height = 12, dpi = 600)





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






