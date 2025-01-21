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

Train_counts <- read.table("/home/workspace/jogrady/ML4TB/data/human/GSE255071_counts.train.tsv")
Test_counts <- read.table("/home/workspace/jogrady/ML4TB/data/human/GSE255073_counts.test.tsv")
Validation_counts <- read.table("/home/workspace/jogrady/ML4TB/data/human/GSE255074_counts.val.tsv")



metadata_train <- read.csv("/home/workspace/jogrady/ML4TB/data/human/metadata.train.tsv", sep = "\t") 
metadata_test <- read.csv("/home/workspace/jogrady/ML4TB/data/human/metadata.test.tsv", sep = "\t")
metadata_val <- read.csv("/home/workspace/jogrady/ML4TB/data/human/metadata.val.tsv", sep = "\t")




#Concatenate Metadata
metadata_file_modif = rbind(metadata_train, metadata_test, metadata_val)
#Recode site to match the country it was collected from across cohorts
metadata_file_modif <- metadata_file_modif %>%
  mutate(site = recode(site, 
                       "Kisenyi" = "Uganda", 
                       "Mulago" = "Uganda"))
metadata_file_modif <- metadata_file_modif %>%
  mutate(site = ifelse(is.na(site) | site == "", "Uganda", site))

#Assign variables for differential abundance analysis
SAMPLE_ID_VAR = "sample_id"
COMP_VAR = "tb"
GROUPS = c("positive","negative")

gene_name_key = "/home/workspace/jogrady/ML4TB/data/human/gencode.biotype.name.key.tsv"

DESIGN = paste0("~",COMP_VAR)
#Read in the count matrices
counts_train = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255071_counts.train.tsv",row.names=1)
counts_test = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255073_counts.test.tsv", row.names = 1)
counts_val = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255074_counts.val.tsv",row.names=1)

##------------------------------------
# Add Gene metadata
annotation = fread(file=gene_name_key)
annotation <- annotation[match(rownames(dds), annotation$gene_id),]



# train
samples <- metadata_train %>% 
  filter(.data[[COMP_VAR]] %in% GROUPS)
rownames(samples) <- samples$sample_id
counts_train <- counts_train[, rownames(samples)]
dds_train <- DESeqDataSetFromMatrix(counts_train,
                              colData = samples,
                              design = ~ 1)

# test
samples <- metadata_test %>% 
  filter(.data[[COMP_VAR]] %in% GROUPS)
rownames(samples) <- samples$sample_id
counts_test <- counts_test[, rownames(samples)]
dds_test <- DESeqDataSetFromMatrix(counts_test,
                                    colData = samples,
                                    design = ~ 1)

#val
samples <- metadata_val %>% 
filter(.data[[COMP_VAR]] %in% GROUPS)
rownames(samples) <- samples$sample_id
counts_val <- counts_val[, rownames(samples)]
dds_val <- DESeqDataSetFromMatrix(counts_val,
                                   colData = samples,
                                   design = ~ 1)

dds_train = DESeq(dds_train)
dds_test = DESeq(dds_test)
dds_val = DESeq(dds_val)


saveRDS(dds_train, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_train.rds")
saveRDS(dds_test, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_test.rds")
saveRDS(dds_val, file = "/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_val.rds")

## -----------------
# Train data normalising
# Training vst and then applying to test and val

vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(dds_train, blind = FALSE)

meanSdPlot(assay(vstNormalizedExpressionDataForTrain))

# Apply the dispersion function on dds test
dispersionFunction(dds_test) <- dispersionFunction(dds_train) # This is related to point above from Mike Love
dispersionFunction(dds_val) <- dispersionFunction(dds_train) # This is related to point above from Mike Love


vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(dds_test, blind = FALSE) # Now perform the normalisation
vstNormalizedExpressionDataForVal <- varianceStabilizingTransformation(dds_val, blind = FALSE) # Now perform the normalisation

# extract the counts form transformed object
train_counts_normalised  <- assay(vstNormalizedExpressionDataForTrain)
test_counts_normalised_train <- assay(vstNormalizedExpressionDataForTest)
val_counts_normalised_train <- assay(vstNormalizedExpressionDataForVal)


all(annotation$gene_id == rownames(test_counts_normalised_train))

rownames(train_counts_normalised) <- annotation$gene_name
rownames(test_counts_normalised_train) <- annotation$gene_name
rownames(val_counts_normalised_train) <- annotation$gene_name


write.table(train_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/human_train_vst_normalised_data.txt", sep = "\t", quote = FALSE)
write.table(val_counts_normalised_train, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_val_vst_normalised_human_train.txt", quote = FALSE, sep = "\t")
write.table(test_counts_normalised_train, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_test_vst_normalised_human_train.txt", quote = FALSE, sep = "\t")

## -----------------
# Now test data
dds_train = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_train.rds")
dds_test = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_test.rds")
dds_val = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_val.rds")


vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(dds_test, blind = FALSE)
meanSdPlot(assay(vstNormalizedExpressionDataForTest))


dispersionFunction(dds_train) <- dispersionFunction(dds_test) # This is related to point above from Mike Love
dispersionFunction(dds_val) <- dispersionFunction(dds_test) # This is related to point above from Mike Love
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(dds_train, blind = FALSE) # Now perform the normalisation
vstNormalizedExpressionDataForVal <- varianceStabilizingTransformation(dds_val, blind = FALSE) # Now perform the normalisation

# extract the counts form transformed object
test_counts_normalised  <- assay(vstNormalizedExpressionDataForTest)
train_counts_normalised_test <- assay(vstNormalizedExpressionDataForTrain)
val_counts_normalised_test <- assay(vstNormalizedExpressionDataForVal)

all(annotation$gene_id == rownames(test_counts_normalised))

rownames(test_counts_normalised) <- annotation$gene_name
rownames(train_counts_normalised_test) <- annotation$gene_name
rownames(val_counts_normalised_test) <- annotation$gene_name

write.table(test_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/human_test_vst_normalised_data.txt", sep = "\t", quote = FALSE)
write.table(val_counts_normalised_test, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_val_vst_normalised_human_test.txt", quote = FALSE, sep = "\t")
write.table(train_counts_normalised_test, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_train_vst_normalised_human_test.txt", quote = FALSE, sep = "\t")

## -----------------
# Now Val
dds_train = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_train.rds")
dds_test = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_test.rds")
dds_val = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_val.rds")




vstNormalizedExpressionDataForVal <- varianceStabilizingTransformation(dds_val, blind = FALSE)
meanSdPlot(assay(vstNormalizedExpressionDataForVal))
dispersionFunction(dds_train) <- dispersionFunction(dds_val) # This is related to point above from Mike Love
dispersionFunction(dds_test) <- dispersionFunction(dds_val) # This is related to point above from Mike Love
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(dds_train, blind = FALSE) # Now perform the normalisation
vstNormalizedExpressionDataForTest <- varianceStabilizingTransformation(dds_test, blind = FALSE) # Now perform the normalisation

# extract the counts form transformed object
val_counts_normalised  <- assay(vstNormalizedExpressionDataForVal)
train_counts_normalised_val <- assay(vstNormalizedExpressionDataForTrain)
test_counts_normalised_val <- assay(vstNormalizedExpressionDataForTest)


all(annotation$gene_id == rownames(val_counts_normalised))

rownames(val_counts_normalised) <- annotation$gene_name
rownames(train_counts_normalised_val) <- annotation$gene_name
rownames(test_counts_normalised_val) <- annotation$gene_name

write.table(val_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/human_val_vst_normalised_data.txt", sep = "\t", quote = FALSE)
write.table(train_counts_normalised_val, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_train_vst_normalised_human_val.txt", quote = FALSE, sep = "\t")
write.table(test_counts_normalised_val, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_combined/human_test_vst_normalised_human_val.txt", quote = FALSE, sep = "\t")




#Read in the count matrices
counts_train = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255071_counts.train.tsv",row.names=1)
counts_test = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255073_counts.test.tsv", row.names = 1)
counts_val = read.delim("/home/workspace/jogrady/ML4TB/data/human/GSE255074_counts.val.tsv",row.names=1)


### ALL
#Concatenate count matrices



counts <- cbind(counts_train, counts_test, counts_val)
samples <- metadata_file_modif %>% 
  filter(.data[[COMP_VAR]] %in% GROUPS)
rownames(samples) <- samples$sample_id
counts <- counts[, rownames(samples)]


DESIGN
##------------------------------------
# Contstruct DESeq Data Set
dds <- DESeqDataSetFromMatrix(counts,
                              colData = samples,
                              design = formula(DESIGN))


##------------------------------------
# Add Gene metadata
annotation = fread(file=gene_name_key)
annotation <- annotation[match(rownames(dds), annotation$gene_id),]
all(rownames(dds) == annotation$ftcount_id)
mcols(dds) <- cbind(mcols(dds), annotation)

dds
##------------------------------------
# Re-factor
dds[[COMP_VAR]] <- factor(dds[[COMP_VAR]], levels = GROUPS)

dds <- DESeq(dds)
vsd <- vst(dds)
colnames(metadata_test)

pcaData <- plotPCA(vsd, intgroup=c("tb", "cohort"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_all_tb <- ggplot(pcaData, aes(PC1, PC2, color=tb)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()

p_all_cohort <-  ggplot(pcaData, aes(PC1, PC2, color=cohort)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() + theme_bw()

library(cowplot)
plot_grid(p_all_tb, p_all_cohort, ncol = 1)

ggsave("/home/workspace/jogrady/ML4TB/work/human/PCA_all_samples_status_cohort.pdf", width = 12, height = 8, dpi = 600)




# Read in the data for each
dds_train = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_train.rds")
dds_test = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_test.rds")
dds_val = readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/individual_deseq2_objects/human_val.rds")

vsd_train <- vst(dds_train)
vsd_test <- vst(dds_test)
vsd_val <- vst(dds_val)

# Train
pcaData <- plotPCA(vsd_train, intgroup=c("tb", "cohort"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_train_tb <- ggplot(pcaData, aes(PC1, PC2, color=tb, shape = cohort)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ylim = c(-30,30), xlim = c(-70,70))+ theme_bw()
p_train_tb


pcaData <- plotPCA(vsd_test, intgroup=c("tb", "cohort"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_test_tb <- ggplot(pcaData, aes(PC1, PC2, color=tb, shape = cohort)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ylim = c(-30,30), xlim = c(-70,70)) + theme_bw()
p_test_tb


pcaData <- plotPCA(vsd_val, intgroup=c("tb", "cohort"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
p_val_tb <- ggplot(pcaData, aes(PC1, PC2, color=tb, shape = cohort)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(ylim = c(-40,40), xlim = c(-70,70)) + theme_bw()
p_val_tb

plot_grid(p_train_tb, p_test_tb, p_val_tb, ncol = 1)
ggsave("/home/workspace/jogrady/ML4TB/work/human/PCA_each_set_status_cohort.pdf", width = 12, height = 8, dpi = 600)


plotPCA(vsd, intgroup = c("cohort"), shape = c("hiv"))


View(metadata_file_modif)
res <- results(dds,alpha=0.05,contrast = c(COMP_VAR,GROUPS[1], GROUPS[2]))
resgene_name <- mcols(dds)gene_name
resgene_type <- mcols(dds)gene_type

summary(res)
