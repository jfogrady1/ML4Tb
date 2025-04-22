# Pairwise DE analysis

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
source("~/ML4TB/bin/Merged_analysis/Funcitons.R")
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
ensemble$gene_name <- gsub(" ", "", ensemble$gene_name)
all(tested_genes$Geneid == ensemble$gene_id)



rownames(data_raw)
data_raw <- data_raw[,-1]
data_raw <- as.matrix(data_raw)
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

ogrady_cov <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt") %>% select(Sample, Batch)
labels <- left_join(labels, ogrady_cov)
labels$Batch
labels$Age
# Read in PCA data

# Ogrady
ogrady_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_filtered_ALL_Pruned.eigenvec") %>% select(1,3,4)

genotype_PCs <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec")
genotype_PCs <- genotype_PCs %>% select(-1)
genotype_PCs$V2 <- gsub("_.*", "", genotype_PCs$V2)
rownames(genotype_PCs) <- genotype_PCs$V2
genotype_PCs <- genotype_PCs %>% select(-1)
genotype_PCs <- genotype_PCs %>% select(1,2)
colnames(genotype_PCs) <- paste0("Genotype_PC", 1:2)

head(genotype_PCs)


colnames(ogrady_eigen_vec) <- c("Sample", "PC1", "PC2")
labels <- left_join(labels, ogrady_eigen_vec)
rownames(labels) <- labels$Sample
labels$Age <- as.numeric(labels$Age)
labels$Age <- scale(labels$Age, center = TRUE)
labels$Batch <- factor(labels$Batch, levels = c("1", "2"))
labels$Batch
labels

labels <- cbind(labels, genotype_PCs)
head(labels)
# DESEQ2 analysis
ddsMat_ogrady <- DESeqDataSetFromMatrix(countData = data_raw,
                                 colData = labels,
                                 design = ~ Batch + Genotype_PC1 + Genotype_PC2 + Age + Condition)





ddsMat_ogrady$Age
# WIARDA
wiarda = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt", sep = "\t") %>% select(-1)
wiarda <- as.matrix(wiarda)
wiarda_labels = fread("/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv", sep = "\t")
rownames(wiarda_labels) <- wiarda_labels$Animal_Code
rownames(wiarda) <- ensemble$gene_name
wiarda
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


wiarda_labels$Time = gsub("^\\d{3}_.+?_W", "", wiarda_labels$Sample)
wiarda_labels[1:5,2] = "Infected"
wiarda_labels[11:13, 2] = "Infected"

wiarda_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_filtered_ALL_Pruned.eigenvec") %>% select(1,3,4)
colnames(wiarda_eigen_vec) <- c("Sample", "PC1", "PC2")
wiarda_labels$Animal_Code <- as.numeric(wiarda_labels$Animal_Code)
wiarda_labels <- left_join(wiarda_labels, wiarda_eigen_vec, by = c("Animal_Code" = "Sample"))


wiarda_labels
wiarda_labels$Time = factor(wiarda_labels$Time, levels = c("0", "4", "10"))
wiarda_labels$Condition = factor(wiarda_labels$Condition, levels = c("Control", "Infected"))
wiarda_labels$Group = paste0(wiarda_labels$Condition, "_", wiarda_labels$Time)                            
table(wiarda_labels$Group)
wiarda_labels$Group <- factor(wiarda_labels$Group, levels = c("Control_0",  "Infected_0",   "Control_4",  "Infected_4",  "Control_10", "Infected_10"))
dds_wiarda = DESeqDataSetFromMatrix(countData = wiarda,
                                     colData = wiarda_labels,
                                     design = ~ PC1 + PC2 +Group)

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

# Ogrady
kirsten_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_filtered_ALL_Pruned.eigenvec") %>% select(1,3,4)
colnames(kirsten_eigen_vec) <- c("Sample", "PC1", "PC2")
kirsten_eigen_vec
kirsten_labels <- left_join(kirsten_labels, kirsten_eigen_vec, by = c("Animal_Code" = "Sample"))

kirsten_labels$Time <- gsub("^A\\d{4}_W", "", kirsten_labels$Sample)
kirsten_labels$Time <- factor(kirsten_labels$Time, levels = c("-1", "1", "2", "6", "10", "12"))
dds_kirsten = DESeqDataSetFromMatrix(countData = Kirsten,
                                    colData = kirsten_labels,
                                    design = ~ PC1 + PC2 + Time)



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
kirsten_pbl_labels$Animal_Code <- gsub("_.*", "", kirsten_pbl_labels$Animal_Code)
rownames(kirsten_pbl_labels)

kirsten_pbl_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_filtered_ALL_Pruned.eigenvec") %>% select(1,3,4)
colnames(kirsten_pbl_eigen_vec) <- c("Sample", "PC1", "PC2")
kirsten_pbl_eigen_vec
head(kirsten_pbl_labels)
kirsten_pbl_labels <- left_join(kirsten_pbl_labels, kirsten_pbl_eigen_vec, by = c("Animal_Code" = "Sample"))


dds_kirsten_pbl = DESeqDataSetFromMatrix(countData = Kirsten_pbl,
                                     colData = kirsten_pbl_labels,
                                     design = ~ PC1 + PC2 + Condition)


# Abdelaal
abdelaal = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/abdelaal/Quantification/abdelaal_count_matrix_clean.txt", sep = "\t") %>% select(-1)
abdelaal <- as.matrix(abdelaal)
rownames(abdelaal) <- ensemble$gene_name
rownames(abdelaal)
abdelaal_labels = fread("/home/workspace/jogrady/ML4TB/data/abdelaal/abdelaal_samples.csv", sep = "\t")
abdelaal_labels
abdelaal_labels <- abdelaal_labels[seq(1,48,2),]
rownames(abdelaal_labels) <- unique(abdelaal_labels$Animal_Code)
row.names(abdelaal) <- ensemble$gene_name
abdelaal_labels$Infection_administration <- "Experimental"
abdelaal_labels$Sex <- "M"
abdelaal_labels$Location <- "US"
abdelaal_labels$Tissue <- "PB"
abdelaal_labels$Animal_Code <- gsub("_8", "", abdelaal_labels$Animal_Code)
abdelaal_labels$Animal_Code <- gsub("_20", "", abdelaal_labels$Animal_Code)
abdelaal_labels$Age <- rep(c(14,11), 12)
abdelaal_labels$Study <- "Abdelaal"
abdelaal_labels

abdelaal_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/abdelaal_filtered_ALL_Pruned.eigenvec") %>% select(1,2,3,4) %>% mutate(V1 = paste0(V1, "_", V2)) %>% select(1,3,4)

colnames(abdelaal_eigen_vec) <- c("Sample", "PC1", "PC2")


abdelaal_labels <- left_join(abdelaal_labels, abdelaal_eigen_vec, by = c("Animal_Code" = "Sample"))



abdelaal_labels
abdelaal_labels$Time <- gsub("W", "", abdelaal_labels$Week)
abdelaal_labels$Time
abdelaal_labels$Group = paste0(abdelaal_labels$Status, "_", abdelaal_labels$Time)
abdelaal_labels$Group <- factor(abdelaal_labels$Group, levels = c("Control_8", "Infected_8", "Control_20", "Infected_20"))
dds_abdelaal = DESeqDataSetFromMatrix(countData = abdelaal,
                                     colData = abdelaal_labels,
                                     design = ~ PC1 + PC2 + Group)

abdelaal_labels$Group <- factor(abdelaal_labels$Group, levels = c("Control_20", "Infected_8", "Control_8", "Infected_20"))
dds_abdelaal_20 = DESeqDataSetFromMatrix(countData = abdelaal,
                                      colData = abdelaal_labels,
                                      design = ~ PC1 + PC2 + Group)



#########################################
# DE Analysis
#########################################

#keep_ogrady
# Filter for low gene counts
keep_ogrady <- rowSums(counts(ddsMat_ogrady) >= 6) >= (ncol(ddsMat_ogrady) *.2) # remove low count genes
keep_ogrady <- keep_ogrady[keep_ogrady == TRUE]

keep_wiarda <- rowSums(counts(dds_wiarda) >= 6) >= (ncol(dds_wiarda) * .2) # remove low count genes
keep_wiarda <- keep_wiarda[keep_wiarda == TRUE]


keep_kirsten <- rowSums(counts(dds_kirsten) >= 6) >= (ncol(dds_kirsten) * .2) # remove low count genes
keep_kirsten <- keep_kirsten[keep_kirsten == TRUE]

keep_kirsten_pbl <- rowSums(counts(dds_kirsten_pbl) >= 6) >= (ncol(dds_kirsten_pbl) *.2) # remove low count genes
keep_kirsten_pbl <- keep_kirsten_pbl[keep_kirsten_pbl == TRUE]

keep_abdelaal <- rowSums(counts(dds_abdelaal) >= 6) >= (ncol(dds_abdelaal) * .2) # remove low count genes
keep_abdelaal <- keep_abdelaal[keep_abdelaal == TRUE]


length(keep_ogrady)
length(keep_wiarda)
expressed_genes <- intersect(names(keep_ogrady), names(keep_wiarda)) 
expressed_genes <- intersect(expressed_genes, names(keep_kirsten))
expressed_genes <- intersect(expressed_genes, names(keep_kirsten_pbl))

expressed_genes_abdelaal <- intersect(expressed_genes, names(keep_abdelaal))


expressed_genes
# filter
ddsMat_ogrady <- ddsMat_ogrady[expressed_genes,] 
dds_kirsten <- dds_kirsten[expressed_genes,]
dds_abdelaal <- dds_abdelaal[expressed_genes_abdelaal,]
dds_abdelaal_20 <- dds_abdelaal_20[expressed_genes_abdelaal,]
dds_kirsten_pbl <- dds_kirsten_pbl[expressed_genes,]
dds_wiarda <- dds_wiarda[expressed_genes,]



### O'Grady
# Pairwise
ddsMat_ogrady
ddsMat_ogrady <- DESeq(ddsMat_ogrady)
dds_kirsten_pbl <- DESeq(dds_kirsten_pbl)
dds_kirsten <- DESeq(dds_kirsten)

# Perform the LF shrink method
res_ogrady <- lfcShrink(ddsMat_ogrady, coef="Condition_Infected_vs_Control", type="apeglm")

res_kirsten_pbl <- lfcShrink(dds_kirsten_pbl, coef = "Condition_Infected_vs_Control", type = "apeglm")

res_kirsten_1_V_0 <- lfcShrink(dds_kirsten, coef = "Time_1_vs_.1", type = "apeglm")
res_kirsten_2_V_0 <- lfcShrink(dds_kirsten, coef = "Time_2_vs_.1", type = "apeglm")
res_kirsten_6_V_0 <- lfcShrink(dds_kirsten, coef = "Time_6_vs_.1", type = "apeglm")
res_kirsten_10_V_0 <- lfcShrink(dds_kirsten, coef = "Time_10_vs_.1", type = "apeglm")
res_kirsten_12_V_0 <- lfcShrink(dds_kirsten, coef = "Time_12_vs_.1", type = "apeglm")


# LRT not powerful enough, will have to just consider timepoint specific effects

res_wiarda = DESeq(dds_wiarda)
res_wiarda_4_V_0 <- lfcShrink(res_wiarda, coef = "Group_Infected_4_vs_Control_0", type = "apeglm")
res_wiarda_10_V_0 <- lfcShrink(res_wiarda, coef = "Group_Infected_10_vs_Control_0", type = "apeglm")





#dds_abdelaal <- DESeq(dds_abdelaal)
#resultsNames(dds_abdelaal)
#dds_abdelaal_20 <- DESeq(dds_abdelaal_20)
#resultsNames(dds_abdelaal_20)


# Infected_8 V Cotrol_ 8
#res_abdelaal_8_v_8 = lfcShrink(dds_abdelaal, coef = "Group_Infected_8_vs_Control_8", type = "apeglm")
#res_abdelaal_20_v_20 = lfcShrink(dds_abdelaal_20, coef = "Group_Infected_20_vs_Control_20", type = "apeglm")


#summary(res_abdelaal_8_v_8)
#summary(res_abdelaal_20_v_20)


# Need to get everything now into a format that is suitable


# Convert to DF
# Natural
res_ogrady_df <- as.data.frame(res_ogrady)

res_kirsten_pbl_df <- as.data.frame(res_kirsten_pbl)

# Experimental
res_kirsten_1_V_0_df <- as.data.frame(res_kirsten_1_V_0)
res_kirsten_2_V_0_df <- as.data.frame(res_kirsten_2_V_0)
res_kirsten_6_V_0_df <- as.data.frame(res_kirsten_6_V_0)
res_kirsten_10_V_0_df <- as.data.frame(res_kirsten_10_V_0)
res_kirsten_12_V_0_df <- as.data.frame(res_kirsten_12_V_0)

res_wiarda_4_V_0_df <- as.data.frame(res_wiarda_4_V_0)
res_wiarda_10_V_0_df <- as.data.frame(res_wiarda_10_V_0)

#res_abdelaal_8_V_8_df <- as.data.frame(res_abdelaal_8_v_8)
#res_abdelaal_20_V_20_df <- as.data.frame(res_abdelaal_20_v_20)


# extract_signif_genes
ogrady_genes = res_ogrady_df %>% filter(padj < 0.05) %>% rownames()
ogrady_genes
kirsten_pbl_genes = res_kirsten_pbl_df %>% filter(padj < 0.05) %>% rownames()

kirsten_1_V_0_genes = res_kirsten_1_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_2_V_0_genes = res_kirsten_2_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_6_V_0_genes = res_kirsten_6_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_10_V_0_genes = res_kirsten_10_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_12_V_0_genes = res_kirsten_12_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_genes = c(kirsten_1_V_0_genes, kirsten_2_V_0_genes, kirsten_6_V_0_genes, kirsten_10_V_0_genes, kirsten_12_V_0_genes)
kirsten_genes = unique(kirsten_genes)

wiarda_4_V_0_genes = res_wiarda_4_V_0_df %>% filter(padj < 0.05) %>% rownames()
wiarda_10_V_0_genes = res_wiarda_10_V_0_df %>% filter(padj < 0.05) %>% rownames()
wiarda_genes = c(wiarda_4_V_0_genes,wiarda_10_V_0_genes)
wiarda_genes = unique(wiarda_genes)


#res_abdelaal_8_V_8_genes <- res_abdelaal_8_V_8_df %>% filter(padj < 0.05) %>% rownames()
#res_abdelaal_20_V_20_genes <- res_abdelaal_20_V_20_df %>% filter(padj < 0.05) %>% rownames()

#abdelaal_genes = c(res_abdelaal_8_V_8_genes,res_abdelaal_20_V_20_genes)
#abdelaal_genes = unique(abdelaal_genes)
#abdelaal_genes



length(ogrady_genes) #  2563
length(kirsten_genes) # 2391
length(wiarda_genes) # 1527
length(kirsten_pbl_genes) # 4101


library(UpSetR)
#save.image(file = "DESEQ2.RData")
#load("DESEQ2.RData")
listInput <- list(OGrady_2025 = ogrady_genes,
                  Mcloughlin_2021 = kirsten_genes,
                  McLoughlin_2014 = kirsten_pbl_genes,
                  Wiarda_2020 = wiarda_genes)

pdf("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/Upset_plot.pdf", width = 15, height = 12)
upset(fromList(listInput), order.by = "freq", nintersects = 25, nsets = 10, sets = names(listInput), query.legend = "top",
      point.size = 4, line.size = 2,  text.scale = c(4, 2.5, 1, 1, 2, 2.5), 
      keep.order = TRUE,sets.x.label = "DE genes", 
      mainbar.y.label = "DE gene intersections", 
      sets.bar.color = scico(n = 4, palette = "glasgow"))
dev.off()

#data_intersection <- data.frame(OGrady_2025 = ogrady_genes,
 #                 Mcloughlin_2021 = kirsten_genes,
  #                McLoughlin_2014 = kirsten_pbl_genes,
   #               Wiarda_2020 = wiarda_genes)

library(VennDiagram)
listInput <- list(ogrady_genes,kirsten_genes,kirsten_pbl_genes,wiarda_genes)
fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}
test = fromList(listInput)

colnames(test) = c("ogrady","Mcloughlin_2021", "Mcloughlin_2014", "Wiarda")
test$sum = rowSums(test)
test
test %>% filter(test$sum == 4)
test
head(res_ogrady_df)
all_ogrady_lfc = res_ogrady_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(OGrady_FC = log2FoldChange, Ogrady_padj = padj) %>% select(OGrady_FC, Ogrady_padj)
all_Mcloughlin_lfc = res_kirsten_pbl_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(Mcloughlin_pbl_FC = log2FoldChange, Mcloughlin_pbl_padj = padj) %>% select(Mcloughlin_pbl_FC, Mcloughlin_pbl_padj)
all_wiarda_W4 = res_wiarda_4_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(wiarda_W4_FC= log2FoldChange, Wiarda_W4_padj = padj) %>% select(wiarda_W4_FC, Wiarda_W4_padj)
all_wiarda_W10 = res_wiarda_10_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(wiarda_W10_FC= log2FoldChange, Wiarda_W10_padj = padj) %>% select(wiarda_W10_FC, Wiarda_W10_padj)

all_kirsten_1_V_0_lfc = res_kirsten_1_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(kirsten_W1_FC= log2FoldChange, kirsten_W1_padj = padj) %>% select(kirsten_W1_FC, kirsten_W1_padj)
all_kirsten_2_V_0_lfc = res_kirsten_2_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(kirsten_W2_FC= log2FoldChange, kirsten_W2_padj = padj) %>% select(kirsten_W2_FC, kirsten_W2_padj)
all_kirsten_6_V_0_lfc = res_kirsten_6_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(kirsten_W6_FC= log2FoldChange, kirsten_W6_padj = padj) %>% select(kirsten_W6_FC, kirsten_W6_padj)
all_kirsten_10_V_0_lfc = res_kirsten_10_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(kirsten_W10_FC= log2FoldChange, kirsten_W10_padj = padj) %>% select(kirsten_W10_FC, kirsten_W10_padj)
all_kirsten_12_V_0_lfc = res_kirsten_12_V_0_df[rownames(test %>% filter(test$sum == 4)),] %>% as.data.frame() %>% select(kirsten_W12_FC= log2FoldChange, kirsten_W12_padj = padj) %>% select(kirsten_W12_FC, kirsten_W12_padj)

all_combined_lfc <- cbind(
  all_ogrady_lfc,
  all_Mcloughlin_lfc,
  all_wiarda_W4,
  all_wiarda_W10,
  all_kirsten_1_V_0_lfc,
  all_kirsten_2_V_0_lfc,
  all_kirsten_6_V_0_lfc,
  all_kirsten_10_V_0_lfc,
  all_kirsten_12_V_0_lfc
)


# Extract LFC values for heatmap
lfc_matrix <- as.matrix(all_combined_lfc[, grepl("_FC$", colnames(all_combined_lfc))])
lfc_matrix
# Create a significance matrix where * marks significant values (adjusted p-value < 0.05)
significance_matrix <- all_combined_lfc[, grepl("_padj$", colnames(all_combined_lfc))] < 0.05
significance_matrix[is.na(significance_matrix)] <- FALSE  # Handle NAs
significance_matrix
# Convert significance matrix to text annotations
annotation_matrix <- ifelse(significance_matrix, "*", "")

length((lfc_matrix))
myColor <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(207)
myBreaks <- c(seq(min(lfc_matrix), 0, length.out=ceiling(length((lfc_matrix))/2) + 1), 
              seq(max(lfc_matrix)/length((lfc_matrix)), max(lfc_matrix), length.out=floor(length((lfc_matrix))/2)))
myBreaks
library(pheatmap)

# Plot heatmap with significance annotation
pdf("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/Heatmap_DE_genes.pdf", width = 18, height = 12)
pheatmap(lfc_matrix, 
         cluster_rows = F, cluster_cols = F,  # Hierarchical clustering
         display_numbers = annotation_matrix,  # Show * for significant values
         color = myColor,
         breaks = myBreaks,# Blue to Red gradient
         fontsize_number = 15,# Adjust * size
         fontsize_row = 12,  # Adjust gene font size
         fontsize_col = 12,  # Adjust column font size
         angle_col = 45,  # Rotate column names,
         legend_position = "bottom",  # Move legend to the top
         labels_row = parse(text = rownames(lfc_matrix)))  # Italicized gene names  # Adjust * size

dev.off()




test= fromList(listInput)
colnames(test) = c("ogrady","Mcloughlin_2021", "Mcloughlin_2014", "Wiarda")
test_ogrady = test %>% filter(ogrady == 1 & Mcloughlin_2014 == 0 & Mcloughlin_2021 == 0 & Wiarda == 0)
test_mcloughlin_pbl = test %>% filter(ogrady == 0 & Mcloughlin_2014 == 1 & Mcloughlin_2021 == 0 & Wiarda == 0)
test_mcloughlin = test %>% filter(ogrady == 0 & Mcloughlin_2014 == 0 & Mcloughlin_2021 == 1 & Wiarda == 0)
test_wiarda = test %>% filter(ogrady == 0 & Mcloughlin_2014 == 0 & Mcloughlin_2021 == 0 & Wiarda == 1)





library(gprofiler2)
rownames(test)
res_ogrady_df
# Individual
input_ogrady_up = res_ogrady_df[rownames(test_ogrady),] %>% arrange(padj) %>% as.data.frame() %>% filter(log2FoldChange > 0) %>% rownames()
input_ogrady_down = res_ogrady_df[rownames(test_ogrady),] %>% arrange(padj) %>% as.data.frame() %>% filter(log2FoldChange < 0) %>% rownames()
results_ogrady_up <- gost(query = input_ogrady_up, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)
results_ogrady_down <- gost(query = input_ogrady_down, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)

input_mcloughlin_pbl_up = res_kirsten_pbl_df[rownames(test_mcloughlin_pbl),] %>% arrange(padj) %>% as.data.frame() %>% filter(log2FoldChange > 0) %>% rownames()
input_mcloughlin_pbl_down = res_kirsten_pbl_df[rownames(test_mcloughlin_pbl),] %>% arrange(padj) %>% as.data.frame() %>% filter(log2FoldChange < 0) %>% rownames()
results_mcloughlin_up <- gost(query = input_mcloughlin_pbl_up, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)
results_mcloughlin_down <- gost(query = input_mcloughlin_pbl_down, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)


res_kirsten_1_V_0_df <- res_kirsten_1_V_0_df %>% mutate(Symbol = rownames(.), Group = "1_V_0")
res_kirsten_2_V_0_df <- res_kirsten_2_V_0_df %>% mutate(Symbol = rownames(.),Group = "2_V_0")
res_kirsten_6_V_0_df <- res_kirsten_6_V_0_df %>% mutate(Symbol = rownames(.),Group = "6_V_0")
res_kirsten_10_V_0_df <- res_kirsten_10_V_0_df %>% mutate(Symbol = rownames(.),Group = "10_V_0")
res_kirsten_12_V_0_df <- res_kirsten_12_V_0_df %>% mutate(Symbol = rownames(.),Group = "12_V_0")
tail(res_kirsten_10_V_0_df)
res_kirsten_all <- rbind(res_kirsten_1_V_0_df,res_kirsten_2_V_0_df,res_kirsten_6_V_0_df,res_kirsten_10_V_0_df,res_kirsten_12_V_0_df)

res_kirsten_input <- res_kirsten_all %>% filter(Symbol %in% rownames(test_mcloughlin)) %>% arrange(padj) %>% filter(!duplicated(Symbol)) %>% rownames()
res_kirsten_input

results_mcloughlin <- gost(query = res_kirsten_input, organism =  "btaurus", correction_method = "fdr", ordered_query = F, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)

head(results_mcloughlin$result)
rownames(res_kirsten_input)


rownames(res_kirsten_all)
# Wiarda


venn.diagram(results_mcloughlin_up$result$term_name, results_mcloughlin_down$result$term_name)

library(ggvenn)

head(results)
results$result$term_name

test_pbl = test %>% filter(ogrady == 0 & Mcloughlin_2014 == 1 & Mcloughlin_2021 == 0 & Wiarda == 1)
results <- gost(query = rownames(test_pbl),organism = "btaurus", correction_method = "fdr", ordered_query = F, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = T)
results$result$term_name
head(test)
