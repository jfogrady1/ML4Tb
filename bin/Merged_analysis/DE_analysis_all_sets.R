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
ogrady_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ogrady_SNPRELATE_eigenvec.txt") %>% select(Sample,PC1,PC2)
ogrady_eigen_vec



colnames(ogrady_eigen_vec) <- c("Sample", "Genotype_PC1", "Genotype_PC2")
labels <- left_join(labels, ogrady_eigen_vec)
rownames(labels) <- labels$Sample
labels$Age <- as.numeric(labels$Age)
labels$Age <- scale(labels$Age, center = TRUE)
labels$Batch <- factor(labels$Batch, levels = c("1", "2"))
labels$Batch
labels

head(labels)
# DESEQ2 analysis
ddsMat_ogrady <- DESeqDataSetFromMatrix(countData = data_raw,
                                 colData = labels,
                                 design = ~ Batch + Genotype_PC1 + Genotype_PC2 + Age + Condition)





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

wiarda_eigen_vec = read.table("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/wiarda_SNPRELATE_eigenvec.txt") %>% select(Sample,PC1,PC2)
wiarda_eigen_vec
colnames(wiarda_eigen_vec) <- c("Sample", "PC1", "PC2")
wiarda_labels$Animal_Code <- as.numeric(wiarda_labels$Animal_Code)
wiarda_labels <- left_join(wiarda_labels, wiarda_eigen_vec, by = c("Animal_Code" = "Sample"))


wiarda_labels
wiarda_labels$Time = factor(wiarda_labels$Time, levels = c("0", "4", "10"))
wiarda_labels$Condition = if_else(wiarda_labels$Condition == "Infected" & wiarda_labels$Time == "0", "Control", wiarda_labels$Condition)

wiarda_labels$Condition = factor(wiarda_labels$Condition, levels = c("Control", "Infected"))
wiarda_labels$Group = paste0(wiarda_labels$Condition, "_", wiarda_labels$Time)                            
table(wiarda_labels$Group)
wiarda_labels$Group <- factor(wiarda_labels$Group, levels = c("Control_0",   "Control_4",  "Infected_4",  "Control_10", "Infected_10"))
wiarda_labels$Group
dds_wiarda = DESeqDataSetFromMatrix(countData = wiarda,
                                     colData = wiarda_labels,
                                     design = ~  PC1 + PC2 + Group)
wiarda_labels
wiarda_labels$Group
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
kirsten_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_SNPRELATE_eigenvec.txt") %>% select(Sample,PC1,PC2)
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

kirsten_pbl_eigen_vec = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/kirsten_pbl_SNPRELATE_eigenvec.txt") %>% select(Sample,PC1,PC2)
colnames(kirsten_pbl_eigen_vec) <- c("Sample", "PC1", "PC2")
kirsten_pbl_eigen_vec
head(kirsten_pbl_labels)
kirsten_pbl_labels <- left_join(kirsten_pbl_labels, kirsten_pbl_eigen_vec, by = c("Animal_Code" = "Sample"))


dds_kirsten_pbl = DESeqDataSetFromMatrix(countData = Kirsten_pbl,
                                     colData = kirsten_pbl_labels,
                                     design = ~ PC1 + PC2 + Condition)




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


length(keep_ogrady)
length(keep_wiarda)
expressed_genes <- intersect(names(keep_ogrady), names(keep_wiarda)) 
expressed_genes <- intersect(expressed_genes, names(keep_kirsten))
expressed_genes <- intersect(expressed_genes, names(keep_kirsten_pbl))


expressed_genes
# filter
ddsMat_ogrady <- ddsMat_ogrady[expressed_genes,] 
dds_kirsten <- dds_kirsten[expressed_genes,]
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

colnames(coef(res_wiarda))





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
kirsten_pbl_genes
kirsten_1_V_0_genes = res_kirsten_1_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_2_V_0_genes = res_kirsten_2_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_6_V_0_genes = res_kirsten_6_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_10_V_0_genes = res_kirsten_10_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_12_V_0_genes = res_kirsten_12_V_0_df %>% filter(padj < 0.05) %>% rownames()
kirsten_genes = c(kirsten_1_V_0_genes, kirsten_2_V_0_genes, kirsten_6_V_0_genes, kirsten_10_V_0_genes, kirsten_12_V_0_genes)
kirsten_genes = unique(kirsten_genes)
kirsten_genes

wiarda_4_V_0_genes = res_wiarda_4_V_0_df %>% filter(padj < 0.05) %>% rownames()
wiarda_10_V_0_genes = res_wiarda_10_V_0_df %>% filter(padj < 0.05) %>% rownames()
wiarda_genes = c(wiarda_4_V_0_genes,wiarda_10_V_0_genes)
wiarda_genes = unique(wiarda_genes)




length(ogrady_genes) #  2979
length(kirsten_genes) # 2228
length(wiarda_genes) # 2486
length(kirsten_pbl_genes) # 4325


library(UpSetR)
#save.image(file = "DESEQ2.RData")
#load("DESEQ2.RData")
listInput <- list(OGrady_2025 = ogrady_genes,
                  Mcloughlin_2021 = kirsten_genes,
                  McLoughlin_2014 = kirsten_pbl_genes,
                  Wiarda_2020 = wiarda_genes)

pdf("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/Upset_plot.pdf", width = 20, height = 15)
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
myColor <- colorRampPalette(c("#2166ac", "white", "#b2182b"))(369)
myBreaks <- c(seq(min(lfc_matrix), 0, length.out=ceiling(length((lfc_matrix))/2) + 1), 
              seq(max(lfc_matrix)/length((lfc_matrix)), max(lfc_matrix), length.out=floor(length((lfc_matrix))/2)))
myBreaks
library(pheatmap)

# Plot heatmap with significance annotation
pdf("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/Heatmap_DE_genes.pdf", width = 20, height =15)
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




# Individual
input_ogrady_up = res_ogrady_df %>% arrange(padj) %>% filter(padj < 0.05) %>% filter(log2FoldChange > 0) %>%  as.data.frame() %>% rownames()
input_ogrady_down = res_ogrady_df %>% arrange(padj) %>% filter (padj < 0.05) %>% filter(log2FoldChange < 0) %>% as.data.frame() %>% rownames()


results_ogrady_up <- gost(query = input_ogrady_up, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC", "WIKI"), evcodes = F)
results_ogrady_down <- gost(query = input_ogrady_down, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "KEGG", "REAC"), evcodes = F)






results_ogrady_down$result <- results_ogrady_down$result %>% mutate(Group = "DE down") %>% mutate(p_value_log10 = -log10(p_value) * -1)
results_ogrady_up$result <- results_ogrady_up$result %>% mutate(Group = "DE up") %>% mutate(p_value_log10 = -log10(p_value))

terms = c("Conavirus disease - COVID-19",
          "aerobic respiration",
          "ribosome biogenesis",
          "Oxidative phosphorylation",
          "Cellular responses to stress",
          "Fcgamma receptor (FCGR) dependent phagocytosis",
          "Macroautophagy",
          "Interleukin-1 signaling",
          "Interleukin-17 signaling",
          "Oxidative phosphorylation",
          "Endocytosis",
          "type I interferon-mediated signaling pathway",
          "IRF3 mediated activation of type 1 IFN",
          "RIG-I-like receptor signaling pathway",
          "defense response to virus",
          "ISG15 antiviral mechanism",
          "rRNA processing in the nucleus and cytosol"
          )


results_ogrady_gprofiler = data.frame(rbind(results_ogrady_up$result, results_ogrady_down$result))


dummy <- as.data.frame(matrix(ncol = ncol(results_ogrady_gprofiler), nrow = 1))
colnames(dummy) <- colnames(results_ogrady_gprofiler)

# Fill in with correct types
dummy[1, "source"] <- "GO:CC"
dummy[1, "p_value"] <- 0  # assuming this column exists and is numeric
# Leave others as NA (will be logical by default), or explicitly set them

# Ensure types match (important for character columns)
dummy <- type.convert(dummy, as.is = TRUE)

results_ogrady_gprofiler <- rbind(results_ogrady_gprofiler, dummy)
results_ogrady_gprofiler <- results_ogrady_gprofiler %>% mutate(label = ifelse(term_name %in% terms, term_name, NA))
results_ogrady_gprofiler <- results_ogrady_gprofiler %>% mutate(alpha = ifelse(term_name %in% terms, 1, 0.4))
results_ogrady_gprofiler <- results_ogrady_gprofiler %>% mutate(Study = "O'Grady et al., 2025")



ggplot(results_ogrady_gprofiler, aes(x = source, y = as.numeric(p_value_log10), col = source, label = label, alpha = alpha)) + geom_jitter(size = 2,position = position_jitter(width = 0.15,seed = 1)) + 
  scale_y_continuous(breaks = seq(from = -21, to = 12, by = 3), limits = c(-21,12)) + 
  geom_hline(yintercept = 1.3, linetype = 2, col = "darkgrey") +
  geom_hline(yintercept = -1.3, linetype = 2, col = "darkgrey") + 
  scale_color_brewer(palette = "Dark2" , name = "Source") +
  theme_bw() +
  geom_text_repel(position = position_jitter(width = 0.15, seed = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 0, color = "black", face = "bold"),
        axis.title.x = element_text(size = 0, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none")

ggsave("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/ogrady_gprofiler.pdf", width = 10, height = 10)

input_mcloughlin_pbl_up = res_kirsten_pbl_df %>% arrange(padj) %>% filter(padj < 0.01) %>%  as.data.frame() %>% filter(log2FoldChange > 0) %>% rownames()
input_mcloughlin_pbl_down = res_kirsten_pbl_df %>% arrange(padj) %>% filter(padj < 0.01) %>% as.data.frame() %>% filter(log2FoldChange < 0) %>% rownames()


results_mcloughlin_up <- gost(query = input_mcloughlin_pbl_up, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = F)
results_mcloughlin_down <- gost(query = input_mcloughlin_pbl_down, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC", "WIKI"), evcodes = F)


terms = c("Deadenylation of mRNA",
          "antigen receptor-mediated signaling pathway",
          "macrophage apoptotic process",
          "Neutrophil degranulation",
          "lysosome",
          "inflammatory response",
          "innate immune response",
          "NOD-like receptor signaling pathway",
          "Immune System",
          "Fc gamma R-mediated phagocytosis"
          )

results_mcloughlin_down$result <- results_mcloughlin_down$result %>% mutate(Group = "DE down") %>% mutate(p_value_log10 = -log10(p_value) * -1)
results_mcloughlin_up$result <- results_mcloughlin_up$result %>% mutate(Group = "DE up") %>% mutate(p_value_log10 = -log10(p_value))
results_mcloughlin_gprofiler = data.frame(rbind(results_mcloughlin_up$result, results_mcloughlin_down$result))
results_mcloughlin_gprofiler <- results_mcloughlin_gprofiler %>% mutate(label = ifelse(term_name %in% terms, term_name, NA))
results_mcloughlin_gprofiler <- results_mcloughlin_gprofiler %>% mutate(alpha = ifelse(term_name %in% terms, 1, 0.4))
results_mcloughlin_gprofiler <- results_mcloughlin_gprofiler %>% mutate(Study = "Mcloughlin et al., 2014")


ggplot(results_mcloughlin_gprofiler, aes(x = source, y = as.numeric(p_value_log10), col = source, label = label, alpha = alpha)) + geom_jitter(size = 2,position = position_jitter(width = 0.15,seed = 1)) + 
  scale_y_continuous(breaks = seq(from = -21, to = 12, by = 3), limits = c(-21,12)) + 
  geom_hline(yintercept = 1.3, linetype = 2, col = "darkgrey") +
  geom_hline(yintercept = -1.3, linetype = 2, col = "darkgrey") + 
  scale_color_brewer(palette = "Dark2" , name = "Source") +
  theme_bw() +
  geom_text_repel(position = position_jitter(width = 0.15, seed = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 0, color = "black", face = "bold"),
        axis.title.x = element_text(size = 0, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none")

ggsave("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/mcloughlin_pbl_gprofiler.pdf", width = 10, height = 10)

# Mcloughlin et al., 2021

res_kirsten_1_V_0_df <- res_kirsten_1_V_0_df %>% mutate(Symbol = rownames(.), Group = "1_V_0")
res_kirsten_2_V_0_df <- res_kirsten_2_V_0_df %>% mutate(Symbol = rownames(.),Group = "2_V_0")
res_kirsten_6_V_0_df <- res_kirsten_6_V_0_df %>% mutate(Symbol = rownames(.),Group = "6_V_0")
res_kirsten_10_V_0_df <- res_kirsten_10_V_0_df %>% mutate(Symbol = rownames(.),Group = "10_V_0")
res_kirsten_12_V_0_df <- res_kirsten_12_V_0_df %>% mutate(Symbol = rownames(.),Group = "12_V_0")


res_kirsten_all <- rbind(res_kirsten_1_V_0_df,res_kirsten_2_V_0_df,res_kirsten_6_V_0_df,res_kirsten_10_V_0_df,res_kirsten_12_V_0_df)

res_kirsten_all <- res_kirsten_all %>% filter(padj < 0.05) 


res_kirsten_input <- res_kirsten_all %>%  group_by(Symbol) %>% arrange(padj) %>% reframe(LFC = log2FoldChange,
                                                                                                                                            FDR = padj,
                                                                                                                                               Group = Group) %>% ungroup() %>% arrange(FDR) %>% select(-c(FDR)) %>%
  pivot_wider(names_from = Group, values_from = LFC) 
 
dim(res_kirsten_input)

consistently_negative <- apply(res_kirsten_input[,2:6], 1, function(row) {
  valid_values <- row[!is.na(row)]
  length(valid_values) > 0 && all(valid_values < 0)
})


consistently_positive <- apply(res_kirsten_input[,2:6], 1, function(row) {
  valid_values <- row[!is.na(row)]
  length(valid_values) > 0 && all(valid_values > 0)
})

consistently_negative
# 3. Subset the original data to show only rows with consistently negative values
negative_rows_kirsten <- res_kirsten_input[consistently_negative, ] 
positive_rows_kirsten <- res_kirsten_input[consistently_positive, ]


dim(negative_rows_kirsten) #  1178 
dim(positive_rows_kirsten) # 1050



results_kirsten_down <- gost(query = negative_rows_kirsten$Symbol, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = F)
results_kirsten_up <- gost(query = positive_rows_kirsten$Symbol, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = F)

results_kirsten_down$result$term_name
results_kirsten_up$result$term_name




terms <- c("mitotic sister chromatid segregation",
           "cell cycle process",
           "leukocyte activation",
           "IL-17 signaling pathway",
           "Regulation of HSF1-mediated heat shock response",
           "mitochondrion",
           "nucleus",
           "Oxidative phosphorylation",
           "mitotic cell cycle",
           "NF-kappa B signaling pathway",
           "T cell differentiation",
           "I-kappaB/NF-kappaB complex",
           "Bcl3/NF-kappaB2 complex",
           "Chemokine receptors bind chemokines",
           "Toll Like Receptor 4 (TLR4) Cascade",
           "Cytokine-cytokine receptor interaction",
           "Cell Cycle",
           "endolysosome"
)


results_kirsten_down$result <- results_kirsten_down$result %>% mutate(Group = "DE down") %>% mutate(p_value_log10 = -log10(p_value) * -1)
results_kirsten_up$result <- results_kirsten_up$result %>% mutate(Group = "DE up") %>% mutate(p_value_log10 = -log10(p_value))

results_kirsten_gprofiler = data.frame(rbind(results_kirsten_up$result, results_kirsten_down$result))


results_kirsten_gprofiler <- results_kirsten_gprofiler %>% mutate(label = ifelse(term_name %in% terms, term_name, NA))
results_kirsten_gprofiler <- results_kirsten_gprofiler %>% mutate(alpha = ifelse(term_name %in% terms, 1, 0.25))
results_kirsten_gprofiler <- results_kirsten_gprofiler %>% mutate(Study = "Mcloughlin et al., 2021")




ggplot(results_kirsten_gprofiler, aes(x = source, y = p_value_log10, col = source, label = label, alpha = alpha)) + geom_jitter(size = 2,position = position_jitter(width = 0.15,seed = 1)) + 
  scale_y_continuous(breaks = seq(from = -21, to = 12, by = 3), limits = c(-21,12)) + 
  geom_hline(yintercept = 1.3, linetype = 2, col = "darkgrey") +
  geom_hline(yintercept = -1.3, linetype = 2, col = "darkgrey") + 
  scale_color_brewer(palette = "Dark2" , name = "Source") +
  theme_bw() +
  geom_text_repel(position = position_jitter(width = 0.15, seed = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 0, color = "black", face = "bold"),
        axis.title.x = element_text(size = 0, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none")

ggsave("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/mcloughlin_gprofiler.pdf", width = 10, height = 10)


# Wiarda


res_wiarda_4_V_0_df$Symbol = rownames(res_wiarda_4_V_0_df)
res_wiarda_10_V_0_df$Symbol = rownames(res_wiarda_10_V_0_df)
res_wiarda_4_V_0_df$Group = "4_V_0"
res_wiarda_10_V_0_df$Group = "10_V_0"

res_wiarda_all <- rbind(res_wiarda_4_V_0_df,res_wiarda_10_V_0_df)

res_wiarda_all <- res_wiarda_all %>% filter(padj < 0.05) 

res_wiarda_all %>% filter(Group == "10_V_0") %>% filter(log2FoldChange > 0) %>% dim()

res_wiarda_input <- res_wiarda_all %>% arrange(padj) %>% group_by(Symbol) %>% arrange(padj) %>% reframe(LFC = log2FoldChange,
                                                                                                        FDR = padj,
                                                                                                        Group = Group) %>% ungroup() %>% arrange(FDR) %>% select(-c(FDR)) %>%
  pivot_wider(names_from = Group, values_from = LFC) 


consistently_negative_wiarda <- apply(res_wiarda_input[,2:3], 1, function(row) {
  valid_values <- row[!is.na(row)]
  length(valid_values) > 0 && all(valid_values < 0)
})


consistently_positive_wiarda <- apply(res_wiarda_input[,2:3], 1, function(row) {
  valid_values <- row[!is.na(row)]
  length(valid_values) > 0 && all(valid_values > 0)
})

negative_rows_wiarda <- res_wiarda_input[consistently_negative_wiarda, ] 
positive_rows_wiarda <- res_wiarda_input[consistently_positive_wiarda, ]



results_wiarda_down <- gost(query = negative_rows_wiarda$Symbol, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = F)
results_wiarda_up <- gost(query = positive_rows_wiarda$Symbol, organism =  "btaurus", correction_method = "fdr", ordered_query = T, custom_bg = expressed_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = F)


results_wiarda_down$result <- results_wiarda_down$result %>% mutate(Group = "DE down") %>% mutate(p_value_log10 = -log10(p_value) * -1)
results_wiarda_up$result <- results_wiarda_up$result %>% mutate(Group = "DE up") %>% mutate(p_value_log10 = -log10(p_value))


negative_rows_wiarda

ograd
# wiarda


results_wiarda_gprofiler = data.frame(rbind(results_wiarda_up$result, results_wiarda_down$result))




terms <- c("defense response",
           "response to cytokine",
           "response to type II interferon",
           "programmed cell death",
           "Complement and coagulation cascades",
           "Integrin cell surface interactions",
           "lytic vacuole",
           "lysosome",
           "Signaling by Interleukin",
           "Cytosolic DNA-sensing pathway",
           "ISG15 antiviral mechanism",
           "Innate Immune System",
           "Influenza A",
           "Nucleotide metabolism",
           "Motor proteins"
           )

results_wiarda_gprofiler <- results_wiarda_gprofiler %>% mutate(label = ifelse(term_name %in% terms, term_name, NA))
results_wiarda_gprofiler <- results_wiarda_gprofiler %>% mutate(alpha = ifelse(term_name %in% terms, 1, 0.25))
results_wiarda_gprofiler <- results_wiarda_gprofiler %>% mutate(Study = "Wiarda et al., 2020")


ggplot(results_wiarda_gprofiler, aes(x = source, y = p_value_log10, col = source, label = label, alpha = alpha)) + geom_jitter(size = 2,position = position_jitter(width = 0.15,seed = 1)) + 
  scale_y_continuous(breaks = seq(from = -21, to = 12, by = 3), limits = c(-21,12)) + 
  geom_hline(yintercept = 1.3, linetype = 2, col = "darkgrey") +
  geom_hline(yintercept = -1.3, linetype = 2, col = "darkgrey") + 
  scale_color_brewer(palette = "Dark2" , name = "Source") +
  theme_bw() +
  geom_text_repel(position = position_jitter(width = 0.15, seed = 1)) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 0, color = "black", face = "bold"),
        axis.title.x = element_text(size = 0, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none")
ggsave("/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/wiarda_gprofiler.pdf", width = 10, height = 10)
dev.off()
colnames(results_wiarda_gprofiler)
colnames(results_ogrady_gprofiler)
colnames(results_kirsten_gprofiler)
colnames(results_mcloughlin_gprofiler)

gprofiler_all = rbind(results_ogrady_gprofiler, results_mcloughlin_gprofiler, results_kirsten_gprofiler,results_wiarda_gprofiler)
gprofiler_all <- apply(gprofiler_all,2,as.character)
write.table(as.data.frame(gprofiler_all), "/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/gprofiler_results.txt", sep = "\t", quote = F, row.names = FALSE)
head(res_ogrady_df %>% mutate(Study = "O'Grady et al., 2025"))


res_ogrady_df <- res_ogrady_df %>% mutate(Symbol = rownames(.), Group = "bTB+_v_bTB-") %>% filter(padj < 0.05)
res_kirsten_pbl_df <- res_kirsten_pbl_df %>% mutate(Symbol = rownames(.), Group = "bTB+_v_bTB-") %>% filter(padj < 0.05)
de_all <- rbind(res_ogrady_df %>% mutate(Study = "O'Grady et al., (2025)"), res_kirsten_pbl_df %>% mutate(Study = "McLoughlin et al., (2014)"), res_wiarda_all %>% mutate(Study = "Wiarda et al., (2020)"), res_kirsten_all %>% mutate(Study = "McLoughlin et al., (2021)"))

head(de_all)
write.table(as.data.frame(de_all), "/home/workspace/jogrady/ML4TB/work/RNA_seq/DE_analysis/all_de_results.txt", sep = "\t", quote = F, row.names = FALSE)
