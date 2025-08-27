# Merging bovine samples
library(tidyverse)
library(gdsfmt)
library(SNPRelate)
library(edgeR)
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
library(magrittr)
library(SNPRelate)
library(dplyr)
library(matrixStats)
args = commandArgs(trailingOnly = TRUE)
source(args[1])
set.seed(42) # For reproducibility
print(args[54])
set_base_url("https://biit.cs.ut.ee/gprofiler_archive3/e112_eg59_p19") # gprofiler
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
data_raw <- fread(args[2])

tested_genes <- data_raw %>% select(Geneid) %>% as.vector()

ensemble <- fread(args[3])
ensemble <- ensemble %>% filter(V3 == "gene")
head(ensemble)
ensemble <- ensemble %>% separate(., V9, into = c("gene_id", "gene_version", "gene_name"), sep = ";")
ensemble$gene_id <- gsub("^gene_id ", "", ensemble$gene_id)
ensemble$gene_id <- gsub('"', '', ensemble$gene_id)
ensemble$gene_name <- gsub("gene_name ", "", ensemble$gene_name)
ensemble$gene_name <- gsub("gene_source ", "", ensemble$gene_name)
ensemble$gene_name <- gsub('"', '', ensemble$gene_name)
ensemble$gene_name <- gsub(" ", "", ensemble$gene_name)
ensemble$gene_name <- if_else(ensemble$gene_name == " ensembl", ensemble$gene_id, ensemble$gene_name)
ensemble$gene_name <- if_else(ensemble$gene_name == " 5S_rRNA", ensemble$gene_id, ensemble$gene_name)
colnames(ensemble)[1] <- "chr"
ensemble <- ensemble %>% dplyr::select(gene_id, gene_name, chr, V4)
colnames(ensemble)[4] <- "pos"
ensemble <- ensemble %>% select(1:2)



ensemble$gene_name <- if_else(duplicated(ensemble$gene_name), ensemble$gene_id, ensemble$gene_name)

print(head(ensemble))

all(tested_genes$Geneid == ensemble$gene_id)

print(head(data_raw))
rownames(data_raw) <- ensemble$gene_name





# Labels and wrangling
labels <- read.table(args[4])
colnames(labels) <- labels[1,]
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
wiarda = fread(args[5], sep = "\t") %>% select(-1)
wiarda_labels = fread(args[6])
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
Kirsten <- fread(args[7], sep = "\t") %>% select(-1)
Kirsten <- as.matrix(Kirsten)
rownames(Kirsten) <- ensemble$gene_name
kirsten_labels <- fread(args[8]) %>% as.data.frame() #%>% select(1,2)
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

print("HERE")

# Kirsten PBL data
Kirsten_pbl <- fread(args[9], sep = "\t") %>% select(-1)
Kirsten_pbl <- as.matrix(Kirsten_pbl)
rownames(Kirsten_pbl) <- ensemble$gene_name
kirsten_pbl_labels <- fread(args[10]) %>% as.data.frame()
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

print("HERE")

wiarda_labels <- wiarda_labels %>% select(colnames(labels))
wiarda_labels <- data.frame(wiarda_labels)
rownames(wiarda_labels) <- wiarda_labels$Sample

kirsten_labels <- kirsten_labels[,colnames(labels)]
kirsten_pbl_labels <- kirsten_pbl_labels[,colnames(labels)]


data_raw <- data_raw %>% select(-1)
data_raw <- as.matrix(data_raw)
rownames(data_raw) <- ensemble$gene_name

dim(Kirsten_pbl)
dim(wiarda)
dim(Kirsten)
dim(data_raw)

data_all <- cbind(data_raw, wiarda, Kirsten, Kirsten_pbl)
metadata_all <- rbind(labels, wiarda_labels, kirsten_labels, kirsten_pbl_labels)

rownames(data_all)
rownames


#####################################################
#####################################################
####
#     Test train split - make sure to ensure time series replicates are not in test and train
####
####################################################
####################################################
####################################################
dim(metadata_all)
dim(data_all)
# Step 1: Identify unique Animal_Codes
unique_animals <- unique(metadata_all$Animal_Code)

unique_animals

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

table(train_set$Animal_Code %in% test_set$Animal_Code) # make sure unique



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

saveRDS(ddsTest,file = args[11])
saveRDS(ddsTrain, file = args[12])


vsdTrain <- vst(ddsTrain)
vsdTest <- vst(ddsTest)

library(scico)


pcaData <- plotPCA(vsdTrain, ntop=750, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Train_PCA <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  xlim(c(-50,50)) + ylim(c(-50, 50))+ theme_bw() + scale_colour_scico_d(palette = "glasgow") +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 


Train_PCA
ggsave(args[13], width = 12, height = 12, dpi = 600)


pcaData <- plotPCA(vsdTest, ntop=750, intgroup=c("Study"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
Test_PCA <- ggplot(pcaData, aes(PC1, PC2, color=Study)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + theme_bw() + scale_colour_scico_d(palette = "glasgow") + xlim(-50,50) + scale_y_reverse(limits = c(50, -50)) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
Test_PCA
ggsave(args[14], width = 12, height = 12, dpi = 600)

# Perform the batch correction

batches= sapply(as.character(train_set$Study), switch, "1_OGrady" = 1, "2_OGrady" = 2, "Wiarda" = 3,  "Mcloughlin" = 4, "Mcloughlin_pbl" = 5, USE.NAMES = FALSE)
groups= sapply(as.character(train_set$Condition), switch, "Control" = 1, "Infected" = 2, USE.NAMES = FALSE)


train_adjusted_counts <- sva::ComBat_seq(as.matrix(train_counts), batch = batches, group = groups)

# Perform the batch correction for test set setting the training data (adjusted) as reference
training_vector = rep("Training", length(train_set$Study))

batches= sapply(as.character(test_set$Study), switch, "1_OGrady" = 1, "2_OGrady" = 2, "Wiarda" = 3, "Mcloughlin" = 4, "Mcloughlin_pbl" = 5, USE.NAMES = FALSE)
groups= sapply(as.character(test_set$Condition), switch, "Control" = 1, "Infected" = 2, USE.NAMES = FALSE)

batches
groups


test_adjusted_counts <- sva::ComBat_seq(as.matrix(test_counts), batch = batches, group = groups)
length(colnames(as.matrix(test_counts)))
length(batches)
test_set$Sample
test_adjusted_counts <- test_adjusted_counts %>% as.data.frame() %>% select(test_set$Sample) %>% as.matrix()
head(test_adjusted_counts)


rownames(train_counts) <- ensemble$gene_name
rownames(test_counts) <- ensemble$gene_name
rownames(train_adjusted_counts) <- ensemble$gene_name
rownames(test_adjusted_counts) <- ensemble$gene_name
rownames(train_adjusted_counts) <- gsub(" ", "", rownames(train_adjusted_counts))
rownames(test_adjusted_counts) <- gsub(" ", "", rownames(test_adjusted_counts))


labels = unique(train_set$Animal_Code)
labels[1:101] = paste0(labels[1:101], "_", labels[1:101])
labels
labels = gsub("_C0.*", "", labels)
labels = gsub("_8.*", "", labels)
labels = gsub("_A.*", "", labels)
labels = gsub("_T0.*", "", labels)

# PCA of genotype data

system(paste0("/home/workspace/jogrady/ML4TB/software/bcftools view -s ", paste(labels, collapse=","), " /home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/ALL_filtered.vcf.gz -Oz -o /home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.vcf.gz && tabix -p vcf /home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.vcf.gz"))

showfile.gds(closeall=TRUE)
snpgdsVCF2GDS("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.vcf.gz", "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.gds", method="biallelic.only")
snpgdsSummary("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.gds")
genofile = snpgdsOpen("/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec2/Train_filtered.gds")
pca <- snpgdsPCA(genofile, num.thread=2)

train_set$Animal_Code

train_set$Animal_Code2 <- gsub("_.*", "", train_set$Animal_Code)
train_set$Animal_Code2
pc.percent <- pca$varprop*100


pc.percent

pca$sample.id


pca$sample.id
# sort out the pca data
# set names
pca$sample.id <- gsub("_.*", "", pca$sample.id)

colnames(pca$eigenvect)[1:ncol(pca$eigenvect)] <- paste0("PC", 1:(ncol(pca$eigenvect)))

pca_df = data.frame(pca$eigenvect)
pca_df$Sample = as.vector(pca$sample.id)
colnames(pca_df)

pca_df <- left_join(pca_df, train_set[!duplicated(train_set$Animal_Code2),], by = c("Sample" = "Animal_Code2"))
dim(pca_df)
pca_df$Study <- gsub("1_", "", pca_df$Study)
pca_df$Study <- gsub("2_", "", pca_df$Study)

pve <- data.frame(PC = 1:20, pve = pca$eigenval[1:20]/sum(pca$eigenval[1:20])*100)
pve

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve$pve
a + ylab("Percentage variance explained") + theme_bw() +
  scale_x_continuous(breaks = c(1:20)) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
sum(pve$pve)
pve$pve

ggsave(args[15], width = 15, height = 12, dpi = 600)
pca_df
pca_df$Study = factor(pca_df$Study, levels = c("Wiarda", "Mcloughlin_pbl", "Mcloughlin", "OGrady"))
pca_df$Condition <- if_else(pca_df$Study == "Wiarda" | pca_df$Study == "Mcloughlin", "Time-series", pca_df$Condition)
ggplot(data = pca_df, aes(x = PC1, y = PC2, fill = Condition, shape = Study)) +
  geom_point(size = 3, stroke = 1) +
  scale_shape_manual(breaks = c("Wiarda", "Mcloughlin_pbl", "Mcloughlin", "OGrady"), values = c(21:24)) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_fill_manual(values = c("darkblue", "darkred", "purple")) +
  labs(colour = "Study", fill = "Condition") + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                  legend.text = element_text(size = 15),
                                                                  axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))
ggsave(args[16], width = 12, height = 12, dpi = 600)

hist(as.numeric(train_set$Age))
train_set$Age
train_set$Age <- cut(as.numeric(train_set$Age), breaks=c(0,6,12,24, 50), labels=c("1-6","6-12","12-24", "24+"))
train_set$Age <- factor(train_set$Age, levels=c("1-6","6-12","12-24", "24+"))
train_set$Sex <- factor(train_set$Sex, levels = c("M", "F"))
train_set$Location <- factor(train_set$Location, levels = c("EU", "US"))
train_set$Infection_administration <- factor(train_set$Infection_administration, levels = c("Natural", "Experimental"))
train_set$Tissue <- factor(train_set$Tissue, levels = c("PB", "PBL"))


# read in genotype data
dim(pca_df)

pca_df$Sample[102:106] <- paste0(pca_df$Sample[102:106],"_CON")
pca_df$Sample[107:112] <- paste0(pca_df$Sample[107:112],"_TB")

pca_df <- pca_df %>% select(PC1, PC2, Sample)
pca_df
train_set <- left_join(train_set, pca_df, by = c("Animal_Code" = "Sample"))
train_set$Condition <- factor(train_set$Condition, levels = c("Control", "Infected"))
colnames(train_set)
train_set






####################################################
####################################################
####################################################
####
#     Differential Expression analysis on training set
####
####################################################
####################################################
####################################################



train_counts <- as.matrix(train_counts)
rownames(train_counts) <- ensemble$gene_name

test_counts <- as.matrix(test_counts)
rownames(train_counts) <- ensemble$gene_name

ddsTrain_adjusted  <- DESeqDataSetFromMatrix(countData = train_counts, 
                                             colData = train_set, 
                                             design = ~ Age + Study + PC1 + PC2 + Condition)

keep_train <- rowSums(counts(ddsTrain_adjusted) >= 6) >= (ncol(ddsTrain_adjusted) *.2) # remove low count genes
ddsTrain_adjusted <- ddsTrain_adjusted[keep_train,]
ddsTrain_adjusted <- DESeq(ddsTrain_adjusted)
assay(ddsTrain_adjusted)
results(ddsTrain_adjusted, alpha = 0.05, name="Condition_Infected_vs_Control")
res <- lfcShrink(ddsTrain_adjusted, coef="Condition_Infected_vs_Control", type="apeglm")

res_df <- as.data.frame(res)


summary(res, alpha =0.05)
res_df <- res_df %>% mutate(Symbol = rownames(res_df),
                            diffexpressed = case_when(
                              log2FoldChange < 0 & padj < 0.05 ~ "#2166ac",
                              log2FoldChange > 0 & padj < 0.05 ~ "#b2182b",
                              FALSE ~ "grey"))


# Print the modified dataframe
res_df <- res_df %>% mutate(Retained = case_when(baseMean > 100 & padj < 0.05 ~ "Retained",
                                                 .default =  "Excluded"))
res_df <- res_df %>% mutate(alpha = if_else(Retained == "Retained", 0.8, 0.15))


symbols_plot <-  res_df %>%
  group_by(diffexpressed, Retained) %>%
  slice_min(n = 10, order_by = padj) %>%
  ungroup() %>%
  mutate(PLOTSYMBOL = TRUE) %>% filter(diffexpressed != "grey" & Retained != "Excluded") %>% select(Symbol, PLOTSYMBOL) %>% mutate(PLOTSYMBOL = Symbol)  # Add a new column marking selected genes

symbols_plot
res_df
symbols_plot$PLOTSYMBOL

res_df <- left_join(res_df, symbols_plot)




write.table(res_df, file = args[17], sep = "\t", quote = FALSE)


# Plotting
ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=PLOTSYMBOL, shape = Retained, alpha = alpha)) +
  geom_point(size = 4) + 
  scale_color_manual("Comparison", values=c("#2166ac", "#b2182b", "grey")) +
  labs(x=bquote(~log[2]~'fold change'),
       y=bquote(~-log[10]~italic(P)[adj]))+
  scale_x_continuous(limits = c(-2,2), breaks = c(-2,-1,0,1,2)) +
  scale_y_continuous(limits = c(0,7)) +
  geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
  theme_bw() +
  geom_text_repel(colour = "black", max.overlaps = 40, size = 3) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 21, color = "black", face = "bold"),
        axis.title.x = element_text(size = 21, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold")) +
  guides(color = guide_legend(override.aes = list(size = 5)))

ggsave(args[18], width = 12, height = 12, dpi = 600)


#out of 16256 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 381, 2.3%
#LFC < 0 (down)     : 393, 2.4%
#outliers [1]       : 16, 0.098%
#low counts [2]     : 1576, 9.7%
#(mean count < 8)
#[1] see 'cooksCutoff' argument of ?results
#[2] see 'independentFiltering' argument of ?results
#[1] see 'cooksCutoff


####################################################
####################################################
#####################################################
####
#     Over representation analysis
####
####################################################
####################################################



ALL_genes <- res_df %>% filter(!is.na(padj))
ALL_genes <- ALL_genes$Symbol # 17787

DE_genes_ORA <- res_df %>% filter(Retained == "Retained") %>% arrange(padj)

DE_genes_ORA <- DE_genes_ORA %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))


DE_genes_ORA_UP <- DE_genes_ORA %>% filter(Direction == "Positive")
DE_genes_ORA_Down <- DE_genes_ORA %>% filter(Direction == "Negative")

ALL_genes
results_ORA_UP <- gost(query = DE_genes_ORA_UP$Symbol,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = ALL_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
results_ORA_UP$result
results_ORA_Down <- gost(query = DE_genes_ORA_Down$Symbol,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = ALL_genes, user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)
results_ORA_Down$result


results_ORA_UP$result$log10 = -log10(results_ORA_UP$result$p_value)
results_ORA_UP$result$Group = "DE Up"
results_ORA_Down$result$log10 = (-log10(results_ORA_Down$result$p_value) * -1)
results_ORA_Down$result$Group = "DE Down"

ORA_results = rbind(results_ORA_UP$result, results_ORA_Down$result)

both_terms <- c("Interleukin-21 signaling",
                "interferon-alpha production",
                "interferon-beta production",
                "NF-kappa B signaling pathway",
                "Antiviral mechanism by IFN-stimulated genes",
                "response to lipopolysaccharide",
                "negative regulation of innate immune response",
                "CD4-positive, alpha-beta T cell differentiation",
                "NOD-like receptor signaling pathway",
                "regulation of macroautophagy",
                "regulation of interleukin-1-mediated signaling pathway",
                "Neutrophil degranulation",
                "innate immune response",
                "cellular response to interleukin-1",
                "regulation of tumor necrosis factor production")



ORA_results = ORA_results %>% mutate(labels = if_else(term_name %in% both_terms, term_name, NA))
ORA_results = ORA_results %>% mutate(alpha = if_else(term_name %in% both_terms, 1, 0.4))



ggplot(ORA_results, aes(x = source, y = as.numeric(log10), col = source, label = labels, alpha = alpha)) + geom_jitter(size = 2,position = position_jitter(width = 0.15,seed = 1)) + 
  scale_y_continuous(breaks = seq(from = -3, to = 12, by = 3), limits = c(-3,12)) + 
  geom_hline(yintercept = 1.3, linetype = 2, col = "darkgrey") +
  geom_hline(yintercept = -1.3, linetype = 2, col = "darkgrey") + 
  scale_colour_manual("Source", values = c("#1b9e77", "#7570b3", "#e7298a", "orange"))+
  theme_bw() +
  geom_text_repel(position = position_jitter(width = 0.15, seed = 1), max.overlaps = 40) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.y = element_text(size = 0, color = "black", face = "bold"),
        axis.title.x = element_text(size = 0, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15, colour = "black", face = "bold"),
        legend.position = "none")

ggsave(args[19], width = 12, height = 12)



####################################################
####################################################
#####################################################
####
#     get CPM couints for GFS method
####
####################################################
####################################################



DE_genes_train <- as.data.frame(res)

DE_genes_train <- DE_genes_train %>% filter(padj < 0.05) %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))

table(DE_genes_train$Direction)

print("HERE")

DE_genes_train_expressed <- DE_genes_train %>% filter(baseMean >= 100 & padj < 0.05)

head(DE_genes_train_expressed)



######
DE_genes_train_expressed = DE_genes_train_expressed %>% arrange(padj) #%>% slice_head(., n = 150)


#Negative Positive 
#344      306 



#Negative Positive 
#69       81 


Pos_genes <- DE_genes_train_expressed %>% filter(Direction == "Positive") %>% rownames()
Neg_genes <- DE_genes_train_expressed %>% filter(Direction == "Negative") %>% rownames()

#train_edge_R_counts <- edgeR::DGEList(counts = train_adjusted_counts, group = train_set$Condition)
#test_edge_R_counts <- edgeR::DGEList(counts = test_adjusted_counts, group = test_set$Condition)


#/ get the normalized counts:
#train_counts_cpms <- edgeR::cpm(train_edge_R_counts, log=FALSE)
#test_counts_cpms <- edgeR::cpm(test_edge_R_counts, log = FALSE)

#write.table(train_counts_cpms, file = "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/train_CPM_filtered_counts.txt", sep = "\t", quote = FALSE)
#write.table(test_counts_cpms, file = "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/test_CPM_filtered_counts.txt", sep = "\t", quote = FALSE)

DE_genes_train_expressed$Symbol <- rownames(DE_genes_train_expressed)
DE_genes_train_expressed_to_subset <- as.vector(DE_genes_train_expressed$Symbol)
length(DE_genes_train_expressed_to_subset)
print(head(DE_genes_train_expressed_to_subset))
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

# Regenerate the DESeq2 object
ddsTrain_adjusted  <- DESeqDataSetFromMatrix(countData = train_counts, 
                                             colData = train_set, 
                                             design = ~ Age + Study + PC1 + PC2 + Condition)


ddsTrain_adjusted <- DESeq(ddsTrain_adjusted)
# Regenerate the DESeq2 object

rownames(test_counts) <- ensemble$gene_name
all(rownames(test_counts) == rownames(train_counts))
ddsTest_adjusted  <- DESeqDataSetFromMatrix(countData = test_counts, 
                                            colData = test_set, 
                                            design = ~ 1) # Note variables arent really important as not doing DE on them, just easy to copy and paste

saveRDS(ddsTrain_adjusted, file = args[20])
write.table(test_counts,file = args[21], row.names = F, quote = F)


ddsTest_adjusted <- DESeq(ddsTest_adjusted)

# variance stabilised transformation on the training data
vstNormalizedExpressionDataForTrain <- varianceStabilizingTransformation(ddsTrain_adjusted, blind = TRUE)



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


set.seed(42)

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

ORA_results <- apply(ORA_results,2,as.character)


train_normalised_filtered_counts <- as.data.frame(train_normalised_filtered_counts)
train_normalised_filtered_counts$Condition <- factor(train_normalised_filtered_counts$Condition, labels = c("Control","Infected"))

write.table(res_df, file =args[22], sep = "\t", row.names = FALSE, quote =F )
write.table(data.frame(ORA_results), file =args[23], sep = "\t", row.names = FALSE, quote = F )
write.table(train_set, file = args[24], sep = "\t", row.names = TRUE, quote =F)
write.table(train_normalised_filtered_counts, file = args[25], sep = "\t", row.names = T)
write.table(train_normalised_filtered_counts[custom_folds[[1]],], file = args[26], sep = "\t", row.names = T)
write.table(train_normalised_filtered_counts[custom_folds[[2]],], file = args[27], sep = "\t", row.names = T)
write.table(train_normalised_filtered_counts[custom_folds[[3]],], file = args[28], sep = "\t", row.names = T)
write.table(train_normalised_filtered_counts[custom_folds[[4]],], file = args[29], sep = "\t", row.names = T)
write.table(train_normalised_filtered_counts[custom_folds[[5]],], file = args[30], sep = "\t", row.names = T)


write.table(test_normalised_filtered_counts, file = args[31], sep = "\t", row.names = T)

length(custom_folds[[1]])

write.table(train_control$index[[1]], file = args[32], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(train_control$index[[2]], file = args[33], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(train_control$index[[3]], file = args[34], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(train_control$index[[4]], file = args[35], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(train_control$index[[5]], file = args[36], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)


test_set
write.table(test_set, file = args[37], sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
dim(train_normalised_filtered_counts)
####################################################
####################################################
####################################################
####
#     Model development --> Simple - logistic regression, penalised (ridge, lasso, elastic net) --> MODELS 1, 2, 3, 4
####
####################################################
####################################################
####################################################



if (!file.exists(args[38])) {
  set.seed(42)
  GLM_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                     method = "glm",
                     preProcess = c("center", "scale"),
                     trControl = train_control,
                     metric = "ROC",
                     family="binomial")
  saveRDS(GLM_model, args[38])
  
  set.seed(42)
  GLMRIDGE_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                          method = "glmnet",
                          tuneGrid = expand.grid(alpha = 0, lambda = seq(0.00001, 1, length=10000)),
                          trControl = train_control,
                          metric = "ROC",
                          verbose = TRUE)
  saveRDS(GLMRIDGE_model, args[39])
  
  
  set.seed(42)
  GLMLASSO_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                          method = "glmnet",
                          tuneGrid = expand.grid(alpha = 1, lambda = seq(0.00001, 1, length=10000)),
                          trControl = train_control,
                          metric = "ROC",
                          verbose = TRUE)
  saveRDS(GLMLASSO_model, args[40])
  
  
  set.seed(42)
  GLMENET_model <- train(Condition ~ ., data = train_normalised_filtered_counts, 
                         method = "glmnet",
                         tuneGrid = expand.grid(alpha = seq(0.05,0.95,0.05), lambda = seq(0.00001, 1, length=10000)),
                         trControl = train_control,
                         metric = "ROC",
                         verbose = TRUE)
  saveRDS(GLMENET_model, args[41])
  
  
  
  
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
  ggsave(args[42], width = 12, height = 12, dpi = 600)
  
  GLMRIDGE_train_ROC <- CV_ROC_plot(GLMRIDGE_model)
  GLMRIDGE_train_ROC
  ggsave(args[43], width = 12, height = 12, dpi = 600)
  
  GLMLASSO_train_ROC <- CV_ROC_plot(GLMLASSO_model)
  GLMLASSO_train_ROC
  ggsave(args[44], width = 12, height = 12, dpi = 600)
  
  GLMENET_train_ROC <- CV_ROC_plot(GLMENET_model)
  GLMENET_train_ROC
  ggsave(args[45], width = 12, height = 12, dpi = 600)
}



GLM_model = readRDS(args[38])
GLMRIDGE_model = readRDS(args[39])
GLMLASSO_model = readRDS(args[40])
GLMENET_model = readRDS(args[41])

GLM_train_ROC <- CV_ROC_plot(GLM_model)
GLM_train_ROC
ggsave(args[42], width = 12, height = 12, dpi = 600)

GLMRIDGE_train_ROC <- CV_ROC_plot(GLMRIDGE_model)
GLMRIDGE_train_ROC
ggsave(args[43], width = 12, height = 12, dpi = 600)

GLMLASSO_train_ROC <- CV_ROC_plot(GLMLASSO_model)
GLMLASSO_train_ROC
ggsave(args[44], width = 12, height = 12, dpi = 600)

GLMENET_train_ROC <- CV_ROC_plot(GLMENET_model)
GLMENET_train_ROC
ggsave(args[45], width = 12, height = 12, dpi = 600)



# Dataframe to save predicitons on test set
Test_AUC <- data.frame(matrix(ncol = 2, nrow = 1))
colnames(Test_AUC) <- c("Model", "AUC")

## Prediction

GLM_predict <- predict(GLM_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLM")
GLM_predict_ROC <- ROC_predict_individual(GLM_predict, test_set$Condition)
GLM_predict_ROC
Test_AUC <- rbind(Test_AUC, c("GLM", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLM_predict[[2]])$auc))
Test_AUC <- Test_AUC[-1,]


GLMRIDGE_predict <- predict(GLMRIDGE_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMRIDGE")
Test_AUC <- rbind(Test_AUC, c("GLMRIDGE", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMRIDGE_predict[[2]])$auc))
GLMRIDGE_ROC <- ROC_predict_individual(GLMRIDGE_predict, test_set$Condition)
GLMRIDGE_ROC

GLMLASSO_predict <- predict(GLMLASSO_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMLASSO")
Test_AUC <- rbind(Test_AUC, c("GLMLASSO", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMLASSO_predict[[2]])$auc))
GLMLASSO_ROC <- ROC_predict_individual(GLMLASSO_predict, test_set$Condition)
GLMLASSO_ROC


GLMENET_predict <- predict(GLMENET_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "GLMENET")
Test_AUC <- rbind(Test_AUC, c("GLMENET", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), GLMENET_predict[[2]])$auc))
GLMENET_ROC <- ROC_predict_individual(GLMENET_predict, test_set$Condition)

GLMRIDGE_ROC
GLMLASSO_ROC
GLMENET_ROC




# save
#ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Train_CV_results/GLM_CV_ROC.pdf", width = 12, height = 12, dpi = 600)


# How to get coefficients for model and filter for those that are non-zero
#test = exp(coef(GLMENET_model$finalModel, GLMENET_model$bestTune$lambda)[,1])
#test
enet_coeff = exp(coef(GLMENET_model$finalModel, GLMENET_model$bestTune$lambda))[,1]
enet_coeff
enet_coeff <- enet_coeff[enet_coeff > 1]
names(enet_coeff)
length(enet_coeff)




####################################################
####################################################
####################################################


# RF model tuning with loops and parallelisation --> MODEL 4

####################################################
####################################################
####################################################


if (!file.exists(args[46])) {
  
  
  library(doParallel) 
  registerDoParallel(20) 
  FOREST <- list()
  
  
  n_features <- length(setdiff(names(train_normalised_filtered_counts), "Condition"))
  n_features
  
  
  sqrt(n_features)
  rf_grid = expand.grid(
    mtry=seq(2, sqrt(n_features), 5),
    min.node.size = seq(2,10),
    splitrule = c("gini")
  )
  colnames(rf_grid)
  
  FOREST <- foreach(num.trees = c(501, 751, 1001, 1501,3001,5001), .packages = c("caret", "ranger"), .combine = 'c') %dopar% {
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
  
  names(FOREST) <- c(501, 751, 1001, 1501,3001,5001)
  
  saveRDS(FOREST, args[46])
  
}


FOREST <- readRDS(args[46])


plot_grid(plot(FOREST[["501"]], metric = "ROC", y = c(0.82,0.95)),
          plot(FOREST[["751"]], metric = "ROC", y = c(0.82,0.95)),
          plot(FOREST[["1001"]], metric = "ROC", y = c(0.82,0.95)),
          plot(FOREST[["1501"]], metric = "ROC", y = c(0.82,0.95)),
          plot(FOREST[["3001"]], metric = "ROC", y = c(0.82,0.95)),
          plot(FOREST[["5001"]], metric = "ROC", y = c(0.82,0.95)))

ggsave(args[47], width = 15, height = 12, dpi = 600)

# Compare models
resamples_rf <- resamples(FOREST[1:6])

summary(resamples_rf)


# Select the best model based on ROC
best_rf_model_name <- names(sort(summary(resamples_rf)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_model_name
RF <- FOREST[[best_rf_model_name]]
RF
best_rf_model_name
RF_train_ROC <- CV_ROC_plot(RF)

RF_train_ROC
ggsave(args[48], width = 12, height = 12, dpi = 600)

RF_predict <- predict(RF, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "RF")
Test_AUC <- rbind(Test_AUC, c("RF", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), RF_predict[[2]])$auc))

rownames(Test_AUC) <- 1:nrow(Test_AUC)

Test_AUC


####################################################
####################################################
####################################################


# RF (EXTRATREES) model tuning with loops and parallelisation --> MODEL 5

####################################################
####################################################
####################################################


if (!file.exists(args[49])) {
  
  
  library(doParallel) 
  registerDoParallel(20) 
  FOREST <- list()
  
  
  n_features <- length(setdiff(names(train_normalised_filtered_counts), "Condition"))
  n_features
  
  
  sqrt(n_features)
  rf_grid = expand.grid(
    mtry=seq(5, 100, 5),
    min.node.size = seq(2,10),
    splitrule = c("extratrees")
  )
  colnames(rf_grid)
  
  EXTRATREES <- foreach(num.trees = c(501, 751, 1001, 1501,3001,5001), .packages = c("caret", "ranger"), .combine = 'c') %dopar% {
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
  
  names(EXTRATREES) <- c(501, 751, 1001, 1501,3001,5001)
  
  saveRDS(EXTRATREES, args[49])
  
}

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
unregister_dopar()

EXTRATREES <- readRDS(args[49])


plot_grid(plot(EXTRATREES[["501"]], metric = "ROC", y = c(0.82,0.95)),
          plot(EXTRATREES[["751"]], metric = "ROC", y = c(0.82,0.95)),
          plot(EXTRATREES[["1001"]], metric = "ROC", y = c(0.82,0.95)),
          plot(EXTRATREES[["1501"]], metric = "ROC", y = c(0.82,0.95)),
          plot(EXTRATREES[["3001"]], metric = "ROC", y = c(0.82,0.95)),
          plot(EXTRATREES[["5001"]], metric = "ROC", y = c(0.82,0.95)))


ggsave(args[50], width = 15, height = 12, dpi = 600)

# Compare models
resamples_rf_et <- resamples(EXTRATREES[1:6])

summary(resamples_rf_et)


# Select the best model based on ROC
best_rf_et_model_name <- names(sort(summary(resamples_rf_et)$statistics$ROC[, "Mean"], decreasing = TRUE)[1])
best_rf_et_model_name
RF_ET <- EXTRATREES[[best_rf_et_model_name]]
RF_ET
best_rf_et_model_name
RF_ET_train_ROC <- CV_ROC_plot(RF_ET)

RF_ET_train_ROC
ggsave(args[51], width = 12, height = 12, dpi = 600)

RF_ET_predict <- predict(RF_ET, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "RF_ET")
Test_AUC <- rbind(Test_AUC, c("RF_ET", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), RF_ET_predict[[2]])$auc))

rownames(Test_AUC) <- 1:nrow(Test_AUC)


Test_AUC

####################################################
####################################################
####################################################


# naive bayes --> MODEL 6

####################################################
####################################################
####################################################


# Define tuning grid
if (!file.exists(args[52])) {
  set.seed(42)
  nb_grid <- expand.grid(usekernel = TRUE,
                         laplace = seq(0,3,0.1), 
                         adjust = seq(0.1,10,0.1))
  
  
  NB_model <- train(Condition ~ ., 
                    data = train_normalised_filtered_counts, 
                    method = "naive_bayes",
                    usepoisson = TRUE,
                    metric = "ROC",
                    trControl = train_control,
                    tuneGrid = nb_grid)
  
  saveRDS(NB_model, args[52])
}

NB_model = readRDS(args[52])

NB_train_ROC <- CV_ROC_plot(NB_model)
NB_train_ROC
ggsave(args[53], width = 12, height = 12, dpi = 600)
NB_predict <- predict(NB_model, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "NB")
Test_AUC <- rbind(Test_AUC, c("NB", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), NB_predict[[2]])$auc))


NB_ROC <- ROC_predict_individual(NB_predict, test_set$Condition)


NB_ROC
tst = roc(test_set$Condition, GLMENET_predict$Infected)
tst$sensitivities
tst$specificities

GLMENET_ROC

####################################################
####################################################
####################################################


# MLP model tuning --> MODEL 7 

# Note - some issues with reproducing the *exact* results on snakemake versus outside in R, results are very similar, however.
# Important to note that optimal model hyperparameters may be different if run on different systems yielding different results
# on the training set and testing set, respectively, however, overall performance of the model is very similar in training and testing.
# Obviously, increasing maxit may improve performance, however, this is a trade-off with time and resources. 
# Feel free to play around as you like.

####################################################
####################################################
####################################################
if (!file.exists(args[54])) {
  
  
  model_weights <- ifelse(train_normalised_filtered_counts$Condition == names(table(train_normalised_filtered_counts$Condition)[1]),
                          (1/table(train_normalised_filtered_counts$Condition)[1]) * 0.5,
                          (1/table(train_normalised_filtered_counts$Condition)[2]) * 0.5)
  model_weights
  
  set.seed(42)
  nnet<- train(Condition ~ ., data = train_normalised_filtered_counts, 
               method = "nnet",           
               trControl = train_control,
               metric = "ROC",
               maxit=200,
               MaxNWts = 150000,
               trace = FALSE,
               tuneGrid=expand.grid(size = c(2,3,4,5,6,7,8,9,10), 
                                    decay = 10^seq(-9,0,by = 1)), 
               na.action = "na.omit",
               allowParallel = FALSE,
               weights = model_weights)
  
  saveRDS(nnet, args[54])
}
nnet = readRDS(args[54])
nnet_train_ROC <- CV_ROC_plot(nnet)
nnet_train_ROC
ggsave(args[55], width = 12, height = 12, dpi = 600)
nnet_predict <- predict(nnet, test_normalised_filtered_counts, type="prob") %>% mutate(Model = "NNET")
nnet_predict
NN_ROC <- ROC_predict_individual(nnet_predict, test_set$Condition)
NN_ROC
Test_AUC <- rbind(Test_AUC, c("MLP", pROC::roc(ifelse(test_set[,"Condition"] == "Infected", 1, 0), nnet_predict[[2]])$auc))





# dotplot of training models
GLMENET_model$resample


Train_dot_plot = rbind(GLM_model$resample,
                       GLMRIDGE_model$resample,
                       GLMLASSO_model$resample,
                       GLMENET_model$resample,
                       NB_model$resample,
                       RF$resample,
                       RF_ET$resample,
                       nnet$resample) %>% mutate(Model = c(rep("Logistic Regression",5),
                                                           rep("Ridge Regression",5),
                                                           rep("Lasso Regression",5),
                                                           rep("Elastic-net Regression",5),
                                                           rep("Naive Bayes",5),
                                                           rep("Random Forest",5),
                                                           rep("Random Forest extra trees",5),
                                                           rep("Multi-layered Perceptron",5)))


model_order <- Train_dot_plot %>%
  group_by(Model) %>%
  summarise(mean_ROC = mean(ROC)) %>%
  arrange(mean_ROC) %>%
  pull(Model)
Train_dot_plot

GLM_model$resample
GLMRIDGE_model$resample
GLMLASSO_model$resample
GLMENET_model$resample
NB_model$resample
RF$resample
RF_ET$resample
nnet$resample


GLM_model
GLMRIDGE_model
GLMLASSO_model
GLMENET_model
NB_model
RF
RF_ET
nnet


Train_dot_plot$Model <- factor(Train_dot_plot$Model, levels = model_order)

ggplot(Train_dot_plot, aes(x = ROC, y = Model, color = Resample, group = Model)) +
  geom_line(linewidth = 2, col = "lightgrey") +
  stat_summary(fun = mean, geom = "point", shape = 21, fill = "black", color = "black", size = 4, stroke = 1.5, show.legend = FALSE) +
  geom_point(size = 4, alpha = 0.7) +
  theme_bw(base_size = 14) +
  scale_colour_npg() +
  labs(x = "AUROC",
       y = "Model") +
  theme(legend.position = c(0.2,0.7))
ggsave(args[56], width = 12, height = 12, dpi = 600)
write.table(Train_dot_plot, file = args[57], sep = "\t", quote = FALSE, row.names = FALSE)

NB_predict

test <- rbind(GLM_predict, GLMRIDGE_predict, GLMLASSO_predict, GLMENET_predict, RF_predict, RF_ET_predict, nnet_predict, NB_predict)
real_test <- ROC_test_combined(test, test_set$Condition)

real_test
ggsave(args[58], width = 12, height = 12, dpi = 600)



####################################################
####################################################
####################################################


# Infected versus experimental performance

####################################################
####################################################
###################################################

test$Sample <- rep(rownames(GLM_predict, 8))
test <- test %>% left_join(., test_set)


combined_predictions = test
combined_predictions$Model <- factor(combined_predictions$Model)

roc_values <- combined_predictions %>%
  group_by(Model, Infection_administration) %>%
  summarise(
    roc_values = as.numeric(pROC::roc(
      response = ifelse(Condition == "Infected", 1, 0),
      predictor = Infected,
      quiet = TRUE
    )$auc),
    AUC_CI_lower = (pROC::ci(pROC::roc(Condition, Infected)))[1],
    AUC_CI_upper = (pROC::ci(pROC::roc(Condition, Infected)))[3])
roc_values <- roc_values %>% arrange(desc(roc_values))

roc_values

write.table(wilcox.test(roc_values ~ factor(Infection_administration), data = roc_values, paired = TRUE)$p.value, file = args[59])
write.table(roc_values, file = args[60], sep = "\t", quote = FALSE, row.names = FALSE)
#data:  roc_values by factor(Infection_administration)
#V = 27, p-value = 0.25
#alternative hypothesis: true location shift is not equal to 0

# Paired differences
diffs <- roc_values[roc_values$Infection_administration == "Experimental",]$roc_values - roc_values[roc_values$Infection_administration == "Natural",]$roc_values 

# Cohen's d for paired samples
d <- mean(diffs) / sd(diffs)
# Example values:
n_per_group <- 8         # number of observations per group
sig_level <- 0.05        # significance level (α)
alternative <- "two.sided"  # or "greater", "less"

library(pwr)
# Calculate power
pwr.t.test(n = n_per_group, d = d, sig.level = sig_level, 
           type = "paired", alternative = alternative)

ggplot(roc_values, aes(x = Infection_administration, y = roc_values,fill = Infection_administration)) + geom_boxplot(alpha = 0.4) + 
  geom_line(aes(group=Model), colour="darkgrey", linetype="11", linewidth = 1.2) + geom_point(aes(col = Model), size = 3) + theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        legend.text = element_text(size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) +
  scale_y_continuous(breaks = c(0.5,0.6,0.7,0.8,0.9,1)) + ylim(0.45,1.05) +
  stat_compare_means(paired = TRUE) + scale_fill_manual(values = c("#FED789FF", "#476F84FF"))
ggsave(args[61], width = 12, height = 12, dpi = 600)


roc_text <- paste(paste0("ROC Resample ", seq_along(roc_values$Model), " = ", roc_values$roc_values, "(95% CI:", round(roc_values$AUC_CI_lower, 2), " - ", round(roc_values$AUC_CI_upper,2), ")"), collapse = "\n")
resample_colors <- ggsci::pal_npg("nrc")(length(roc_values$roc_values))


################################################
# Heatmap based on thresholds
################################################

GLM_predict <- cbind(GLM_predict, test_set) 
GLMRIDGE_predict <- cbind(GLMRIDGE_predict, test_set) 
GLMLASSO_predict<- cbind(GLMLASSO_predict, test_set) 
GLMENET_predict<- cbind(GLMENET_predict, test_set) 
RF_predict<- cbind(RF_predict, test_set) 
RF_ET_predict<- cbind(RF_ET_predict, test_set) 
nnet_predict<- cbind(nnet_predict, test_set) 
NB_predict<- cbind(NB_predict, test_set) 


GLM_ROC <- pROC::roc(GLM_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMRIDGE_ROC <- pROC::roc(GLMRIDGE_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMLASSO_ROC <- pROC::roc(GLMLASSO_predict, response = Condition, predictor = Infected, quiet = TRUE)
GLMENET_ROC <- pROC::roc(GLMENET_predict, response = Condition, predictor = Infected, quiet = TRUE)
RF_ROC <- pROC::roc(RF_predict, response = Condition, predictor = Infected, quiet = TRUE)
RF_ET_ROC <- pROC::roc(RF_ET_predict, response = Condition, predictor = Infected, quiet = TRUE)
nnet_ROC <- pROC::roc(nnet_predict, response = Condition, predictor = Infected, quiet = TRUE)
NB_ROC <- pROC::roc(NB_predict, response = Condition, predictor = Infected, quiet = TRUE)


GLM_thresholds = cbind(GLM_ROC$thresholds, GLM_ROC$specificities, GLM_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLM_thresholds$Combined = rowSums(GLM_thresholds[,c("Spec", "Sens")])

GLMRIDGE_thresholds = cbind(GLMRIDGE_ROC$thresholds, GLMRIDGE_ROC$specificities, GLMRIDGE_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMRIDGE_thresholds$Combined = rowSums(GLMRIDGE_thresholds[,c("Spec", "Sens")])

GLMLASSO_thresholds = cbind(GLMLASSO_ROC$thresholds, GLMLASSO_ROC$specificities, GLMLASSO_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMLASSO_thresholds$Combined = rowSums(GLMLASSO_thresholds[,c("Spec", "Sens")])

GLMENET_thresholds = cbind(GLMENET_ROC$thresholds, GLMENET_ROC$specificities, GLMENET_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GLMENET_thresholds$Combined = rowSums(GLMENET_thresholds[,c("Spec", "Sens")])


RF_thresholds = cbind(RF_ROC$thresholds, RF_ROC$specificities, RF_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
RF_thresholds$Combined = rowSums(RF_thresholds[,c("Spec", "Sens")])


RF_ET_thresholds = cbind(RF_ET_ROC$thresholds, RF_ET_ROC$specificities, RF_ET_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
RF_ET_thresholds$Combined = rowSums(RF_ET_thresholds[,c("Spec", "Sens")])


nnet_thresholds = cbind(nnet_ROC$thresholds, nnet_ROC$specificities, nnet_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
nnet_thresholds$Combined = rowSums(nnet_thresholds[,c("Spec", "Sens")])


NB_thresholds = cbind(NB_ROC$thresholds, NB_ROC$specificities, NB_ROC$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
NB_thresholds$Combined = rowSums(NB_thresholds[,c("Spec", "Sens")])

GLM_final_threshold = GLM_thresholds %>% filter(Sens > 0.7) %>% filter(Spec > 0.1) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMRIDGE_final_threshold = GLMRIDGE_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMLASSO_final_threshold = GLMLASSO_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
GLMENET_final_threshold = GLMENET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
RF_final_threshold = RF_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
RF_ET_final_threshold = RF_ET_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
nnet_final_thresholds = nnet_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()
NB_final_thresholds = NB_thresholds %>% filter(Sens > 0.9) %>% filter(Combined == max(Combined)) %>% slice_max(Sens) %>% select(Threshold) %>% as.numeric()


GLM_thresholds_predict = predict(GLM_model, test_normalised_filtered_counts, type="prob")[,"Infected"]
GLM_thresholds_predict <- ifelse(GLM_thresholds_predict > GLM_final_threshold, "Infected", "Control")

GLMENET_thresholds_predict = predict(GLMENET_model, test_normalised_filtered_counts, type="prob")[,"Infected"]
GLMENET_thresholds_predict <- ifelse(GLMENET_thresholds_predict > GLMENET_final_threshold, "Infected", "Control")

GLMRIDGE_thresholds_predict = predict(GLMRIDGE_model, test_normalised_filtered_counts, type="prob")[,"Infected"]
GLMRIDGE_thresholds_predict <- ifelse(GLMRIDGE_thresholds_predict > GLMRIDGE_final_threshold, "Infected", "Control")

GLMLASSO_thresholds_predict = predict(GLMLASSO_model, test_normalised_filtered_counts, type="prob")[,"Infected"]
GLMLASSO_thresholds_predict <- ifelse(GLMLASSO_thresholds_predict > GLMLASSO_final_threshold, "Infected", "Control")

RF_thresholds_predict = predict(RF, test_normalised_filtered_counts, type="prob")[,"Infected"]
RF_thresholds_predict <- ifelse(RF_thresholds_predict > RF_final_threshold, "Infected", "Control")

RF_ET_thresholds_predict = predict(RF_ET, test_normalised_filtered_counts, type="prob")[,"Infected"]
RF_ET_thresholds_predict <- ifelse(RF_ET_thresholds_predict > RF_ET_final_threshold, "Infected", "Control")

nnet_thresholds_predict = predict(nnet, test_normalised_filtered_counts, type="prob")[,"Infected"]
nnet_thresholds_predict <- ifelse(nnet_thresholds_predict > nnet_final_thresholds, "Infected", "Control")

NB_thresholds_predict = predict(NB_model, test_normalised_filtered_counts, type="prob")[,"Infected"]
NB_thresholds_predict <- ifelse(NB_thresholds_predict > NB_final_thresholds, "Infected", "Control")

NB_final_thresholds

thresholds_df <- as.matrix(cbind(GLMENET_thresholds_predict, GLMLASSO_thresholds_predict, nnet_thresholds_predict, RF_ET_thresholds_predict, RF_thresholds_predict, GLMRIDGE_thresholds_predict, NB_thresholds_predict, GLM_thresholds_predict))
thresholds_df <- ifelse(thresholds_df == "Infected", 1, 0)
mode(thresholds_df) <- "numeric"
rownames(thresholds_df) <- test_set$Sample

test_set <- test_set %>% mutate(Infection_administration = if_else(Study == "Mcloughlin_pbl", "Natural", Infection_administration))
test_set$Study <- if_else(test_set$Study == "1_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "2_OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Wiarda", "Wiarda et al., (2020)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin", "McLoughlin et al., (2021)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin_pbl", "McLoughlin et al., (2014)", test_set$Study)

test_set$Study <- factor(test_set$Study, levels = c("O'Grady et al., (2025)", "McLoughlin et al., (2014)", "Wiarda et al., (2020)", "McLoughlin et al., (2021)"))
test_set <- test_set %>% arrange(Study)
test_set

thresholds_df <- thresholds_df[rownames(test_set),]
thresholds_df

thresholds_df
colnames(thresholds_df) <- c("ENET", "LASSO", "MLP", "RF(ET)", "RF", "RIDGE", "NB", "GLM")

ann_colors <- list(
  Study = c("O'Grady et al., (2025)" = "#351338",
            "Wiarda et al., (2020)" = "#dad2ff",
            "McLoughlin et al., (2021)" = "#734500",
            "McLoughlin et al., (2014)" = "#69a3a5"),
  Condition = c("Control" = "#2166ac",
                "Infected" = "#b2182b"),
  Infection_administration = c(Natural = "#708238",
                               Experimental = "beige"))

annotation_columns = test_set %>% select(c(Condition, Study, Infection_administration))

write.table(thresholds_df, file = args[62], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
library(pheatmap)
pdf(args[63], width = 15, height = 20)
pheatmap(t(thresholds_df),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_row = c(seq(2:8)),
         annotation_col = annotation_columns,
         annotation_colors = ann_colors,
         border_color = "black",
         show_colnames = F,
         color = c("lightblue", "orange"))