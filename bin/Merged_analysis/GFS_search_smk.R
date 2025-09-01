# ROC feed forward evaluation

library(tidyverse)
library(DESeq2)

VST_counts <- readRDS("/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/Train_adjusted_DESEQ2_matrix.rds")
VST_counts <- varianceStabilizingTransformation(VST_counts, blind = TRUE)


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



# Labensemble# Labels and wrangling
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




wiarda_labels <- wiarda_labels %>% select(colnames(labels))
wiarda_labels <- data.frame(wiarda_labels)
rownames(wiarda_labels) <- wiarda_labels$Sample

kirsten_labels <- kirsten_labels[,colnames(labels)]
kirsten_pbl_labels <- kirsten_pbl_labels[,colnames(labels)]

metadata_all <- rbind(labels, wiarda_labels, kirsten_labels, kirsten_pbl_labels)

metadata_all



DE_genes <- fread("/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/DE_results_integrated.txt") %>% filter(padj < 0.05) %>% filter(baseMean > 100) %>% select(Symbol)
DE_genes <- as.vector(DE_genes$Symbol)
length(DE_genes)
DE_genes
train_normalised_filtered_counts <- assay(VST_counts)[DE_genes,]

#metadata <- read.table("/home/workspace/jogrady/ML4TB/work/")


source("~/ML4TB/bin/Merged_analysis/Funcitons.R")


train_set <- read.table(file = "/home/workspace/jogrady/ML4TB/snktest/work/merged/Temp_files/train_data_manuscript.txt")

rownames(train_set) = train_set$Sample

DE_results = fread("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/DE_results_integrated.txt") %>% filter(padj < 0.05) %>% filter(baseMean > 100) %>% select(-V1) %>% mutate(Direction = if_else(log2FoldChange < 0, "Negative", "Positive"))
rownames(DE_results) <- DE_results$Symbol
dim(DE_results)
metadata_train <- metadata_all %>% filter(rownames(.) %in% colnames(train_normalised_filtered_counts))
metadata_train$Study <- gsub("1_", "", metadata_train$Study)
metadata_train$Study <- gsub("2_", "", metadata_train$Study)
train_set

GFS_train_counts = t(train_normalised_filtered_counts[,-c(length(colnames(train_normalised_filtered_counts)))])
GFS_train_counts <- t(GFS_train_counts)
head(GFS_train_counts)
head(GFS_train_counts)

hist(GFS_train_counts["UBASH3A",])


DE_genes <- DE_genes[DE_genes %in% rownames(GFS_train_counts[1:10,])]

ROC = greedy_forward_search(gene.list = DE_genes, cpm_matrix = GFS_train_counts[1:10,], metadata = train_set, de_results = DE_results[DE_results$Symbol %in% rownames(GFS_train_counts[1:10,]),])

View(ROC)
View(DE_results)

DE_genes_train_expressed_to_subset
DE_genes
GFS_train_counts

colnames(GFS_train_counts) == rownames(train_set)
write.table(ROC, file = "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC1_genes.txt", sep = "/t", row.names = F, quote = FALSE)
GFS_1_genes <- ROC$combination[14]
GFS_1_genes <- str_split_1(GFS_1_genes, pattern = "_")


DE_results1 = DE_results %>% filter(!(Symbol %in% GFS_1_genes))
DE_results
ROC2 = greedy_forward_search(gene.list = DE_results1$Symbol, cpm_matrix = GFS_train_counts[!rownames(GFS_train_counts) %in% GFS_1_genes, ], metadata = train_set, de_results = DE_results1)
write.table(ROC2, file = "/home/workspace/jogrady/ML4TB/work/merged/Temp_files/ROC2_genes.txt", sep = "/t", row.names = F, quote = FALSE)
GFS_2_genes <- ROC2$combination[18]
GFS_2_genes <- str_split_1(GFS_2_genes, pattern = "_")



Pos_genes1 = DE_results %>% filter(Symbol %in% GFS_1_genes & log2FoldChange > 0) %>% select(Symbol)
Neg_genes1 = DE_results %>% filter(Symbol %in% GFS_1_genes & log2FoldChange < 0) %>% select(Symbol)

Pos_genes2 = DE_results %>% filter(Symbol %in% GFS_2_genes & log2FoldChange > 0) %>% select(Symbol)
Neg_genes2 = DE_results %>% filter(Symbol %in% GFS_2_genes & log2FoldChange < 0) %>% select(Symbol)

Pos_genes3 = DE_results %>% filter(Symbol %in% c(GFS_1_genes,GFS_2_genes) & log2FoldChange > 0) %>% select(Symbol)
Neg_genes3 = DE_results %>% filter(Symbol %in% c(GFS_1_genes,GFS_2_genes) & log2FoldChange < 0) %>% select(Symbol)



if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =GFS_train_counts[rownames(GFS_train_counts) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =GFS_train_counts[rownames(GFS_train_counts) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}




Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result

Scores_for_plotting1$Score <- Scores_for_plotting1$pos_Score

# Calculate Score
Scores_for_plotting2$Score <- Scores_for_plotting2$pos_Score - Scores_for_plotting2$neg_Score
#Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))

Scores_for_plotting3$Score <- Scores_for_plotting3$pos_Score - Scores_for_plotting3$neg_Score
#Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))

Scores_for_plotting1
Scores_for_plotting2
Scores_for_plotting3

fold_assignments

Scores_for_plotting1$Sample <- rownames(Scores_for_plotting1)
Scores_for_plotting2$Sample <- rownames(Scores_for_plotting2)
Scores_for_plotting3$Sample <- rownames(Scores_for_plotting3)


Scores_for_plotting1 <- left_join(Scores_for_plotting1, train_set)
Scores_for_plotting2 <- left_join(Scores_for_plotting2, train_set)
Scores_for_plotting3 <- left_join(Scores_for_plotting3, train_set)


Scores_for_plotting1$Study <- gsub("1_", "", Scores_for_plotting1$Study)
Scores_for_plotting1$Study <- gsub("2_", "", Scores_for_plotting1$Study)


Scores_for_plotting2$Study <- gsub("1_", "", Scores_for_plotting2$Study)
Scores_for_plotting2$Study <- gsub("2_", "", Scores_for_plotting2$Study)


Scores_for_plotting3$Study <- gsub("1_", "", Scores_for_plotting3$Study)
Scores_for_plotting3$Study <- gsub("2_", "", Scores_for_plotting3$Study)



Scores_for_plotting1$Study <- factor(Scores_for_plotting1$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Scores_for_plotting2$Study <- factor(Scores_for_plotting2$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Scores_for_plotting3$Study <- factor(Scores_for_plotting3$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))



Scores_for_plotting1$Condition <- factor(Scores_for_plotting1$Condition, levels = c(0,1), labels = c("Control", "Infected"))
Scores_for_plotting2$Condition <- factor(Scores_for_plotting2$Condition, levels = c(0,1), labels = c("Control", "Infected"))
Scores_for_plotting3$Condition <- factor(Scores_for_plotting3$Condition, levels = c(0,1), labels = c("Control", "Infected"))

Scores_for_plotting1$Set = "Set1"
Scores_for_plotting2$Set = "Set2"
Scores_for_plotting3$Set = "Set3"

Scores_for_plotting1$Score <- as.numeric(scale(Scores_for_plotting1$Score))
Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))
Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))


Scores_for_plotting2 <- Scores_for_plotting2[,colnames(Scores_for_plotting1)]
Scores_for_plotting3 <- Scores_for_plotting3[,colnames(Scores_for_plotting1)]



Scores_for_plotting_combined <- rbind(Scores_for_plotting1, Scores_for_plotting2, Scores_for_plotting3)
Scores_for_plotting_combined$Set <- factor(Scores_for_plotting_combined$Set, levels = c("Set1", "Set2", "Set3"))

ggplot(Scores_for_plotting_combined, aes(y = Score, x = factor(Study, levels = c("OGrady","Mcloughlin_pbl", "Wiarda", "Mcloughlin")), fill = Set, shape = Condition, col = Condition)) + 
  geom_boxplot(position = position_dodge2(width = 2, padding = 1.9), outlier.colour = NA, alpha = 1) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5, size = 3) + 
  scale_colour_manual(values = c("#2166ac", "#b2182b")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  xlab("Study") +
  ylab("Score") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        #legend.position = c(.9, .25),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        #legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14)) #change legend title font size

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/TB_scoe_training.pdf", width = 12, height = 12, dpi = 600)
Scores_for_plotting_combined %>%
  group_by(Study, Set) %>%
  wilcox_test(Score ~ Condition) %>%
  adjust_pvalue(method = "BH") %>%  # Optional: adjust p-values
  add_significance() %>% arrange(desc(Study))   

# A tibble: 12 × 11
#Study          Set   .y.   group1  group2      n1    n2 statistic        p    p.adj p.adj.signif
#<fct>          <fct> <chr> <chr>   <chr>    <int> <int>     <dbl>    <dbl>    <dbl> <chr>       
#  1 Mcloughlin_pbl Set1  Score Control Infected     5     6         0 4.33e- 3 4.33e- 3 **          
#  2 Mcloughlin_pbl Set2  Score Control Infected     5     6         0 4.33e- 3 4.33e- 3 **          
#  3 Mcloughlin_pbl Set3  Score Control Infected     5     6         0 4.33e- 3 4.33e- 3 **          
#  4 Mcloughlin     Set1  Score Control Infected     6    28         0 1.49e- 6 1.99e- 6 ****        
#  5 Mcloughlin     Set2  Score Control Infected     6    28         0 1.49e- 6 1.99e- 6 ****        
#  6 Mcloughlin     Set3  Score Control Infected     6    28         0 1.49e- 6 1.99e- 6 ****        
#  7 Wiarda         Set1  Score Control Infected    17    16         3 1.20e- 8 2.40e- 8 ****        
#  8 Wiarda         Set2  Score Control Infected    17    16         0 1.71e- 9 4.10e- 9 ****        
#  9 Wiarda         Set3  Score Control Infected    17    16         0 1.71e- 9 4.10e- 9 ****        
#  10 OGrady         Set1  Score Control Infected    41    43        39 2.17e-19 8.68e-19 ****        
#  11 OGrady         Set2  Score Control Infected    41    43        20 3.31e-21 3.30e-20 ****        
#  12 OGrady         Set3  Score Control Infected    41    43        22 5.50e-21 3.30e-20 ****  



GFS_thresholds_Set1 = roc(Condition ~ Score, data = Scores_for_plotting_combined[Scores_for_plotting_combined$Set == "Set1",]) 
GFS_thresholds_Set2 = roc(Condition ~ Score, data = Scores_for_plotting_combined[Scores_for_plotting_combined$Set == "Set2",])
GFS_thresholds_Set3 = roc(Condition ~ Score, data = Scores_for_plotting_combined[Scores_for_plotting_combined$Set == "Set3",]) 

GFS_thresholds_Set1 <- cbind(GFS_thresholds_Set1$thresholds, GFS_thresholds_Set2$sensitivities, GFS_thresholds_Set3$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set1")
GFS_thresholds_Set1

View(cbind(GFS_thresholds_Set2$thresholds, GFS_thresholds_Set2$sensitivities, GFS_thresholds_Set2$specificities))
GFS_train_roc_AUC1 <- Scores_for_plotting1 %>% group_by(Folds) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                             AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                             AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3]) %>% mutate(set = "Set1")


GFS_train_roc_AUC2 <-  Scores_for_plotting2 %>% group_by(Folds) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                              AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                              AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3]) %>% mutate(set = "Set2")


GFS_train_roc_AUC3 <-  Scores_for_plotting3 %>% group_by(Folds) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                              AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                              AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3]) %>% mutate(set = "Set3")






GFS_train_combined = rbind(GFS_train_roc_AUC1,GFS_train_roc_AUC2,GFS_train_roc_AUC3)

GFS_train_combined


ggplot(GFS_train_combined, aes(x = Folds, y = AUC, col = set, group = factor(set))) + geom_pointrange(aes(ymin = AUC_CI_lower, ymax = AUC_CI_upper, color = set),
                                                                                                      position = position_dodge(width = 0.2)) + ylim(0.95,1) + 
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = c(.9, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14)) #change legend title font size

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/AUC_0.95_3Sets_GFS.pdf", width = 12, height = 12, dpi = 600)



GFS_train_roc_AUC3





Scores_for_plotting_combined$Study <- as.character(Scores_for_plotting_combined$Study)
Scores_for_plotting_combined$Study <- c(Scores_for_plotting_combined$Study, rep("Overall", 972))
Scores_for_plotting_long$Study <- factor(Scores_for_plotting_long$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Overall"))
Scores_for_plotting_wide = Scores_for_plotting_long %>% pivot_wider(names_from = "Study", values_from = "Score")

GFS_thresholds_overall = roc(Condition ~ Overall, data = Scores_for_plotting_wide)
GFS_thresholds_overall = cbind(GFS_thresholds_overall$thresholds, GFS_thresholds_overall$specificities, GFS_thresholds_overall$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GFS_thresholds_overall$Combined = rowSums(GFS_thresholds_overall[,c("Spec", "Sens")])

max(GFS_thresholds_overall$Combined)
GFS_thresholds_overall
#0.95652174 0.95698925 1.913511
GFS_train_roc_AUC <- Scores_for_plotting_long %>% group_by(Study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                                AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                                AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3])


GFS_train_score = ggplot(Scores_for_plotting_long, aes(m=Score, d=factor(Condition, levels = c("Control", "Infected")), colour = Study)) + # Note do not need to filter unless saving all predictions
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
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + 
  scale_colour_manual(breaks = c("OGrady", "Mcloughlin", "Mcloughlin_pbl", "Wiarda", "Overall"), values = c("#351338","#734500","#5F927D","#DAD2FF", "black"))

resample_colors <-  c("#5F927D","#734500","black","#351338","#DAD2FF")
GFS_train_roc_AUC <- GFS_train_roc_AUC %>% arrange(desc(AUC))
for (i in seq_along(GFS_train_roc_AUC$AUC)) {
  GFS_train_score <- GFS_train_score + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("AUROC ", GFS_train_roc_AUC$Study[i], " = ", round(GFS_train_roc_AUC$AUC[i],3), " (", round(GFS_train_roc_AUC$AUC_CI_lower[i],2), " - ", round(GFS_train_roc_AUC$AUC_CI_upper[i],2), ")"), 
             hjust = .25, 
             size = 5,
             col = resample_colors[i])
}

GFS_train_score
# Extract CI data (All Score only)
ci_vals <- ci.se(roc(Condition ~ All_Score, data = Scores_for_plotting_wide, legacy.axes = TRUE), specificities = seq(0, 1, 0.01))
df <- data.frame(
  specificity = as.numeric(rownames(ci_vals)),
  sensitivity = ci_vals[, "50%"],
  lower = ci_vals[, "2.5%"],
  upper = ci_vals[, "97.5%"],
  model = "Overall")

GFS_train_score + 
  geom_ribbon(data = df, inherit.aes = FALSE, aes(x = 1 - specificity, ymin = lower, ymax = upper), alpha = 0.2, color = "grey") +
  annotate("point", x = 1 -0.95652174, y = 0.95698925, colour = "darkred", size = 8, shape = 18)

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/GFS_Train_ROC_curves.pdf", width = 12, height = 12, dpi = 600)



# Correlation between Score and IGRA
IGRA = fread("/home/workspace/jogrady/ML4TB/data/ogrady/Ogrady_IGRA_results.txt")
Scores_for_plotting <- left_join(Scores_for_plotting, IGRA, by = c("Sample" = "Sample ID"))
ggplot(Scores_for_plotting, aes(x = Score, y = `IFGA_bovine - IFGA_avian`)) + geom_point(size = 3, aes(col = Condition)) +
  scale_color_manual(values = c("#2166ac", "#b2182b")) +
  theme_bw() +
  ylab("IFGA bovine - IFGA avian") +
  xlab("Score") +
  geom_smooth(method = "lm", col = "black") +
  geom_hline(yintercept = 80, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept =  5.307279, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position = c(.9, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=2)) +
  annotate("text", x = 5.1, y = 2500,
           label = "paste(italic(ρ), \" = 0.593\")", parse = TRUE) +
  annotate("text", x = 5.1, y = 2200,
           label = "paste(italic(P), ' = 2.79 x 10^-9')", parse = TRUE)

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/IGRA_TB_Score_correlation.pdf", width = 12, height = 12, dpi = 600)


# Test set


Test_data = as.data.frame(t(read.table("~/ML4TB/work/merged/Temp_files/test_normalised_filtered_counts.txt", sep = "\t", check.names=F)))


if (length(Pos_genes1$Symbol)==1){
  Pos_counts_ROC1 = Test_data[rownames(Test_data) %in% Pos_genes1$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC1 =Test_data[rownames(Test_data) %in% Pos_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC1=NA
}


if (length(Neg_genes1$Symbol)==1){
  Neg_counts_ROC1 =Test_data[rownames(Test_data) %in% Neg_genes1$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes)>1) {
  Neg_counts_ROC1 =Test_data[rownames(Test_data) %in% Neg_genes1$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC1 = 0
}


# Set 2
if (length(Pos_genes2$Symbol)==1){
  Pos_counts_ROC2 = Test_data[rownames(Test_data) %in% Pos_genes2$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes1$Symbol)>1) {
  Pos_counts_ROC2 =Test_data[rownames(Test_data) %in% Pos_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC2=NA
}


if (length(Neg_genes2$Symbol)==1){
  Neg_counts_ROC2 =Test_data[rownames(Test_data) %in% Neg_genes2$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes2$Symbol)>1) {
  Neg_counts_ROC2 =Test_data[rownames(Test_data) %in% Neg_genes2$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC2 = 0
}


# Set 3
if (length(Pos_genes3$Symbol)==1){
  Pos_counts_ROC3 = Test_data[rownames(Test_data) %in% Pos_genes3$Symbol,] %>% as.data.frame() %>% mutate_all(as.numeric) %>% set_colnames("pos_Score")
} else if (length(Pos_genes2$Symbol)>1) {
  Pos_counts_ROC3 =Test_data[rownames(Test_data) %in% Pos_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("pos_Score")    
} else{
  Pos_counts_ROC3=NA
}


if (length(Neg_genes3$Symbol)==1){
  Neg_counts_ROC3 =Test_data[rownames(Test_data) %in% Neg_genes3$Symbol,]%>% as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")
} else if (length(Neg_genes3$Symbol)>1) {
  Neg_counts_ROC3 =Test_data[rownames(Test_data) %in% Neg_genes3$Symbol,] %>% 
    as.data.frame() %>% mutate_all(as.numeric) %>% 
    summarise_all(mean) %>% t() %>% set_colnames("neg_Score")    
} else{
  Neg_counts_ROC3 = 0
}


Test_Scores_for_plotting1 = as.data.frame(data.frame(Pos_counts_ROC1)) # No neg genes in this GFS result
Test_Scores_for_plotting2 = data.frame(cbind(Pos_counts_ROC2, Neg_counts_ROC2)) # No neg genes in this GFS result
Test_Scores_for_plotting3 = as.data.frame(cbind(Pos_counts_ROC3, Neg_counts_ROC3)) # No neg genes in this GFS result


Test_Scores_for_plotting1$Score <- Test_Scores_for_plotting1$pos_Score

# Calculate Score
Test_Scores_for_plotting2$Score <- Test_Scores_for_plotting2$pos_Score - Test_Scores_for_plotting2$neg_Score
#Scores_for_plotting2$Score <- as.numeric(scale(Scores_for_plotting2$Score))

Test_Scores_for_plotting3$Score <- Test_Scores_for_plotting3$pos_Score - Test_Scores_for_plotting3$neg_Score
#Scores_for_plotting3$Score <- as.numeric(scale(Scores_for_plotting3$Score))

Test_Scores_for_plotting1
Test_Scores_for_plotting2
Test_Scores_for_plotting3



Test_Scores_for_plotting1$Sample <- rownames(Test_Scores_for_plotting1)
Test_Scores_for_plotting2$Sample <- rownames(Test_Scores_for_plotting2)
Test_Scores_for_plotting3$Sample <- rownames(Test_Scores_for_plotting3)


Test_Scores_for_plotting1 <- left_join(Test_Scores_for_plotting1, metadata_all)
Test_Scores_for_plotting2 <- left_join(Test_Scores_for_plotting2, metadata_all)
Test_Scores_for_plotting3 <- left_join(Test_Scores_for_plotting3, metadata_all)


Test_Scores_for_plotting1$Study <- gsub("1_", "", Test_Scores_for_plotting1$Study)
Test_Scores_for_plotting1$Study <- gsub("2_", "", Test_Scores_for_plotting1$Study)


Test_Scores_for_plotting2$Study <- gsub("1_", "", Test_Scores_for_plotting2$Study)
Test_Scores_for_plotting2$Study <- gsub("2_", "", Test_Scores_for_plotting2$Study)


Test_Scores_for_plotting3$Study <- gsub("1_", "", Test_Scores_for_plotting3$Study)
Test_Scores_for_plotting3$Study <- gsub("2_", "", Test_Scores_for_plotting3$Study)


Test_Scores_for_plotting1$Study <- factor(Test_Scores_for_plotting1$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Test_Scores_for_plotting2$Study <- factor(Test_Scores_for_plotting2$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Test_Scores_for_plotting3$Study <- factor(Test_Scores_for_plotting3$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))



Test_Scores_for_plotting1$Condition <- factor(Test_Scores_for_plotting1$Condition, levels = c("Control", "Infected"))
Test_Scores_for_plotting2$Condition <- factor(Test_Scores_for_plotting2$Condition, levels = c("Control", "Infected"))
Test_Scores_for_plotting3$Condition <- factor(Test_Scores_for_plotting3$Condition, levels = c("Control", "Infected"))

Test_Scores_for_plotting1$Set = "Set1"
Test_Scores_for_plotting2$Set = "Set2"
Test_Scores_for_plotting3$Set = "Set3"

Test_Scores_for_plotting1$Score <- as.numeric(scale(Test_Scores_for_plotting1$Score))
Test_Scores_for_plotting2$Score <- as.numeric(scale(Test_Scores_for_plotting2$Score))
Test_Scores_for_plotting3$Score <- as.numeric(scale(Test_Scores_for_plotting3$Score))


Test_Scores_for_plotting2 <- Test_Scores_for_plotting2[,colnames(Test_Scores_for_plotting1)]
Test_Scores_for_plotting3 <- Test_Scores_for_plotting3[,colnames(Test_Scores_for_plotting1)]



Test_Scores_for_plotting_combined <- rbind(Test_Scores_for_plotting1, Test_Scores_for_plotting2, Test_Scores_for_plotting3)
Test_Scores_for_plotting_combined$Set <- factor(Test_Scores_for_plotting_combined$Set, levels = c("Set1", "Set2", "Set3"))


ggplot(Test_Scores_for_plotting_combined, aes(y = Score, x = factor(Study, levels = c("OGrady","Mcloughlin_pbl", "Wiarda", "Mcloughlin")), fill = Set, shape = Condition, col = Condition)) + 
  geom_boxplot(position = position_dodge2(width = 2, padding = 1.9), outlier.colour = NA, alpha = 1) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.5, size = 3) + 
  scale_colour_manual(values = c("#2166ac", "#b2182b")) +
  scale_fill_brewer(palette = "Set2") +
  theme_bw() +
  xlab("Study") +
  ylab("Score") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        #legend.position = c(.9, .25),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        #legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14)) #change legend title font size

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/TB_scoe_testing.pdf", width = 12, height = 12, dpi = 600)
Test_Scores_for_plotting_combined %>% filter(Study != "Wiarda") %>%
  group_by(Study, Set) %>%
  wilcox_test(Score ~ Condition) %>%
  adjust_pvalue(method = "BH") %>%  # Optional: adjust p-values
  add_significance() %>% arrange(desc(Study))   


#1 Mcloughlin_pbl Set1  Score Control Infected     3     2         0 0.2       0.2      ns          
#2 Mcloughlin_pbl Set2  Score Control Infected     3     2         0 0.2       0.2      ns          
#3 Mcloughlin_pbl Set3  Score Control Infected     3     2         0 0.2       0.2      ns          
#4 Mcloughlin     Set1  Score Control Infected     3    15         3 0.0172    0.0258   *           
#  5 Mcloughlin     Set2  Score Control Infected     3    15         0 0.00245   0.00441  **          
#  6 Mcloughlin     Set3  Score Control Infected     3    15         0 0.00245   0.00441  **          
#  7 OGrady         Set1  Score Control Infected    22    17        54 0.0000728 0.000218 ***         
#  8 OGrady         Set2  Score Control Infected    22    17        51 0.0000464 0.000209 ***         
#  9 OGrady         Set3  Score Control Infected    22    17        47 0.0000247 0.000209 ***

Test_GFS_thresholds_Set1 = roc(Condition ~ Score, data = Test_Scores_for_plotting_combined[Test_Scores_for_plotting_combined$Set == "Set1",]) 
Test_GFS_thresholds_Set2 = roc(Condition ~ Score, data = Test_Scores_for_plotting_combined[Test_Scores_for_plotting_combined$Set == "Set2",])
Test_GFS_thresholds_Set3 = roc(Condition ~ Score, data = Test_Scores_for_plotting_combined[Test_Scores_for_plotting_combined$Set == "Set3",]) 


Test_GFS_thresholds_Set1 <- cbind(Test_GFS_thresholds_Set1$thresholds, Test_GFS_thresholds_Set1$sensitivities, Test_GFS_thresholds_Set1$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set1")
Test_GFS_thresholds_Set2 <- cbind(Test_GFS_thresholds_Set2$thresholds, Test_GFS_thresholds_Set2$sensitivities, Test_GFS_thresholds_Set2$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set2")
Test_GFS_thresholds_Set3 <- cbind(Test_GFS_thresholds_Set3$thresholds, Test_GFS_thresholds_Set3$sensitivities, Test_GFS_thresholds_Set3$specificities) %>% as.data.frame() %>% set_colnames(c("Threshold", "Sensitivity", "Specificity")) %>% mutate(Set = "Set3")

Test_GFS_thresholds_Set1$Combined = rowSums(Test_GFS_thresholds_Set1[,c("Sensitivity", "Specificity")])
Test_GFS_thresholds_Set2$Combined = rowSums(Test_GFS_thresholds_Set2[,c("Sensitivity", "Specificity")])
Test_GFS_thresholds_Set3$Combined = rowSums(Test_GFS_thresholds_Set3[,c("Sensitivity", "Specificity")])

Test_GFS_thresholds_Set1

Test_GFS_thresholds_Set1 %>% filter(Sensitivity > 0.9) %>% filter(Combined == max(Combined)) # 1 -0.3459629   0.9117647   0.8529412 Set1 1.764706
Test_GFS_thresholds_Set2 %>% filter(Sensitivity > 0.9) %>% filter(Combined == max(Combined)) # -0.2379047   0.9411765   0.7647059 Set2 1.705882
Test_GFS_thresholds_Set3 %>% filter(Sensitivity > 0.9) %>% filter(Combined == max(Combined))#  -0.3740885   0.9411765   0.7647059 Set3 1.705882

threshold_set1 = -0.3459629
threshold_set2 = -0.2379047 
threshold_set3 = -0.3740885

thresholds_df

Test_Scores_for_plotting1 = Test_Scores_for_plotting1 %>% mutate(Prediction = if_else(Score > threshold_set1, "Infected", "Control"))
Test_Scores_for_plotting2 = Test_Scores_for_plotting2 %>% mutate(Prediction = if_else(Score > threshold_set2, "Infected", "Control"))
Test_Scores_for_plotting3 = Test_Scores_for_plotting3 %>% mutate(Prediction = if_else(Score > threshold_set3, "Infected", "Control"))


confusionMatrix(factor(Test_Scores_for_plotting1$Prediction, levels = c("Control", "Infected")), Test_Scores_for_plotting1$Condition, positive = "Infected")


thresholds_df <- as.matrix(cbind(Test_Scores_for_plotting1$Prediction, Test_Scores_for_plotting2$Prediction, Test_Scores_for_plotting3$Prediction))
thresholds_df
thresholds_df <- ifelse(thresholds_df == "Infected", 1, 0)
mode(thresholds_df) <- "numeric"
rownames(thresholds_df) <- Test_Scores_for_plotting1$Sample

test_set <- test_set %>% mutate(Infection_administration = if_else(Study == "Mcloughlin_pbl", "Natural", Infection_administration))
test_set$Study <- if_else(test_set$Study == "OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "OGrady", "O'Grady et al., (2025)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Wiarda", "Wiarda et al., (2020)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin", "McLoughlin et al., (2021)", test_set$Study)
test_set$Study <- if_else(test_set$Study == "Mcloughlin_pbl", "McLoughlin et al., (2014)", test_set$Study)

test_set$Study <- factor(test_set$Study, levels = c("O'Grady et al., (2025)", "McLoughlin et al., (2014)", "Wiarda et al., (2020)", "McLoughlin et al., (2021)"))
test_set <- test_set %>% arrange(Study)
test_set

thresholds_df <- thresholds_df[rownames(test_set),]

colnames(thresholds_df)
colnames(thresholds_df) <- c("Set1", "Set2", "Set3")
thresholds_df

ROC2

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



pdf("/home/workspace/jogrady/ML4TB/work/merged/figures/Heatmap_3_models.pdf", width = 15, height = 20)
pheatmap(t(thresholds_df),
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         gaps_row = c(seq(2:3)),
         annotation_col = annotation_columns,
         annotation_colors = ann_colors,
         border_color = "black",
         show_colnames = F,
         color = c("lightblue", "orange"))


dev.off()
thresholds_df
plot(1:5,1:5)
Test_Scores_for_plotting_combined
thr
GFS_test_ROCs = Test_Scores_for_plotting_combined  %>% group_by(Set) %>%
  summarise(roc_values = as.numeric(pROC::roc(
    response = ifelse(Condition == "Infected", 1, 0),
    predictor = Score,
    quiet = TRUE
  )$auc),
  AUC_CI_lower = (pROC::ci(pROC::roc(Condition, Score)))[1],
  AUC_CI_upper = (pROC::ci(pROC::roc(Condition, Score)))[3])

GFS_test_ROCs
roc_values <- GFS_test_ROCs %>% arrange(desc(roc_values))
roc_values
roc_text <- paste(paste0("ROC ", seq_along(roc_values$Set), " = ", roc_values$roc_values, "(95% CI:", round(roc_values$AUC_CI_lower, 2), " - ", round(roc_values$AUC_CI_upper,2), ")"), collapse = "\n")
roc_text
resample_colors <- ggsci::pal_npg("nrc")(length(roc_values$roc_values))

GFS_test_score = ggplot(Test_Scores_for_plotting_combined, aes(m=Score, d=factor(Condition, levels = c("Control", "Infected")), col = factor(Set, levels = c(roc_values$Set)))) + # Note do not need to filter unless saving all predictions
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
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"))



# Add annotations for each resample
for (i in seq_along(roc_values$Set)) {
  GFS_test_score  <- GFS_test_score  + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("ROC ", roc_values$Set[i], " = ", round(roc_values$roc_values[i],3)," (95% CI: ", round(roc_values$AUC_CI_lower[i], 2), " - ", round(roc_values$AUC_CI_upper[i],2), ")"), 
             hjust = 0, 
             size = 5, 
             color = resample_colors[i])
}


GFS_test_score
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/GFS_Test_ROC_curves.pdf", width = 12, height = 12, dpi = 600)

ROC <- ROC %>% filter(!is.na(Fold1)) %>% mutate(Iteration = 1:n(),
                                                Group = "Set1")
ROC2 <- ROC2 %>% filter(!is.na(Fold1)) %>% mutate(Iteration = 1:n(),
                                                  Group = "Set2")

ROC2



ROC_line = rbind(ROC, ROC2)
ROC_line$Group <- factor(ROC_line$Group)

ggplot(ROC_line, aes(
  x = as.numeric(Iteration),
  y = round(as.numeric(Average), 3),
  group = Group,
  color = Group
)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  theme_bw() +
  labs(
    x = "Iteration",
    y = "Average AUC",
    color = "Gene Set"
  ) + theme_bw() +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        #legend.position = c(.9, .25),
        #legend.justification = c("right", "top"),
        #legend.box.just = "right",
        #legend.margin = margin(6, 6, 6, 6),
        #legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14)) +
  scale_x_continuous(breaks = seq(1,18,2))#change legend title font size
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/ROC_iteration_curve.pdf", width = 15, height = 12, dpi = 600)


#########################################
# Lasso genes
#########################################


GLMLASSO_model <- readRDS("/home/workspace/jogrady/ML4TB/work/merged/Models/GLMLASSO_model.rds")
lasso_coeff = exp(coef(GLMLASSO_model$finalModel, GLMLASSO_model$bestTune$lambda))[,1]
lasso_coeff
lasso_coeff <- lasso_coeff[lasso_coeff > 1]
names(lasso_coeff)



Pos_genes = DE_results %>% filter(Symbol %in% names(lasso_coeff) & log2FoldChange > 0)
Neg_genes = DE_results %>% filter(Symbol %in% names(lasso_coeff) & log2FoldChange < 0)

Scores_for_plotting = ScoreGenesMtx(CPM_counts, Pos_genes$Symbol, Neg_genes$Symbol)
colnames(Scores_for_plotting) <- "Score"
Scores_for_plotting <- as.data.frame(Scores_for_plotting)
Scores_for_plotting$Sample <- rownames(Scores_for_plotting)
Scores_for_plotting <- left_join(Scores_for_plotting, metadata_all)
Scores_for_plotting$Study <- gsub("1_", "", Scores_for_plotting$Study)
Scores_for_plotting$Study <- gsub("2_", "", Scores_for_plotting$Study)
Scores_for_plotting$Study <- factor(Scores_for_plotting$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl"))
Scores_for_plotting$Study
Scores_for_plotting$Condition <- factor(Scores_for_plotting$Condition, levels = c("Control", "Infected"))
Scores_for_plotting$Condition
Scores_for_plotting$All_Score <- Scores_for_plotting$Score

class(Scores_for_plotting)
ggplot(Scores_for_plotting, aes(y = Score, x = factor(Study, levels = c("OGrady","Mcloughlin_pbl", "Wiarda", "Mcloughlin")), fill = Condition, shape = Condition)) + geom_boxplot(outlier.colour = NA, alpha = 0.5) + 
  geom_point(position = position_jitterdodge(jitter.width = 0.25), size = 3, aes(col = Condition)) + 
  scale_fill_manual(values = c("#2166ac", "#b2182b")) +
  scale_colour_manual(values = c("#2166ac", "#b2182b")) +
  theme_bw() +
  xlab("Study") +
  ylab("Score") +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.position = c(.9, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=2),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12)) + geom_pwc(lable = "p.adj", method = "wilcoxon", p.adjust.method = "bonferroni")
ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Lasso_score_training.pdf", width = 12, height = 12, dpi = 600)
Scores_for_plotting %>%
  group_by(Study) %>%
  wilcox_test(Score ~ Condition) %>%
  adjust_pvalue(method = "bonferroni") %>%  # Optional: adjust p-values
  add_significance() %>% arrange(desc(Study))   

GFS_train_roc_AUC <- Scores_for_plotting %>% group_by(Study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                           AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                           AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3])
Scores_for_plotting$Study <- as.character(Scores_for_plotting$Study)
Scores_for_plotting_long <- rbind(Scores_for_plotting, Scores_for_plotting)
Scores_for_plotting_long$Study <- c(Scores_for_plotting$Study, rep("Overall", 162))
Scores_for_plotting_long$Study <- factor(Scores_for_plotting_long$Study, levels = c("OGrady", "Wiarda", "Mcloughlin", "Mcloughlin_pbl", "Overall"))
Scores_for_plotting_wide = Scores_for_plotting_long %>% pivot_wider(names_from = "Study", values_from = "Score")

GFS_thresholds_overall = roc(Condition ~ Overall, data = Scores_for_plotting_wide)
GFS_thresholds_overall = cbind(GFS_thresholds_overall$thresholds, GFS_thresholds_overall$specificities, GFS_thresholds_overall$sensitivities) %>% data.frame() %>% set_colnames(., c("Threshold", "Spec", "Sens"))
GFS_thresholds_overall$Combined = rowSums(GFS_thresholds_overall[,c("Spec", "Sens")])

max(GFS_thresholds_overall$Combined)
#-0.018848162 0.98550725 (Specificity) 0.86021505 (Sensitivity) 1.845722

GFS_train_roc_AUC <- Scores_for_plotting_long %>% group_by(Study) %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Score)))),
                                                                                AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[1],
                                                                                AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Score))))[3])


GFS_train_score = ggplot(Scores_for_plotting_long, aes(m=Score, d=factor(Condition, levels = c("Control", "Infected")), colour = Study)) + # Note do not need to filter unless saving all predictions
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
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) + 
  scale_colour_manual(breaks = c("OGrady", "Mcloughlin", "Mcloughlin_pbl", "Wiarda", "Overall"), values = c("#351338","#734500","#5F927D","#DAD2FF", "black"))

GFS_train_score
resample_colors <-  c("#5F927D","#DAD2FF","black","#351338","#734500")
GFS_train_roc_AUC <- GFS_train_roc_AUC %>% arrange(desc(AUC))
GFS_train_roc_AUC
for (i in seq_along(GFS_train_roc_AUC$AUC)) {
  GFS_train_score <- GFS_train_score + 
    annotate("text", 
             x = 0.75, 
             y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
             label = paste0("AUROC ", GFS_train_roc_AUC$Study[i], " = ", round(GFS_train_roc_AUC$AUC[i],3), " (", round(GFS_train_roc_AUC$AUC_CI_lower[i],2), " - ", round(GFS_train_roc_AUC$AUC_CI_upper[i],2), ")"), 
             hjust = .25, 
             size = 5,
             col = resample_colors[i])
}


GFS_train_score

# Extract CI data (All Score only)
ci_vals <- ci.se(roc(Condition ~ All_Score, data = Scores_for_plotting_wide, legacy.axes = TRUE), specificities = seq(0, 1, 0.01))
df <- data.frame(
  specificity = as.numeric(rownames(ci_vals)),
  sensitivity = ci_vals[, "50%"],
  lower = ci_vals[, "2.5%"],
  upper = ci_vals[, "97.5%"],
  model = name)

GFS_train_score + 
  geom_ribbon(data = df, inherit.aes = FALSE, aes(x = 1 - specificity, ymin = lower, ymax = upper), alpha = 0.2, color = "grey") +
  annotate("point", x = 1 - 0.98550725, y = 0.86021505, colour = "darkred", size = 8, shape = 18)

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/Lasso_Train_ROC_curves.pdf", width = 12, height = 12, dpi = 600)

IGRA = fread("/home/workspace/jogrady/ML4TB/data/ogrady/Ogrady_IGRA_results.txt")
Scores_for_plotting <- left_join(Scores_for_plotting, IGRA, by = c("Sample" = "Sample ID"))

cor.test(Scores_for_plotting$Score, Scores_for_plotting$`IFGA_bovine - IFGA_avian`, method = "spearman", exact = FALSE)

ggplot(Scores_for_plotting, aes(x = Score, y = `IFGA_bovine - IFGA_avian`)) + geom_point(size = 3, aes(col = Condition)) +
  scale_color_manual(values = c("#2166ac", "#b2182b")) +
  theme_bw() +
  ylab("IFGA bovine - IFGA avian") +
  xlab("Score") +
  geom_smooth(method = "lm", col = "black") +
  geom_hline(yintercept = 80, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept =-0.018848162, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position = c(.9, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=2)) +
  annotate("text", x = -2, y = 2500,
           label = "paste(italic(ρ), \" = 0.724\")", parse = TRUE) +
  annotate("text", x = -2, y = 2200,
           label = "paste(italic(P), ' = 7.36 x 10^-15')", parse = TRUE)

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/IGRA_Lasso_TB_Score_correlation.pdf", width = 12, height = 12, dpi = 600)
Test_data = read.table("~/ML4TB/work/merged/Temp_files/test_CPM_filtered_counts.txt", sep = "\t", check.names=F)
Test_data_score = ScoreGenesMtx(GeneMtx = Test_data, pos.genes = Pos_genes$Symbol, neg.genes = Neg_genes$Symbol)
Test_data_score <- data.frame(Test_data_score)
Test_labels <- read.table("/home/workspace/jogrady/ML4TB/work/merged/Temp_files/Test_set_labels.txt")
Test_labels$Sample <- colnames(Test_data)
colnames(Test_labels)[1] <- "Condition"
colnames(Test_data_score)[1] <- "Score"
Test_data_score$Sample <- colnames(Test_data)
Test_data_score <- left_join(Test_data_score, Test_labels)
Test_data_score
Test_data_score

GFS_test_AUC = Test_data_score %>% summarize(AUC = as.numeric(pROC::auc(pROC::roc(Condition, as.numeric(Test_data_score$Score)))),
                                             AUC_CI_lower = (pROC::ci(pROC::roc(Condition, as.numeric(Test_data_score$Score))))[1],
                                             AUC_CI_upper = (pROC::ci(pROC::roc(Condition, as.numeric(Test_data_score$Score))))[3])

GFS_test_AUC
GFS_test_score = ggplot(Test_data_score, aes(m=Score, d=factor(Condition, levels = c("Control", "Infected")))) + # Note do not need to filter unless saving all predictions
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
        legend.position = "none",
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black")) +
  annotate("text", 
           x = 0.75, 
           y = 0.25 - (i - 1) * 0.025, # Adjust y position for each line
           label = paste0("AUROC",  " = ", round(GFS_test_AUC$AUC[1],2)," (", round(GFS_test_AUC$AUC_CI_lower[1],2), " - ", round(GFS_test_AUC$AUC_CI_upper[1],2), ")"), 
           hjust = .25, 
           size = 5)
ci_vals <- ci.se(roc(Condition ~ Score, data = Test_data_score, legacy.axes = TRUE), specificities = seq(0, 1, 0.01))
df <- data.frame(
  specificity = as.numeric(rownames(ci_vals)),
  sensitivity = ci_vals[, "50%"],
  lower = ci_vals[, "2.5%"],
  upper = ci_vals[, "97.5%"],
  model = name)

GFS_test_score = GFS_test_score + geom_ribbon(data = df, inherit.aes = FALSE, aes(x = 1 - specificity, ymin = lower, ymax = upper), alpha = 0.2, color = "grey")
GFS_test_score

ggsave("/home/workspace/jogrady/ML4TB/work/merged/figures/LASSO_Score_Test_ROC_curves.pdf", width = 12, height = 12, dpi = 600)


IGRA = fread("/home/workspace/jogrady/ML4TB/data/ogrady/Ogrady_IGRA_results.txt")
Test_data_score <- left_join(Test_data_score, IGRA, by = c("Sample" = "Sample ID"))

cor.test(Test_data_score$Score, Test_data_score$`IFGA_bovine - IFGA_avian`, method = "spearman", exact = FALSE)

ggplot(Test_data_score, aes(x = Score, y = `IFGA_bovine - IFGA_avian`)) + geom_point(size = 3, aes(col = Condition)) +
  scale_color_manual(values = c("#2166ac", "#b2182b")) +
  theme_bw() +
  ylab("IFGA bovine - IFGA avian") +
  xlab("Score") +
  geom_smooth(method = "lm", col = "black") +
  geom_hline(yintercept = 80, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept =-0.018848162, col = "darkgrey", linetype = "dashed", linewidth = 1) +
  theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
        axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        legend.key.size = unit(1.2, 'cm'), #change legend key size
        legend.key.height = unit(1.2, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=12),
        legend.position = c(.9, .25),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.box.background = element_rect(color="black", size=2)) +
  annotate("text", x = -1.5, y = 2000,
           label = "paste(italic(ρ), \" = 0.329\")", parse = TRUE) +
  annotate("text", x = -1.5, y = 1800,
           label = "paste(italic(P), ' = 0.04')", parse = TRUE)
