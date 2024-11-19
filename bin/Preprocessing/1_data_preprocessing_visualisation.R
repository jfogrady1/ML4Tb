# This script reads in data, generates dds objects, performs vst transformation, visualiation using PCA, and writes, vst normalised counts and dds objects.
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(data.table)
library(ggthemes)

 

theme_Publication <- function(base_size=12, base_family = "sans") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(
                                         size = 16, hjust = 0.5, face = "bold"),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = "black", linewidth = 1.5),
               axis.title = element_text(face = "bold",size = 14),
               axis.title.y = element_text(angle=90,vjust =0, size = 14),
               axis.title.x = element_text(vjust = -0.2, size = 14),
               axis.text = element_text(size = 12), 
               axis.line = element_line(colour="white"),
               axis.ticks = element_line(),
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "right",
               legend.direction = "vertical",
               legend.title = element_text(size = 14, face = "bold"),
               legend.text = element_text(size = 14, face = "bold"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
          ))
      
}


ogrady_data_raw <- fread("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/RAW/count_matrix_clean.txt")
wiarda_data_raw <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/wiarda/Quantification/wiarda_count_matrix_clean.txt")
kirsten_data_raw <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten/Quantification/kirsten_count_matrix_clean.txt")
kirsten_pbl_raw <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/kirsten_pbl/Quantification/kirsten_pbl_count_matrix_clean.txt")
abdelaal_data_raw <- fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/abdelaal/Quantification/abdelaal_count_matrix_clean.txt")

# Covariate data
# Labels and wrangling
ogrady_covariate <- read.table("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/data/RNA_seq/covariate_RNA_seq.txt", header = TRUE)
wiarda_covariate <- fread("/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv")
kirsten_covariate <- fread("/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_samples.csv")
kirsten_covariate <- kirsten_covariate[!duplicated(kirsten_covariate$Animal_Code),] 
kirsten_covariate <- kirsten_covariate %>% select(3,4,5)
kirsten_covariate



kirsten_pbl_covariate <- fread("/home/workspace/jogrady/ML4TB/data/kirsten_pbl/kirsten_pbl_samples.csv")
kirsten_pbl_covariate



abdelaal_covariate <- fread("/home/workspace/jogrady/ML4TB/data/abdelaal/abdelaal_samples.csv", sep = "\t") 
abdelaal_covariate <- abdelaal_covariate[!duplicated(abdelaal_covariate$Animal_Code),] 
abdelaal_covariate <- abdelaal_covariate %>% select(3,4,5)
abdelaal_covariate





### Here we will focus on Wiarda
### INfer sex

row.names(wiarda_covariate) <- unique(wiarda_covariate$Animal_Code)
row.names(wiarda_data_raw) <- wiarda_data_raw$Geneid
wiarda_data_raw <- wiarda_data_raw %>% select(-1)



dds_wiarda_full <- DESeqDataSetFromMatrix(countData = wiarda_data_raw, 
                              colData = wiarda_covariate, 
                              design = ~ Status)  # Adjust design based on the outcome variable


dds_wiarda_full_vst <- vst(dds_wiarda_full)

saveRDS(dds_wiarda_full_vst, file = "/home/workspace/jogrady/ML4TB/work/normalisation/wiarda_dds_vst.rds")

dds_wiarda_full_vst[,dds_wiarda_full_vst$Week == "W0" | dds_wiarda_full_vst$Week == "W4"]



dds_wiarda_full_test <- readRDS("/home/workspace/jogrady/ML4TB/work/normalisation/wiarda_dds_vst.rds")



dds_wiarda_full_pca_data <- plotPCA(dds_wiarda_full_vst, intgroup=c("Status", "Week"), returnData=TRUE)
percentVar <- round(100 * attr(dds_wiarda_full_pca_data, "percentVar"))
percentVar
wiarda_all_samples_ggplot <- ggplot(dds_wiarda_full_pca_data, aes(PC1, PC2, color=Status, shape=Week)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("All samples")



dds_wiarda_full_pca_data_W0 <- plotPCA(dds_wiarda_full_vst[,dds_wiarda_full_vst$Week == "W0"], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_wiarda_full_pca_data_W0, "percentVar"))
percentVar
wiarda_W0_samples_ggplot <- ggplot(dds_wiarda_full_pca_data_W0, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 0")



dds_wiarda_full_pca_data_W4 <- plotPCA(dds_wiarda_full_vst[,dds_wiarda_full_vst$Week == "W4"], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_wiarda_full_pca_data_W4, "percentVar"))
percentVar
wiarda_W4_samples_ggplot <- ggplot(dds_wiarda_full_pca_data_W4, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 4")


dds_wiarda_full_pca_data_W10 <- plotPCA(dds_wiarda_full_vst[,dds_wiarda_full_vst$Week == "W10"], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_wiarda_full_pca_data_W10, "percentVar"))
percentVar
wiarda_W10_samples_ggplot <- ggplot(dds_wiarda_full_pca_data_W10, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))  + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) +scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 10")


library(patchwork)
library(cowplot)



wiarda_row <- plot_grid(wiarda_all_samples_ggplot + theme(legend.position="none"), wiarda_W0_samples_ggplot + theme(legend.position="none"), 
                        wiarda_W4_samples_ggplot + theme(legend.position="none"), wiarda_W10_samples_ggplot + theme(legend.position="none"), labels = c('A', 'B', 'C', 'D'), label_size = 20, nrow = 2, ncol = 2)
legend <- get_legend(
  # create some space to the left of the legend
  wiarda_all_samples_ggplot + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(wiarda_row, legend, rel_widths = c(1, .15))

ggsave2("/home/workspace/jogrady/ML4TB/work/normalisation/Figures/wiarda_PCA.pdf", width = 12, height = 12, dpi = 600)



### Here we will focus on Abdelaal
row.names(abdelaal_covariate) <- unique(abdelaal_covariate$Animal_Code)
row.names(abdelaal_data_raw) <- abdelaal_data_raw$Geneid
abdelaal_data_raw <- abdelaal_data_raw %>% select(-1)
abdelaal_covariate$Week
dds_abdelaal_full <- DESeqDataSetFromMatrix(countData = abdelaal_data_raw, 
                              colData = abdelaal_covariate, 
                              design = ~ Status)  # Adjust design based on the outcome variable


dds_abdelaal_full_vst <- vst(dds_abdelaal_full)
saveRDS(dds_wiarda_full_vst, file = "/home/workspace/jogrady/ML4TB/work/normalisation/abdelaal_dds_vst.rds")
dds_abdelaal_full_vst[,dds_abdelaal_full_vst$Week == "W0" | dds_abdelaal_full_vst$Week == "W4"]

dds_abdelaal_full




dds_abdelaal_full_pca_data <- plotPCA(dds_abdelaal_full_vst, intgroup=c("Status", "Week"), returnData=TRUE)
dds_abdelaal_full_pca_data
percentVar <- round(100 * attr(dds_abdelaal_full_pca_data, "percentVar"))
percentVar
abdelaal_all_samples_ggplot <- ggplot(dds_abdelaal_full_pca_data, aes(PC1, PC2, color=Status, shape=Week)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("All samples")
abdelaal_all_samples_ggplot



dds_abdelaal_full_pca_data_W8 <- plotPCA(dds_abdelaal_full_vst[,dds_abdelaal_full_vst$Week == "W8"], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_abdelaal_full_pca_data_W8, "percentVar"))
percentVar
abdelaal_W8_samples_ggplot <- ggplot(dds_abdelaal_full_pca_data_W8, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 8")


dds_abdelaal_full_pca_data_W20 <- plotPCA(dds_abdelaal_full_vst[,dds_abdelaal_full_vst$Week == "W20"], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_abdelaal_full_pca_data_W20, "percentVar"))
percentVar
abdelaal_W20_samples_ggplot <- ggplot(dds_abdelaal_full_pca_data_W20, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))  + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) +scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 20")




abdelaal_row <- plot_grid(abdelaal_all_samples_ggplot + theme(legend.position="none"), abdelaal_W8_samples_ggplot + theme(legend.position="none"), 
                        abdelaal_W20_samples_ggplot + theme(legend.position="none"), labels = c('A', 'B', 'C'), label_size = 20, nrow = 2, ncol = 2)
legend <- get_legend(
  # create some space to the left of the legend
  abdelaal_all_samples_ggplot + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(abdelaal_row, legend, rel_widths = c(1, .15))

ggsave2("/home/workspace/jogrady/ML4TB/work/normalisation/Figures/abdelaal_PCA.pdf", width = 12, height = 12, dpi = 600)





# here is kirsten_raw
### Here we will focus on Abdelaal
row.names(kirsten_covariate) <- unique(kirsten_covariate$Animal_Code)
row.names(kirsten_data_raw) <- kirsten_data_raw$Geneid
kirsten_data_raw <- kirsten_data_raw %>% select(-1)
kirsten_covariate$Week
dds_kirsten_full <- DESeqDataSetFromMatrix(countData = kirsten_data_raw, 
                              colData = kirsten_covariate, 
                              design = ~ Status)  # Adjust design based on the outcome variable


dds_kirsten_full_vst <- vst(dds_kirsten_full)
saveRDS(dds_wiarda_full_vst, file = "/home/workspace/jogrady/ML4TB/work/normalisation/kirsten_dds_vst.rds")
dds_kirsten_full




dds_kirsten_full_pca_data <- plotPCA(dds_kirsten_full_vst, intgroup=c("Status", "Week"), returnData=TRUE)
dds_kirsten_full_pca_data
percentVar <- round(100 * attr(dds_kirsten_full_pca_data, "percentVar"))
percentVar
kirsten_all_samples_ggplot <- ggplot(dds_kirsten_full_pca_data, aes(PC1, PC2, color=Status, shape=Week)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("All samples")



dds_kirsten_full_pca_data_W1 <- plotPCA(dds_kirsten_full_vst[,dds_kirsten_full_vst$Week %in% c("W-1","W1")], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_full_pca_data_W1, "percentVar"))
percentVar
kirsten_W1_samples_ggplot <- ggplot(dds_kirsten_full_pca_data_W1, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 1")
kirsten_W1_samples_ggplot

dds_kirsten_full_pca_data_W2 <- plotPCA(dds_kirsten_full_vst[,dds_kirsten_full_vst$Week %in% c("W-1","W2")], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_full_pca_data_W2, "percentVar"))
percentVar
kirsten_W2_samples_ggplot <- ggplot(dds_kirsten_full_pca_data_W2, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 2")


dds_kirsten_full_pca_data_W6 <- plotPCA(dds_kirsten_full_vst[,dds_kirsten_full_vst$Week %in% c("W-1","W6")], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_full_pca_data_W6, "percentVar"))
percentVar
kirsten_W6_samples_ggplot <- ggplot(dds_kirsten_full_pca_data_W6, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 6")




dds_kirsten_full_pca_data_W10 <- plotPCA(dds_kirsten_full_vst[,dds_kirsten_full_vst$Week %in% c("W-1","W10")], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_full_pca_data_W10, "percentVar"))
percentVar
kirsten_W10_samples_ggplot <- ggplot(dds_kirsten_full_pca_data_W10, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 10")




dds_kirsten_full_pca_data_W12 <- plotPCA(dds_kirsten_full_vst[,dds_kirsten_full_vst$Week %in% c("W-1","W12")], intgroup=c("Status"), returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_full_pca_data_W12, "percentVar"))
percentVar
kirsten_W12_samples_ggplot <- ggplot(dds_kirsten_full_pca_data_W12, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))  + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-25,25)) +scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("Week 12")




kirsten_row <- plot_grid(kirsten_all_samples_ggplot + theme(legend.position="none"), kirsten_W1_samples_ggplot + theme(legend.position="none"), kirsten_W2_samples_ggplot + theme(legend.position="none"),
                        kirsten_W6_samples_ggplot + theme(legend.position="none"), kirsten_W10_samples_ggplot + theme(legend.position="none"), kirsten_W12_samples_ggplot + theme(legend.position="none"), 
                        labels = c('A', 'B', 'C', 'D', 'E', 'F'), label_size = 20, nrow = 3, ncol = 2)
legend <- get_legend(
  # create some space to the left of the legend
  kirsten_all_samples_ggplot + theme(legend.box.margin = margin(0, 0, 0, 0))
)
plot_grid(kirsten_row, legend, rel_widths = c(1, .15))

ggsave2("/home/workspace/jogrady/ML4TB/work/normalisation/Figures/kirsten_PCA.pdf", width = 12, height = 12, dpi = 600)









### Here we will focus on kirsten_pbl
row.names(kirsten_pbl_covariate) <- unique(kirsten_pbl_covariate$Animal_Code)
row.names(kirsten_pbl_raw) <- kirsten_pbl_raw$Geneid
kirsten_pbl_raw <- kirsten_pbl_raw %>% select(-1)
kirsten_pbl_raw
dds_kirsten_pbl_full <- DESeqDataSetFromMatrix(countData = kirsten_pbl_raw, 
                              colData = kirsten_pbl_covariate, 
                              design = ~ Status)  # Adjust design based on the outcome variable


dds_kirsten_pbl_full_vst <- vst(dds_kirsten_pbl_full)
saveRDS(dds_wiarda_full_vst, file = "/home/workspace/jogrady/ML4TB/work/normalisation/kirsten_pbl_dds_vst.rds")
dds_kirsten_pbl_full




dds_kirsten_pbl_full_pca_data <- plotPCA(dds_kirsten_pbl_full_vst, intgroup=c("Status"), ntop = 500, returnData=TRUE)
percentVar <- round(100 * attr(dds_kirsten_pbl_full_pca_data, "percentVar"))
percentVar
kirsten_pbl_all_samples_ggplot <- ggplot(dds_kirsten_pbl_full_pca_data, aes(PC1, PC2, color=Status)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed(xlim = c(-30, 30), ylim = c(-20,20)) + scale_colour_manual(labels = c("Control", "Infected"), values = c("#386cb0", "#662506")) + theme_Publication() + ggtitle("All samples")
kirsten_pbl_all_samples_ggplot

ggsave2("/home/workspace/jogrady/ML4TB/work/normalisation/Figures/kirsten_pbl_PCA.pdf", width = 12, height = 12, dpi = 600)


##################################
##################################


# Write VST files
wiarda_counts_normalised  <- assay(dds_wiarda_full_vst)
kirsten_counts_normalised  <- assay(dds_kirsten_full_vst)
kirsten_pbl_counts_normalised  <- assay(dds_kirsten_pbl_full_vst)
abdelaal_counts_normalised  <- assay(dds_abdelaal_full_vst)




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

rownames(wiarda_counts_normalised)
tested_genes <- data.frame(rownames(wiarda_data_raw))
colnames(tested_genes) <- "gene_id"
head(tested_genes$gene_id,20)
head(ensemble$gene_name,20)

all(tested_genes$gene_id == ensemble$gene_id)
ensemble$gene_name <- gsub(' ', '', ensemble$gene_name)
ensemble$gene_id <- gsub(' ', '', ensemble$gene_id)
head(tested_genes$gene_id,20)
head(ensemble$gene_id,20)

tested_genes <- left_join(tested_genes, ensemble, by = c("gene_id" = "gene_id"))

head(wiarda_counts_normalised)
tested_genes$gene_id
rownames(wiarda_counts_normalised) <- tested_genes$gene_name
rownames(kirsten_counts_normalised) <- tested_genes$gene_name
rownames(kirsten_pbl_counts_normalised) <- tested_genes$gene_name
rownames(abdelaal_counts_normalised) <- tested_genes$gene_name


write.table(wiarda_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/wiarda_vst_normalised_data.txt", quote = FALSE, sep = "\t")
write.table(kirsten_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/kirsten_vst_normalised_data.txt", quote = FALSE, sep = "\t")
write.table(kirsten_pbl_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/kirsten_pbl_vst_normalised_data.txt", quote = FALSE, sep = "\t")
write.table(abdelaal_counts_normalised, "/home/workspace/jogrady/ML4TB/work/normalisation/vst_individual/abdelaal_vst_normalised_data.txt", quote = FALSE, sep = "\t")