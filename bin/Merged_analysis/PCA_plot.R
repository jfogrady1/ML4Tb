# load tidyverse package
library(tidyverse)
library(pophelper)
library(SNPRelate)
library(vcfR)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
args = commandArgs(trailingOnly = TRUE)


###################################
# Wiarda
###################################

#args[1] <- "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv"
#args[2] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_PVE.pdf" #output
#args[3] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_PCA.pdf" # output
#args[4] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_filtered_ALL_Pruned.vcf.gz"
#args[5] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda2.gds" # output

#args[6] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_PVE.pdf" # output
#args[7] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_PCA.pdf" #output
#args[8] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_filtered_ALL_Pruned.vcf.gz"
#args[9] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl2.gds" # output


#args[10] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_PVE.pdf" # output
#args[11] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_PCA.pdf" #output
#args[12] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_filtered_ALL_Pruned.vcf.gz"
#args[13] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten.gds" #output


#args[14] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_PVE.pdf" #output
#args[15] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_PCA.pdf" #output
#args[16] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_filtered_ALL_Pruned.vcf.gz"
#args[17] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady.gds" #output
#args[18] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec"
#args[19] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenval"
#args[20] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/PC1_RNA_SNP_correlation.pdf" #output - FigS02
#args[21] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/PC2_RNA_SNP_correlation.pdf" #output -FigS02
#args[23] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/PC1_RNA_SNP_wilcox.pdf" #output -FigS02
#args[24] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/PC2_RNA_SNP_wilcox.pdf" #output -FigS02


info = fread(args[1])
info <- info[order(info$Animal_ID),]
info <- info %>% filter(duplicated(Animal_ID))
info <- info %>% filter(duplicated(Animal_ID))


snpgdsVCF2GDS(args[4], args[5], method="biallelic.only")
snpgdsSummary(args[5])
genofile = snpgdsOpen(args[5])
pca <- snpgdsPCA(genofile, num.thread=2)



pc.percent <- pca$varprop*100


pc.percent

# abdelaal
pca$sample.id


pca$sample.id
# sort out the pca data
# set names
pca$sample.id <- gsub("_.*", "", pca$sample.id)

colnames(pca$eigenvect)[1:ncol(pca$eigenvect)] <- paste0("PC", 1:(ncol(pca$eigenvect)))
pca_wiarda = data.frame(pca$eigenvect)
pca_wiarda$Sample = as.vector(pca$sample.id)
pca$varprop
info <- info[match(pca$sample.id, info$Animal_ID), ]
sum(pca$varprop)
pca$Status <- info$Status
pca$sample.id

pve <- data.frame(PC = 1:(ncol(pca$eigenvect)), pve = pca$eigenval/sum(pca$eigenval)*100)
pve
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
pve$pve
a + ylab("Percentage variance explained") + theme_bw() +
  scale_x_continuous(breaks = c(1:13)) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
sum(pve$pve)
pve$pve
ggsave(args[2], width = 12, height = 12, dpi = 600)

pca$eigenvect
ggplot(data = data.frame(pca$eigenvect), aes(x = PC1, y = PC2, colour = pca$Status, label = pca$sample.id)) +
  geom_point(size = 3) +
  geom_text_repel() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_colour_manual(values = c("darkblue", "darkred")) +
  labs(colour = "Infection Status", shape = "Group", size = "Pedigree Holstein %") + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                        axis.title.x = element_text(colour = "black", size = 18),
                                                                                                        axis.title.y = element_text(colour = "black", size = 18),
                                                                                                        axis.text.x = element_text(colour = "black", size = 15),
                                                                                                        legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                        legend.text = element_text(size = 15),
                                                                                                        axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))
ggsave(args[3], width = 12, height = 12, dpi = 600)





###################################
# Kirsten pbl
###################################
snpgdsVCF2GDS(args[8], args[9], method="biallelic.only")
snpgdsSummary(args[9])
genofile = snpgdsOpen(args[9])
pca <- snpgdsPCA(genofile, num.thread=2)


pca$sample.id <- gsub("_.*", "", pca$sample.id)


colnames(pca$eigenvect)[1:ncol(pca$eigenvect)] <- paste0("PC", 1:(ncol(pca$eigenvect)))
pca_kirsten_pbl = data.frame(pca$eigenvect)
pca_kirsten_pbl$Sample = as.vector(pca$sample.id)
pca$Status <- c(rep("Control", 8), rep("Infected", 8))
pca$Status
pve <- data.frame(PC = 1:(ncol(pca$eigenvect)), pve = pca$eigenval/sum(pca$eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  scale_x_continuous(breaks = c(1:16)) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve
ggsave(args[6], width = 12, height = 12, dpi = 600)


ggplot(data = data.frame(pca$eigenvect), aes(x = PC1, y = PC2, colour = pca$Status, label = pca$sample.id)) +
  geom_label_repel() +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_colour_manual(values = c("darkblue", "darkred")) +
  labs(colour = "Infection Status", shape = "Group", size = "Pedigree Holstein %") + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                        axis.title.x = element_text(colour = "black", size = 18),
                                                                                                        axis.title.y = element_text(colour = "black", size = 18),
                                                                                                        axis.text.x = element_text(colour = "black", size = 15),
                                                                                                        legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                        legend.text = element_text(size = 15),
                                                                                                        axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))
ggsave(args[7], width = 12, height = 12, dpi = 600)



###################################
# Kirsten time series
###################################
snpgdsVCF2GDS(args[12], args[13], method="biallelic.only")
snpgdsSummary(args[13])
genofile = snpgdsOpen(args[13])
pca <- snpgdsPCA(genofile, num.thread=2)



pca$sample.id <- gsub("_.*", "", pca$sample.id)

colnames(pca$eigenvect)[1:ncol(pca$eigenvect)] <- paste0("PC", 1:(ncol(pca$eigenvect)))
pca_kirsten = data.frame(pca$eigenvect)
pca_kirsten$Sample = as.vector(pca$sample.id)
pca$eigenvect
pve <- data.frame(PC = 1:(ncol(pca$eigenvect)), pve = pca$eigenval/sum(pca$eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  scale_x_continuous(breaks = c(1:9)) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve
ggsave(args[10], width = 12, height = 12, dpi = 600)


ggplot(data = data.frame(pca$eigenvect), aes(x = PC1, y = PC2, col = pca$sample.id, label = pca$sample.id)) +
  geom_label_repel() +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  labs(colour = "Animal ID", shape = "Group", size = "Pedigree Holstein %") + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                        axis.title.x = element_text(colour = "black", size = 18),
                                                                                                        axis.title.y = element_text(colour = "black", size = 18),
                                                                                                        axis.text.x = element_text(colour = "black", size = 15),
                                                                                                        legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                        legend.text = element_text(size = 15),
                                                                                                        axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))
ggsave(args[11], width = 12, height = 12, dpi = 600)



###################################
# O'Grady
###################################



snpgdsVCF2GDS(args[16], args[17], method="biallelic.only")
snpgdsSummary(args[17])
genofile = snpgdsOpen(args[17])
pca <- snpgdsPCA(genofile, num.thread=2)


pca$sample.id <- gsub("_.*", "", pca$sample.id)

pca$sample.id

colnames(pca$eigenvect)[1:ncol(pca$eigenvect)] <- paste0("PC", 1:(ncol(pca$eigenvect)))
pca_ogrady = data.frame(pca$eigenvect)
pca_ogrady$Sample = as.vector(pca$sample.id)
pca$Status <- c(rep("Control", 63), rep("Infected",60))

pca
pve <- data.frame(PC = 1:20, pve = pca$eigenval[1:20]/sum(pca$eigenval[1:20])*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  scale_x_continuous(breaks = c(seq(1,20,2))) +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 

ggsave(args[14], width = 12, height = 12, dpi = 600)

ggplot(data = data.frame(pca$eigenvect), aes(x = PC1 * -1, y = PC2, colour = pca$Status)) +
  geom_point(size = 3) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  scale_colour_manual(values = c("darkblue", "darkred")) +
  labs(colour = "Infection Status", shape = "Group", size = "Pedigree Holstein %") + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                        axis.title.x = element_text(colour = "black", size = 18),
                                                                                                        axis.title.y = element_text(colour = "black", size = 18),
                                                                                                        axis.text.x = element_text(colour = "black", size = 15),
                                                                                                        legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                        legend.text = element_text(size = 15),
                                                                                                        axis.text.y = element_text(colour = "black", size = 15)) +
  guides(shape = guide_legend(override.aes = list(size = 4)))

ggsave(args[15], width = 12, height = 12, dpi = 600)



pca_WGS <- read_table(args[18], col_names = FALSE)
eigenval <- scan(args[19])
eigenval
pca_WGS
# sort out the pca data
# set names
names(pca_WGS)[1] <- "ind"
pca_WGS$ind <- gsub("_.*", "", pca_WGS$X2)
pca_WGS
names(pca_WGS)[3:ncol(pca_WGS)] <- paste0("PC", 1:(ncol(pca_WGS)-2))
colnames(pca_WGS)[2] <- "Status"
pca_WGS$Status <- c(rep("Control", 63), rep("Infected",60))
pca_WGS





#################################
# Correlation between RNA called SNPs and SNP-array derived SNPs
#################################

PCA_comparison = data.frame(cbind(pca_WGS$ind, pca_WGS$PC1, pca_WGS$PC2, as.numeric(data.frame(pca$eigenvect)$PC1), as.numeric(data.frame(pca$eigenvect)$PC2)))


colnames(PCA_comparison) <- c("Sample", "WGS_PC1", "WGS_PC2", "RNA_PC1", "RNA_PC2")

PCA_comparison$RNA_PC1 <- as.numeric(PCA_comparison$RNA_PC1)
PCA_comparison$RNA_PC2 <- as.numeric(PCA_comparison$RNA_PC2)

PC1 = cor.test(as.numeric(PCA_comparison$WGS_PC1), (PCA_comparison$RNA_PC1* -1), method = "spearman", exact = FALSE)
PC2 = cor.test(as.numeric(PCA_comparison$WGS_PC2), (PCA_comparison$RNA_PC2), method = "spearman", exact = FALSE)
PC1$estimate


PC2$p.value
PC2$estimate

ggplot(PCA_comparison, aes(x = as.numeric(WGS_PC1), y = (RNA_PC1 * -1))) + labs(x = "Eigenvector 1 derived from SNP-array data",
                                                                                y = "Eigenvector 1 derived from RNA-called variants") + geom_point(size = 3) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                  legend.text = element_text(size = 15),
                                                                                                  axis.text.y = element_text(colour = "black", size = 15)) +
  annotate("text", label = PC1$p.value, x = 0, y = 0) +
  annotate("text", label = PC1$estimate, x = 0, y = 0)
ggsave(args[20], width = 12, height = 12, dpi = 600)


ggplot(PCA_comparison, aes(x = as.numeric(WGS_PC2), y = (RNA_PC2))) + labs(x = "Eigenvector 2 derived from SNP-array data",
                                                                                y = "Eigenvector 2 derived from RNA-called variants") + geom_point(size = 3) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                  legend.text = element_text(size = 15),
                                                                                                  axis.text.y = element_text(colour = "black", size = 15)) +
  annotate("text", label = PC2$p.value, x = 0, y = 0) +
  annotate("text", label = PC2$estimate, x = 0, y = 0)
ggsave(args[21], width = 12, height = 12, dpi = 600)


PCA_comparison$WGS_PC1 <- as.numeric(PCA_comparison$WGS_PC1)
PCA_comparison$WGS_PC2 <- as.numeric(PCA_comparison$WGS_PC2)
PCA_comparison$RNA_PC1 <- PCA_comparison$RNA_PC1 * -1
PCA_long <- PCA_comparison %>%
  pivot_longer(
    cols = -Sample,
    names_to = c("Type", "PC"),
    names_sep = "_",
    values_to = "Value"
  )


PCA_long %>% filter(PC == "PC1") %>% ggplot(., aes(x = Type, y = Value,fill = Type)) + geom_boxplot(outlier.colour = NA, alpha = 0.5) + geom_jitter(width = 0.12, size = 3, col = "black", alpha = 0.5) + stat_compare_means(method = "wilcox", paired = TRUE) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                                                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                                                                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                                                                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                                                                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                                                                                                                  legend.text = element_text(size = 15),
                                                                                                                                                                                                  axis.text.y = element_text(colour = "black", size = 15),
                                                                                                                                                                                                  legend.position = "none") +
  labs(x = "Group", y = "PC1 coordinates")
ggsave(args[22], width = 12, height = 12, dpi = 600)
PCA_long %>% filter(PC == "PC2") %>% ggplot(., aes(x = Type, y = Value, fill = Type)) + geom_boxplot(outlier.colour = NA, alpha = 0.5) + geom_jitter(width = 0.12,size = 3, col = "black", alpha = 0.5) + stat_compare_means(method = "wilcox", paired = TRUE) + theme_bw() + theme(axis.title  = element_text(colour = "black"),
                                                                                                                                                                                                  axis.title.x = element_text(colour = "black", size = 18),
                                                                                                                                                                                                  axis.title.y = element_text(colour = "black", size = 18),
                                                                                                                                                                                                  axis.text.x = element_text(colour = "black", size = 15),
                                                                                                                                                                                                  legend.title = element_text(size = 15, color = "black", face = "bold"),
                                                                                                                                                                                                  legend.text = element_text(size = 15),
                                                                                                                                                                                                  axis.text.y = element_text(colour = "black", size = 15),
                                                                                                                                                                                                  legend.position = "none") +
  labs(x = "Group", y = "PC2 coordinates")
ggsave(args[23], width = 12, height = 12, dpi = 600)

head(data.frame(pca_wiarda))
head(pca_kirsten)
head(pca_kirsten_pbl)
head(pca_ogrady)
print("HERE")

write.table(pca_wiarda, file = args[24], sep = "\t", quote = FALSE)
write.table(pca_kirsten, file = args[25], sep = "\t", quote = FALSE)
write.table(pca_kirsten_pbl, file = args[26], sep = "\t", quote = FALSE)
write.table(pca_ogrady, file = args[27], sep = "\t", quote = FALSE)

showfile.gds(closeall=TRUE)
