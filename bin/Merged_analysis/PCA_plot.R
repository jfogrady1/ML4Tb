# load tidyverse package
library(tidyverse)
library(pophelper)

args = commandArgs(trailingOnly = TRUE)


###################### Read in PCA data


# read in data
#args[1] and #args[2]
args[1] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/abdelaal_filtered_ALL_Pruned.eigenvec"
args[2] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/abdelaal_filtered_ALL_Pruned.eigenval"
args[3] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/abdelaal_PVE.pdf"
args[4] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/abdelaal_PCA.pdf"



# abdelaal
pca <- read_table(args[1], col_names = FALSE)
eigenval <- scan(args[2])
eigenval
# sort out the pca data
# set names
names(pca)[1] <- "ind"
pca$ind <- paste0(pca$ind, "_", pca$X2)
pca$ind
pca
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
colnames(pca)[2] <- "Status"
pca$Status <- rep(c(rep("Infected", 6), rep("Control", 6)))
pca
pve <- data.frame(PC = 1:(ncol(pca)-2), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve

ggsave(args[3], width = 15, height = 10, dpi = 600)


ggplot(data = pca, aes(x = PC1, y = PC2, colour = Status)) +
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
ggsave(args[4], width = 15, height = 10, dpi = 600)




args[5] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_filtered_ALL_Pruned.eigenvec"
args[6] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_filtered_ALL_Pruned.eigenval"
args[7] <- "/home/workspace/jogrady/ML4TB/data/wiarda/wiarda_samples.csv"
args[8] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_PVE.pdf"
args[9] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/wiarda_PCA.pdf"


# abdelaal
pca <- read_table(args[5], col_names = FALSE)
eigenval <- scan(args[6])
eigenval
# sort out the pca data
# set names
names(pca)[1] <- "ind"
colnames(pca)[2] <- "Status"
info = fread(args[7])
info <- info[order(info$Animal_ID),]
info <- info %>% filter(duplicated(Animal_ID))
info <- info %>% filter(duplicated(Animal_ID))


pca <- pca[order(pca$ind),]
pca$ind == info$Animal_ID
pca$Status <- info$Status
pca$ind <- paste0(pca$ind, "_", pca$Status)


pca
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))

pve <- data.frame(PC = 1:(ncol(pca)-2), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve

ggsave(args[8], width = 15, height = 10, dpi = 600)


ggplot(data = pca, aes(x = PC1, y = PC2, colour = Status)) +
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
ggsave(args[9], width = 15, height = 10, dpi = 600)



# Kirsten

# read in data
#args[1] and #args[2]
args[10] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_filtered_ALL_Pruned.eigenvec"
args[11] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_filtered_ALL_Pruned.eigenval"
args[12] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_PVE.pdf"
args[13] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_pbl_PCA.pdf"



# abdelaal
pca <- read_table(args[10], col_names = FALSE)
eigenval <- scan(args[11])
eigenval
# sort out the pca data
# set names
names(pca)[1] <- "ind"
pca$ind <- paste0(pca$ind, "_", pca$X2)
pca$ind
pca


names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
colnames(pca)[2] <- "Status"

pca$Status <- rep(c(rep("Control", 8), rep("Infected", 8)))
pca
pve <- data.frame(PC = 1:(ncol(pca)-2), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve
ggsave(args[12], width = 15, height = 10, dpi = 600)


ggplot(data = pca, aes(x = PC1, y = PC2, colour = Status)) +
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
ggsave(args[13], width = 15, height = 10, dpi = 600)



args[14] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_filtered_ALL_Pruned.eigenvec"
args[15] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_filtered_ALL_Pruned.eigenval"
args[16] <- "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt"
args[17] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_PVE.pdf"
args[18] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/kirsten_PCA.pdf"



# abdelaal
pca <- read_table(args[14], col_names = FALSE)
eigenval <- scan(args[15])
eigenval
pca
# sort out the pca data
# set names
names(pca)[1] <- "ind"
pca$ind <- paste0(pca$ind, "_", pca$X2)
pca$ind
pca
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
colnames(pca)[2] <- "Status"

info = fread(args[16])
info
info$Sample <- gsub("_W.*", "", info$Sample)
info
info <- info[order(info$Sample),]
info <- info[1:9,]


pca
pve <- data.frame(PC = 1:(ncol(pca)-2), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve
ggsave(args[17], width = 15, height = 10, dpi = 600)


ggplot(data = pca, aes(x = PC1, y = PC2, colour = Status)) +
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
ggsave(args[18], width = 15, height = 10, dpi = 600)



# Ogrady
args[19] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_filtered_ALL_Pruned.eigenvec"
args[20] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_filtered_ALL_Pruned.eigenval"
#args[16] <- "/home/workspace/jogrady/ML4TB/data/kirsten/kirsten_covariate.txt"
args[21] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_PVE.pdf"
args[22] <- "/home/workspace/jogrady/ML4TB/work/RNA_seq/vcf_isec/ogrady_PCA.pdf"



# abdelaal
pca <- read_table(args[19], col_names = FALSE)
eigenval <- scan(args[20])
eigenval
pca
# sort out the pca data
# set names
names(pca)[1] <- "ind"
pca$ind <- paste0(pca$ind, "_", pca$X2)
pca$ind
pca
names(pca)[3:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-2))
colnames(pca)[2] <- "Status"
pca$Status <- c(rep("Control", 63), rep("Infected",60))

pca
pve <- data.frame(PC = 1:(ncol(pca)-2), pve = eigenval/sum(eigenval)*100)
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_bw() +
  theme(axis.title  = element_text(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour = "black", size = 18),
        axis.text.x = element_text(colour = "black", size = 15),
        legend.title = element_text(size = 15, color = "black", face = "bold"),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(colour = "black", size = 15)) 
pve
ggsave(args[20], width = 15, height = 10, dpi = 600)

ggplot(data = pca, aes(x = PC1, y = PC2, colour = Status)) +
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

ggsave(args[21], width = 15, height = 10, dpi = 600)


# Ogrady
args[22] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenvec"
args[23] <- "/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/ADMIXTURE/SNP_Pruned.eigenval"

pca_WGS <- read_table(args[22], col_names = FALSE)
eigenval <- scan(args[23])
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

cor(pca_WGS$PC1, pca$PC1, method = "spearman")

cor(pca_WGS$PC2, pca$PC2, method = "spearman")



plot(x = pca$PC1, y = pca_WGS$PC1)
plot(x = pca$PC2, y = pca_WGS$PC2)


wilcox.test(pca$PC1, pca_WGS$PC1, paired = TRUE)
#Wilcoxon signed rank test with continuity correction

#data:  pca$PC1 and pca_WGS$PC1
#V = 4016, p-value = 0.6093
#alternative hypothesis: true location shift is not equal to 0

wilcox.test(pca$PC2, pca_WGS$PC2, paired = TRUE)

#Wilcoxon signed rank test with continuity correction

#data:  pca$PC2 and pca_WGS$PC2
#V = 3929, p-value = 0.7706
#alternative hypothesis: true location shift is not equal to 0

shapiro.test(pca_WGS$PC2)
pca
