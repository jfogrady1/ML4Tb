# Script to assess concordance at SNP sites

library(data.table)
library(tidyverse)
DNA = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/DNA_ref_alt.txt") %>% rename(., CHROM = V1,  POS = V2, REF = V3, ALT = V4) %>% mutate(LOC = paste0(CHROM, ":", POS))
RNA = fread("/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/RNA_ref_alt.txt") %>% rename(., CHROM_RNA = V1,  POS_RNA = V2, REF_RNA = V3, ALT_RNA = V4) %>% mutate(LOC = paste0(CHROM_RNA, ":", POS_RNA))

Combined = inner_join(DNA, RNA, by = c("LOC" = "LOC"))
# dim =  277302      9



Combined <- Combined %>% mutate(Concordance = case_when(
  REF == REF_RNA & ALT == ALT_RNA ~ "Concordant",
  .default = "Discordant"
))


Combined_no_triallelic <- Combined %>% filter(ALT_RNA == "A" | ALT_RNA == "T" | ALT_RNA == "C" | ALT_RNA == "G")
Combined_no_triallelic$LOC <- paste0(Combined_no_triallelic$CHROM, ":", Combined_no_triallelic$POS)

table(Combined_no_triallelic$Concordance)

get_concordance_proportions <- function(data) {
  data %>%
    count(Concordance) %>%
    mutate(Percentage = round(((n / sum(n))*100),2))
}

# Get proportions for both datasets
concordance_with_triallelic <- get_concordance_proportions(Combined)
concordance_no_triallelic <- get_concordance_proportions(Combined_no_triallelic)

concordance_with_triallelic
#Concordance      n Percentage
#1:  Concordant 275039      99.18
#2:  Discordant   2263       0.82


concordance_no_triallelic
#Concordance      n Percentage
#1:  Concordant 275039      99.98
#2:  Discordant     54       0.02

plot_concordance <- function(data, title) {
  ggplot(data, aes(x = factor(Concordance), y = Percentage, fill = Concordance)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Percentage), vjust = -0.5, size = 5) +  # Add labels above bars
    scale_y_continuous(limits = c(0, 100)) +  # Set y-axis limits
    labs(x = "Concordance Status", y = "Percentage (%)") +
    theme_bw() +
    scale_fill_manual(values = c("Concordant" = "steelblue", "Discordant" = "firebrick"))
}


# Generate both plots
plot1 <- plot_concordance(concordance_with_triallelic) + theme(legend.position = "none")
plot2 <- plot_concordance(concordance_no_triallelic) + theme(legend.position = c(0.75,0.85))

# Display plots side by side
cowplot::plot_grid(plot1, plot2)


library(vcfR)

DNA = read.vcfR("/home/workspace/jogrady/eqtl_study/eqtl_nextflow/results/SNP_data/imputation/SNP_Data_Final.vcf.gz")
RNA = read.vcfR("/home/workspace/jogrady/ML4TB/work/RNA_seq/ogrady/variant_call/ogrady_RNA_merged.vcf.gz")


DNA@fix
DNA <- cbind(DNA@fix, DNA@gt)
RNA <- cbind(RNA@fix, RNA@gt)

DNA <- as.data.frame(DNA)
RNA <- as.data.frame(RNA)
head(DNA)

DNA$ID <- paste0(DNA$CHROM, ":", DNA$POS)
RNA$ID <- paste0(RNA$CHROM, ":", RNA$POS)

DNA_filter <- DNA %>% filter(ID %in% Combined_no_triallelic$LOC)
RNA_filter <- RNA %>% filter(ID %in% Combined_no_triallelic$LOC)




rm(DNA)
rm(RNA)

get_gt_concordance_proportions <- function(data) {
  data %>%
    count(Match) %>%
    mutate(Percentage = round(((n / sum(n))*100),2)) %>% filter(Match == "MATCH")
}



get_gt_concordance_proportions <- function(data) {
  data %>%
    count(Match) %>%
    mutate(Percentage = round(((n / sum(n))*100),2)) %>% filter(Match == "MATCH") %>% select(Percentage) %>% as.vector()
}


final_df <- data.frame(matrix(ncol = 3, nrow = 1))

colnames(final_df) <- c("Sample", "N_snps", "Percentage")

final_df
# DNA
for (i in 10:ncol(DNA_filter)) {
  #print(i)

  temp <- DNA_filter[,c(1:9,i)]
  temp[,10] <- gsub(":.*", "", temp[,10])
  temp[,10] <- if_else(temp[,10] == "1/0", "1", temp[,10])
  temp[,10] <- if_else(temp[,10] == "0/1", "1", temp[,10])
  temp[,10] <- if_else(temp[,10] == "0/0", "0", temp[,10])
  temp[,10] <- if_else(temp[,10] == "1/1", "2", temp[,10])
  temp <- data.frame(temp)
  # RNA
  temp_RNA <- RNA_filter[,c(1:9,i)]
  temp_RNA[,10] <- gsub(":.*", "", temp_RNA[,10])
  temp_RNA[,10] <- if_else(temp_RNA[,10] == "./.", NA, temp_RNA[,10])
  temp_RNA[,10] <- if_else(temp_RNA[,10] == "1/0", "1", temp_RNA[,10])
  temp_RNA[,10] <- if_else(temp_RNA[,10] == "0/1", "1", temp_RNA[,10])
  temp_RNA[,10] <- if_else(temp_RNA[,10] == "0/0", "0", temp_RNA[,10])
  temp_RNA[,10] <- if_else(temp_RNA[,10] == "1/1", "2", temp_RNA[,10])
  
  
  
  temp_RNA <- data.frame(temp_RNA[!is.na(temp_RNA[,10]),])
  temp <- temp %>% filter(ID %in% temp_RNA$ID)
  temp_RNA <- temp_RNA %>% filter(ID %in% temp$ID)
  temp <- temp %>% arrange(ID) 
  temp_RNA <- temp_RNA %>% arrange(ID)
  all(temp$ID == temp_RNA$ID)
  
  RNA_final <- temp_RNA[,10]
  
  temp <- cbind(temp, RNA_final)
  
  temp$Match <- if_else(temp[,10] == temp[,11], "MATCH", "NOMATCH")
  
  print(table(temp[,10]))
  N = dim(temp_RNA)[1]
  Sample <- colnames(temp)[10]
  Percentage <- get_gt_concordance_proportions(temp)
  vec = c(Sample, N, Percentage$Percentage)
  names(vec) <- colnames(final_df)
  
  final_df <- rbind(final_df, vec)
  final_df
}
final_df


T064 = read.vcfR("~/eqtl_study/work/RNA-SNP-Calling/CONCORDANCE/selected.T064.bam.filtered.vcf.gz")

T064 <- as.data.frame(cbind(T064@fix, T064@gt))
# RNA
T064_RNA <- T064[,c(1:9,10)]
T064_RNA
T064_RNA[,10] <- gsub(":.*", "", T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "./.", NA, T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "1/0", "1", T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "0/1", "1", T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "0/0", "0", T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "1/1", "2", T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "1/2", NA, T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "2/2", NA, T064_RNA[,10])
T064_RNA[,10] <- if_else(T064_RNA[,10] == "0/2", NA, T064_RNA[,10])

table(T064_RNA[,10])


T064_RNA <- T064_RNA %>% mutate(ID = paste0(T064_RNA$CHROM, ":", T064_RNA$POS))

head(T064_RNA)
head(temp_RNA)
T064_RNA <- data.frame(T064_RNA[!is.na(T064_RNA[,10]),])

temp <- temp %>% filter(ID %in% T064_RNA$ID)
T064_RNA <- T064_RNA %>% filter(ID %in% temp$ID)


temp <- temp %>% arrange(ID) 
T064_RNA <- T064_RNA %>% arrange(ID)
all(T064_RNA$ID == temp$ID)

table(temp_RNA$T064)
table(temp$T064)
table(T064_RNA[,10])

dim(temp)
dim(T064_RNA)


T064_RNA <- T064_RNA
table(T064_RNA[,10])
temp
head(temp_RNA, 30)

