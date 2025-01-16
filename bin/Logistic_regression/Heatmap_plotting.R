library(tidyverse)





# Abdellal, Wiarda, Mcloughlin 2021, Mcloughlin 2014, O'Grady test

gene_based = matrix(nrow = 5, ncol = 5)

# abdelaal
gene_based[1,2] = 0.720
gene_based[1,3] = 0.307
gene_based[1,4] = 0.078
gene_based[1,5] = 0.448

# Wiarda
gene_based[2,1] = 0.674
gene_based[2,3] = 0.721
gene_based[2,4] = 0.344
gene_based[2,5] = 0.540

# Kirsten
gene_based[3,1] = 0.389
gene_based[3,2] = 1.00
gene_based[3,4] = 1.00
gene_based[3,5] = 0.559


# Kirsten_pbl
gene_based[4,1] = 0.420
gene_based[4,2] = 1.00
gene_based[4,3] = 0.811
gene_based[4,5] = 0.420

gene_based[5,1] = 0.576
gene_based[5,2] = 0.861
gene_based[5,3] = 0.757
gene_based[5,4] = 0.875



####################################

## PCA
####################################

# Abdellal, Wiarda, Mcloughlin 2021, Mcloughlin 2014, O'Grady test
PCA_based = matrix(nrow = 5, ncol = 5)
# abdelaal
PCA_based[1,2] = 0.755
PCA_based[1,3] = 0.261
PCA_based[1,4] = 0.031
PCA_based[1,5] = 0.321

#wiarda
PCA_based[2,1] = 0.847
PCA_based[2,3] = 0.429
PCA_based[2,4] = 0.062
PCA_based[2,5] = 0.377

#Kirsten
PCA_based[3,1] = 0.097
PCA_based[3,2] = 0.972
PCA_based[3,4] = 1.00
PCA_based[3,5] = 0.707

# Kirsten_pbl
PCA_based[4,1] = 0.347
PCA_based[4,2] = 1.00
PCA_based[4,3] = 0.747
PCA_based[4,5] = 0.725

# O'grady
PCA_based[5,1] = 0.333
PCA_based[5,2] = 0.940
PCA_based[5,3] = 0.615
PCA_based[5,4] = 0.766


####################################
## ICA
####################################
# Abdellal, Wiarda, Mcloughlin 2021, Mcloughlin 2014, O'Grady test
ica_based = matrix(nrow = 5, ncol = 5)

# abdelaal
ica_based[1,2] = 0.742
ica_based[1,3] = 0.261
ica_based[1,4] = 0.016
ica_based[1,5] = 0.333

# wiarda
ica_based[2,1] = 0.847
ica_based[2,3] = 0.426
ica_based[2,4] = 0.062
ica_based[2,5] = 0.377

ica_based[3,1] = 0.111
ica_based[3,2] = 0.966
ica_based[3,4] = 1.00
ica_based[3,5] = 0.722

# Kirsten_pbl
ica_based[4,1] = 0.347
ica_based[4,2] = 1.00
ica_based[4,3] = 0.747
ica_based[4,5] = 0.725

# O'grady
ica_based[5,1] = 0.347
ica_based[5,2] = 0.908
ica_based[5,3] = 0.620
ica_based[5,4] = 0.656


####################################

## NMF
####################################


nmf_based = matrix(nrow = 5, ncol = 5)

# abdelaal
nmf_based[1,2] = 0.856
nmf_based[1,3] = 0.349
nmf_based[1,4] = 0.047
nmf_based[1,5] = 0.315


# wiarda
nmf_based[2,1] = 0.840
nmf_based[2,3] = 0.426
nmf_based[2,4] = 0.062
nmf_based[2,5] = 0.377

# Kirsten
nmf_based[3,1] = 0.118
nmf_based[3,2] = 0.969
nmf_based[3,4] = 1.00
nmf_based[3,5] = 0.719


# Kirsten_pbl
nmf_based[4,1] = 0.340 
nmf_based[4,2] = 1.00
nmf_based[4,3] = 0.796
nmf_based[4,5] = 0.731

# O'Grady
nmf_based[5,1] = 0.431
nmf_based[5,2] = 0.889
nmf_based[5,3] = 0.486
nmf_based[5,4] = 0.781



# Load required library
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
library(pheatmap)


library(RColorBrewer)
# Assign row and column names
rownames(gene_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))
colnames(gene_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))

rownames(PCA_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))
colnames(PCA_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))

rownames(ica_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))
colnames(ica_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))


rownames(nmf_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))
colnames(nmf_based) <- factor(c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"), levels = c("abdelaal", "wiarda", "Kirsten", "Kirsten_pbl", "O'Grady"))



save_pheatmap_pdf <- function(x, filename, width=12, height=12) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}


breaksList = seq(0, 1, by = 0.01)

NMF = pheatmap(nmf_based,              # Scale rows for better compariso,  # Color gradient
         cluster_rows = FALSE,                 # Disable row clustering
         cluster_cols = FALSE,                 # Disable column clustering
         display_numbers = T, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         show_colnames = F, show_rownames = F,# Optionally display numbers in cells           # Title of the heatmap
         fontsize_row = 10,                    # Font size for row labels
         fontsize_col = 10,                    # Font size for column labels
         angle_col = 0,                        # Labels are horizontal at the top
         border_color = "black")  # Title of the heatmap
save_pheatmap_pdf(NMF, "NMF_heatmap.pdf")

Gene = pheatmap(gene_based,              # Scale rows for better compariso,  # Color gradient
         cluster_rows = FALSE,                 # Disable row clustering
         cluster_cols = FALSE,                 # Disable column clustering
         display_numbers = T, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         show_colnames = F, show_rownames = F,# Optionally display numbers in cells           # Title of the heatmap
         fontsize_row = 10,                    # Font size for row labels
         fontsize_col = 10,                    # Font size for column labels
         angle_col = 0,                        # Labels are horizontal at the top
         border_color = "black")  # Title of the heatmap
save_pheatmap_pdf(Gene, "Gene_heatmap.pdf")

PCA = pheatmap(PCA_based,              # Scale rows for better compariso,  # Color gradient
         cluster_rows = FALSE,                 # Disable row clustering
         cluster_cols = FALSE,                 # Disable column clustering
         display_numbers = T, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         show_colnames = F, show_rownames = F,# Optionally display numbers in cells           # Title of the heatmap
         fontsize_row = 10,                    # Font size for row labels
         fontsize_col = 10,                    # Font size for column labels
         angle_col = 0,                        # Labels are horizontal at the top
         border_color = "black")  # Title of the heatmap
save_pheatmap_pdf(PCA, "PCA_heatmap.pdf")


ICA = pheatmap(ica_based,              # Scale rows for better compariso,  # Color gradient
         cluster_rows = FALSE,                 # Disable row clustering
         cluster_cols = FALSE,                 # Disable column clustering
         display_numbers = T, 
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)),
         show_colnames = F, show_rownames = F,# Optionally display numbers in cells           # Title of the heatmap
         fontsize_row = 10,                    # Font size for row labels
         fontsize_col = 10,                    # Font size for column labels
         angle_col = 0,                        # Labels are horizontal at the top
         border_color = "black")  # Title of the heatmap
save_pheatmap_pdf(ICA, "ICA_heatmap.pdf")
