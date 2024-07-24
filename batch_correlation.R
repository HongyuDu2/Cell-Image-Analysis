library(harmony)
library(umap)
library(ggplot2)

# Import dataset
data1 <- read.csv("C:/Users/duho/Desktop/cell_table_size_normalized.csv")
data2 <- read.csv("C:/Users/duho/Desktop/cell_table_size_normalized (1).csv")
data3 <- read.csv("C:/Users/duho/Desktop/cell_table_size_normalized (2).csv")

# Add batch indicators
data1$batch <- 'Fuentes'
data2$batch <- 'SLE39'
data3$batch <- 'SLE40'

new_names <- c("B-Tubulin", "C3d29","CC3", "CD11b",
               "CD11c", "CD14", "CD16", "CD163", "CD20",
               "CD21", "CD3", "CD31", "CD38", "CD4", "CD45",
               "CD45RO", "CD56", "CD68", "CD8", "CXCR5",
               "Col III", "DC-SIGN", "Factor H", "Foxp3",
               "Granzyme B", "HLA Class 1", "HLADR", "IgA",
               "IgG-PE", "IgM", "Ki67", "MPO", "PAN-KERATIN",
               "PD-1", "PE", "Podoplanin", "SMA", "Vimentin", "dsDNA")

# change name
colnames(data1)[2:40] <- new_names
colnames(data2)[2:40] <- new_names
colnames(data3)[2:40] <- new_names

combined <-  rbind(data1, data2, data3)

# only with markers
combined_markers <- combined[2:40]

# Run PCA and Harmony for batch correction
combined_markers <- scale(combined_markers)

# Perform PCA
pca_result <- prcomp(combined_markers, center = TRUE, scale. = TRUE)
npcs <- 30 # set the number of PCA used
pca_scores <- pca_result$x[, 1:npcs]
pca_loadings <- pca_result$rotation[,1:npcs] # rotation matrix

# Prepare data for Harmony
combined_batches <- combined$batch

# Apply Harmony for batch correction (HarmonyMatrix function is used here for demonstration)
# HarmonyMatrix requires the `harmony` package, make sure it's installed
harmony_result <- HarmonyMatrix(data = pca_scores, meta_data = combined_batches, vars_use = "batch")

# Go back to original vector space
batch_corrected_original_space <- harmony_result %*% t(pca_loadings)

# Convert back to dataframe
batch_corrected_df <- as.data.frame(batch_corrected_original_space)
batch_corrected_df$batch <- combined$batch
write.csv(batch_corrected_df, "C:/Users/duho/Desktop/Umap_data.csv")