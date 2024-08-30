library(dplyr)
library(ggplot2)
library(edgeR)
library(ggrepel)

# Define the path to the folder containing normalized count values
data_path <- "./adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define group labels for each tissue
group_labels_list <- list(
  "Brain-Cortex" = c(rep("20-49", 9), rep("50-59", 18), rep("60-79", 55)),
  "Heart-AtrialAppendage" = c(rep("20-49", 9), rep("50-59", 20), rep("60-79", 49)), 
  "Heart-LeftVentricle" = c(rep("20-49", 10), rep("50-59", 14), rep("60-79", 38)),
  "Kidney-Cortex" = c(rep("20-49", 4), rep("50-59", 3), rep("60-79", 15)),
  "Liver" = c(rep("20-49", 4), rep("50-59", 8), rep("60-79", 21)),
  "Lung" = c(rep("20-49", 12), rep("50-59", 24), rep("60-79", 59)),
  "Muscle-Skeletal" = c(rep("20-49", 16), rep("50-59", 28), rep("60-79", 89)),
  "Skin-NotSunExposed" = c(rep("20-49", 11), rep("50-59", 20), rep("60-79", 76)),
  "Skin-SunExposed" = c(rep("20-49", 15), rep("50-59", 26), rep("60-79", 79))
)

# Define gene names manually
gene_names <- c("MTOR", "FOXO3", "CYC1", "SIRT1", "PTEN", "SIRT3", "RRAS2", "CAT", 
                "IGF1", "KLOTHO", "IGF1R", "TP53", "SIRT6", "APOE", "RRAS", "SOD1")

# Loop through each tissue file
for (tissue_file in tissue_files) {
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Check if the tissue is in the group labels list
  if (!tissue_name %in% names(group_labels_list)) next
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Ensure the number of gene names matches the number of rows
  if (nrow(normalized_counts) != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  # Replace rownames with the manually defined gene names
  rownames(normalized_counts) <- gene_names
  
  # Handle negative counts by replacing them with a small positive value
  normalized_counts[normalized_counts < 0] <- 1e-9
  
  # Convert to data frame and transpose for PCA
  data_for_pca <- as.data.frame(t(normalized_counts))
  
  # Z-score standardize the data
  data_for_pca <- scale(data_for_pca)
  
  # Perform PCA
  pca_result <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)
  
  # Extract PCA results
  pca_data <- as.data.frame(pca_result$x)
  pca_data$Group <- group_labels_list[[tissue_name]]

  
  # Calculate percentage of variance explained
  pca_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  pc1_var <- round(pca_variance[1], 2)
  pc2_var <- round(pca_variance[2], 2)
  
  # Plot PCA
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 0.7, size = 3) + # Adjust size for better visibility
    scale_color_manual(values = c("black", "red", "blue", "green", "purple", "orange")) + # Adjust colors as needed
    labs(
      title = paste("PCA of", tissue_name, "\nPC1:", pc1_var, "%, PC2:", pc2_var, "%"),
      x = paste("Principal Component 1 (", pc1_var, "%)", sep = ""),
      y = paste("Principal Component 2 (", pc2_var, "%)", sep = ""),
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"), # Set background to white
      plot.background = element_rect(fill = "white"), # Set plot background to white
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_PCA", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Save the PCA plot as PDF and PNG in the respective tissue directory
  ggsave(file.path(result_directory, paste0('pca_plot_', tissue_name, '.pdf')), pca_plot, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("pca_plot_", tissue_name, ".png")), plot = pca_plot, units = "cm", width = 18, height = 18, dpi = 300)
}


rm(list = ls())
