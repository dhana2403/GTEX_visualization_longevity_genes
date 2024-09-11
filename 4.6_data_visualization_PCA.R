library(dplyr)
library(ggplot2)
library(edgeR)
library(ggrepel)

# path to the folder containing normalized count values
data_path <- "./adjusted_sva"

tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

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

gene_names <- c("MTOR", "FOXO3", "CYC1", "SIRT1", "PTEN", "SIRT3", "RRAS2", "CAT", 
                "IGF1", "KLOTHO", "IGF1R", "TP53", "SIRT6", "APOE", "RRAS", "SOD1")

for (tissue_file in tissue_files) {

  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  if (!tissue_name %in% names(group_labels_list)) next
  
  normalized_counts <- readRDS(tissue_file)
  
  if (nrow(normalized_counts) != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  rownames(normalized_counts) <- gene_names
  
  normalized_counts[normalized_counts < 0] <- 1e-9
  
  data_for_pca <- as.data.frame(t(normalized_counts))
  
  data_for_pca <- scale(data_for_pca)
  
  pca_result <- prcomp(data_for_pca, center = TRUE, scale. = TRUE)
  
  pca_data <- as.data.frame(pca_result$x)
  pca_data$Group <- group_labels_list[[tissue_name]]

  
  pca_variance <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100
  pc1_var <- round(pca_variance[1], 2)
  pc2_var <- round(pca_variance[2], 2)
  
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(alpha = 0.7, size = 3) + 
    scale_color_manual(values = c("black", "red", "blue", "green", "purple", "orange")) + 
    labs(
      title = paste("PCA of", tissue_name, "\nPC1:", pc1_var, "%, PC2:", pc2_var, "%"),
      x = paste("Principal Component 1 (", pc1_var, "%)", sep = ""),
      y = paste("Principal Component 2 (", pc2_var, "%)", sep = ""),
      color = "Group"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"), 
      plot.background = element_rect(fill = "white"), 
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    )
  
  result_directory <- file.path("results_after_sva_PCA", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(result_directory, paste0('pca_plot_', tissue_name, '.pdf')), pca_plot, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("pca_plot_", tissue_name, ".png")), plot = pca_plot, units = "cm", width = 18, height = 18, dpi = 300)
}


rm(list = ls())
