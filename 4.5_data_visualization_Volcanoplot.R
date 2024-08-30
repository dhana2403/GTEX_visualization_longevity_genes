library(dplyr)
library(ggplot2)
library(ggrepel)
library(edgeR)  # Ensure this is loaded for DGEList, calcNormFactors, etc.

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
  
  # Create DGEList object for edgeR analysis
  dge <- DGEList(counts = normalized_counts)
  group_labels <- group_labels_list[[tissue_name]]
  dge$samples$group <- factor(group_labels)
  
  # Perform normalization and differential expression analysis
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  fit <- glmQLFit(dge)
  result <- glmQLFTest(fit)
  
  # Extract results
  top_tags <- topTags(result, n = Inf)
  df <- as.data.frame(top_tags)
  
  # Handle potential edge cases for FDR or PValue
  df$FDR[df$FDR > 1] <- 1
  df$neg_log10_pvalue <- -log10(df$PValue)
  
  # Create a data frame for the volcano plot
  volcano_data <- data.frame(
    Gene = rownames(df),
    logFC = df$logFC,
    neg_log10_pvalue = df$neg_log10_pvalue,
    FDR = df$FDR
  )
  
  # Ensure there are no negative or infinite values in logFC or neg_log10_pvalue
  volcano_data <- na.omit(volcano_data)  # Remove NAs if any
  volcano_data <- volcano_data[!is.infinite(volcano_data$neg_log10_pvalue),]
  
  # Identify significant genes
  volcano_data$Significant <- volcano_data$FDR < 0.05
  significant_genes <- subset(volcano_data, Significant)
  
  # Plot
  p <- ggplot(volcano_data, aes(x = logFC, y = neg_log10_pvalue)) +
    geom_point(aes(color = Significant), alpha = 0.7) +
    scale_color_manual(values = c("black", "red"), labels = c("Not Significant", "Significant")) +
    labs(
      title = paste("Volcano Plot -", tissue_name),
      x = "Log Fold Change (logFC)",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),  # Set panel background to white
      plot.background = element_rect(fill = "white"),  # Set plot background to white
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 10)
    ) +
    geom_text_repel(
      data = significant_genes,
      aes(label = Gene),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = 'grey50'
    )
  
  # Print the plot
  print(p)
  
  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_volcano_plot", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Save the plot as PDF and PNG in the respective tissue directory
  ggsave(file.path(result_directory, paste0('volcano_plot_', tissue_name, '.pdf')), p, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("volcano_plot_", tissue_name, ".png")), plot = p, units = "cm", width = 18, height = 18, dpi = 300)
}

rm(list = ls())

