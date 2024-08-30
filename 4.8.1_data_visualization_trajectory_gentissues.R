library(dplyr)
library(ggplot2)
library(ggsignif)

# Define the path to the folder containing normalized count values
data_path <- "./adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Extract tissue names from the file names
tissue_names <- gsub(".rds$", "", basename(tissue_files))

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
  if (!tissue_name %in% names(group_labels_list)) {
    next  # Skip this tissue if group labels are not defined
  }
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Check the number of rows and columns in the dataframe
  num_rows <- nrow(normalized_counts)
  num_cols <- ncol(normalized_counts)
  print(paste("Number of rows:", num_rows))
  print(paste("Number of columns:", num_cols))
  
  # Ensure the number of gene names matches the number of rows
  if (num_rows != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  # Replace rownames with the manually defined gene names
  rownames(normalized_counts) <- gene_names
  
  # Define a small constant for log transformation
  small_constant <- 1e-9
  
  # Transform the data
  normalized_counts <- as.data.frame(normalized_counts)
  data_transformed_log10 <- as.data.frame(lapply(normalized_counts, function(x) log10(x + small_constant)))
  rownames(data_transformed_log10) <- gene_names
  
  # Define the group labels for the current tissue
  group_labels <- group_labels_list[[tissue_name]]
  
  # Ensure length of group_labels matches the number of columns in the transformed data
  if (length(group_labels) != num_cols) {
    stop(paste("Error: The length of group_labels does not match the number of columns in", tissue_name))
  }
  
  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_trajectory", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Combine the data for all genes into one data frame
  df_combined <- data.frame()
  for (gene in gene_names) {
    # Extract the values from the current gene
    gene_values <- data_transformed_log10[gene, ]
    gene_values <- t(gene_values)
    
    # Create a new data frame with these values
    df <- as.data.frame(gene_values)
    df$age_groups <- group_labels
    df$gene <- gene  # Add gene name as a column
    
    # Rename columns for consistency
    colnames(df) <- c("Expression", "age_groups", "gene")
    
    # Append to the combined data frame
    df_combined <- rbind(df_combined, df)
  }
  
  # Ensure correct column names in df_combined
  colnames(df_combined) <- c("Expression", "age_groups", "gene")
  df_combined$age_groups <- factor(df_combined$age_groups, levels = unique(group_labels))
  
  # Create the plot
  p1 <- ggplot(df_combined, aes(x = age_groups, y = Expression, color = gene, group = gene)) +
    geom_smooth(method = "loess", se = FALSE) +  # Smooth lines for each gene
    geom_point(size = 2) +  # Add points
    labs(title = paste("Expression Trajectories of Longevity Genes in", tissue_name),
         x = "Age Group",
         y = "Log10(Expression + Small Constant)",
         color = "Gene") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      plot.title = element_text(size = 16, face = "bold"),
      plot.background = element_rect(fill = "white"),  # White background for the legend
    )
  
  # Save the plot as PDF and PNG in the respective tissue directory
  ggsave(file.path(result_directory, paste0('trajectory_plot_all_genes.pdf')), p1, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("trajectory_plot_all_genes.png")), plot = p1, units = "cm", width = 18, height = 18, dpi = 300)
}
