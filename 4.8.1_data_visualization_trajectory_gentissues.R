library(dplyr)
library(ggplot2)
library(ggsignif)

#path to the folder containing normalized count values
data_path <- "./adjusted_sva"

tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

tissue_names <- gsub(".rds$", "", basename(tissue_files))


# --- Hidden Code (Available after publication) ---
# The group labels and gene names are hidden and will be available after the results are published.


for (tissue_file in tissue_files) {

  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  if (!tissue_name %in% names(group_labels_list)) {
    next  
  }
  
  normalized_counts <- readRDS(tissue_file)
  
  num_rows <- nrow(normalized_counts)
  num_cols <- ncol(normalized_counts)
  print(paste("Number of rows:", num_rows))
  print(paste("Number of columns:", num_cols))
  
  if (num_rows != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  rownames(normalized_counts) <- gene_names
  
  small_constant <- 1e-9
  
  normalized_counts <- as.data.frame(normalized_counts)
  data_transformed_log10 <- as.data.frame(lapply(normalized_counts, function(x) log10(x + small_constant)))
  rownames(data_transformed_log10) <- gene_names
  
  group_labels <- group_labels_list[[tissue_name]]
  
  if (length(group_labels) != num_cols) {
    stop(paste("Error: The length of group_labels does not match the number of columns in", tissue_name))
  }
  
  result_directory <- file.path("results_after_sva_trajectory", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  df_combined <- data.frame()
  for (gene in gene_names) {

    gene_values <- data_transformed_log10[gene, ]
    gene_values <- t(gene_values)
    
    df <- as.data.frame(gene_values)
    df$age_groups <- group_labels
    df$gene <- gene  
    
    colnames(df) <- c("Expression", "age_groups", "gene")
    
    df_combined <- rbind(df_combined, df)
  }
  
  colnames(df_combined) <- c("Expression", "age_groups", "gene")
  df_combined$age_groups <- factor(df_combined$age_groups, levels = unique(group_labels))
  
  p1 <- ggplot(df_combined, aes(x = age_groups, y = Expression, color = gene, group = gene)) +
    geom_smooth(method = "loess", se = FALSE) +  
    geom_point(size = 2) + 
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
      plot.background = element_rect(fill = "white"),  
    )
  
  ggsave(file.path(result_directory, paste0('trajectory_plot_all_genes.pdf')), p1, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("trajectory_plot_all_genes.png")), plot = p1, units = "cm", width = 18, height = 18, dpi = 300)
}
