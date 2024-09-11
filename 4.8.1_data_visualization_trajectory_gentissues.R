library(dplyr)
library(ggplot2)
library(ggsignif)

#path to the folder containing normalized count values
data_path <- "./adjusted_sva"

tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

tissue_names <- gsub(".rds$", "", basename(tissue_files))

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
