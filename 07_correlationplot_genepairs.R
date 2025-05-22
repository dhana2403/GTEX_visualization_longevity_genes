
########################################individual correlation plot################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(tools)
library(tidyr)

# Define the path to the folder containing normalized count values
data_path <- "./adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define the genes you want to correlate
gene_names <- c("MTOR", "FOXO3", "CYC1", "SIRT1", "PTEN", "SIRT3", "RRAS2", "CAT",
                "IGF1", "KLOTHO", "IGF1R", "TP53", "SIRT6", "APOE", "RRAS", "SOD1")

# Generate all gene pairs
gene_pairs <- combn(gene_names, 2, simplify = FALSE)

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

# Main results directory
main_results_directory <- "results_after_sva_correlation_age_genes_separatelines"

# Initialize a data frame to store results
results_df <- data.frame(
  Tissue = character(),
  Gene1 = character(),
  Gene2 = character(),
  Age_Group = character(),
  Spearman_Correlation = numeric(),
  P_Value = numeric(),
  FDR_P_Value = numeric(),  # Add column for FDR-adjusted p-values
  Significance = character(),
  Direction = character(),
  stringsAsFactors = FALSE
)

# Loop through each tissue file and perform correlation analysis
for (tissue_file in tissue_files) {
  
  # Extract tissue name from file path
  tissue_name <- file_path_sans_ext(basename(tissue_file))
  
  # Create the tissue-specific directory
  tissue_directory <- file.path(main_results_directory, tissue_name)
  dir.create(tissue_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Ensure rownames are correct
  rownames(normalized_counts) <- gene_names
  
  # Check if all genes are present in the data
  missing_genes <- setdiff(gene_names, rownames(normalized_counts))
  if (length(missing_genes) > 0) {
    warning(paste("The following genes are not found in", tissue_name, ":", paste(missing_genes, collapse = ", ")))
    next
  }
  
  # Get the age group labels for the current tissue
  age_groups <- group_labels_list[[tissue_name]]
  
  # Ensure that the length of age_groups matches the number of columns in normalized_counts
  if (length(age_groups) != ncol(normalized_counts)) {
    stop(paste("The number of age group labels does not match the number of samples in", tissue_name))
  }
  
  # Add age groups as a new row to the normalized_counts data
  normalized_counts <- rbind(normalized_counts, Age_Group = age_groups)
  
  # Initialize a list to store p-values for FDR adjustment
  p_values_list <- c()
  
  # Loop through each gene pair
  for (pair in gene_pairs) {
    gene1 <- pair[1]
    gene2 <- pair[2]
    
    # Extract the expression values for the two genes
    gene1_expression <- as.numeric(normalized_counts[gene1, ])
    gene2_expression <- as.numeric(normalized_counts[gene2, ])
    
    # Check for and handle NA values
    if (any(is.na(gene1_expression)) || any(is.na(gene2_expression))) {
      valid_indices <- complete.cases(gene1_expression, gene2_expression)
      gene1_expression <- gene1_expression[valid_indices]
      gene2_expression <- gene2_expression[valid_indices]
      age_groups <- age_groups[valid_indices]
    }
    
    # Prepare the data for plotting
    df_plot <- data.frame(
      Age_Group = factor(age_groups, levels = c("20-49", "50-59", "60-79")),
      Gene1_Expression = gene1_expression,
      Gene2_Expression = gene2_expression
    )
    
    # Ensure no negative values for log transformation
    df_plot$Gene1_Expression[df_plot$Gene1_Expression < 0] <- NA
    df_plot$Gene2_Expression[df_plot$Gene2_Expression < 0] <- NA
    
    # Apply log transformation to the expression data for better visualization
    df_plot$Gene1_Expression_log <- log10(df_plot$Gene1_Expression + 1)
    df_plot$Gene2_Expression_log <- log10(df_plot$Gene2_Expression + 1)
    
    # Handle cases where log transformation resulted in NA values
    df_plot <- df_plot %>% drop_na(Gene1_Expression_log, Gene2_Expression_log)
    
    # Calculate Spearman correlation and p-value for each age group
    correlation_by_age_group <- df_plot %>%
      group_by(Age_Group) %>%
      summarize(
        Spearman_Correlation = if (n() > 1) cor(Gene1_Expression_log, Gene2_Expression_log, method = "spearman") else NA,
        P_Value = if (n() > 1) cor.test(Gene1_Expression_log, Gene2_Expression_log, method = "spearman")$p.value else NA,
        .groups = 'drop'
      )
    
    # Collect p-values for FDR adjustment
    p_values_list <- c(p_values_list, correlation_by_age_group$P_Value)
    
    # Determine significance and direction
    correlation_by_age_group <- correlation_by_age_group %>%
      mutate(
        Significance = ifelse(P_Value < 0.05, "Significant", "Not Significant"),
        Direction = ifelse(Spearman_Correlation > 0, "Positive", "Negative")
      )
    
    # Create a plot with separate regression lines and annotations
    p <- ggplot(df_plot, aes(x = Gene1_Expression_log, y = Gene2_Expression_log, color = Age_Group)) +
      geom_point() +
      geom_smooth(method = "lm", aes(group = Age_Group), se = TRUE) +  # Separate regression lines
      labs(
        title = paste("Correlation between", gene1, "and", gene2, "Expression in", tissue_name),
        x = paste(gene1, "Expression (log10)"),
        y = paste(gene2, "Expression (log10)")
      ) +
      theme_minimal() +
      theme(
        plot.background = element_blank(),  # Remove the plot outline
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
      ) +
      annotate("text", x = Inf, y = Inf,
               label = paste("Spearman Correlation:\n",
                             paste(sapply(1:nrow(correlation_by_age_group), function(i) {
                               paste(
                                 correlation_by_age_group$Age_Group[i], ": ",
                                 round(correlation_by_age_group$Spearman_Correlation[i], 2),
                                 " (p=", sprintf("%.2e", correlation_by_age_group$P_Value[i]), ")",
                                 sep = ""
                               )
                             }), collapse = "\n")),
               hjust = 1.1, vjust = 1.1, size = 5, color = c("black"), family = "mono")
    
    # Save the plot as PNG
    ggsave(
      filename = file.path(tissue_directory, paste0("correlation_plot_", gene1, "_vs_", gene2, ".png")),
      plot = p, units = "cm", width = 18, height = 18, dpi = 300
    )
    
    # Append results to results_df
    results_df <- results_df %>%
      bind_rows(
        data.frame(
          Tissue = tissue_name,
          Gene1 = gene1,
          Gene2 = gene2,
          Age_Group = correlation_by_age_group$Age_Group,
          Spearman_Correlation = correlation_by_age_group$Spearman_Correlation,
          P_Value = correlation_by_age_group$P_Value,
          FDR_P_Value = NA,  # Initialize FDR column
          Significance = correlation_by_age_group$Significance,
          Direction = correlation_by_age_group$Direction,
          stringsAsFactors = FALSE
        )
      )
  }
  
  # Calculate FDR-adjusted p-values
  results_df$FDR_P_Value <- p.adjust(results_df$P_Value, method = "BH")
  
  # Save the results to a CSV file
  write.csv(results_df, file.path(main_results_directory, "correlation_results.csv"), row.names = FALSE)
  
  # Print a message indicating the process is complete
  cat("Correlation analysis complete. Results saved to", file.path(main_results_directory, "correlation_results.csv"), "\n")
}
