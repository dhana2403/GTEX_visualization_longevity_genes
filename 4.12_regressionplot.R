
library(ggplot2)

data_path <- "./adjusted_sva"

tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define gene names manually
gene_names <- #INPUT gene names of interest

selected_genes <- # SELECT YOUR GENE OF INTEREST

tissues_of_interest <- # SELECT YOUR TISSUE OF INTEREST

# Initialize a data frame to store expression values
expression_data <- data.frame()

# Loop through each tissue file
for (tissue_file in tissue_files) {
  
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Check if the tissue is one of the tissues of interest
  if (!tissue_name %in% tissues_of_interest) {
    next  # Skip this tissue if it's not in the list
  }
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Ensure the number of gene names matches the number of rows
  if (nrow(normalized_counts) != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  # Replace rownames with the manually defined gene names
  rownames(normalized_counts) <- gene_names
  
  # Subset the data for the selected genes (GENE1 and GENE2)
  tissue_data <- normalized_counts[selected_genes, , drop = FALSE]
  
  # Convert to data frame and transpose
  tissue_data <- as.data.frame(t(tissue_data))
  
  # Add a tissue label
  tissue_data$Tissue <- tissue_name
  
  # Combine expression data
  expression_data <- rbind(expression_data, tissue_data)
}

# Log-transform the expression data
expression_data$GENE1_log <- log10(expression_data$GENE1 + 1e-9)  
expression_data$GENE2_log <- log10(expression_data$GENE2 + 1e-9) 

# Perform linear regression analysis
# Fit a linear regression model: GENE1_log ~ GENE2_log
model <- lm(GENE1_log ~ GENE2_log, data = expression_data)

model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]  # p-value for the slope
r_squared <- model_summary$r.squared          # R-squared value

# Plot the regression results
p <- ggplot(expression_data, aes(x = GENE1_log, y = GENE2_log)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Regression Analysis: Log-transformed gene1 vs gene2 Expression in tissue of interest",
       x = "Log10(GENE1 Expression + 1e-9)",
       y = "Log10(GENE2 Expression + 1e-9)") +
  theme_bw() +  
  
  # Annotate the plot with p-value and R-squared
  annotate("text", x = Inf, y = Inf,
           label = paste("R² =", round(r_squared, 3), "\np-value =", formatC(p_value, format = "e", digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 4, color = "black")

# Save the plot
result_directory <- "results_regression_analysis"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(result_directory, 'regression_plot_gene1_vs_gene2_beforesva.pdf'), p, width = 10, height = 8)
ggsave(file.path(result_directory, 'regression_plot_gene1_vs_gene2_beforesva.png'), p, width = 10, height = 8, dpi = 300)

cat("Regression plot saved to:", result_directory, "\n")
