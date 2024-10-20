# Load necessary libraries
library(ggplot2)

# Define the path to the folder containing normalized count values
data_path <- "./data/processed/expression/readcounts_tmm_specific"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define gene names manually
gene_names <- c("MTOR", "FOXO3", "CYC1", "SIRT1", "PTEN", "SIRT3", "RRAS2", "CAT",
                "IGF1", "KLOTHO", "IGF1R", "TP53", "SIRT6", "APOE", "RRAS", "SOD1")

# Define the selected genes (CYC1 and FOXO3)
selected_genes <- c("CYC1", "FOXO3")

# Define the tissue of interest (Muscle-Skeletal)
tissues_of_interest <- c("Muscle-Skeletal")

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
  
  # Subset the data for the selected genes (CYC1 and FOXO3)
  tissue_data <- normalized_counts[selected_genes, , drop = FALSE]
  
  # Convert to data frame and transpose
  tissue_data <- as.data.frame(t(tissue_data))
  
  # Add a tissue label
  tissue_data$Tissue <- tissue_name
  
  # Combine expression data
  expression_data <- rbind(expression_data, tissue_data)
}

# Check if expression_data is populated
if (nrow(expression_data) == 0) {
  stop("No data found for the selected genes in the specified tissues.")
}

# Log-transform the expression data
expression_data$CYC1_log <- log10(expression_data$CYC1 + 1e-9)  # Log-transform CYC1 expression
expression_data$FOXO3_log <- log10(expression_data$FOXO3 + 1e-9)  # Log-transform FOXO3 expression

# Perform linear regression analysis
# Fit a linear regression model: CYC1_log ~ FOXO3_log
model <- lm(CYC1_log ~ FOXO3_log, data = expression_data)

# Summary of the model to get p-value and R-squared
model_summary <- summary(model)

# Extract p-value and R-squared
p_value <- model_summary$coefficients[2, 4]  # p-value for the slope
r_squared <- model_summary$r.squared          # R-squared value

# Plot the regression results
p <- ggplot(expression_data, aes(x = FOXO3_log, y = CYC1_log)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "blue") +
  labs(title = "Regression Analysis: Log-transformed CYC1 vs FOXO3 Expression in Muscle-Skeletal",
       x = "Log10(FOXO3 Expression + 1e-9)",
       y = "Log10(CYC1 Expression + 1e-9)") +
  theme_bw() +  # Change to a white background theme
  
  # Annotate the plot with p-value and R-squared
  annotate("text", x = Inf, y = Inf,
           label = paste("R² =", round(r_squared, 3), "\np-value =", formatC(p_value, format = "e", digits = 2)),
           hjust = 1.1, vjust = 1.1, size = 4, color = "black")

# Save the plot
result_directory <- "results_regression_analysis"
dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)

# Save the plot as PDF
ggsave(file.path(result_directory, 'regression_plot_cyc1_vs_foxo3_beforesva.pdf'), p, width = 10, height = 8)

# Save the plot as PNG
ggsave(file.path(result_directory, 'regression_plot_cyc1_vs_foxo3_beforesva.png'), p, width = 10, height = 8, dpi = 300)

# Print message after saving
cat("Regression plot saved to:", result_directory, "\n")
