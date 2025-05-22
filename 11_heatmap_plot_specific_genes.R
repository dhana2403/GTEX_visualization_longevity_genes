
# Load necessary libraries
library(ComplexHeatmap)
library(circlize)  # For colorRamp2 function

# Define the path to the folder containing normalized count values
data_path <- "./adjusted_sva"

# List all .rds files under the data path
tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# Define group labels for the Muscle-Skeletal tissue
group_labels_brain_cortex <- c(rep("20-49", 16), rep("50-59", 28), rep("60-79", 89))

# Define gene names of interest
genes_of_interest <- c("CYC1", "FOXO3")

# Loop through each tissue file
for (tissue_file in tissue_files) {
  # Extract the tissue name
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  # Process only the Muscle-Skeletal tissue
  if (tissue_name != "Muscle-Skeletal") {
    next  # Skip other tissues
  }
  
  # Load the normalized read counts data from the .rds file
  normalized_counts <- readRDS(tissue_file)
  
  # Replace rownames with the manually defined gene names
  rownames(normalized_counts) <- gene_names
  
  # Convert to data frame and log-transform
  normalized_counts <- as.data.frame(normalized_counts)
  data_log10 <- log10(normalized_counts + 1e-9)  # Add a small constant to avoid log(0)
  
  # Z-score normalization
  data_zscore <- as.data.frame(t(scale(t(data_log10))))
  
  # Filter to include only APOE and SOD1
  data_zscore_filtered <- data_zscore[rownames(data_zscore) %in% genes_of_interest, ]
  
  # Check for any NA values and remove rows/columns that are fully NA
  data_zscore_filtered <- na.omit(data_zscore_filtered)  # Removes any rows with NA values
  
  # Ensure that the length of group_labels matches the number of columns in the filtered data
  if (length(group_labels_brain_cortex) != ncol(data_zscore_filtered)) {
    stop(paste("Error: The length of group_labels does not match the number of columns in Muscle-Skeletal"))
  }
  
  # Convert filtered data_zscore to a data frame and transpose it for horizontal heatmap
  data_zscore_transposed <- as.data.frame(t(data_zscore_filtered))
  
  # Convert 'group_labels' to a factor with explicit levels
  group_labels_factor <- factor(group_labels_brain_cortex, levels = c("20-49", "50-59", "60-79"))
  
  # Create the annotation data frame for rows (genes)
  annotation_row <- data.frame(Group = group_labels_factor)
  
  # Define colors for row annotations
  annotation_colors <- list(
    Group = c(
      "20-49" = "lightblue",
      "50-59" = "mediumblue",
      "60-79" = "violet"
    )
  )
  
  # Create the heatmap annotation for rows
  ha_row <- rowAnnotation(
    Group = annotation_row$Group,
    col = annotation_colors
  )
  
  data_zscore_reordered <- as.matrix(data_zscore_transposed)
  
  # Create a directory to save the plots for the current tissue
  result_directory <- file.path("results_after_sva_cluster", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  # Define a custom color gradient using circlize::colorRamp2
  custom_colors <- colorRamp2(c(-2, 0, 2), c("pink", "white", "skyblue"))  # Corrected this line
  
  
  # Generate the horizontal heatmap with adjusted parameters
  p <- Heatmap(
    data_zscore_reordered,
    name = "Expression",
    left_annotation = ha_row,  # Apply row annotation on the left
    cluster_rows = FALSE,       # Do not cluster rows
    cluster_columns = FALSE,     # Do not cluster columns to preserve order
    show_column_names = TRUE,    # Show column names
    show_row_names = FALSE,       # Show row names
    col = custom_colors,      # Color gradient
    width = unit(15, "cm"),       # Adjust the overall heatmap width to reduce cell size
    height = unit(15, "cm"),      # Adjust the overall heatmap height to reduce cell size
    column_names_gp = gpar(fontsize = 8),  # Reduce column name text size
    row_names_gp = gpar(fontsize = 8)      # Reduce row name text size if row names are shown
  )
  
  # Save the heatmap as PDF and PNG with larger dimensions
  pdf(file.path(result_directory, paste0('heatmap_CYC1_FOXO3_Muscle_Skeletal.pdf')), width = 20, height = 16)
  draw(p)  # Use draw() to render the heatmap
  dev.off()
  
  png(file.path(result_directory, paste0('heatmap_CYC1_FOXO3_Muscle_Skeletal.png')), width = 6000, height = 8000, res = 300)
  draw(p)  # Use draw() to render the heatmap
  dev.off()
}
