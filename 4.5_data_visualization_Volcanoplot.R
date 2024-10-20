library(dplyr)
library(ggplot2)
library(ggrepel)
library(edgeR)  

# path to the folder containing normalized count values
data_path <- "./adjusted_sva"

tissue_files <- list.files(path = data_path, pattern = "*.rds", full.names = TRUE)

# --- Hidden Code (Available after publication) ---
# The group labels and gene names will be available after the results are published.

for (tissue_file in tissue_files) {
  
  tissue_name <- gsub(".rds$", "", basename(tissue_file))
  
  if (!tissue_name %in% names(group_labels_list)) next
  
  normalized_counts <- readRDS(tissue_file)
  
  if (nrow(normalized_counts) != length(gene_names)) {
    stop(paste("Error: The number of gene names does not match the number of rows in", tissue_name))
  }
  
  rownames(normalized_counts) <- gene_names
  
  normalized_counts[normalized_counts < 0] <- 1e-9
  
  dge <- DGEList(counts = normalized_counts)
  group_labels <- group_labels_list[[tissue_name]]
  dge$samples$group <- factor(group_labels)
  
  dge <- calcNormFactors(dge)
  dge <- estimateDisp(dge)
  fit <- glmQLFit(dge)
  result <- glmQLFTest(fit)
  
  top_tags <- topTags(result, n = Inf)
  df <- as.data.frame(top_tags)
  
  df$FDR[df$FDR > 1] <- 1
  df$neg_log10_pvalue <- -log10(df$PValue)
  
  volcano_data <- data.frame(
    Gene = rownames(df),
    logFC = df$logFC,
    neg_log10_pvalue = df$neg_log10_pvalue,
    FDR = df$FDR
  )
  
  volcano_data <- na.omit(volcano_data)  
  volcano_data <- volcano_data[!is.infinite(volcano_data$neg_log10_pvalue),]
  
  volcano_data$Significant <- volcano_data$FDR < 0.05
  significant_genes <- subset(volcano_data, Significant)
  
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
      panel.background = element_rect(fill = "white"),  
      plot.background = element_rect(fill = "white"),  
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
  
  print(p)
  
  result_directory <- file.path("results_after_sva_volcano_plot", tissue_name)
  dir.create(result_directory, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(file.path(result_directory, paste0('volcano_plot_', tissue_name, '.pdf')), p, units = 'cm', width = 18, height = 18, useDingbats = FALSE)
  ggsave(filename = file.path(result_directory, paste0("volcano_plot_", tissue_name, ".png")), plot = p, units = "cm", width = 18, height = 18, dpi = 300)
}

rm(list = ls())

