#!/usr/bin/env Rscript
#
# SparCC Results Processing
# This script processes SparCC correlation results and generates network files

library(dplyr)
library(reshape2)
library(igraph)

cat("=== Processing SparCC Results ===\n")

# Function to process a single group
process_sparcc_group <- function(group) {
  cat("Processing group:", group, "\n")
  
  # File paths
  cor_file <- paste0(group, "/cor_sparcc.norm_rarefy_", group, ".txt")
  pval_file <- paste0("pvals/", group, "/pvals.", group, ".txt")
  
  # Check if files exist
  if (!file.exists(cor_file)) {
    cat("  Correlation file not found:", cor_file, "\n")
    return(NULL)
  }
  
  if (!file.exists(pval_file)) {
    cat("  P-value file not found:", pval_file, "\n")
    return(NULL)
  }
  
  # Read correlation matrix
  cor_matrix <- as.matrix(read.table(cor_file, header = TRUE, 
                                    row.names = 1, sep = "\t"))
  
  # Read p-value matrix
  pval_matrix <- as.matrix(read.table(pval_file, header = TRUE, 
                                     row.names = 1, sep = "\t"))
  
  cat("  Correlation matrix dimensions:", dim(cor_matrix), "\n")
  cat("  P-value matrix dimensions:", dim(pval_matrix), "\n")
  
  # Ensure matrices have the same dimensions
  common_asvs <- intersect(rownames(cor_matrix), rownames(pval_matrix))
  
  if (length(common_asvs) < 10) {
    cat("  Too few common ASVs:", length(common_asvs), "\n")
    return(NULL)
  }
  
  cor_matrix <- cor_matrix[common_asvs, common_asvs]
  pval_matrix <- pval_matrix[common_asvs, common_asvs]
  
  # Create combined results data frame
  results_df <- data.frame(
    ASV1 = character(),
    ASV2 = character(),
    Correlation = numeric(),
    P_value = numeric(),
    Significant = logical(),
    stringsAsFactors = FALSE
  )
  
  # Fill the data frame (upper triangle only)
  for (i in 1:(nrow(cor_matrix) - 1)) {
    for (j in (i + 1):ncol(cor_matrix)) {
      asv1 <- rownames(cor_matrix)[i]
      asv2 <- colnames(cor_matrix)[j]
      
      correlation <- cor_matrix[i, j]
      p_value <- pval_matrix[i, j]
      significant <- p_value < 0.05 & abs(correlation) > 0.3
      
      results_df <- rbind(results_df, data.frame(
        ASV1 = asv1,
        ASV2 = asv2,
        Correlation = correlation,
        P_value = p_value,
        Significant = significant,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Calculate network statistics
  n_edges <- nrow(results_df)
  n_significant <- sum(results_df$Significant)
  positive_cor <- sum(results_df$Significant & results_df$Correlation > 0)
  negative_cor <- sum(results_df$Significant & results_df$Correlation < 0)
  
  cat("  Total edges:", n_edges, "\n")
  cat("  Significant edges:", n_significant, "\n")
  cat("  Positive correlations:", positive_cor, "\n")
  cat("  Negative correlations:", negative_cor, "\n")
  
  # Save detailed results
  output_file <- paste0("results/sparcc_network_", group, ".txt")
  write.table(results_df, output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Results saved to:", output_file, "\n")
  
  # Create and save significant network only
  significant_edges <- results_df[results_df$Significant, ]
  sig_output_file <- paste0("results/sparcc_significant_network_", group, ".txt")
  write.table(significant_edges, sig_output_file, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("  Significant network saved to:", sig_output_file, "\n")
  
  # Generate network summary
  summary_stats <- data.frame(
    Group = group,
    Total_ASVs = length(common_asvs),
    Total_Edges = n_edges,
    Significant_Edges = n_significant,
    Positive_Correlations = positive_cor,
    Negative_Correlations = negative_cor,
    Significance_Ratio = round(n_significant / n_edges * 100, 2),
    Mean_Correlation = round(mean(results_df$Correlation), 4),
    Mean_Abs_Correlation = round(mean(abs(results_df$Correlation)), 4),
    stringsAsFactors = FALSE
  )
  
  return(summary_stats)
}

# Process all groups
groups <- c("none", "mild", "moderate", "severe")
all_summaries <- list()

for (group in groups) {
  summary_stats <- process_sparcc_group(group)
  if (!is.null(summary_stats)) {
    all_summaries[[group]] <- summary_stats
  }
}

# Combine all summaries
if (length(all_summaries) > 0) {
  combined_summary <- do.call(rbind, all_summaries)
  
  # Save combined summary
  write.table(combined_summary, "results/sparcc_network_summary.txt",
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat("\n=== SparCC Network Summary ===\n")
  print(combined_summary)
  cat("\nSummary saved to: results/sparcc_network_summary.txt\n")
  
  # Create visualization of network statistics
  if (require(ggplot2)) {
    p <- ggplot(combined_summary, aes(x = Group, y = Significant_Edges, fill = Group)) +
      geom_col() +
      labs(title = "Significant Network Edges by Group",
           x = "Group", y = "Number of Significant Edges") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave("results/sparcc_network_stats.png", p, width = 8, height = 6, dpi = 300)
    cat("Network statistics plot saved to: results/sparcc_network_stats.png\n")
  }
} else {
  cat("No valid SparCC results found for any group\n")
}

cat("=== SparCC Results Processing Completed ===\n")
