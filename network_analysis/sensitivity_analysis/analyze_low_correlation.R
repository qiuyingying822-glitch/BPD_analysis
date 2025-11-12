#!/usr/bin/env Rscript
#
# Sensitivity Analysis: Distribution Analysis
# This script analyzes correlation distributions and method consistency

library(ggplot2)
library(reshape2)
library(dplyr)

cat("=== Microbial Network Method Comparison Analysis ===\n")

groups <- c("none", "mild", "moderate", "severe")

analyze_distributions_fixed <- function(group) {
  cat("\n=== Analyzing", group, "group ===\n")
  
  # Read results
  sparcc_file <- paste0(group, "/cor_sparcc.norm_rarefy_", group, ".txt")
  spearman_file <- paste0("sensitivity_analysis/spearman_cor_", group, ".txt")
  
  sparcc_cor <- as.matrix(read.table(sparcc_file, header = TRUE, 
                                    row.names = 1, sep = "\t"))
  spearman_cor <- as.matrix(read.table(spearman_file, header = TRUE, 
                                      row.names = 1, sep = "\t"))
  
  # Ensure common ASVs
  common_asvs <- intersect(rownames(sparcc_cor), rownames(spearman_cor))
  cat("  Common ASVs:", length(common_asvs), "\n")
  
  if (length(common_asvs) < 10) {
    cat("  Too few common ASVs, skipping\n")
    return(NULL)
  }
  
  sparcc_cor <- sparcc_cor[common_asvs, common_asvs]
  spearman_cor <- spearman_cor[common_asvs, common_asvs]
  
  # Convert to vectors (remove diagonal)
  sparcc_vec <- sparcc_cor[upper.tri(sparcc_cor)]
  spearman_vec <- spearman_cor[upper.tri(spearman_cor)]
  
  # Basic statistics
  cat("SparCC Statistics:\n")
  cat("  Range:", round(range(sparcc_vec), 3), "\n")
  cat("  Mean:", round(mean(sparcc_vec), 4), "\n")
  cat("  SD:", round(sd(sparcc_vec), 4), "\n")
  cat("  |r| > 0.3 proportion:", round(sum(abs(sparcc_vec) > 0.3) / length(sparcc_vec) * 100, 2), "%\n")
  
  cat("Spearman Statistics:\n")
  cat("  Range:", round(range(spearman_vec), 3), "\n")
  cat("  Mean:", round(mean(spearman_vec), 4), "\n")
  cat("  SD:", round(sd(spearman_vec), 4), "\n")
  cat("  |r| > 0.3 proportion:", round(sum(abs(spearman_vec) > 0.3) / length(spearman_vec) * 100, 2), "%\n")
  
  # Analyze strong correlation edge consistency
  n_edges <- length(sparcc_vec)
  strong_sparcc <- abs(sparcc_vec) > 0.3
  strong_spearman <- abs(spearman_vec) > 0.3
  
  cat("Strong Correlation Edge Analysis (|r| > 0.3):\n")
  cat("  SparCC strong edges:", sum(strong_sparcc), "\n")
  cat("  Spearman strong edges:", sum(strong_spearman), "\n")
  cat("  Common strong edges:", sum(strong_sparcc & strong_spearman), "\n")
  
  total_strong <- sum(strong_sparcc | strong_spearman)
  if (total_strong > 0) {
    jaccard <- sum(strong_sparcc & strong_spearman) / total_strong
    cat("  Strong edge Jaccard index:", round(jaccard, 4), "\n")
  } else {
    cat("  Strong edge Jaccard index: 0 (no strong edges)\n")
  }
  
  # Plot distributions (separate plots)
  df_sparcc <- data.frame(
    Method = "SparCC",
    Correlation = sparcc_vec
  )
  
  df_spearman <- data.frame(
    Method = "Spearman",
    Correlation = spearman_vec
  )
  
  df <- rbind(df_sparcc, df_spearman)
  
  p <- ggplot(df, aes(x = Correlation, fill = Method)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ Method, ncol = 2, scales = "free_y") +
    labs(title = paste0(group, " Group: Correlation Coefficient Distribution"),
         x = "Correlation Coefficient",
         y = "Density") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  ggsave(paste0("sensitivity_analysis/distribution_", group, "_english.png"), p,
         width = 10, height = 5, dpi = 300)
  
  # Return key findings
  return(list(
    sparcc_strong = sum(strong_sparcc),
    spearman_strong = sum(strong_spearman),
    common_strong = sum(strong_sparcc & strong_spearman),
    jaccard = ifelse(total_strong > 0, jaccard, 0)
  ))
}

# Execute analysis and summarize
results <- list()
for(group in groups) {
  results[[group]] <- analyze_distributions_fixed(group)
}

# Generate summary table
cat("\n=== Key Findings Summary ===\n")
summary_table <- data.frame(
  Group = groups,
  SparCC_Strong_Edges = sapply(results, function(x) if(!is.null(x)) x$sparcc_strong else NA),
  Spearman_Strong_Edges = sapply(results, function(x) if(!is.null(x)) x$spearman_strong else NA),
  Common_Strong_Edges = sapply(results, function(x) if(!is.null(x)) x$common_strong else NA),
  Jaccard_Index = sapply(results, function(x) if(!is.null(x)) round(x$jaccard, 4) else NA)
)

print(summary_table)

# Save results
write.table(summary_table, "sensitivity_analysis/strong_edges_comparison.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\nResults saved to: sensitivity_analysis/strong_edges_comparison.txt\n")
