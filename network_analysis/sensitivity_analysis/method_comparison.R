#!/usr/bin/env Rscript
#
# Method Comparison Analysis
# This script compares SparCC and Spearman correlation methods in detail

library(ggplot2)
library(reshape2)
library(dplyr)

cat("=== Method Comparison: SparCC vs Spearman ===\n")

groups <- c("none", "mild", "moderate", "severe")

compare_methods_detailed <- function(group) {
  cat("Comparing", group, "group...\n")
  
  # Read SparCC results
  sparcc_file <- paste0(group, "/cor_sparcc.norm_rarefy_", group, ".txt")
  if (!file.exists(sparcc_file)) {
    cat(" SparCC file does not exist\n")
    return(NULL)
  }
  
  # Read Spearman results
  spearman_file <- paste0("sensitivity_analysis/spearman_cor_", group, ".txt")
  if (!file.exists(spearman_file)) {
    cat(" Spearman file does not exist\n")
    return(NULL)
  }
  
  sparcc_cor <- as.matrix(read.table(sparcc_file, header = TRUE, 
                                    row.names = 1, sep = "\t"))
  spearman_cor <- as.matrix(read.table(spearman_file, header = TRUE, 
                                      row.names = 1, sep = "\t"))
  
  cat(" SparCC dimensions:", dim(sparcc_cor), "\n")
  cat(" Spearman dimensions:", dim(spearman_cor), "\n")
  
  # Ensure matrix dimensions are consistent
  common_asvs <- intersect(rownames(sparcc_cor), rownames(spearman_cor))
  cat(" Common ASVs:", length(common_asvs), "\n")
  
  if (length(common_asvs) < 10) {
    cat(" Too few common ASVs, skipping\n")
    return(NULL)
  }
  
  sparcc_cor <- sparcc_cor[common_asvs, common_asvs]
  spearman_cor <- spearman_cor[common_asvs, common_asvs]
  
  # Convert to long format
  sparcc_melt <- melt(sparcc_cor)
  spearman_melt <- melt(spearman_cor)
  
  combined <- data.frame(
    ASV1 = sparcc_melt$Var1,
    ASV2 = sparcc_melt$Var2,
    SparCC = sparcc_melt$value,
    Spearman = spearman_melt$value
  )
  
  # Remove diagonal
  combined <- combined[combined$ASV1 != combined$ASV2, ]
  cat(" Valid edges:", nrow(combined), "\n")
  
  # Calculate correlation
  cor_test <- cor.test(combined$SparCC, combined$Spearman, method = "spearman")
  
  # Create scatter plot
  p <- ggplot(combined, aes(x = SparCC, y = Spearman)) +
    geom_hex(bins = 50) +
    scale_fill_viridis_c(trans = "log10") +
    geom_smooth(method = "lm", color = "red") +
    labs(title = paste0(group, " Group: SparCC vs Spearman"),
         subtitle = paste0("Spearman ρ = ", round(cor_test$estimate, 4),
                          ", p = ", format.pval(cor_test$p.value, digits = 4)),
         x = "SparCC Correlation Coefficient",
         y = "Spearman Correlation Coefficient") +
    theme_bw()
  
  ggsave(paste0("sensitivity_analysis/comparison_", group, "_hex.png"), p,
         width = 8, height = 6, dpi = 300)
  
  # Save detailed results
  result <- list(
    cor_coef = cor_test$estimate,
    p_value = cor_test$p.value,
    n_edges = nrow(combined),
    data = combined
  )
  
  return(result)
}

# Execute comparison
results <- list()
for(group in groups) {
  results[[group]] <- compare_methods_detailed(group)
}

# Generate summary
summary_df <- data.frame(
  Group = groups,
  Correlation = sapply(results, function(x) if(!is.null(x)) round(x$cor_coef, 4) else NA),
  P_value = sapply(results, function(x) if(!is.null(x)) format.pval(x$p_value, digits = 4) else NA),
  Number_of_edges = sapply(results, function(x) if(!is.null(x)) x$n_edges else NA)
)

write.table(summary_df, "sensitivity_analysis/method_comparison_summary_new.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("\n=== Method Comparison Summary ===\n")
print(summary_df)
cat("\nComparison plots saved to sensitivity_analysis/\n")
