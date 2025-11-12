
#!/usr/bin/env Rscript
#
# Keystone Taxa Sensitivity Analysis
# This script analyzes the sensitivity of keystone taxa identification across different correlation thresholds

library(dplyr)
library(ggplot2)
library(tidyr)
library(igraph)

cat("=== Keystone Taxa Sensitivity Analysis ===\n")

# Read node centrality data
if (file.exists("network/node_centrality_properties.txt")) {
  node_data <- read.table("network/node_centrality_properties.txt", 
                         header = TRUE, sep = "\t")
  
  # Use multiple thresholds for sensitivity analysis
  thresholds <- c(0.2, 0.3, 0.4)
  
  # Enhanced keystone taxa identification function
  identify_keystone_taxa_enhanced <- function(group_data, network_cor, threshold, top_percent = 0.1) {
    
    # Build network based on current threshold
    adj_matrix <- ifelse(abs(network_cor) > threshold, 1, 0)
    diag(adj_matrix) <- 0
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    
    # Remove isolated nodes
    isolated_nodes <- which(degree(g) == 0)
    if (length(isolated_nodes) > 0) {
      g <- delete.vertices(g, isolated_nodes)
    }
    
    # Calculate transitivity (clustering coefficient) - node level
    transitivity_local <- transitivity(g, type = "local")
    names(transitivity_local) <- V(g)$name
    
    # Calculate eigenvector centrality
    eigen_centrality <- eigen_centrality(g)$vector
    names(eigen_centrality) <- V(g)$name
    
    # Merge all topological features
    enhanced_data <- group_data %>%
      left_join(
        data.frame(
          node = names(transitivity_local),
          transitivity = transitivity_local,
          eigen_centrality = eigen_centrality[names(transitivity_local)]
        ),
        by = "node"
      )
    
    # Handle NA values
    enhanced_data$transitivity[is.na(enhanced_data$transitivity)] <- 0
    
    # Calculate standardized scores
    enhanced_data <- enhanced_data %>%
      mutate(
        degree_z = as.numeric(scale(degree)),
        closeness_z = as.numeric(scale(closeness)),
        transitivity_z = as.numeric(scale(transitivity)),
        eigen_z = as.numeric(scale(eigen_centrality)),
        betweenness_z = -as.numeric(scale(betweenness))
      )
    
    # Handle NA values
    enhanced_data[is.na(enhanced_data)] <- 0
    
    # Calculate comprehensive keystone score
    enhanced_data <- enhanced_data %>%
      mutate(
        keystone_score = (2 * degree_z) +
          (1.5 * closeness_z) +
          (1.5 * transitivity_z) +
          (1 * eigen_z) +
          (-1 * betweenness_z)
      )
    
    # Identify keystone taxa (top percent)
    n_keystone <- max(1, round(nrow(enhanced_data) * top_percent))
    keystone_taxa <- enhanced_data %>%
      arrange(desc(keystone_score)) %>%
      slice_head(n = n_keystone) %>%
      mutate(keystone_rank = 1:n(),
             threshold = threshold)
    
    return(keystone_taxa)
  }
  
  # Identify keystone taxa for each group and each threshold
  all_keystone_results <- list()
  
  for (threshold in thresholds) {
    cat("=== Analyzing threshold:", threshold, "===\n")
    node_data_threshold <- node_data %>% filter(threshold == !!threshold)
    
    if (nrow(node_data_threshold) > 0) {
      for (group in unique(node_data_threshold$group)) {
        group_data <- node_data_threshold %>% filter(group == !!group)
        
        # Read corresponding SparCC correlation matrix
        sparcc_file <- paste0(group, "/cor_sparcc.norm_rarefy_", group, ".txt")
        
        if (file.exists(sparcc_file) && nrow(group_data) >= 10) {
          sparcc_cor <- as.matrix(read.table(sparcc_file, header = TRUE, 
                                            row.names = 1, sep = "\t"))
          
          # Ensure common ASVs
          common_asvs <- intersect(rownames(sparcc_cor), group_data$node)
          
          if (length(common_asvs) >= 10) {
            sparcc_cor <- sparcc_cor[common_asvs, common_asvs]
            group_data <- group_data %>% filter(node %in% common_asvs)
            
            keystone_taxa <- identify_keystone_taxa_enhanced(group_data, sparcc_cor, threshold, top_percent = 0.1)
            keystone_taxa$group <- group
            all_keystone_results[[paste(group, threshold, sep = "_")]] <- keystone_taxa
            
            cat("  Group", group, ": identified", nrow(keystone_taxa), "keystone taxa\n")
          }
        }
      }
    }
  }
  
  # Combine all results
  if (length(all_keystone_results) > 0) {
    all_keystone <- do.call(rbind, all_keystone_results)
    
    # Save keystone taxa results
    write.table(all_keystone, "network/keystone_taxa_sensitivity_analysis.txt",
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Generate keystone taxa summary table
    keystone_summary <- all_keystone %>%
      group_by(group, threshold) %>%
      summarise(
        n_keystone = n(),
        avg_degree = mean(degree),
        avg_betweenness = mean(betweenness),
        avg_closeness = mean(closeness),
        avg_transitivity = mean(transitivity),
        avg_keystone_score = mean(keystone_score),
        .groups = 'drop'
      )
    
    write.table(keystone_summary, "network/keystone_taxa_sensitivity_summary.txt",
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Visualize sensitivity analysis results
    if (nrow(all_keystone) > 0) {
      # 1. Number of keystone taxa across different thresholds
      p1 <- ggplot(keystone_summary, aes(x = factor(threshold), y = n_keystone, fill = group)) +
        geom_col(position = "dodge") +
        labs(title = "Number of Keystone Taxa Across Different Thresholds",
             x = "Correlation Threshold", y = "Number of Keystone Taxa") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave("network/keystone_sensitivity_count.pdf", p1, width = 10, height = 6)
      
      # 2. Keystone score distribution across thresholds
      p2 <- ggplot(all_keystone, aes(x = factor(threshold), y = keystone_score, fill = group)) +
        geom_violin(alpha = 0.7) +
        geom_boxplot(width = 0.2, alpha = 0.8, position = position_dodge(0.9)) +
        facet_wrap(~ group, ncol = 2) +
        labs(title = "Keystone Score Distribution Across Thresholds",
             x = "Correlation Threshold", y = "Keystone Score") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
      
      ggsave("network/keystone_sensitivity_scores.pdf", p2, width = 12, height = 8)
      
      # 3. Keystone taxa overlap analysis
      overlap_analysis <- all_keystone %>%
        group_by(group) %>%
        summarise(
          overlap_0.2_0.3 = length(intersect(
            node[threshold == 0.2],
            node[threshold == 0.3]
          )),
          overlap_0.3_0.4 = length(intersect(
            node[threshold == 0.3],
            node[threshold == 0.4]
          )),
          overlap_0.2_0.4 = length(intersect(
            node[threshold == 0.2],
            node[threshold == 0.4]
          )),
          total_unique_0.2 = length(unique(node[threshold == 0.2])),
          total_unique_0.3 = length(unique(node[threshold == 0.3])),
          total_unique_0.4 = length(unique(node[threshold == 0.4]))
        ) %>%
        mutate(
          jaccard_0.2_0.3 = overlap_0.2_0.3 / (total_unique_0.2 + total_unique_0.3 - overlap_0.2_0.3),
          jaccard_0.3_0.4 = overlap_0.3_0.4 / (total_unique_0.3 + total_unique_0.4 - overlap_0.3_0.4)
        )
      
      write.table(overlap_analysis, "network/keystone_overlap_analysis.txt",
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # 4. Network topology changes across thresholds
      topology_summary <- all_keystone %>%
        group_by(group, threshold) %>%
        summarise(
          mean_degree = mean(degree),
          mean_betweenness = mean(betweenness),
          mean_closeness = mean(closeness),
          mean_transitivity = mean(transitivity),
          .groups = 'drop'
        )
      
      p3 <- ggplot(topology_summary %>%
                   pivot_longer(cols = c(mean_degree, mean_betweenness, mean_closeness, mean_transitivity),
                               names_to = "metric", values_to = "value"),
                 aes(x = factor(threshold), y = value, color = metric, group = metric)) +
        geom_line(size = 1) +
        geom_point(size = 2) +
        facet_wrap(~ group, scales = "free_y") +
        labs(title = "Network Topology Changes Across Thresholds",
             x = "Correlation Threshold", y = "Metric Value") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5),
              legend.position = "bottom")
      
      ggsave("network/topology_sensitivity.pdf", p3, width = 12, height = 8)
      
      cat("=== Keystone Taxa Sensitivity Analysis Completed ===\n")
      cat("Analyzed thresholds:", paste(thresholds, collapse = ", "), "\n")
      cat("\nKeystone Taxa Statistics:\n")
      print(keystone_summary)
      cat("\nKeystone Taxa Overlap Analysis:\n")
      print(overlap_analysis)
      cat("\nSensitivity analysis results saved to: network/keystone_taxa_sensitivity_analysis.txt\n")
      
      # Display robust keystone taxa (identified in multiple thresholds)
      robust_keystone <- all_keystone %>%
        group_by(node, group) %>%
        summarise(
          n_thresholds = n_distinct(threshold),
          thresholds = paste(sort(unique(threshold)), collapse = ","),
          .groups = 'drop'
        ) %>%
        filter(n_thresholds >= 2) %>%
        arrange(group, desc(n_thresholds))
      
      if (nrow(robust_keystone) > 0) {
        cat("\n=== Robust Keystone Taxa (identified in ≥2 thresholds) ===\n")
        print(robust_keystone)
        write.table(robust_keystone, "network/robust_keystone_taxa.txt",
                    sep = "\t", row.names = FALSE, quote = FALSE)
      }
    }
  } else {
    cat("Cannot compute keystone taxa, possibly insufficient data\n")
  }
} else {
  cat("Node centrality file does not exist\n")
}
