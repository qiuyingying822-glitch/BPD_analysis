#!/usr/bin/env Rscript
#
# Enhanced Keystone Taxa Identification
# This script identifies keystone taxa in microbial networks using multiple topological features

library(dplyr)
library(ggplot2)
library(tidyr)
library(igraph)

# Set working directory (adjust as needed)
# setwd("path/to/your/project")

cat("=== Enhanced Keystone Taxa Identification ===\n")

# Read node centrality data
if (file.exists("network/node_centrality_properties.txt")) {
  node_data <- read.table("network/node_centrality_properties.txt", 
                         header = TRUE, sep = "\t")
  
  # Use threshold 0.3 data for keystone identification
  node_data_0.3 <- node_data %>% filter(threshold == 0.3)
  
  if (nrow(node_data_0.3) > 0) {
    
    # Enhanced keystone taxa identification function
    identify_keystone_taxa_enhanced <- function(group_data, network_cor, top_percent = 0.1) {
      
      # Build network and calculate additional topological features
      adj_matrix <- ifelse(abs(network_cor) > 0.3, 1, 0)
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
      
      # Handle NA values (for nodes with degree 1, local transitivity might be NA)
      enhanced_data$transitivity[is.na(enhanced_data$transitivity)] <- 0
      
      # Calculate standardized scores - adjust weights based on observed patterns
      enhanced_data <- enhanced_data %>%
        mutate(
          # Positive features: high degree, high closeness, high transitivity
          degree_z = as.numeric(scale(degree)),
          closeness_z = as.numeric(scale(closeness)),
          transitivity_z = as.numeric(scale(transitivity)),
          eigen_z = as.numeric(scale(eigen_centrality)),
          # Negative feature: low betweenness (negative value)
          betweenness_z = -as.numeric(scale(betweenness))  # Negative sign indicates we want low betweenness
        )
      
      # Handle NA values
      enhanced_data[is.na(enhanced_data)] <- 0
      
      # Calculate comprehensive keystone score - adjust weights based on observed patterns
      enhanced_data <- enhanced_data %>%
        mutate(
          # Emphasize high degree, high closeness, high transitivity
          # Reduce weight for betweenness centrality (can even be negative weight)
          keystone_score = (2 * degree_z) +           # High degree - double weight
            (1.5 * closeness_z) +                     # High closeness - 1.5x weight
            (1.5 * transitivity_z) +                  # High transitivity - 1.5x weight
            (1 * eigen_z) +                           # Eigenvector centrality
            (-1 * betweenness_z)                      # Low betweenness (negative weight)
        )
      
      # Identify keystone taxa (top percent)
      n_keystone <- max(1, round(nrow(enhanced_data) * top_percent))
      keystone_taxa <- enhanced_data %>%
        arrange(desc(keystone_score)) %>%
        slice_head(n = n_keystone) %>%
        mutate(keystone_rank = 1:n())
      
      return(keystone_taxa)
    }
    
    # Identify keystone taxa for each group
    keystone_results <- list()
    
    for (group in unique(node_data_0.3$group)) {
      group_data <- node_data_0.3 %>% filter(group == !!group)
      
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
          
          keystone_taxa <- identify_keystone_taxa_enhanced(group_data, sparcc_cor, top_percent = 0.1)
          keystone_results[[group]] <- keystone_taxa
        }
      }
    }
    
    # Combine keystone taxa from all groups
    if (length(keystone_results) > 0) {
      all_keystone <- do.call(rbind, keystone_results)
      
      # Save keystone taxa results
      write.table(all_keystone, "network/keystone_taxa_enhanced.txt",
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Generate keystone taxa summary table
      keystone_summary <- all_keystone %>%
        group_by(group) %>%
        summarise(
          n_keystone = n(),
          avg_degree = mean(degree),
          avg_betweenness = mean(betweenness),
          avg_closeness = mean(closeness),
          avg_transitivity = mean(transitivity),
          avg_keystone_score = mean(keystone_score)
        )
      
      write.table(keystone_summary, "network/keystone_taxa_enhanced_summary.txt",
                  sep = "\t", row.names = FALSE, quote = FALSE)
      
      # Visualize keystone taxa features
      if (nrow(all_keystone) > 0) {
        # Feature distribution of keystone taxa
        feature_long <- all_keystone %>%
          select(group, node, degree, betweenness, closeness, transitivity, keystone_score) %>%
          pivot_longer(cols = c(degree, betweenness, closeness, transitivity, keystone_score),
                      names_to = "feature", values_to = "value")
        
        p1 <- ggplot(feature_long, aes(x = value, fill = group)) +
          geom_density(alpha = 0.5) +
          facet_grid(feature ~ group, scales = "free") +
          labs(title = "Topological Features of Enhanced Keystone Taxa",
               x = "Feature Value", y = "Density") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave("network/keystone_enhanced_features.pdf", p1, width = 12, height = 10)
        
        # Keystone score comparison
        p2 <- ggplot(all_keystone, aes(x = group, y = keystone_score, fill = group)) +
          geom_violin(alpha = 0.7) +
          geom_boxplot(width = 0.1, alpha = 0.8) +
          labs(title = "Enhanced Keystone Scores by Group",
               x = "Group", y = "Enhanced Keystone Score") +
          theme_bw() +
          theme(plot.title = element_text(hjust = 0.5))
        
        ggsave("network/keystone_enhanced_scores.pdf", p2, width = 10, height = 6)
        
        # Feature correlation heatmap
        feature_cor <- all_keystone %>%
          group_by(group) %>%
          summarise(
            deg_bet_cor = cor(degree, betweenness, method = "spearman"),
            deg_close_cor = cor(degree, closeness, method = "spearman"),
            deg_trans_cor = cor(degree, transitivity, method = "spearman"),
            bet_close_cor = cor(betweenness, closeness, method = "spearman")
          )
        
        write.table(feature_cor, "network/keystone_enhanced_correlations.txt",
                    sep = "\t", row.names = FALSE, quote = FALSE)
        
        cat("=== Enhanced Keystone Taxa Identification Completed ===\n")
        cat("Keystone Taxa Statistics:\n")
        print(keystone_summary)
        cat("\nNumber of Keystone Taxa per Group:\n")
        print(table(all_keystone$group))
        cat("\nEnhanced keystone taxa details saved to: network/keystone_taxa_enhanced.txt\n")
        
        # Display typical features of keystone taxa
        cat("\nTypical Features of Keystone Taxa (averages):\n")
        typical_features <- all_keystone %>%
          summarise(
            avg_degree = mean(degree),
            avg_betweenness = mean(betweenness),
            avg_closeness = mean(closeness),
            avg_transitivity = mean(transitivity)
          )
        print(typical_features)
      }
    } else {
      cat("Cannot compute enhanced keystone taxa, possibly insufficient data\n")
    }
  } else {
    cat("No node data found for threshold 0.3\n")
  }
} else {
  cat("Node centrality file does not exist\n")
}
