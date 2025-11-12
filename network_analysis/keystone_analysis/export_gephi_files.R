
#!/usr/bin/env Rscript
#
# Export Network Files for Gephi Visualization
# This script exports network data in formats compatible with Gephi

library(dplyr)
library(igraph)

cat("=== Exporting Network Files for Gephi ===\n")

# Create Gephi output directory
gephi_dir <- "gephi"
if (!dir.exists(gephi_dir)) {
  dir.create(gephi_dir)
  cat("Created Gephi output directory:", gephi_dir, "\n")
}

# Read taxonomy data
taxonomy_data <- read.csv("asv_taxa.csv", header = TRUE, 
                         stringsAsFactors = FALSE, fileEncoding = "UTF-8")

if ("ASVID" %in% colnames(taxonomy_data)) {
  colnames(taxonomy_data)[colnames(taxonomy_data) == "ASVID"] <- "ASV"
  cat("Renamed ASVID column to ASV\n")
}

# Read node centrality data
if (file.exists("network/node_centrality_properties.txt")) {
  node_data <- read.table("network/node_centrality_properties.txt", 
                         header = TRUE, sep = "\t")
  
  # Use threshold 0.3 data
  node_data_0.3 <- node_data %>% filter(threshold == 0.3)
  
  # Generate Gephi files for each group
  groups <- c("none", "mild", "moderate", "severe")
  
  for(group in groups) {
    cat("\nProcessing group:", group, "\n")
    
    # Read SparCC correlation matrix
    sparcc_file <- paste0(group, "/cor_sparcc.norm_rarefy_", group, ".txt")
    
    if (!file.exists(sparcc_file)) {
      cat("  File does not exist:", sparcc_file, "\n")
      next
    }
    
    sparcc_cor <- as.matrix(read.table(sparcc_file, header = TRUE, 
                                      row.names = 1, sep = "\t"))
    
    # Get group-specific node data
    group_nodes <- node_data_0.3 %>% filter(group == !!group)
    
    if (nrow(group_nodes) == 0) {
      cat("  No node data for group:", group, "\n")
      next
    }
    
    # Ensure common ASVs
    common_asvs <- intersect(rownames(sparcc_cor), group_nodes$node)
    
    if (length(common_asvs) < 10) {
      cat("  Too few common ASVs:", length(common_asvs), "\n")
      next
    }
    
    # Filter matrices
    sparcc_cor <- sparcc_cor[common_asvs, common_asvs]
    group_nodes <- group_nodes %>% filter(node %in% common_asvs)
    
    # Create adjacency matrix with threshold 0.3
    adj_matrix <- ifelse(abs(sparcc_cor) > 0.3, 1, 0)
    diag(adj_matrix) <- 0
    
    # Create graph object
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
    
    # Remove isolated nodes
    isolated_nodes <- which(degree(g) == 0)
    if (length(isolated_nodes) > 0) {
      g <- delete.vertices(g, isolated_nodes)
    }
    
    if (vcount(g) == 0) {
      cat("  No edges after filtering isolated nodes\n")
      next
    }
    
    # Create nodes table for Gephi
    nodes_df <- data.frame(
      Id = V(g)$name,
      Label = V(g)$name,
      Degree = degree(g),
      Betweenness = betweenness(g),
      Closeness = closeness(g),
      Eigenvector = eigen_centrality(g)$vector,
      Group = group,
      stringsAsFactors = FALSE
    )
    
    # Add taxonomy information if available
    if (exists("taxonomy_data") && nrow(taxonomy_data) > 0) {
      nodes_df <- nodes_df %>%
        left_join(taxonomy_data, by = c("Id" = "ASV"))
    }
    
    # Create edges table for Gephi
    edges_df <- as_data_frame(g, what = "edges")
    if (nrow(edges_df) > 0) {
      edges_df <- edges_df %>%
        mutate(
          Type = "Undirected",
          Weight = 1,
          Correlation = sparcc_cor[cbind(from, to)],
          Source = from,
          Target = to
        ) %>%
        select(Source, Target, Type, Weight, Correlation)
    }
    
    # Export files
    nodes_file <- paste0(gephi_dir, "/", group, "_nodes.csv")
    edges_file <- paste0(gephi_dir, "/", group, "_edges.csv")
    
    write.csv(nodes_df, nodes_file, row.names = FALSE, fileEncoding = "UTF-8")
    write.csv(edges_df, edges_file, row.names = FALSE, fileEncoding = "UTF-8")
    
    cat("  Exported:", nodes_file, "\n")
    cat("  Exported:", edges_file, "\n")
    cat("  Network statistics - Nodes:", vcount(g), "Edges:", ecount(g), "\n")
  }
  
  cat("\n=== Gephi Export Completed ===\n")
  cat("All network files exported to:", gephi_dir, "\n")
  
} else {
  cat("Node centrality file does not exist\n")
  cat("Please ensure network analysis scripts have been run\n")
}
