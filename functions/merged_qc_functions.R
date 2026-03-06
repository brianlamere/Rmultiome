
#' Check PC correlations with technical covariates
#' @param seurat_obj Merged Seurat object with PCA reduction
#' @param n_pcs Number of PCs to test (default 10)
#' @param reduction Name of reduction to use (default "pca")
#' @return ggplot object showing correlation heatmap
check_pc_technical_bias <- function(seurat_obj, n_pcs = 10, reduction = "pca") {
  require(ggplot2)
  require(reshape2)
  
  # Technical covariates to test
  covariates <- c("percent.mt", "nCount_RNA", "nFeature_RNA", 
                  "TSS.enrichment", "nCount_ATAC", "nFeature_ATAC")
  
  # Check which covariates exist
  available_covariates <- covariates[covariates %in% colnames(seurat_obj@meta.data)]
  
  # Get PC embeddings
  pcs <- Embeddings(seurat_obj, reduction = reduction)[, 1:n_pcs]
  
  # Compute correlations
  cor_matrix <- matrix(NA, nrow = n_pcs, ncol = length(available_covariates))
  rownames(cor_matrix) <- paste0("PC", 1:n_pcs)
  colnames(cor_matrix) <- available_covariates
  
  for (i in 1:n_pcs) {
    for (j in seq_along(available_covariates)) {
      cov_name <- available_covariates[j]
      cor_matrix[i, j] <- cor(pcs[, i], seurat_obj@meta.data[[cov_name]], 
                              use = "complete.obs")
    }
  }
  
  # Melt for ggplot
  cor_df <- melt(cor_matrix)
  colnames(cor_df) <- c("PC", "Covariate", "Correlation")
  
  # Create heatmap
  p1 <- ggplot(cor_df, aes(x = Covariate, y = PC, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("PC Correlation with Technical Covariates -", 
                       seurat_obj@project.name),
         x = "Technical Covariate", y = "Principal Component")
  
  # Create line plot (absolute correlations)
  cor_df$AbsCorrelation <- abs(cor_df$Correlation)
  p2 <- ggplot(cor_df, aes(x = PC, y = AbsCorrelation, 
                           color = Covariate, group = Covariate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Absolute Correlation by PC (Technical Bias Decay)",
         x = "Principal Component", y = "Absolute Correlation",
         color = "Covariate") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
    annotate("text", x = n_pcs * 0.8, y = 0.12, 
             label = "Threshold (|r| = 0.1)", color = "gray50")
  
  # Print summary stats
  cat("\n=== PC1 Correlations (Technical Bias Check) ===\n")
  pc1_cors <- cor_matrix[1, ]
  print(sort(abs(pc1_cors), decreasing = TRUE))
  
  cat("\n=== PC2 Correlations (Should be lower) ===\n")
  pc2_cors <- cor_matrix[2, ]
  print(sort(abs(pc2_cors), decreasing = TRUE))
  
  cat("\n=== Maximum Absolute Correlation by PC ===\n")
  max_cors <- apply(abs(cor_matrix), 1, max)
  print(max_cors)
  
  # Return both plots
  list(heatmap = p1, lineplot = p2, cor_matrix = cor_matrix)
}

#' Compute clustering quality metrics
#' @param seurat_obj Seurat object with clustering
#' @param reduction Which reduction to use for silhouette/distance
#' @return List of metric values
compute_cluster_metrics <- function(seurat_obj, reduction = "pca",
                                   dims = 2:30) {
  require(cluster)
  require(igraph)

  clusters <- Idents(seurat_obj)
  embeddings <- Embeddings(seurat_obj, reduction = reduction)[, dims]

  # Silhouette score (range -1 to 1, higher = better separated clusters)
  # BUT: Can favor over-clustering
  dist_matrix <- dist(embeddings)
  sil <- silhouette(as.numeric(clusters), dist_matrix)
  silhouette_avg <- mean(sil[, 3])

  # Modularity (from graph, higher = better community structure)
  # More biologically relevant for scRNA-seq
  graph <- seurat_obj@graphs$wsnn
  if (!is.null(graph)) {
    ig <- graph_from_adjacency_matrix(as.matrix(graph), mode = "undirected",
                                     weighted = TRUE)
    modularity_score <- modularity(ig, as.numeric(clusters))
  } else {
    modularity_score <- NA
  }

  # Cluster size distribution (detect over-fragmentation)
  cluster_sizes <- table(clusters)
  n_singletons <- sum(cluster_sizes == 1)
  size_cv <- sd(cluster_sizes) / mean(cluster_sizes)  # coefficient of variation

  # Number of clusters
  n_clusters <- length(unique(clusters))

  list(
    silhouette = silhouette_avg,
    modularity = modularity_score,
    n_clusters = n_clusters,
    n_singletons = n_singletons,
    size_cv = size_cv,
    median_cluster_size = median(cluster_sizes),
    min_cluster_size = min(cluster_sizes),
    max_cluster_size = max(cluster_sizes)
  )
}

#' Test clustering stability across multiple runs
#' @param seurat_obj Seurat object
#' @param params Single parameter set (dims, knn, res)
#' @param n_runs Number of times to re-run clustering
#' @return Adjusted Rand Index (ARI) between runs (higher = more stable)
test_clustering_stability <- function(seurat_obj, params, n_runs = 5) {
  require(mclust)  # for adjustedRandIndex

  clustering_results <- list()

  for (i in 1:n_runs) {
    # Re-run with different seed
    obj <- FindClusters(seurat_obj, graph.name = "wsnn",
                       resolution = params$res,
                       random.seed = i * 1000)
    clustering_results[[i]] <- Idents(obj)
  }

  # Compute pairwise ARI between all runs
  ari_scores <- c()
  for (i in 1:(n_runs - 1)) {
    for (j in (i + 1):n_runs) {
      ari <- adjustedRandIndex(clustering_results[[i]],
                              clustering_results[[j]])
      ari_scores <- c(ari_scores, ari)
    }
  }

  mean(ari_scores)  # Average stability across runs
}
