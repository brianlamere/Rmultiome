
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

#' Define parameter sweep combinations
#' @param dims_range List of dimension ranges (e.g., list(c(2:30), c(2:40)))
#' @param knn_values Vector of k-nearest neighbors values
#' @param res_values Vector of resolution values
#' @return Data frame with all parameter combinations
define_parameter_sweep <- function(dims_range, knn_values, res_values) {
  # Create all combinations
  params <- expand.grid(
    dims_idx = seq_along(dims_range),
    knn = knn_values,
    res = res_values,
    stringsAsFactors = FALSE
  )

  # Add the actual dims vectors
  params$dims <- lapply(params$dims_idx, function(i) dims_range[[i]])

  # Create readable dims string for labeling
  params$dims_str <- sapply(params$dims, function(d) {
    sprintf("%d:%d", min(d), max(d))
  })

  # Remove temporary index column
  params$dims_idx <- NULL

  # Add unique ID for each parameter set
  params$param_id <- seq_len(nrow(params))

  return(params)
}

#' Run parameter sweep with metrics
#' @param seurat_obj Harmonized Seurat object (post-harmony, pre-FMMN)
#' @param dims_range List of dimension ranges
#' @param knn_values Vector of k values
#' @param res_values Vector of resolution values
#' @param save_dir Directory to save intermediate objects (optional)
#' @param compute_stability Whether to test clustering stability (slow)
#' @return Data frame with metrics for each parameter combination
run_parameter_sweep_with_metrics <- function(seurat_obj, dims_range, knn_values,
                                             res_values, save_dir = NULL,
                                             compute_stability = FALSE) {
  # Define parameter combinations
  params <- define_parameter_sweep(dims_range, knn_values, res_values)

  cat(sprintf("\n=== Running parameter sweep: %d combinations ===\n", nrow(params)))

  results <- list()

  for (i in seq_len(nrow(params))) {
    param_set <- params[i, ]
    cat(sprintf("\n[%d/%d] Testing: dims=%s, knn=%d, res=%.3f\n",
               i, nrow(params), param_set$dims_str, param_set$knn, param_set$res))

    # Run FMMN with this parameter set
    obj <- FindMultiModalNeighbors(
      object = seurat_obj,
      reduction.list = list("pca", "harmony"),
      dims.list = list(param_set$dims[[1]], param_set$dims[[1]]),
      k.nn = param_set$knn,
      knn.graph.name = "wknn",
      snn.graph.name = "wsnn",
      weighted.nn.name = "weighted.nn"
    )

    # Cluster with this resolution
    obj <- FindClusters(
      obj,
      graph.name = "wsnn",
      algorithm = 3,
      resolution = param_set$res,
      random.seed = 1984
    )

    # Compute UMAP for visualization
    obj <- RunUMAP(obj, reduction = "harmony", dims = param_set$dims[[1]])

    # Compute metrics
    metrics <- compute_cluster_metrics(obj, reduction = "harmony",
                                      dims = param_set$dims[[1]])

    # Optional: Test stability (slow, only for finalists)
    if (compute_stability) {
      stability <- test_clustering_stability(obj, param_set, n_runs = 3)
      metrics$stability_ari <- stability
    }

    # Combine parameters and metrics
    result <- c(as.list(param_set), metrics)
    results[[i]] <- result

    # Optionally save intermediate object
    if (!is.null(save_dir)) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
      filename <- sprintf("sweep_dims%s_k%d_r%.3f.rds",
                         gsub(":", "_", param_set$dims_str),
                         param_set$knn, param_set$res)
      saveRDS(obj, file.path(save_dir, filename))
    }

    cat(sprintf("  → %d clusters, modularity=%.3f, silhouette=%.3f\n",
               metrics$n_clusters, metrics$modularity, metrics$silhouette))
  }

  # Convert results to data frame
  results_df <- do.call(rbind, lapply(results, function(x) {
    as.data.frame(x, stringsAsFactors = FALSE)
  }))

  return(results_df)
}

#' Load a specific parameter sweep result
#' @param sweep_dir Directory where sweep objects were saved
#' @param dims_str Dimension string (e.g., "2:40")
#' @param knn k-nearest neighbors value
#' @param res Resolution value
#' @return Seurat object with those parameters
load_sweep_result <- function(sweep_dir, dims_str, knn, res) {
  filename <- sprintf("sweep_dims%s_k%d_r%.3f.rds",
                     gsub(":", "_", dims_str), knn, res)
  filepath <- file.path(sweep_dir, filename)

  if (!file.exists(filepath)) {
    stop(sprintf("Sweep result not found: %s", filepath))
  }

  readRDS(filepath)
}
