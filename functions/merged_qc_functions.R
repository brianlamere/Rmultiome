
#' Check PC correlations with technical covariates
#' @param seurat_obj Merged Seurat object with PCA reduction
#' @param n_pcs Number of PCs to test (default 10)
#' @param reduction Name of reduction to use (default "pca")
#' @return ggplot object showing correlation heatmap
check_pc_technical_bias <- function(seurat_obj, n_pcs = 10, reduction = "pca") {
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

#' Run parameter sweep and display plots for visual inspection
#' @param seurat_obj Harmonized Seurat object
#' @param dims_range List of dimension ranges
#' @param knn_values Vector of k values
#' @param res_values Vector of resolution values
#' @param alg Clustering algorithm
#' @param cluster_seed Random seed
#' @return Data frame with basic results (for reference only)
run_parameter_sweep_plots <- function(seurat_obj, dims_range, knn_values,
                                     res_values, alg, cluster_seed) {

  n_combos <- length(dims_range) * length(knn_values) * length(res_values)
  cat(sprintf("\n=== Parameter Sweep: %d combinations ===\n", n_combos))
  cat("Plots will be displayed on workspace 9 (Hyprland)\n\n")

  results <- list()
  counter <- 1

  # OUTER LOOP: dims + knn (run FMMN once per combination)
  for (dims_idx in seq_along(dims_range)) {
    dims <- dims_range[[dims_idx]]
    dims_min <- min(dims)
    dims_max <- max(dims)
    dims_str <- sprintf("%d-%d", dims_min, dims_max)

      for (knn in knn_values) {
        cat(sprintf("\n=== FMMN: dims=%s, knn=%d ===\n",
                   gsub("-", ":", dims_str), knn))

        # NOW passing dims to FMMN_task!
        obj_fmmn <- FMMN_task(seurat_obj, knn = knn, dims = dims)

      # INNER LOOP: resolutions (reuse FMMN result)
      for (res in res_values) {
        cat(sprintf("[%d/%d] dims=%s, knn=%d, res=%.3f ... ",
                   counter, n_combos, dims_str, knn, res))

        # Copy FMMN result and cluster
        obj_clustered <- obj_fmmn
        obj_clustered <- cluster_data(
          obj_clustered,
          alg = alg,
          res = res,
          cluster_seed = cluster_seed,
          singleton_handling = "keep",
          run_umap = TRUE  # Need UMAP for plotting
        )

        # Get cluster info
        cluster_assignments <- obj_clustered@meta.data$seurat_clusters
        singleton_count <- sum(cluster_assignments == "singleton")
        cluster_count <- length(setdiff(unique(cluster_assignments), "singleton"))

        cat(sprintf("%d clusters (%d singletons)\n", cluster_count, singleton_count))

        # Store results
        results[[counter]] <- list(
          dims_min = dims_min,
          dims_max = dims_max,
          dims_str = dims_str,
          knn = knn,
          resolution = res,
          n_clusters = cluster_count,
          n_singletons = singleton_count
        )

        # Display plot on workspace 9
        maybe_new_device_workspace(width = 8, height = 6, workspace = 9,
                                   title = "R_ParamSweep")

        plot_title <- sprintf("dims=%s, knn=%d, res=%.3f\n%d clusters (%d singletons)",
                             dims_str, knn, res, cluster_count, singleton_count)

        p <- DimPlot(obj_clustered, reduction = "wnn.umap",
                    group.by = "seurat_clusters", label = TRUE, raster = FALSE) +
          ggtitle(plot_title)

        print(p)

        # Clean up
        rm(obj_clustered, p)
        gc(verbose = FALSE)

        counter <- counter + 1
      }

      # Clean up FMMN result after all resolutions
      rm(obj_fmmn)
      gc(verbose = FALSE, full = TRUE)
    }
  }

  # Convert results to data frame for reference
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      dims_min = x$dims_min,
      dims_max = x$dims_max,
      knn = x$knn,
      resolution = x$resolution,
      n_clusters = x$n_clusters,
      n_singletons = x$n_singletons,
      stringsAsFactors = FALSE
    )
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

#' Check system memory availability
#' @return Available memory in GB
check_available_memory <- function() {
  # Read /proc/meminfo
  meminfo <- readLines("/proc/meminfo")

  # Extract MemAvailable
  avail_line <- grep("^MemAvailable:", meminfo, value = TRUE)
  avail_kb <- as.numeric(gsub("^MemAvailable:\\s+(\\d+).*", "\\1", avail_line))

  avail_gb <- avail_kb / 1024 / 1024
  return(avail_gb)
}

#' Warn if memory is getting low
check_memory_threshold <- function(threshold_gb = 50) {
  avail <- check_available_memory()
  if (avail < threshold_gb) {
    warning(sprintf("Low memory warning: Only %.1f GB available (threshold: %d GB)",
                   avail, threshold_gb))
  }
  return(avail)
}
