parameter_sweep <- function(merged_obj, param_grid, resolutions, output_dir, plot_dir) {
  results <- list()
  result_counter <- 1
  init_df <- data.frame(
    dims_min = numeric(0),
    dims_max = numeric(0),
    knn = numeric(0),
    resolution = numeric(0),
    variance_total = numeric(0),
    n_clusters = numeric(0),
    n_singletons = numeric(0)
  )
  
  # Ensure directories exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
  
  csv_file <- file.path(output_dir, "results_debug.csv")
  
  # If the CSV exists, rename with timestamp
  if (file.exists(csv_file)) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    file.rename(csv_file, file.path(output_dir, paste0("results_debug_", timestamp, ".csv")))
  }
  
  # Write header to new CSV
  write.table(init_df, file=csv_file, sep=",", row.names=FALSE, col.names=TRUE, append=FALSE)
  
  for (i in 1:nrow(param_grid)) {
    dims_min <- param_grid$dims_min[i]
    dims_max <- param_grid$dims_max[i]
    dims <- dims_min:dims_max
    knn <- param_grid$knn[i]
    
    obj_neighbors <- FMMN_task(
      merged_obj,
      dims_pca = dims,
      dims_harmony = dims,
      knn = knn
    )
    DefaultAssay(obj_neighbors) <- "RNA"
    obj_clustered <- RunPCA(obj_neighbors, assay = "RNA",
                            features = VariableFeatures(obj_neighbors))
    
    for (resolution in resolutions) {
      obj_clustered <- obj_neighbors
      
      obj_clustered <- cluster_data(obj_clustered, alg = 3, res = resolution,
                                    cluster_dims = dims)
      # Save DimPlot
      plot_file <- file.path(
        plot_dir,
        paste0("DimPlot_dims", dims_min, "-", dims_max, "_knn", knn, "_res", resolution, ".png")
      )
      png(plot_file)
      print(DimPlot(obj_clustered, reduction = "wnn.umap",
                    group.by = "seurat_clusters", raster = FALSE))
      dev.off()
      
      print(DimPlot(obj_clustered, reduction = "wnn.umap",
                    group.by = "seurat_clusters", raster = FALSE)) +
        ggtitle(plot_file)
      
      variance <- Stdev(obj_clustered[["pca"]])^2
      jackstraw <- obj_clustered[["pca"]]@jackstraw$overall.p.values
      cluster_assignments <- obj_clustered@meta.data$seurat_clusters
      singleton_count <- sum(cluster_assignments == "singleton")
      cluster_count <- length(setdiff(unique(cluster_assignments), "singleton"))
      
      results[[result_counter]] <- list(
        params = data.frame(
          dims_min = dims_min,
          dims_max = dims_max,
          knn = knn,
          resolution = resolution
        ),
        variance = variance,
        jackstraw = jackstraw,
        n_clusters = cluster_count,
        n_singletons = singleton_count
      )
      row_df <- data.frame(
        dims_min = dims_min,
        dims_max = dims_max,
        knn = knn,
        resolution = resolution,
        variance_total = sum(variance[dims_min:dims_max]),
        n_clusters = cluster_count,
        n_singletons = singleton_count
      )
      
      # Append row to CSV (header omitted on appends)
      write.table(row_df, file = csv_file,
                  append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE)
      result_counter <- result_counter + 1
    }
  }
  invisible(results)
}