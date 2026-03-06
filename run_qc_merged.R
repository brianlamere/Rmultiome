source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

# Load harmonized data
harmony_obj <- readRDS(file.path(rdsdir, "harmonized.rds"))

# === STEP 1: Diagnostic checks ===
# Check for technical bias in PCs
pc_check <- check_pc_technical_bias(harmony_obj, n_pcs = 10)
maybe_new_device(width = 10, height = 8)
print(pc_check$heatmap)

# Find elbow (informative but not prescriptive)
elbow <- findElbow(harmony_obj)
cat(sprintf("Elbow detected at PC%d\n", elbow))

# === STEP 2: Define parameter sweep ranges ===
sweep_params <- define_parameter_sweep(
  dims_range = list(c(2:30), c(2:40), c(2:50)),
  knn_values = c(20, 30, 40, 50),
  res_values = c(0.03, 0.04, 0.05, 0.1, 0.2),
  exclude_pc1 = TRUE  # Based on PC bias check
)

# === STEP 3: Run parameter sweep ===
sweep_results <- run_parameter_sweep(
  seurat_obj = harmony_obj,
  params = sweep_params,
  save_dir = file.path(rdsdir, "param_sweep"),
  metrics = c("silhouette", "modularity", "stability")
)

# === STEP 4: Evaluate results ===
best_params <- evaluate_sweep_results(sweep_results)
print(best_params)



results <- list()
result_counter <- 1
# Create an empty data.frame with the desired columns
init_df <- data.frame(
  dims_min = numeric(0),
  dims_max = numeric(0),
  knn = numeric(0),
  resolution = numeric(0),
  variance_total = numeric(0),
  n_clusters = numeric(0),
  n_singletons = numeric(0)
)

# Write the header to the CSV (overwrites any existing file)
write.table(init_df, file="results_debug.csv", sep=",",
            row.names=FALSE, col.names=TRUE)

param_grid <- expand.grid(
  dims_min = c(2),
  dims_max = c(40),
  knn = c(40)
)

resolutions <- c(0.06, 0.1)

for (i in 1:nrow(param_grid)) {
  # Set parameters for FMMN
  dims_min <- param_grid$dims_min[i]
  dims_max <- param_grid$dims_max[i]
  dims <- dims_min:dims_max
  knn <- param_grid$knn[i]
  
  # Start with a fresh copy for neighbors
  obj_neighbors <- FMMN_task(merged_obj,
                             dims_pca = dims,
                             dims_harmony = dims,
                             knn = knn)
  DefaultAssay(obj_clustered) <- "RNA"
  obj_clustered <- RunPCA(obj_clustered, assay = "RNA",
                          features = VariableFeatures(obj_clustered))
  
  for (resolution in resolutions) {
    # Optionally: clone the object to ensure clustering doesn't interfere
    obj_clustered <- obj_neighbors # you can use the same object

    obj_clustered <- cluster_data(obj_clustered, alg = 3, res = resolution,
                                  cluster_dims = dims)
    # Save DimPlot
    plot_file <- paste0("/projects/opioid/parameter_sweep/DimPlot_dims",
                        dims_min, "-", dims_max, "_knn",
                        knn, "_res", resolution, ".png")
    png(plot_file)
    print(DimPlot(obj_clustered, reduction = "wnn.umap",
                  group.by = "seurat_clusters", raster = FALSE))
    dev.off()
    
    print(DimPlot(obj_clustered, reduction = "wnn.umap",
                  group.by = "seurat_clusters", raster = FALSE)) +
      ggtitle(plot_file)
    
    # Collect metrics
    variance <- Stdev(obj_clustered[["pca"]])^2
    jackstraw <- obj_clustered[["pca"]]@jackstraw$overall.p.values
    cluster_assignments <- obj_clustered@meta.data$seurat_clusters
    singleton_count <- sum(cluster_assignments == "singleton")
    cluster_count <- length(setdiff(unique(cluster_assignments), "singleton"))
    
    # Store results uniquely for each cluster run
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

    # Write to CSV (append mode)
    write.table(row_df, file="/projects/opioid/parameter_sweep/results_debug.csv",
                append=TRUE, sep=",", row.names=FALSE, col.names=FALSE)
    result_counter <- result_counter + 1
  }
}

# Save log file
write.csv(do.call(rbind, lapply(results, function(x) 
  data.frame(x$params, n_clusters=x$n_clusters, n_singletons=x$n_singletons))),
  "parameter_sweep_results1.csv")

# Summarize all results
summary_df <- summarize_results(results)

# Objective recommendation: best variance & most significant PCs
best_by_variance <- summary_df[which.max(summary_df$variance_total), ]
best_by_jackstraw <- summary_df[which.max(summary_df$jackstraw_sig), ]
# Optionally, filter out excessive singletons
summary_df_clean <- subset(summary_df, n_singletons < 0.1 * n_clusters)
best_balanced <- summary_df_clean[which.max(summary_df_clean$variance_total), ]

# Print recommendations
print("Best by variance:")
print(best_by_variance)
print("Best by jackstraw significance:")
print(best_by_jackstraw)
print("Best balanced (few singletons):")
print(best_balanced)

