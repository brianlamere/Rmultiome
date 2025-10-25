#these are not functions.  They should be, and will be functions, but as of now...


#testing PC1 and PC2 for technical data, to set minimum dim for neighbors and UMAP
#used same for random other PCs and found low correlatoin similar to PC2 for our data

if (FALSE) {
pc1 <- Embeddings(merged_data, "pca")[, 1]
cor.test(pc1, merged_data$percent.mt)
cor.test(pc1, merged_data$nCount_RNA)
cor.test(pc1, merged_data$nFeature_RNA)
cor.test(pc1, merged_data$TSS.enrichment)
cor.test(pc1, merged_data$nCount_ATAC)
cor.test(pc1, merged_data$nFeature_ATAC)

#same for pc2
pc2 <- Embeddings(merged_data, "pca")[, 2]
cor.test(pc2, merged_data$percent.mt)
cor.test(pc2, merged_data$nCount_RNA)
cor.test(pc2, merged_data$nFeature_RNA)
cor.test(pc2, merged_data$TSS.enrichment)
cor.test(pc2, merged_data$nCount_ATAC)
cor.test(pc2, merged_data$nFeature_ATAC)
}

#read in the post-pipeline1, post-merge, post-RNA and ATAC processing, data.
merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
merged_obj <- harmonize_both(merged_data, harmony_max_iter = 50)


#theory: find elbow to set dimensionality for Neighbors.  But mormalized data 
#will throw this off and for this project the elbow made for ugly unclear graphs.
#it is an unclear metric, and I prefer parameter sweeps and objective metrics 
#taken from that, but here's your elbow - you can set max dims to the elbow if you
#want to not worry about how much of a massive difference these parameters make
findElbow(merged_data)

# Here we'll try a bunch of different combinations of a couple settings, to see
# what works better.  The elbow above can have whatever impact you'd like

merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
merged_obj <- harmonize_both(merged_data, harmony_max_iter = 50)


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

