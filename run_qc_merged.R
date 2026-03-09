source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#load object created at end of run_pipeline1.R
merged_obj <- readRDS(file.path(rdsdir,"merged_preharmony.Rds"))

# ============================================================================
# PHASE 1: Clustering Parameter Selection
# ============================================================================

# === STEP 1: Check for technical bias in PCs ===
pc_check <- check_pc_technical_bias(merged_obj, n_pcs = 30)

maybe_new_device(width = 10, height = 8)
print(pc_check$heatmap)

maybe_new_device(width = 12, height = 6)
print(pc_check$lineplot)

# Decision: Exclude PC1 based on high technical correlation

# === STEP 2: Decide harmony settings ===
# Based on PC bias check, decide which PCs to use

harmony_settings <- list(
  dims_use = 2:50,
  random_seed = 1984,
  max_iter = 50,
  project_dim = FALSE
)

saveRDS(harmony_settings, harmony_settings_file))

# === STEP 3: Run Harmony for parameter sweep ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = harmony_settings$random_seed,
  harmony_max_iter = harmony_settings$max_iter,
  harmony_project.dim = harmony_settings$project_dim,
  harmony_dims = harmony_settings$dims_use
)

saveRDS(harmony_obj, file.path(rdsdir, "harmonized.rds"))
#harmony_obj <- readRDS(file.path(rdsdir, "harmonized.rds"))

# === STEP 4: Find elbow (informative but not prescriptive) ===
# Note: Elbow point is informative but not necessarily optimal for clustering.
maybe_new_device(width = 10, height = 6)
print(findElbow(harmony_obj))

# === STEP 5: Run parameter sweep with metrics ===
sweep_results <- run_parameter_sweep_with_metrics(
  seurat_obj = harmony_obj,
  dims_range = list(c(2:30), c(2:40), c(2:50)),
  knn_values = c(20, 30, 40, 50),
  res_values = c(0.03, 0.05, 0.1, 0.2),
  alg = 3,                   # SLM algorithm (required parameter)
  cluster_seed = 1984,       # Reproducibility (required parameter)
  save_dir = sweep_dir,
  compute_stability = FALSE  # Skip during initial sweep (too slow)
)

# Save full sweep results table
write.csv(sweep_results, file.path(tmpfiledir, "param_sweep_results_full.csv"),
         row.names = FALSE)

# === STEP 6: Filter based on metrics ===
filtered_results <- sweep_results %>%
  filter(n_singletons < 10) %>%
  filter(n_clusters >= 5 & n_clusters <= 20) %>%
  filter(modularity > 0.3) %>%
  arrange(desc(modularity)) %>%
  head(20)

cat(sprintf("\nFiltered to %d candidates (from %d total)\n",
           nrow(filtered_results), nrow(sweep_results)))

# === STEP 7: Visual inspection of finalists ===
for (i in seq_len(nrow(filtered_results))) {
  param <- filtered_results[i, ]

  obj <- load_sweep_result(sweep_dir, param$dims_str, param$knn, param$res)

  maybe_new_device_workspace(width = 8, height = 6, workspace = 9,
                            title = "R_ParamSweep")

  # Create labeled UMAP (removed silhouette from title)
  p <- DimPlot(obj, reduction = "wnn.umap", label = TRUE, label.size = 3) +
    ggtitle(sprintf("Rank #%d: dims=%s, knn=%d, res=%.2f\nn_clusters=%d, mod=%.3f",
                   i, param$dims_str, param$knn, param$res,
                   param$n_clusters, param$modularity)) +
    theme(plot.title = element_text(size = 10))

  print(p)
}

# === STEP 8: Manual selection ===
chosen_rank <- 3
chosen_params <- filtered_results[chosen_rank, ]

cat("\n=== Selected Parameters ===\n")
cat(sprintf("dims: %s\n", chosen_params$dims_str))
cat(sprintf("knn: %d\n", chosen_params$knn))
cat(sprintf("resolution: %.3f\n", chosen_params$res))
cat(sprintf("n_clusters: %d\n", chosen_params$n_clusters))
cat(sprintf("modularity: %.3f\n", chosen_params$modularity))

# === STEP 9: Save cluster settings ===
cluster_settings <- data.frame(
  dims_min = min(chosen_params$dims[[1]]),
  dims_max = max(chosen_params$dims[[1]]),
  knn = chosen_params$knn,
  resolution = chosen_params$res,
  algorithm = 3,  # SLM algorithm
  random_seed = 1984  # For reproducibility
)

saveRDS(cluster_settings, cluster_settings_file)


# ============================================================================
# PHASE 2: Cell Type Annotation
# ============================================================================

# === STEP 10: Find marker genes for all clusters ===
DefaultAssay(chosen_obj) <- "RNA"

all_markers <- FindAllMarkers(chosen_obj, 
                             only.pos = TRUE,
                             min.pct = 0.25,
                             logfc.threshold = 0.25)

# Save all markers for reference
write.csv(all_markers, 
          file.path(tmpfiledir, "all_cluster_markers.csv"),
          row.names = FALSE)

cat(sprintf("All markers saved to: %s\n", 
           file.path(tmpfiledir, "all_cluster_markers.csv")))

# Show top markers per cluster
cat("\nTop 5 markers per cluster:\n")
top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

print(top_markers)

# === STEP 11: Visualize marker genes ===
cat("\nStep 9: Visualizing marker expression...\n")

# Example: Plot known brain cell type markers
brain_markers <- c(
  "SLC17A7",  # Excitatory neurons
  "GAD1",     # Inhibitory neurons
  "MBP",      # Oligodendrocytes
  "GFAP",     # Astrocytes
  "AIF1",     # Microglia
  "PDGFRA"    # OPCs
)

maybe_new_device(width = 12, height = 10)
print(DotPlot(chosen_obj, features = brain_markers) + 
      coord_flip() +
      ggtitle("Known Brain Cell Type Markers"))

maybe_new_device(width = 14, height = 10)
print(FeaturePlot(chosen_obj, features = brain_markers, ncol = 3))

# === STEP 12: Manual cell type assignment ===

# USER INPUT: Edit this mapping based on marker genes
# This is an EXAMPLE - adjust based on your actual markers
celltype_mapping <- data.frame(
  cluster = 0:8,  # Adjust based on actual number of clusters
  celltype = c(
    "Excitatory_Neurons",    # Cluster 0
    "Inhibitory_Neurons",    # Cluster 1
    "Oligodendrocytes",      # Cluster 2
    "Astrocytes",            # Cluster 3
    "Microglia",             # Cluster 4
    "OPCs",                  # Cluster 5
    "Endothelial",           # Cluster 6
    "Pericytes",             # Cluster 7
    "Ambiguous"              # Cluster 8 - mark for removal
  ),
  action = c(
    rep("keep", 8),
    "remove"  # Remove ambiguous cluster
  ),
  stringsAsFactors = FALSE
)

#Current cell type mapping:
print(celltype_mapping)

# === STEP 13: Save cell type mapping ===
saveRDS(celltype_mapping, celltype_settings_file)


# === STEP 14: Preview cell type assignment ===

# Apply cell types to preview
cluster_ids <- as.numeric(as.character(Idents(chosen_obj)))
chosen_obj$celltypes <- celltype_mapping$celltype[match(cluster_ids, 
                                                        celltype_mapping$cluster)]

# Plot with cell type labels
maybe_new_device(width = 10, height = 8)
Idents(chosen_obj) <- chosen_obj$celltypes
print(DimPlot(chosen_obj, reduction = "umap", label = TRUE) +
      ggtitle("Cell Type Assignments (Preview)"))

# Show cluster sizes by cell type
cat("\nCells per cell type:\n")
print(table(chosen_obj$celltypes))

# ============================================================================
# Summary
# ============================================================================

cat("\n\n=== QC MERGED COMPLETE ===\n")
cat("Outputs saved:\n")
cat(sprintf("  1. Cluster settings: %s\n", 
           file.path(project_export, "cluster_settings.csv")))
cat(sprintf("  2. Cell type mapping: %s\n", 
           file.path(project_export, "celltype_mapping.csv")))
cat("\nNext step: Run run_pipeline2.R to apply these settings and perform DE/DA\n")

cat("\nTop candidates:\n")
print(filtered_results[, c("dims_str", "knn", "res", "n_clusters", "modularity")])
