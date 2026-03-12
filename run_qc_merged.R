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

saveRDS(harmony_settings, harmony_settings_file)

# === STEP 3: Run Harmony for parameter sweep ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = harmony_settings$random_seed,
  harmony_max_iter = harmony_settings$max_iter,
  harmony_project.dim = harmony_settings$project_dim,
  harmony_dims = harmony_settings$dims_use
)

saveRDS(harmony_obj, file.path(rdsdir, "harmonized.rds"))
#next line is only present due to debugging the code and restarting here
#harmony_obj <- readRDS(file.path(rdsdir, "harmonized.rds"))

# === STEP 4: Find elbow (informative but not prescriptive) ===
# Note: Elbow point is informative but not necessarily optimal for clustering.
maybe_new_device(width = 10, height = 6)
print(findElbow(harmony_obj))

#BELOW IS VERY INTERESTING, DETERMINE WHAT TOP CELLS IN THE TRANSITIONING GROUPS ARE
# === STEP 5: Run parameter sweep with metrics ===
sweep_results <- run_parameter_sweep_plots(
  seurat_obj = harmony_obj,
  dims_range = list(c(1:43),c(1:44),c(1:45)),
  knn_values = c(40),
  res_values = c(0.16,0.17,0.18),
  alg = 3,                   # SLM algorithm (required parameter)
  cluster_seed = 1984       # Reproducibility (required parameter)
)

# Save full sweep results table
write.csv(sweep_results, file.path(tmpfiledir, "param_sweep_results.csv"),
         row.names = FALSE)

cat("\n=== PARAMETER SWEEP COMPLETE ===\n")
cat("Results table saved to:", file.path(tmpfiledir, "param_sweep_results.csv"), "\n")
cat("\nAll plots displayed on workspace 9. Review them and choose your parameters.\n\n")

# Print results sorted by cluster count
cat("Results sorted by cluster count:\n")
print(sweep_results[order(sweep_results$n_clusters), ])

# === STEP 6: USER INPUT - Set chosen parameters ===
cat("\n=== SET YOUR CHOSEN PARAMETERS BELOW ===\n")

# After reviewing plots, set these to your chosen values:
chosen_dims_min <- 1      # CHANGE THIS
chosen_dims_max <- 44     # CHANGE THIS
chosen_knn <- 40          # CHANGE THIS
chosen_resolution <- 0.16 # CHANGE THIS

# === STEP 7: Save cluster settings ===
cluster_settings <- data.frame(
  dims_min = chosen_dims_min,
  dims_max = chosen_dims_max,
  knn = chosen_knn,
  resolution = chosen_resolution,
  algorithm = 3,
  random_seed = 1984
)

saveRDS(cluster_settings, cluster_settings_file)

# hack code keeping briefly
# fmmn_obj <- FMMN_task(harmony_obj, knn = 40, dims = 1:44)
# clustered_obj <- cluster_data(fmmn_obj, alg = 3, res = 0.16,
#                            cluster_seed = 1984,
# 			   singleton_handling = "keep")
# DimPlot(premap_obj,label=T, raster=FALSE)

# === STEP 8: Create final clustered object from settings ===
cat("\n=== STEP 8: Applying cluster settings ===\n")

# Read the settings you just saved
dims <- cluster_settings$dims_min:cluster_settings$dims_max

cat(sprintf("Clustering with: dims=%d:%d, knn=%d, res=%.3f\n",
           cluster_settings$dims_min, cluster_settings$dims_max,
           cluster_settings$knn, cluster_settings$resolution))

# Apply FMMN
cat("Running FindMultiModalNeighbors...\n")
chosen_obj <- FMMN_task(harmony_obj,
                       knn = cluster_settings$knn,
                       dims = dims)

# Apply clustering
cat("Clustering...\n")
chosen_obj <- cluster_data(
  chosen_obj,
  alg = cluster_settings$algorithm,
  res = cluster_settings$resolution,
  cluster_seed = cluster_settings$random_seed,
  singleton_handling = "keep",
  run_umap = TRUE
)

# Quick summary
n_clusters <- length(unique(Idents(chosen_obj))) -
              (if("singleton" %in% levels(Idents(chosen_obj))) 1 else 0)
n_singletons <- sum(Idents(chosen_obj) == "singleton")

cat(sprintf("Result: %d clusters, %d singletons (%.2f%%)\n",
           n_clusters, n_singletons,
           100 * n_singletons / ncol(chosen_obj)))

# Quick viz
maybe_new_device(width = 10, height = 8)
print(DimPlot(chosen_obj, reduction = "wnn.umap", label = TRUE, raster = FALSE))

# Save
saveRDS(chosen_obj, file.path(rdsdir, "chosen_clustered_obj.rds"))
cat(sprintf("Saved to: %s\n", file.path(rdsdir, "chosen_clustered_obj.rds")))

# ============================================================================
# PHASE 2: Cell Type Annotation
# ============================================================================

# === STEP 9: Define and save marker panels ===
cat("\n=== STEP 9: Setting up marker panels ===\n")

# Define markers for YOUR tissue (stressed brain nuclei)
marker_panels <- data.frame(
  celltype = c(
    rep("Oligodendrocytes", 6),
    rep("OPCs", 4),
    rep("Astrocytes", 6),
    rep("Microglia", 5),
    rep("Excitatory_Neurons", 4),
    rep("Inhibitory_Neurons", 6),
    rep("Endothelial", 4),
    rep("Stress_markers", 7)
  ),
  gene = c(
    # Oligodendrocytes (from your manual_mapping2.R)
    "MBP", "MOBP", "PLP1", "MOG", "MAG", "CNP",
    # OPCs
    "VCAN", "PDGFRA", "PCDH15", "OLIG1",
    # Astrocytes
    "GFAP", "ALDH1L1", "GLUL", "AQP4", "SLC1A2", "SLC4A4",
    # Microglia
    "CSF1R", "APBB1IP", "P2RY12", "AIF1", "CX3CR1",
    # Excitatory neurons
    "SLC17A7", "CAMK2A", "SATB2", "TBR1",
    # Inhibitory neurons (from your manual_mapping2.R)
    "GAD1", "GAD2", "SLC32A1", "VIP", "SST", "PVALB",
    # Endothelial
    "FLT1", "CLDN5", "KDR", "VWF",
    # Stress markers (important for YOUR study!)
    "FOS", "JUN", "JUNB", "HSPA1A", "HSPA1B", "ATF3", "DDIT3"
  ),
  notes = c(
    rep("Core oligodendrocyte markers", 6),
    rep("Oligodendrocyte precursors", 4),
    rep("Core astrocyte markers", 6),
    rep("Core microglia markers", 5),
    rep("Excitatory neuron markers", 4),
    rep("Inhibitory neuron markers", 6),
    rep("Endothelial markers", 4),
    rep("Stress/activation markers - critical for this study", 7)
  ),
  stringsAsFactors = FALSE
)

# Save as project data
write.csv(marker_panels,
         file.path(project_outdir, "marker_panels.csv"),
         row.names = FALSE)

cat(sprintf("Marker panels saved to: %s\n",
           file.path(project_outdir, "marker_panels.csv")))
cat("\nWARNING: These markers are optimized for stressed brain nuclei.\n")
cat("For other tissue types, edit marker_panels.csv before running.\n")

# Also save as RDS for easy loading
saveRDS(marker_panels, file.path(project_outdir, "marker_panels.rds"))

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
