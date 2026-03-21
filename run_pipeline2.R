source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

pipeline1_settings <- read_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

# === Load merged data and all settings ===
merged_obj <- readRDS(file.path(rdsdir, "merged_preharmony.Rds"))

EnsDbAnnos <- loadannotations()

# === STEP 1: Run Harmony with saved settings (reproducible) ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = harmony_settings$random_seed,
  harmony_max_iter = harmony_settings$max_iter,
  harmony_project.dim = harmony_settings$project_dim,
  harmony_dims = harmony_settings$dims_use
)

# === STEP 2: Apply clustering settings ===
dims <- cluster_settings$dims_min:cluster_settings$dims_max

clustered_obj <- FMMN_task(harmony_obj,
                          dims = dims,
                          knn = cluster_settings$knn)

clustered_obj <- cluster_data(clustered_obj,
                             alg = cluster_settings$algorithm,
                             res = cluster_settings$resolution,
                             cluster_seed = cluster_settings$random_seed,
			     singleton_handling = "discard",
			     run_umap = TRUE)

# === STEP 3: Apply cell type labels ===
# Apply labels and remove flagged cell types in one call
clustered_obj <- apply_celltype_labels(
  clustered_obj,
  celltype_settings = celltype_settings,
  remove_flagged = TRUE,     # Remove clusters marked as "remove"
  verbose = TRUE
)

obj_assigned <- clustered_obj

# === Visualization ===
#DimPlot(obj_assigned, label = TRUE, raster = FALSE)
#DimPlot(obj_assigned, group.by = "celltype", label = TRUE, raster = FALSE)

# === STEP 4: Add experimental group metadata ===
obj_assigned$group <- NA
obj_assigned$group[obj_assigned$orig.ident %in% c("LG300", "LG301")] <- "No_HIV"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG22", "LG25", "LG38")] <- "Low"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG05", "LG26", "LG31")] <- "Acute"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG08", "LG23", "LG33")] <- "Chronic"

# Save the final annotated object
saveRDS(obj_assigned, file.path(rdsdir, "annotated_obj_final.rds"))
# obj_assigned <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

# === STEP 5: DE/DA Analysis ===
# CRITICAL: Join layers before differential expression testing
# Joining RNA layers
DefaultAssay(obj_assigned) <- "RNA"
obj_assigned <- JoinLayers(obj_assigned)

# Define cell types and comparisons for DE/DA
celltypes_list <- c("Oligodendrocytes", "Microglia", "Astrocytes")
comparisons_list <- list(
  c("Low", "No_HIV"),
  c("Acute", "Low"),
  c("Chronic", "Low")
)

# Run Differential Expression (PSEUDO-BULK)
for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    cat(sprintf("DE: %s - %s vs %s\n", celltype, comparison[1], comparison[2]))
    run_pseudobulk_DE_and_export(
      seurat_obj = obj_assigned,
      celltype_col = "celltype",
      celltype_value = celltype,
      sample_col = "orig.ident",
      group_col = "group",
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = file.path(project_outdir, "DiffExpress_results"),
      min_cells_per_sample = 50
    )
  }
}

# Run Differential Accessibility (PSEUDO-BULK)
cat("\n=== Running Differential Accessibility (Pseudo-bulk) ===\n")
for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    cat(sprintf("DA: %s - %s vs %s\n", celltype, comparison[1], comparison[2]))

    run_pseudobulk_DA_and_export(
      seurat_obj = obj_assigned,
      celltype_col = "celltype",
      celltype_value = celltype,
      sample_col = "orig.ident",
      group_col = "group",
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = file.path(project_outdir, "DiffAccess_results"),
      min_cells_per_sample = 50
    )
  }
}
