source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

trimming_settings <- read_trimming_settings(trimming_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

# === Load merged data and all settings ===
merged_obj <- readRDS(file.path(rdsdir, "merged_preharmony.Rds"))

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

# That's it! Labels applied, bad cells removed, ready for analysis
obj_assigned <- clustered_obj

# === Visualization ===
DimPlot(clustered_obj, group.by = "celltypes", label = TRUE, raster = FALSE)

# Filter out unassigned clusters (cells with empty celltype, if any remain)
assigned_cells <- WhichCells(clustered_obj, expression = celltypes != "")
obj_assigned <- subset(clustered_obj, cells = assigned_cells)

# Verify
DimPlot(obj_assigned, label = TRUE, raster = FALSE)
DimPlot(obj_assigned, group.by = "celltypes", label = TRUE, raster = FALSE)

# === STEP 4: Add experimental group metadata ===
obj_assigned$group <- NA
obj_assigned$group[obj_assigned$orig.ident %in% c("LG300", "LG301")] <- "No_HIV"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG22", "LG25", "LG38")] <- "Low"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG05", "LG26", "LG31")] <- "Acute"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG08", "LG23", "LG33")] <- "Chronic"

# Save the final annotated object
saveRDS(obj_assigned, file.path(rdsdir, "annotated_obj_final.rds"))

# === STEP 5: DE/DA Analysis ===

# CRITICAL: Join layers before differential expression testing
cat("\nJoining layers for DE/DA testing...\n")
obj_assigned <- JoinLayers(obj_assigned)

# Define cell types and comparisons for DE/DA
celltypes_list <- c("Oligodendrocytes", "Microglia", "Astrocytes")  # UPDATE to match your actual cell type names!
comparisons_list <- list(
  c("No_HIV", "Low"),
  c("Low", "Acute"),
  c("Low", "Chronic")
)

# Run Differential Expression
cat("\n=== Running Differential Expression ===\n")
for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    cat(sprintf("DE: %s - %s vs %s\n", celltype, comparison[1], comparison[2]))
    run_DiffExpress_and_export(
      seurat_obj = obj_assigned,        # Use obj_assigned consistently
      celltype_col = "celltypes",
      celltype = celltype,
      group_col = "group",
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = file.path(project_outdir, "DiffExpress_results")
    )
  }
}

# Run Differential Accessibility
cat("\n=== Running Differential Accessibility ===\n")
for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    cat(sprintf("DA: %s - %s vs %s\n", celltype, comparison[1], comparison[2]))
    run_DiffAccess_and_export(
      seurat_obj = obj_assigned,        # Use obj_assigned consistently
      celltype_col = "celltypes",
      celltype = celltype,
      group_col = "group",
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = file.path(project_outdir, "DiffAccess_results")
    )
  }
}

#================old code
# === STEP 3: Apply cell type labels ===
cluster_ids <- as.numeric(as.character(Idents(clustered_obj)))
clustered_obj$celltypes <- celltype_mapping$celltype[
  match(cluster_ids, celltype_mapping$cluster)
]

remove_celltypes <- celltype_mapping$celltype[celltype_mapping$action == "remove"]
if (length(remove_celltypes) > 0) {
  clustered_obj <- subset(clustered_obj, subset = celltypes %in% remove_celltypes,
                         invert = TRUE)
}

Idents(clustered_obj) <- clustered_obj$celltypes

# === STEP 4: DE/DA ===

#major restructuring occuring, below down is no longer valid!  
p <- DimPlot(labeled_obj, group.by = "celltypes", label = TRUE, raster=FALSE)

p + coord_cartesian(xlim = c(-17, 10))

#the below filters based on unassigned clusters.  Comment this out to not do that.
assigned_cells <- WhichCells(labeled_obj, expression = celltypes != "")
obj_assigned <- subset(labeled_obj, cells = assigned_cells)

#to check
DimPlot(obj_assigned,label=T, raster=FALSE)
DimPlot(obj_assigned, group.by = "celltypes", label = TRUE, raster=FALSE)

obj_assigned$group <- NA
obj_assigned$group[obj_assigned$orig.ident %in% c("LG300", "LG301")] <- "No_HIV"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG22", "LG25", "LG38")] <- "Low"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG05", "LG26", "LG31")] <- "Acute"
obj_assigned$group[obj_assigned$orig.ident %in% c("LG08", "LG23", "LG33")] <- "Chronic"

saveRDS(obj_assigned, "/projects/opioid/vault96/tagged_dim240_knn40_res0.04_seed1984.rds")
tagged_obj <- obj_assigned


celltypes_list <- c("Oligodendrocyte", "Microglia", "Astrocyte")
comparisons_list <- list(
  c("No_HIV", "Low"),
  c("Low", "Acute"),
  c("Low", "Chronic")
)


for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    #print(paste("For cell type: ", celltype, "first is ", comparison[1], " and second is ", comparison[2]))
    run_DiffExpress_and_export(
      seurat_obj = obj_assigned,
      celltype_col = "celltypes",    # Your cell type column name
      celltype = celltype,
      group_col = "group",           # Update if your grouping column is named differently
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = "DiffExpress_results"
    )
  }
}


for (celltype in celltypes_list) {
  for (comparison in comparisons_list) {
    #print(paste("For cell type: ", celltype, "first is ", comparison[1], " and second is ", comparison[2]))
    run_DiffAccess_and_export(
      seurat_obj = tagged_obj,
      celltype_col = "celltypes",    # Your cell type column name
      celltype = celltype,
      group_col = "group",           # Update if your grouping column is named differently
      ident.1 = comparison[1],
      ident.2 = comparison[2],
      output_prefix = "DiffAccess_results"
    )
  }
}

##from here below it gets even more string-of-consciousness as things were done
# in multiple tabs.  Above takes you from raw data to a fully integrated object,
# then you do DE and DA on the object

# see if problem areas follow a particular sample.
DimPlot(obj_assigned, split.by = "orig.ident", raster = FALSE)

markers_oligo_low_acute[order(-markers_oligo_low_acute$p_val_adj), ][1:10, ]

cells_of_type <- WhichCells(labeled_obj, expression = celltypes == "Oligodendrocyte")
expr_matrix <- GetAssayData(labeled_obj, slot = "data")[, cells_of_type]

# Count genes with any expression in this cell type
sum(rowSums(expr_matrix > 0) > 0) # Number of genes expressed in any Oligodendrocyte cell

cells_group1 <- WhichCells(labeled_obj, expression = group == "Low" & celltypes == "Oligodendrocyte")
cells_group2 <- WhichCells(labeled_obj, expression = group == "Acute" & celltypes == "Oligodendrocyte")

expr_matrix_group1 <- GetAssayData(labeled_obj, slot = "data")[, cells_group1]
expr_matrix_group2 <- GetAssayData(labeled_obj, slot = "data")[, cells_group2]

pct_expressing_group1 <- rowMeans(expr_matrix_group1 > 0)
pct_expressing_group2 <- rowMeans(expr_matrix_group2 > 0)

genes_passing <- sum(pct_expressing_group1 > 0.01 | pct_expressing_group2 > 0.01)
genes_passing


# Assuming your integrated Seurat object is called 'obj_assigned'
# and your classification is stored in a column named 'customclassif' or similar

# Choose/rename columns to match hers
meta_df <- obj_assigned@meta.data[, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                                      "percent.mt","expt", "dose",
                                      "integrated_snn_res.0.03", "seurat_clusters",
                                      "customclassif")]

# Drop the "group" column from metadata before export
meta_df <- obj_assigned@meta.data[, !(colnames(obj_assigned@meta.data) %in% "group")]

# Write to CSV with cell barcodes as the first column (row names)
write.csv(meta_df, "integrated_opioid.csv", row.names = TRUE)

