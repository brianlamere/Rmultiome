source("/projects/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
source(file.path(project_outdir, "project_settings.R"))

#At this point you stop being production, and move to doing a parameter sweep to
#find the right settings to use.  You will do most of the rest of these tasks
#many times if you do it correctly, since you will be using subjective and objective
#assessments of the different settings, including with cell typing, before moving
#back to differential analysis.  Don't return here until you know what to put on
#the next lines for dims, knn, and resolution

#harmony_obj <- readRDS("/projects/opioid/vault/harmonized.rds")

harmony_obj250_k40 <- FMMN_task(harmony_obj, dims_pca = 2:50, dims_harmony = 2:50, knn = 40)
premap_obj <- cluster_data(harmony_obj250_k40, alg = 3, res = 0.1,
                           cluster_dims = 2:50, cluster_seed = 1984)
DimPlot(premap_obj,label=T, raster=FALSE)

#adding file name info that corresponds to the resolution used for FindClusters
saveRDS(premap_obj, "/projects/opioid/vault/pre_mapping_dim240_knn40_res0.04_seed1984.rds")
premap_obj <- readRDS("/projects/opioid/vault96/pre_mapping_dim240_knn40_res0.04_seed1984.rds")

DimPlot(premap_obj,label=T, raster=FALSE)

labeled_obj <- premap_obj
#Note: these assignments are ONLY TRUE if 2:40/40/0.05 is used and only for this data!
cluster_to_celltype <- c(
  "0" = "Oligodendrocyte",
  "1" = "Glutamatergic neurons",
  "2" = "Astrocyte",
  "3" = "Microglia",
  "4" = "GABAergic neurons",
  "5" = "GABAergic neurons",
  "6" = "Oligodendrocyte precursor cells",
  "7" = "Endothelial",
  "8" = "Mature Neurons",
  "9" = "Unknown", # Myelinating Schwann Cells"
  "10" = "Glutamatergic neurons",
  "11" = "Dopaminergic neurons",
  "12" = "" #Endothelial cells per sctype
)

labeled_obj$celltypes <- cluster_to_celltype[as.character(labeled_obj$seurat_clusters)]

DimPlot(labeled_obj, group.by = "celltypes", label = TRUE, raster=FALSE)

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

