source("/projects/opioid/Rmultiome/Rmultiome-main.R")
source("/projects/opioid/Rmultiome/settings.R")

#The intent of this script is to be able to be run after you've selected all
#settings, done all QC, and this script will re-create your DE/DA/etc results
#from the raw output of cellranger-atac.  Step1 is not this file!  Step1 is to
#use run_base_QC.R and establish your pipeline1 settings, as stored in settings.R
#If you use run_all.R multiple times with settings tailored for your data, and
#the end results change, then you either have a data stability problem you haven't
#resolved yet, your 1D/KDE trimming is not correct, or you need to select different
#dims/knn/resolution during FMMN and clustering.

# If your biological goal is data that is "good enough" to inform another process,
# and the observed changes between runs fall within a threshold acceptable for your
# research context, then some small degree of variability may be tolerable.
# However, if your intent is a fully deterministic pipeline, any change in results
# between runs is a strong indicator of unresolved issues in data stability,
# preprocessing, or parameter selection.

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))

missing <- setdiff(samplelist, trimming_settings$sample)
if (length(missing) > 0) {
  stop(
    "These samples are in samplelist but missing in trimming_settings: ",
    paste(missing, collapse = ", ")
  )
}

EnsDbAnnos <- loadannotations()

for (sample in samplelist) {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  
  # STEP 1: RAW (import)
  base_path <- get_rds_path(sample, "base")
  if (!file.exists(base_path)) {
    base_obj <- base_object(sample)
    print("Adding chromosome mapping information to ATAC assay.")
    base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
    base_obj <- update_provenance(base_obj, "raw_import")
    saveRDS(base_obj, base_path)
  } else {
    print("Reading previous base file from vault\n")
    base_obj <- readRDS(base_path)
  }
  
  # STEP 2: 1D TRIM
  trimmed_path <- get_rds_path(sample, "trimmed")
  if (!file.exists(trimmed_path)) {
    trim_obj <- base_obj
    print("Calculating NucleosomeSignal and TSSEnrichment for ATAC data.")
    DefaultAssay(trim_obj) <- "ATAC"
    trim_obj <- NucleosomeSignal(trim_obj)
    trim_obj <- TSSEnrichment(trim_obj)
    print("Trimming based on QC")
    trim_obj <- trimSample(trim_obj)
    trim_obj <- update_provenance(trim_obj, "trim_data")
    saveRDS(trim_obj, trimmed_path)
  } else {
    print("Reading previous trim file from vault\n")
    trim_obj <- readRDS(trimmed_path)
  }
  
  #For filtering, "union" keeps cells if either atac or rna passes, where
  # "intersection" keeps only those cells that pass both.
  # STEP 3: 2D TRIM
  kde_path <- get_rds_path(sample, "kdetrim")
  if (!file.exists(kde_path)) {
    kde_obj <- trim_obj
    print("Doing n-Dimensional KDE trimming.")
    cat("Cells before trimming:", nrow(kde_obj@meta.data), "\n")
    #percent filters are percent being kept from each
    #union means all that are kept in at least one of the checks
    #intersection means only those that are kept in both checks
    kde_obj <- kdeTrimSample(kde_obj,
                             atac_percentile = 0.95,
                             rna_percentile = 0.95,
                             combine_method = "intersection"
                             )
    cat("Cells after KDE trimming:", nrow(kde_obj@meta.data), "\n")
    saveRDS(kde_obj, kde_path)
  } else {
    print("Reading previous kdetrim file from vault\n")
    kde_obj <- readRDS(kde_path)
  }
    
  # STEP 4: preRNA
  preRNA_path <- get_rds_path(sample, "preRNA")
  if (!file.exists(preRNA_path)) {
    obj <- kde_obj
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2500)
    obj <- update_provenance(obj, "pre-merge_rna")
    saveRDS(obj, preRNA_path)
  } else {
    print("Reading previous preRNA file from vault\n")
    obj <- readRDS(preRNA_path)
  }
  
  # STEP 5: preATAC/pipeline1
  pipeline1_path <- get_rds_path(sample, "pipeline1")
  if (!file.exists(pipeline1_path)) {
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj)
    obj <- update_provenance(obj, "pre-merge_atac")
    saveRDS(obj, pipeline1_path)
  } else {
    print("Files already present for this sample for pipeline1\n")
  }
  #rm(x_atac, y_atac, x_rna, y_rna, dens_atac, dens_rna, pass_atac, pass_rna)
  gc() #R is obnoxious
  cat(sprintf("Sample %s completed successfully.\n", sample))
}

merged_data <- merge_sample_objects(samplelist)

saveRDS(merged_data, "/projects/opioid/vault/merged_samples95.rds")
#merged_data <- readRDS("/projects/opioid/vault/merged_samples95.rds")

#post-merge RNA modality
merged_data <- post_merge_rna(merged_data)

#post-merge ATAC modality
merged_data <- post_merge_atac(merged_data)

#leaving this here in case you want to restart pre-harmony
#merged_data <- readRDS("/projects/opioid/vault/postATAC.rds")
DefaultAssay(merged_data) <- "RNA"

#time for harmony and FindMultiModalNeighbors
harmony_obj <- harmonize_both(merged_data, harmony_max_iter = 50)
saveRDS(harmony_obj, "/projects/opioid/vault/harmonized.rds")

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

