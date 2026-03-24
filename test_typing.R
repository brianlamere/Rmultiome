source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

init_project()
pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

# Load settings
cortex_markers <- readRDS(file.path(Rmultiome_path, "Cortex_Consolidated_Markers.rds"))

# Use the ALREADY HARMONIZED object from early pipeline2
merged_obj <- readRDS(file.path(rdsdir, "merged_preharmony.Rds"))
full_obj <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

samples <- unique(merged_obj$orig.ident)
cat(sprintf("  Samples: %d\n", length(samples)))
cat(sprintf("  Total cells: %d\n\n", ncol(merged_obj)))

celltypes_list <- c("Oligodendrocytes", "Microglia_Macrophages", "Astrocytes")
comparisons_list <- list(c("Low", "No_HIV"),c("Acute", "Low"),c("Chronic", "Low"))

# Extract markers from existing settings
celltype_markers <- extract_celltype_markers(celltype_settings, celltypes_list)

validation_dir <- file.path(project_export, "LOO_validation_plots")
dir.create(validation_dir, showWarnings = FALSE)

DEDA_results <- DEDA_results_compose(celltypes_list = celltypes_list)
full_da_results <- DEDA_results$da_results
full_de_results <- DEDA_results$de_results

loo_results <- list()
overall_start <- Sys.time()

for (i in seq_along(samples)) {
  sample_to_exclude <- samples[i]
  #Subset
  loo_obj <- subset(merged_obj, subset = orig.ident != sample_to_exclude)
  #Harmony, FMMN, clustering
  clustered_obj <- harmony_FMMN_cluster_task(merged_obj)
  cat("join layers, assign celltypes\n")
  DefaultAssay(clustered_obj) <- "RNA"
  loo_obj <- JoinLayers(clustered_obj)
  typing_results <- assign_celltype_from_dotplot(loo_obj, cortex_markers$markers_lists)
  cluster_to_celltype <- typing_results$assignments %>%
    filter(!is.na(cluster), !is.na(celltype)) %>%
    select(cluster, celltype, confidence, score)
  loo_obj$celltype <- cluster_to_celltype$celltype[
    match(as.character(loo_obj$seurat_clusters),
        as.character(cluster_to_celltype$cluster))]
  # Subset to cell types of interest
  loo_obj <- subset(loo_obj, subset = celltype %in% celltypes_list)
  # Good to this line?
  # TODO: WIP, new DE/DA (below is just pulled from run_pipeline2.R, needs mods.
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
        output_prefix = file.path(project_export, "DiffExpress_results"),
        min_cells_per_sample = 50
      )
    }
  }

  # Run Differential Accessibility (PSEUDO-BULK)
  cat("\n=== Running Differential Accessibility (Pseudo-bulk) ===\n")
  for (celltype in celltypes_list) {
    for (comparison in comparisons_list) {
      cat(sprintf("DA: %s - %s vs %s\n", celltype, comparison[1], comparison[2]))
      run_pseudobulk_DA_and_export(seurat_obj = obj_assigned, celltype_col = "celltype",
        celltype_value = celltype, sample_col = "orig.ident", group_col = "group",
        ident.1 = comparison[1], ident.2 = comparison[2], min_cells_per_sample = 50,
        output_prefix = file.path(project_export, "DiffAccess_results"))
  }}

