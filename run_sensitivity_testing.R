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

# Original as labeled in full data set, where more involved typing found microglia had 
# macrophage markers as well
celltypes_list_original <- c("Oligodendrocytes", "Microglia_Macrophages", "Astrocytes")
# Revised, with what the quicker sensitivity assignment will find
celltypes_list <- c("Oligodendrocytes", "Microglia", "Astrocytes")
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
  sample_group <- unique(merged_obj$group[merged_obj$orig.ident == sample_to_exclude])
  cat(paste(rep("=", 40), collapse = ""), "\n")
  cat(sprintf("==Excluding: %s at Time: %s\n", sample_to_exclude, format(Sys.time(), "%X")))
  cat(" INFO: Subsetting, running harmony, FMMN, and clustering\n")
  loo_obj <- subset(merged_obj, subset = orig.ident != sample_to_exclude)
  clustered_obj <- harmony_FMMN_cluster_task(loo_obj)
  cat(" INFO: Joining layers, assigning celltypes\n")
  DefaultAssay(clustered_obj) <- "RNA"
  #loo_obj <- JoinLayers(clustered_obj)
  typing_results <- assign_celltype_from_dotplot(seurat_obj = clustered_obj,
    cortex_markers$marker_lists)
  loo_obj <- JoinLayers(clustered_obj)
  cluster_to_celltype <- typing_results$assignments %>%
    filter(!is.na(cluster), !is.na(celltype)) %>%
    dplyr::select(cluster, celltype, score)
  loo_obj$celltype <- cluster_to_celltype$celltype[
    match(as.character(loo_obj$seurat_clusters),
        as.character(cluster_to_celltype$cluster))]
  cat(" INFO: Subsetting to cell types of interest.\n")
  loo_obj <- subset(loo_obj, subset = celltype %in% celltypes_list)
  # Good to this line?
  # TODO: WIP, new DE/DA (below is just pulled from run_pipeline2.R, needs mods.

  loo_de <- run_LOO_DEDA(loo_obj, celltypes_list, comparisons_list, assay = "DE")
  loo_da <- run_LOO_DEDA(loo_obj, celltypes_list, comparisons_list, assay = "DA")
  
  de_comparison <- compare_LOO_to_full(loo_de, full_de_results, 
                                        sample_to_exclude, sample_group)
  da_comparison <- compare_LOO_to_full(loo_da, full_da_results,
                                        sample_to_exclude, sample_group)
  
  loo_results[[sample_to_exclude]] <- list(de = de_comparison, da = da_comparison)
}
summary_de <- do.call(rbind, lapply(loo_results, `[[`, "de"))
summary_da <- do.call(rbind, lapply(loo_results, `[[`, "da"))
