source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
init_project()

source(file.path(Rmultiome_path, "sarah_func.R"))

pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

EnsDbAnnos <- loadannotations()

obj_assigned <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

DefaultAssay(obj_assigned) <- "RNA"
obj_assigned <- JoinLayers(obj_assigned)

pb <- make_pseudobulk_RNA_for_scores(
  seurat_obj = obj_assigned,
  celltype_value = "Microglia_Macrophages",
  sample_col = "orig.ident",
  group_col = "group",
  min_cells_per_sample = 50
)


#pdf(file.path(project_export, "UMAP_celltype_labeled_font5.pdf"), width = 12, height = 10)
#saveRDS(obj_assigned, file.path(rdsdir, "annotated_obj_final.rds"))
saveRDS(pb, "./newDE/microglia_pseudobulk_for_viral_load_scores.rds")
write.csv(pb$vst_mat, "./newDE/microglia_pseudobulk_vst_matrix.csv")
write.csv(pb$norm_counts, "./newDE/microglia_pseudobulk_normalized_counts.csv")
write.csv(pb$sample_meta, "./newDE/microglia_pseudobulk_sample_meta.csv", row.names = FALSE)
