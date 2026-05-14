source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
init_project()

pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

EnsDbAnnos <- loadannotations()

obj_assigned <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

make_pseudobulk_RNA_for_scores <- function(
    seurat_obj,
    celltype_col = "celltype",
    celltype_value = "Microglia_Macrophages",
    sample_col = "orig.ident",
    group_col = "group",
    min_cells_per_sample = 50
) {
  cell_subset <- subset(
    seurat_obj,
    subset = !!as.name(celltype_col) == celltype_value
  )
 
  DefaultAssay(cell_subset) <- "RNA"
 
  counts <- tryCatch(
    GetAssayData(cell_subset, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(cell_subset, assay = "RNA", slot = "counts")
  )
 
  sample_meta <- cell_subset@meta.data %>%
    group_by(!!sym(sample_col)) %>%
    summarize(
      group = dplyr::first(!!sym(group_col)),
      n_cells = n(),
      .groups = "drop"
    ) %>%
    filter(n_cells >= min_cells_per_sample)
 
  pseudobulk_counts <- matrix(
    0,
    nrow = nrow(counts),
    ncol = nrow(sample_meta),
    dimnames = list(rownames(counts), sample_meta[[sample_col]])
  )
 
  for (sample in sample_meta[[sample_col]]) {
    sample_cells <- colnames(cell_subset)[
      cell_subset@meta.data[[sample_col]] == sample
    ]
    pseudobulk_counts[, sample] <- Matrix::rowSums(
      counts[, sample_cells, drop = FALSE]
    )
  }
 
  # Keep expressed genes
  keep <- rowSums(pseudobulk_counts > 0) > 0
  pseudobulk_counts <- pseudobulk_counts[keep, ]
 
  coldata <- data.frame(
    sample = sample_meta[[sample_col]],
    group = factor(sample_meta$group),
    n_cells = sample_meta$n_cells,
    row.names = sample_meta[[sample_col]]
  )
 
  dds <- DESeqDataSetFromMatrix(
    countData = round(pseudobulk_counts),
    colData = coldata,
    design = ~ group
  )
 
  dds <- estimateSizeFactors(dds)
 
  norm_counts <- counts(dds, normalized = TRUE)
  vst_mat <- assay(vst(dds, blind = TRUE))
 
  list(
    pseudobulk_counts = pseudobulk_counts,
    dds = dds,
    norm_counts = norm_counts,
    vst_mat = vst_mat,
    sample_meta = sample_meta
  )
}

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
