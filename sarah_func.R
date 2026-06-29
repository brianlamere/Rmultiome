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
 
  counts <- GetAssayData(cell_subset, assay = "RNA", layer = "counts")
 
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

