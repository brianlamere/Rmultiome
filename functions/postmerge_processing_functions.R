# Contains rna and ATAC modality functions, and harmony/batch effects reductions

#' @param seurat_obj Seurat object.
#' @return Seurat object with FindVariableFeatures, ScaleData, and RunPCA
post_merge_rna <- function(pm_rna_obj) {
  DefaultAssay(pm_rna_obj) <- "RNA"
  pm_rna_obj <- FindVariableFeatures(pm_rna_obj, assay = "RNA")
  pm_rna_obj <- ScaleData(pm_rna_obj, assay = "RNA",
                         features = VariableFeatures(pm_rna_obj))
  pm_rna_obj <- RunPCA(pm_rna_obj, assay = "RNA",
                      features = VariableFeatures(pm_rna_obj))
  return(pm_rna_obj)
}

#' @param seurat_obj Seurat object.
#' @return Chromatin Assay with RunTFIDF, FindTopFeatures, and RunSVD
post_merge_atac <- function(pm_atac_obj) {
  DefaultAssay(pm_atac_obj) <- "ATAC"
  pm_atac_obj <- RunTFIDF(pm_atac_obj)
  pm_atac_obj <- FindTopFeatures(pm_atac_obj, min.cutoff = 'q0')
  pm_atac_obj <- RunSVD(pm_atac_obj)
}

harmonize_both <- function(harmony_obj, harmony_max_iter = 50,
                         harmony_project.dim = FALSE,
                         harmony_dims = NULL, random_seed = NULL,
			 plot_convergence = FALSE) { 
  DefaultAssay(harmony_obj) <- "RNA"
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }  
  harmony_obj <- RunHarmony(
    harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "pca",
    plot_convergence = plot_convergence,
    max_iter = harmony_max_iter,
    reduction.save = reduction.save.RNA,
    project.dim = harmony_project.dim,
    dims.use = harmony_dims
  )

  DefaultAssay(harmony_obj) <- "ATAC"
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }  
  harmony_obj <- RunHarmony(
    object = harmony_obj,
    group.by.vars = "orig.ident",
    reduction.use = "lsi",
    project.dim = harmony_project.dim,
    max_iter = harmony_max_iter,
    reduction.save = reduction.save.ATAC,
    dims.use = harmony_dims
  )
  return(harmony_obj)
}

FMMN_task <- function(FMMN_obj, knn, dims) {
  FMMN_obj <- FindMultiModalNeighbors(
    object = FMMN_obj,
    reduction.list = list(reduction.save.RNA, reduction.save.ATAC),
    dims.list = list(dims, dims),  # Use all harmony dims
    k.nn = knn,
    knn.graph.name = "wknn",
    snn.graph.name = "wsnn",
    weighted.nn.name = "weighted.nn"
  )
  return(FMMN_obj)
}

cluster_data <- function(harmony_obj, alg, res, run_umap = FALSE, cluster_seed,
                         singleton_handling = c("discard", "merge", "keep")) {
  singleton_handling <- match.arg(singleton_handling)

  # The DefaultAssay is being set for consistent behavior, not because we're doing
  # assay-specific actions; I don't want a minor unintended change to occur just
  # because the assay was ATAC or something else, instead of RNA.
  DefaultAssay(harmony_obj) <- "RNA"

  # Determine group.singletons argument for FindClusters
  group_singletons <- (singleton_handling == "merge")

  # Clustering step
  harmony_obj <- FindClusters(
    harmony_obj,
    graph.name = "wsnn",
    algorithm = alg,
    resolution = res,
    group.singletons = group_singletons,
    # I do not like needing the below, and will work on the data until it isn't
    # needed.  Stable data doesn't change with new random seeds.  Hardcoding the
    # seed is cheating.  All but one cluster/cell type are very very stable, so
    # setting this allows parts of the project to move forward without it
    random.seed = cluster_seed
  )

  # Run UMAP BEFORE discarding singletons
  # This ensures weighted.nn neighbor graph is still present
  if (run_umap) {
    harmony_obj <- RunUMAP(
      harmony_obj,
      nn.name = "weighted.nn",            # This matches weighted.nn.name from FMMN_task
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_"
    )
  }

  # NOW discard singletons (after UMAP is calculated)
  # The UMAP coordinates for retained cells will be preserved during subsetting
  if (singleton_handling == "discard") {
    if ("singleton" %in% levels(harmony_obj$seurat_clusters)) {
      singleton_cells <- WhichCells(harmony_obj, idents = "singleton")
      harmony_obj <- subset(harmony_obj, cells = setdiff(colnames(harmony_obj), singleton_cells))
      # Drop unused cluster level
      harmony_obj$seurat_clusters <- droplevels(harmony_obj$seurat_clusters)

      cat(sprintf("Discarded %d singleton cells\n", length(singleton_cells)))
    }
  }

  return(harmony_obj)
}

target_markers <- function(harmony_obj, numMarks = 5) {
  DefaultAssay(harmony_obj) <- "RNA"
  cluster_counts <- table(harmony_obj$seurat_clusters)
  non_singleton_clusters <- names(cluster_counts[cluster_counts > 1])
  
  # Subset object to exclude singleton clusters
  seurat_nonsingleton <- subset(harmony_obj, idents = non_singleton_clusters)
  
  # Find markers for the non-singleton clusters only
  markers_nonsingleton <- FindAllMarkers(seurat_nonsingleton, assay = "RNA")
  
  # For top 5 markers per cluster (adjust 'n' as needed)
  top_markers_nonsingleton <- markers_nonsingleton %>% group_by(cluster) %>% 
    top_n(n = 5, wt = avg_log2FC)
  #or with few/no singletons...
  #markers <- FindAllMarkers(merged_data, assay = "RNA", only.pos = TRUE)
  #top_markers <- markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  #print(top_markers)
  return(harmony_obj)
}

transfer_labels <- function(loo_obj,
                           celltype_markers,
                           celltypes_of_interest,
                           output_dir = NULL,
                           sample_name = NULL) {

  # Step 1: Find markers for each LOO cluster
  cat("    Finding cluster markers...\n")
  loo_markers <- FindAllMarkers(loo_obj, only.pos = TRUE,
                                min.pct = 0.25, logfc.threshold = 0.25,
                                verbose = FALSE)

  # Step 2: Assign each cluster to best-matching celltype
  cluster_assignments <- data.frame()

  for (cluster_id in unique(loo_obj$seurat_clusters)) {
    # Top markers for this cluster
    cluster_markers <- loo_markers %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = 50) %>%
      pull(gene)

    # Calculate Jaccard similarity with each celltype
    best_match <- NA
    best_score <- 0

    for (celltype in celltypes_of_interest) {
      reference_markers <- celltype_markers[[celltype]]

      overlap <- length(intersect(cluster_markers, reference_markers))
      union_size <- length(union(cluster_markers, reference_markers))
      jaccard <- overlap / union_size

      if (jaccard > best_score) {
        best_score <- jaccard
        best_match <- celltype
      }
    }

    cluster_assignments <- rbind(cluster_assignments, data.frame(
      cluster = cluster_id,
      assigned_celltype = best_match,
      jaccard_score = best_score,
      n_cells = sum(loo_obj$seurat_clusters == cluster_id),
      stringsAsFactors = FALSE
    ))
  }

  # FIX: Reset row names to avoid duplicate row.names error
  rownames(cluster_assignments) <- NULL

  cat("    Cluster assignments:\n")
  print(cluster_assignments[order(cluster_assignments$cluster),
                           c("cluster", "assigned_celltype", "jaccard_score", "n_cells")])

  # Step 3: Assign labels to cells
  loo_obj$celltype <- cluster_assignments$assigned_celltype[
    match(loo_obj$seurat_clusters, cluster_assignments$cluster)
  ]

  # Step 4: Generate validation dotplots
  if (!is.null(output_dir) && !is.null(sample_name)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    for (celltype in celltypes_of_interest) {
      markers_to_plot <- celltype_markers[[celltype]]

      Idents(loo_obj) <- loo_obj$celltype

      p <- DotPlot(loo_obj, features = markers_to_plot) +
        ggtitle(paste0("LOO: Exclude ", sample_name, " - ", celltype, " markers")) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))

      ggsave(filename = file.path(output_dir,
                                  paste0("loo_", sample_name, "_", celltype, "_markers.pdf")),
             plot = p, width = 8, height = 5)
    }

    # Save cluster assignment table
    write.csv(cluster_assignments,
             file.path(output_dir, paste0("loo_", sample_name, "_cluster_assignments.csv")),
             row.names = FALSE)

    cat(sprintf("    Validation plots saved to: %s\n", output_dir))
  }

  # Step 5: Subset to only cell types of interest (smaller object for DE/DA)
  loo_obj <- subset(loo_obj, subset = celltype %in% celltypes_of_interest)

  cat(sprintf("    Kept %d cells across %d cell types\n",
             ncol(loo_obj), length(unique(loo_obj$celltype))))

  return(list(
    obj = loo_obj,
    assignments = cluster_assignments
  ))
}
