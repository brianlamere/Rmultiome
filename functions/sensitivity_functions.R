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
