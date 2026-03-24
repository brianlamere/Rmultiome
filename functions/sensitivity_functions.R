transfer_labels <- function(loo_obj,
                           celltype_markers,
                           celltypes_of_interest,
                           output_dir = NULL,
                           sample_name = NULL) {

  # Step 1: Find markers for each LOO cluster
  cat("    Finding cluster markers...\n")
  loo_markers <- FindAllMarkers(loo_obj, only.pos = TRUE,
                                min.pct = 0.10, logfc.threshold = 0.25,
                                verbose = FALSE)

  # Step 2: Assign each cluster to best-matching celltype
  # FIXED: Use list collection instead of rbind loop
  assignment_list <- list()
  
  cluster_ids <- unique(loo_obj$seurat_clusters)
  
  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    
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

    # Store in list (avoids row name conflicts)
    assignment_list[[i]] <- data.frame(
      cluster = cluster_id,
      assigned_celltype = best_match,
      jaccard_score = best_score,
      n_cells = sum(loo_obj$seurat_clusters == cluster_id),
      stringsAsFactors = FALSE
    )
  }

  # Bind all at once and set clean row names
  cluster_assignments <- do.call(rbind, assignment_list)
  rownames(cluster_assignments) <- as.character(1:nrow(cluster_assignments))

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

# instead of being in the test, putting here.  Just compiles lists from the previous output
DEDA_results_compose <- function(celltypes_list) {
  full_de_results <- list()
  full_da_results <- list()
  for (celltype in celltypes_list) {
    for (comparison in comparisons_list) {
      comp_name <- sprintf("%s_%s_vs_%s", celltype, comparison[1], comparison[2])

      # Load DE
      de_file <- file.path(project_export,
                        sprintf("DiffExpress_results_pseudobulk_%s_%s_vs_%s.csv",
                               celltype, comparison[1], comparison[2]))
      if (file.exists(de_file)) {
        full_de_results[[comp_name]] <- read.csv(de_file)
        cat(sprintf("  ✓ DE: %s (%d genes)\n",
	comp_name, nrow(full_de_results[[comp_name]])))
      } else {
        warning(sprintf("  ✗ Missing DE: %s\n", basename(de_file)))
      }

      # Load DA
      da_file <- file.path(project_export,
                        sprintf("DiffAccess_results_pseudobulk_%s_%s_vs_%s.csv",
                               celltype, comparison[1], comparison[2]))
      if (file.exists(da_file)) {
        full_da_results[[comp_name]] <- read.csv(da_file)
        cat(sprintf("  ✓ DA: %s (%d peaks)\n", comp_name,
	    nrow(full_da_results[[comp_name]])))
      } else {
        warning(sprintf("  ✗ Missing DA: %s\n", basename(da_file)))
      }
    }
  }
  return(list(da_results = full_da_results, de_results = full_de_results))
}

#pulling this directly from run_pipeline2, though as soon as I do that it is out of sync.
# given these all use the settings files after the objects, I should make the parameters 
# all be default values, so I can call the functions with just the objects
harmony_FMMN_cluster_task <- function(preharmony_obj) {
	# doing the reassignment so the lines are exactly the same as run_pipeline2.R
	merged_obj <- preharmony_obj
	harmony_obj <- harmonize_both(
	  merged_obj,
	  random_seed = harmony_settings$random_seed,
	  harmony_max_iter = harmony_settings$max_iter,
	  harmony_project.dim = harmony_settings$project_dim,
	  harmony_dims = harmony_settings$dims_use
	)

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
	return(clustered_obj)
}
