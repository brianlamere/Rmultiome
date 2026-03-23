#I want to be able to see oft-repeated flows in a short chunk, only exposing
# those things which might change.  Since by definition these must be done like
# pipeline2, doing them that way.


# === STEP 1: Run Harmony with saved settings (reproducible) ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = harmony_settings$random_seed,
  harmony_max_iter = harmony_settings$max_iter,
  harmony_project.dim = harmony_settings$project_dim,
  harmony_dims = harmony_settings$dims_use
)

# === STEP 2: Apply clustering settings ===
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
