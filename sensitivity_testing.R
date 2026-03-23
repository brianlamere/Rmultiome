source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

init_project()

pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)

# Loading datasets
merged_obj <- readRDS(file.path(rdsdir, "merged_preharmony.Rds"))
full_obj <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

samples <- unique(merged_obj$orig.ident)
cat(sprintf("  Samples: %d\n", length(samples)))
cat(sprintf("  Total cells: %d\n\n", ncol(merged_obj)))

# Celltypes of interest - unique to your task
celltypes_list <- c("Oligodendrocytes", "Microglia_Macrophages", "Astrocytes")

# === Define comparisons (DOSE-RESPONSE CONVENTION) ===
# Positive log2FC = gene/peak INCREASES with higher exposure
# Negative log2FC = gene/peak DECREASES with higher exposure

comparisons_list <- list(
  c("Low", "No_HIV"),
  c("Acute", "Low"),
  c("Chronic", "Low")
)

# Extract markers from existing settings
celltype_markers <- extract_celltype_markers(celltype_settings, celltypes_list)
cat("Cell type markers:\n")
for (ct in names(celltype_markers)) {
  cat(sprintf("  %s: %s\n", ct, paste(celltype_markers[[ct]], collapse = ", ")))
}
cat("\n")

# Validation output directory
validation_dir <- file.path(project_export, "LOO_validation_plots")
dir.create(validation_dir, showWarnings = FALSE)

# Suppress warnings, necessary because they persisted even with warning resolved
options(future.rng.onMisuse = "ignore")

# === Load full dataset results ===
cat("Loading full dataset DE/DA results...\n")
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
      cat(sprintf("  ✓ DE: %s (%d genes)\n", comp_name, nrow(full_de_results[[comp_name]])))
    } else {
      warning(sprintf("  ✗ Missing DE: %s\n", basename(de_file)))
    }

    # Load DA
    da_file <- file.path(project_export,
                        sprintf("DiffAccess_results_pseudobulk_%s_%s_vs_%s.csv",
                               celltype, comparison[1], comparison[2]))
    if (file.exists(da_file)) {
      full_da_results[[comp_name]] <- read.csv(da_file)
      cat(sprintf("  ✓ DA: %s (%d peaks)\n", comp_name, nrow(full_da_results[[comp_name]])))
    } else {
      warning(sprintf("  ✗ Missing DA: %s\n", basename(da_file)))
    }
  }
}

cat(sprintf("\nLoaded %d DE and %d DA comparisons\n\n",
           length(full_de_results), length(full_da_results)))

# === Storage for results ===
loo_results <- list()
overall_start <- Sys.time()

# ============================================================================
# === LOO LOOP (Sequential iterations, parallelized DA within each) ===
# ============================================================================

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("Starting LOO iterations...\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

for (i in seq_along(samples)) {
  sample_to_exclude <- samples[i]
  iter_start <- Sys.time()

  cat(sprintf("\n[%d/%d] ========== Excluding: %s ==========\n",
             i, length(samples), sample_to_exclude))
  cat(sprintf("Time: %s\n", format(Sys.time(), "%X")))

  # === STEP 1: Subset ===
  cat("  Creating LOO subset...\n")
  loo_obj <- subset(merged_obj, subset = orig.ident != sample_to_exclude)
  cat(sprintf("    %d → %d cells\n", ncol(merged_obj), ncol(loo_obj)))

  # === STEP 2: Harmony ===
  cat("  Running Harmony integration...\n")
  loo_obj <- harmonize_both(loo_obj, random_seed = harmony_settings$random_seed,
    harmony_max_iter = harmony_settings$max_iter, harmony_dims = harmony_settings$dims_use,
    harmony_project.dim = harmony_settings$project_dim, plot_convergence = FALSE)

  # === STEP 3: FMMN ===
  cat("  Running FMMN + Clustering...\n")
  dims <- cluster_settings$dims_min:cluster_settings$dims_max
  loo_obj <- FMMN_task(loo_obj, knn = cluster_settings$knn, dims = dims)
  
  # === STEP 4: Clustering ===
  loo_obj <- cluster_data(loo_obj, alg = cluster_settings$algorithm,
    res = cluster_settings$resolution, cluster_seed = cluster_settings$random_seed,
    singleton_handling = "discard", run_umap = TRUE)

  # === STEP 5: Join layers for DE/DA ===
  cat("  Joining RNA layers...\n")
  DefaultAssay(loo_obj) <- "RNA"
  loo_obj <- JoinLayers(loo_obj)

  # STEP 6: Transfer labels
  transfer_result <- transfer_labels(loo_obj = loo_obj, celltype_markers = celltype_markers,
    celltypes_of_interest = celltypes_list, output_dir = validation_dir,
    sample_name = sample_to_exclude)
  loo_obj <- transfer_result$obj  # Subsetted to 3 cell types

  # Report cell types
  cat("    Cell types transferred:\n")
  print(table(loo_obj$celltype))

  # === STEP 7: Differential Expression (sequential, ~5-10 min) ===
  cat("  Running pseudo-bulk DE (sequential)...\n")
  de_start <- Sys.time()
  de_results <- list()
  for (celltype in celltypes_list) {
    for (comparison in comparisons_list) {
      comp_name <- sprintf("%s_%s_vs_%s", celltype, comparison[1], comparison[2])
      tryCatch({
        de <- run_pseudobulk_DE(seurat_obj = loo_obj, celltype_col = "celltype",
          celltype_value = celltype, sample_col = "orig.ident", group_col = "group",
          ident.1 = comparison[1], ident.2 = comparison[2], min_cells_per_sample = 50)
        if (nrow(de) > 0) {
          de_results[[comp_name]] <- de
          cat(sprintf("    ✓ %s: %d genes, %d sig\n",
                     comp_name, nrow(de), sum(de$padj < 0.05, na.rm = TRUE)))
        }
      }, error = function(e) {
        cat(sprintf("    ✗ %s: ERROR - %s\n", comp_name, e$message))
      })
    }
  }

  de_time <- difftime(Sys.time(), de_start, units = "mins")
  cat(sprintf("  DE complete: %.1f min (%d comparisons)\n",
             de_time, length(de_results)))

  # === STEP 8: Differential Accessibility (PARALLELIZED, ~80-120 min) ===
  cat("  Running pseudo-bulk DA (PARALLELIZED across 9 comparisons)...\n")
  da_start <- Sys.time()

  # Create job list (9 jobs: 3 cell types × 3 comparisons)
  da_jobs <- expand.grid(celltype = celltypes_list,
    comparison_idx = 1:length(comparisons_list), stringsAsFactors = FALSE)

  n_workers <- min(9, parallel::detectCores() - 1)
  plan(multicore, workers = n_workers)

  # Run DA in parallel
  da_results_list <- future_lapply(1:nrow(da_jobs), function(j) {
    job <- da_jobs[j, ]
    comparison <- comparisons_list[[job$comparison_idx]]
    comp_name <- sprintf("%s_%s_vs_%s", job$celltype, comparison[1], comparison[2])
    tryCatch({
      da <- run_pseudobulk_DA(seurat_obj = loo_obj, celltype_col = "celltype",
        celltype_value = job$celltype, sample_col = "orig.ident", group_col = "group",
        ident.1 = comparison[1], ident.2 = comparison[2], min_cells_per_sample = 50)
      if (nrow(da) > 0) {
        return(list(name = comp_name, da = da,
                   n_sig = sum(da$padj < 0.05, na.rm = TRUE)))
      } else {
        return(NULL)
      }
    }, error = function(e) {
      return(list(name = comp_name, error = e$message))
    })
  }, future.seed = TRUE)
  plan(sequential) # Reset plan
  # Convert to named list and report results
  da_results <- list()
  for (result in da_results_list) {
    if (!is.null(result) && !is.null(result$da)) {
      da_results[[result$name]] <- result$da
    }
  }
  da_time <- difftime(Sys.time(), da_start, units = "mins")
  cat(sprintf("  DA complete: %.1f min (%d comparisons)\n",
             da_time, length(da_results)))
  # === STEP 9: Store results (lightweight) ===
  loo_results[[sample_to_exclude]] <- list(sample = sample_to_exclude,
    n_cells = ncol(loo_obj), celltype_props = loo_props, prop_correlation = prop_cor,
    de_results = de_results, da_results = da_results,
    timing = list( de_minutes = as.numeric(de_time), da_minutes = as.numeric(da_time)))
  # Save individual LOO result
  saveRDS(loo_results[[sample_to_exclude]],
         file.path(rdsdir, sprintf("loo_result_%s.rds", sample_to_exclude)))
  iter_time <- difftime(Sys.time(), iter_start, units = "mins")
  cat(sprintf("  ✓ Iteration complete: %.1f min (DE: %.1f, DA: %.1f)\n",
             iter_time, de_time, da_time))
  # Cleanup
  rm(loo_obj, de_results, da_results, da_results_list)
  gc()
  # Estimate remaining time
  if (i < length(samples)) {
    avg_time <- difftime(Sys.time(), overall_start, units = "mins") / i
    remaining <- avg_time * (length(samples) - i)
    eta <- Sys.time() + as.difftime(remaining, units = "mins")
    cat(sprintf("  ETA for completion: %s (%.1f hours remaining)\n",
               format(eta, "%Y-%m-%d %H:%M"), remaining/60))
  }
}

overall_time <- difftime(Sys.time(), overall_start, units = "mins")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("All LOO iterations complete: %.1f minutes (%.1f hours)\n",
           overall_time, overall_time/60))
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ============================================================================
# === ANALYSIS: Compare LOO to full dataset ===
# ============================================================================

cat("=== Analyzing Results ===\n\n")

# === 1. Cell type proportion correlations ===
prop_cors <- sapply(loo_results, function(x) x$prop_correlation)

cat("=== Cell Type Proportions ===\n")
cat(sprintf("Mean correlation: %.4f ± %.4f\n", mean(prop_cors), sd(prop_cors)))
cat(sprintf("Range: %.4f - %.4f\n", min(prop_cors), max(prop_cors)))
cat(sprintf("Lowest: %s (r = %.4f)\n",
           names(which.min(prop_cors)), min(prop_cors)))
cat("\n")

# === 2. DE effect size correlations ===
cat("=== Differential Expression Stability ===\n")
de_correlations <- data.frame()

for (sample in samples) {
  for (comp_name in names(full_de_results)) {
    full_de <- full_de_results[[comp_name]]
    loo_de <- loo_results[[sample]]$de_results[[comp_name]]

    if (is.null(loo_de) || nrow(loo_de) == 0) next

    # Match genes
    common_genes <- intersect(full_de$gene, loo_de$gene)
    if (length(common_genes) < 10) next

    # Effect size correlation
    full_sub <- full_de[match(common_genes, full_de$gene), ]
    loo_sub <- loo_de[match(common_genes, loo_de$gene), ]

    fc_cor <- cor(full_sub$log2FoldChange, loo_sub$log2FoldChange)

    # Significant gene overlap
    full_sig <- full_de$gene[full_de$padj < 0.05 & !is.na(full_de$padj)]
    loo_sig <- loo_de$gene[loo_de$padj < 0.05 & !is.na(loo_de$padj)]
    jaccard <- length(intersect(full_sig, loo_sig)) /
              length(union(full_sig, loo_sig))

    de_correlations <- rbind(de_correlations, data.frame(
      sample_excluded = sample,
      comparison = comp_name,
      fc_correlation = fc_cor,
      jaccard = jaccard,
      n_genes_tested = length(common_genes),
      n_sig_full = length(full_sig),
      n_sig_loo = length(loo_sig)
    ))
  }
}

cat(sprintf("Mean log2FC correlation: %.4f ± %.4f\n",
           mean(de_correlations$fc_correlation, na.rm = TRUE),
           sd(de_correlations$fc_correlation, na.rm = TRUE)))
cat(sprintf("Mean Jaccard (sig genes): %.4f ± %.4f\n",
           mean(de_correlations$jaccard, na.rm = TRUE),
           sd(de_correlations$jaccard, na.rm = TRUE)))

# Check for unstable comparisons
low_cor <- de_correlations[de_correlations$fc_correlation < 0.8, ]
if (nrow(low_cor) > 0) {
  cat("\nComparisons with low correlation (<0.8) when excluding:\n")
  print(low_cor[, c("sample_excluded", "comparison", "fc_correlation", "jaccard")])
}
cat("\n")

# === 3. DA effect size correlations ===
cat("=== Differential Accessibility Stability ===\n")
da_correlations <- data.frame()

for (sample in samples) {
  for (comp_name in names(full_da_results)) {
    full_da <- full_da_results[[comp_name]]
    loo_da <- loo_results[[sample]]$da_results[[comp_name]]

    if (is.null(loo_da) || nrow(loo_da) == 0) next

    # Match peaks
    common_peaks <- intersect(full_da$peak, loo_da$peak)
    if (length(common_peaks) < 10) next

    # Effect size correlation
    full_sub <- full_da[match(common_peaks, full_da$peak), ]
    loo_sub <- loo_da[match(common_peaks, loo_da$peak), ]

    fc_cor <- cor(full_sub$log2FoldChange, loo_sub$log2FoldChange)

    # Significant peak overlap
    full_sig <- full_da$peak[full_da$padj < 0.05 & !is.na(full_da$padj)]
    loo_sig <- loo_da$peak[loo_da$padj < 0.05 & !is.na(loo_da$padj)]
    jaccard <- length(intersect(full_sig, loo_sig)) /
              length(union(full_sig, loo_sig))

    da_correlations <- rbind(da_correlations, data.frame(
      sample_excluded = sample,
      comparison = comp_name,
      fc_correlation = fc_cor,
      jaccard = jaccard,
      n_peaks_tested = length(common_peaks),
      n_sig_full = length(full_sig),
      n_sig_loo = length(loo_sig)
    ))
  }
}

cat(sprintf("Mean log2FC correlation: %.4f ± %.4f\n",
           mean(da_correlations$fc_correlation, na.rm = TRUE),
           sd(da_correlations$fc_correlation, na.rm = TRUE)))
cat(sprintf("Mean Jaccard (sig peaks): %.4f ± %.4f\n",
           mean(da_correlations$jaccard, na.rm = TRUE),
           sd(da_correlations$jaccard, na.rm = TRUE)))

# Check for unstable comparisons
low_cor_da <- da_correlations[da_correlations$fc_correlation < 0.8, ]
if (nrow(low_cor_da) > 0) {
  cat("\nComparisons with low correlation (<0.8) when excluding:\n")
  print(low_cor_da[, c("sample_excluded", "comparison", "fc_correlation", "jaccard")])
}
cat("\n")

# === 4. Sample QC metrics ===
cat("=== Sample Characteristics ===\n")
sample_qc <- data.frame()

for (sample in samples) {
  sample_cells <- full_obj$orig.ident == sample

  sample_qc <- rbind(sample_qc, data.frame(
    sample = sample,
    group = unique(full_obj$group[sample_cells])[1],
    n_cells = sum(sample_cells),
    median_genes = median(full_obj$nFeature_RNA[sample_cells]),
    median_UMI = median(full_obj$nCount_RNA[sample_cells]),
    pct_mt_median = median(full_obj$percent.mt[sample_cells])
  ))
}

print(sample_qc)
cat("\n")

# === 5. Combined summary ===
summary_stats <- data.frame(
  sample = samples,
  group = sample_qc$group,
  n_cells = sample_qc$n_cells,
  prop_correlation = prop_cors,
  mean_de_fc_cor = tapply(de_correlations$fc_correlation,
                         de_correlations$sample_excluded,
                         mean, na.rm = TRUE)[samples],
  mean_de_jaccard = tapply(de_correlations$jaccard,
                          de_correlations$sample_excluded,
                          mean, na.rm = TRUE)[samples],
  mean_da_fc_cor = tapply(da_correlations$fc_correlation,
                         da_correlations$sample_excluded,
                         mean, na.rm = TRUE)[samples],
  mean_da_jaccard = tapply(da_correlations$jaccard,
                          da_correlations$sample_excluded,
                          mean, na.rm = TRUE)[samples],
  median_genes = sample_qc$median_genes,
  median_UMI = sample_qc$median_UMI,
  pct_mt = sample_qc$pct_mt_median
)

# === 6. Save results ===
cat("=== Saving Results ===\n")

write.csv(summary_stats,
         file.path(project_export, "loo_summary.csv"),
         row.names = FALSE)
cat(sprintf("✓ Saved: %s\n", file.path(project_export, "loo_summary.csv")))

write.csv(de_correlations,
         file.path(project_export, "loo_de_correlations.csv"),
         row.names = FALSE)
cat(sprintf("✓ Saved: %s\n", file.path(project_export, "loo_de_correlations.csv")))

write.csv(da_correlations,
         file.path(project_export, "loo_da_correlations.csv"),
         row.names = FALSE)
cat(sprintf("✓ Saved: %s\n", file.path(project_export, "loo_da_correlations.csv")))

saveRDS(loo_results, file.path(rdsdir, "loo_comprehensive_results.rds"))
cat(sprintf("✓ Saved: %s\n", file.path(rdsdir, "loo_comprehensive_results.rds")))

# === 7. Final summary table ===
cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== Summary Table ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

print(summary_stats[, c("sample", "group", "n_cells", "prop_correlation",
                       "mean_de_fc_cor", "mean_da_fc_cor")])

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("=== Sensitivity Analysis Complete ===\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Total time: %.1f hours\n", overall_time/60))
cat(sprintf("Average time per LOO: %.1f minutes\n", overall_time/length(samples)))

# Average timing breakdown
avg_de_time <- mean(sapply(loo_results, function(x) x$timing$de_minutes))
avg_da_time <- mean(sapply(loo_results, function(x) x$timing$da_minutes))
cat(sprintf("  Average DE time: %.1f min\n", avg_de_time))
cat(sprintf("  Average DA time: %.1f min (parallelized)\n", avg_da_time))

cat("\nKey findings:\n")
cat(sprintf("  - Cell type proportions: r = %.3f ± %.3f (stable if > 0.9)\n",
           mean(prop_cors), sd(prop_cors)))
cat(sprintf("  - DE effect sizes: r = %.3f ± %.3f (robust if > 0.8)\n",
           mean(de_correlations$fc_correlation, na.rm = TRUE),
           sd(de_correlations$fc_correlation, na.rm = TRUE)))
cat(sprintf("  - DA effect sizes: r = %.3f ± %.3f (robust if > 0.8)\n",
           mean(da_correlations$fc_correlation, na.rm = TRUE),
           sd(da_correlations$fc_correlation, na.rm = TRUE)))

cat("\n")

