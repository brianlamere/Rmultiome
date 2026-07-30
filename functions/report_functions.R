
#' Generate per-sample pipeline1 report
#'
#' Creates a text report summarizing QC steps and settings
#'
#' @param sample Sample name
#' @param base_obj Base Seurat object (after import)
#' @param trim_obj After 1D trimming
#' @param kde_obj After KDE trimming
#' @param final_obj Final object (after doublet removal if applicable)
#' @param pipeline1_path Path where final object was saved
#' @param pipeline1_settings Dataframe with all settings
#' @param use_scdblfinder Whether doublet removal was performed
#' @param doublet_stats Optional list with doublet statistics (if performed)
#' @param report_dir Directory to save report (default: project_export)
generate_sample_report <- function(
    sample,
    base_obj,
    trim_obj,
    kde_obj,
    final_obj,
    pipeline1_path,
    pipeline1_settings,
    use_scdblfinder = FALSE,
    doublet_stats = NULL,
    report_dir = project_export
) {
  
  report_path <- file.path(report_dir, paste0(sample, "_pipeline1_report.txt"))
  
  # Get settings for this sample
  params <- pipeline1_settings[pipeline1_settings$sample == sample, ]
  
  # Cell counts at each stage
  n_base <- nrow(base_obj@meta.data)
  n_trim <- nrow(trim_obj@meta.data)
  n_kde <- nrow(kde_obj@meta.data)
  n_final <- nrow(final_obj@meta.data)
  
  # Feature counts
  rna_base <- nrow(base_obj[["RNA"]])
  atac_base <- nrow(base_obj[["ATAC"]])
  rna_trim <- nrow(trim_obj[["RNA"]])
  atac_trim <- nrow(trim_obj[["ATAC"]])
  rna_kde <- nrow(kde_obj[["RNA"]])
  atac_kde <- nrow(kde_obj[["ATAC"]])
  rna_final <- nrow(final_obj[["RNA"]])
  atac_final <- nrow(final_obj[["ATAC"]])
  
  # Build trimming settings text
  trim_fields <- c("min_nCount_ATAC", "max_nCount_ATAC", "min_nCount_RNA", "max_nCount_RNA",
                   "min_nss", "max_nss", "max_percentMT", "min_TSS", "max_TSS")
  trim_lines <- sapply(trim_fields, function(f) {
    sprintf("  %s = %s", f, as.character(params[[f]]))
  })
  
  # Build KDE settings text
  kde_fields <- c("atac_percentile", "rna_percentile", "combine_method")
  kde_lines <- sapply(kde_fields, function(f) {
    sprintf("  %s = %s", f, as.character(params[[f]]))
  })
  
  # Start writing report
  report_lines <- c(
    "=====================================",
    sprintf("Pipeline1 Sample Report for: %s", sample),
    sprintf("Date: %s", Sys.time()),
    "",
    sprintf("Step 1: Cells at base (import): %d", n_base),
    "Step 1: Features at base (import):",
    sprintf("  RNA features (base): %d", rna_base),
    sprintf("  ATAC features (base): %d", atac_base),
    "",
    sprintf("Step 2: Cells after 1D trim: %d", n_trim),
    sprintf("  Removed: %d cells (%.1f%%)", n_base - n_trim, 100*(n_base - n_trim)/n_base),
    "Step 2: Features after 1D trim:",
    sprintf("  RNA features (after 1D trim): %d", rna_trim),
    sprintf("  ATAC features (after 1D trim): %d", atac_trim),
    "",
    sprintf("Step 3: Cells after KDE trim: %d", n_kde),
    sprintf("  Removed: %d cells (%.1f%%)", n_trim - n_kde, 100*(n_trim - n_kde)/n_trim),
    "Step 3: Features remaining after KDE trim:",
    sprintf("  RNA features (after KDE trim): %d", rna_kde),
    sprintf("  ATAC features (after KDE trim): %d", atac_kde),
    ""
  )
  
  # Add doublet detection section if performed
  if (use_scdblfinder && !is.null(doublet_stats)) {
    doublet_section <- c(
      sprintf("Step 3.5: Doublet Detection (scDblFinder):"),
      sprintf("  Cells before doublet removal: %d", doublet_stats$n_cells_before),
      sprintf("  Expected doublets: %.1f (%.2f%%)", 
              doublet_stats$expected_dbr, 
              100 * doublet_stats$expected_dbr / doublet_stats$n_cells_before),
      sprintf("  Detected doublets: %d (%.2f%%)", 
              doublet_stats$n_doublets, 
              doublet_stats$pct_doublets),
      sprintf("  Detected singlets: %d (%.2f%%)", 
              doublet_stats$n_singlets,
              100 * doublet_stats$n_singlets / doublet_stats$n_cells_before),
      sprintf("  Difference from expected: %+.1f doublets (%+.2f%%)",
              doublet_stats$n_doublets - doublet_stats$expected_dbr,
              doublet_stats$pct_doublets - 100 * doublet_stats$expected_dbr / doublet_stats$n_cells_before),
      "",
      "  Score statistics:",
      sprintf("    Threshold: %.4f", doublet_stats$threshold),
      sprintf("    Singlet scores (median): %.4f", doublet_stats$singlet_score_median),
      sprintf("    Doublet scores (median): %.4f", doublet_stats$doublet_score_median),
      sprintf("    Separation: %.4f", doublet_stats$doublet_score_median - doublet_stats$singlet_score_median),
      "",
      sprintf("Step 4: Cells after doublet removal: %d", n_final),
      sprintf("  Removed: %d doublet cells", doublet_stats$n_doublets),
      ""
    )
    report_lines <- c(report_lines, doublet_section)
  } else if (use_scdblfinder && is.null(doublet_stats)) {
    # Flag if doublets should have been run but stats missing
    report_lines <- c(report_lines, 
                     "Step 3.5: Doublet detection enabled but no statistics available",
                     "")
  }
  
  # Add settings used
  report_lines <- c(report_lines,
    "The following initial 1D trimming settings were used:",
    trim_lines,
    "",
    "KDE trim settings:",
    kde_lines,
    ""
  )
  
  # Add final summary
  report_lines <- c(report_lines,
    "=== Final Summary ===",
    sprintf("Total cells retained: %d / %d (%.1f%%)", 
            n_final, n_base, 100*n_final/n_base),
    sprintf("Total cells removed: %d (%.1f%%)", 
            n_base - n_final, 100*(n_base - n_final)/n_base),
    "",
    sprintf("Final object saved at: %s", pipeline1_path),
    "====================================="
  )
  
  # Write to file
  writeLines(report_lines, con = report_path)
  
  cat(sprintf("Sample report saved: %s\n", report_path))
}

#' Compare CellBender merge reports across samples
#'
#' @param mode "qc" or "pipeline" (determines report directory)
#' @param samplelist Sample names to compare (if NULL, reads from 
#'   pipeline1_settings_file)
compare_cellbender_reports <- function(
  mode = c("qc", "pipeline"), 
  samplelist = NULL
) {
  mode <- match.arg(mode)
  
  # If samplelist not provided, read from trimming settings
  if (is.null(samplelist)) {
    pipeline1_settings <- read_pipeline1_settings(pipeline1_settings_file)
    samplelist <- pipeline1_settings$sample
    cat(sprintf(
      "Using samples from pipeline1_settings_file: %s\n", 
      pipeline1_settings_file
    ))
  }
  
  # Determine report directory based on mode
  if (mode == "qc") {
    report_dir <- file.path(tmpfiledir, "cellbender_merge_reports")
  } else {
    report_dir <- cellbender_report_dir
  }
  
  # Read per-sample reports
  reports_list <- list()
  missing_samples <- c()
  
  for (sample in samplelist) {
    report_path <- file.path(report_dir, paste0(sample, ".csv"))
    
    if (file.exists(report_path)) {
      reports_list[[sample]] <- read.csv(
        report_path, 
        stringsAsFactors = FALSE
      )
    } else {
      missing_samples <- c(missing_samples, sample)
    }
  }
  
  # Handle missing reports
  if (length(missing_samples) > 0) {
    warning(sprintf(
      "Missing CellBender reports for samples: %s\nRun base_object() 
with cb_report enabled for these samples.",
      paste(missing_samples, collapse = ", ")
    ))
  }
  
  if (length(reports_list) == 0) {
    stop("No CellBender reports found. Check report directory and 
sample names.")
  }
  
  # Combine all reports
  combined <- do.call(rbind, reports_list)
  
  # Sort by key metrics (prop_atac_covered asc, pearson asc, 
  # weighted_removed desc)
  combined <- combined[order(
    combined$prop_atac_covered, 
    combined$pearson, 
    -combined$weighted_removed
  ), ]
  
  # Reset row names
  rownames(combined) <- NULL
  
  # Display
  cat(sprintf("\n=== CellBender Merge Report Comparison (%s mode) ===\n", mode ))
  cat(sprintf("Samples: %d\n", nrow(combined)))
  cat(sprintf("Report directory: %s\n\n", report_dir))
  
  # Write aggregated report to tmp for easy viewing
  agg_path <- file.path(report_dir, "aggregated_comparison.csv")
  write.csv(combined, agg_path, row.names = FALSE)

  print(combined, row.names = FALSE)

  if (interactive() && Sys.getenv("RSTUDIO") == "1") {
    utils::View(combined, title = paste("CellBender Comparison:", mode))
  }
  
  return(invisible(combined))
}

#' Capture session info for reproducibility
#' @param output_file Path to save session info
capture_session_info <- function(output_file = NULL) {
  session <- sessionInfo()

  # Key packages for multiome
  key_packages <- c("Seurat", "Signac", "harmony", "SeuratObject",
                   "GenomicRanges", "EnsDb.Hsapiens.v86")

  # Extract versions
  versions <- sapply(key_packages, function(pkg) {
    if (pkg %in% rownames(installed.packages())) {
      as.character(packageVersion(pkg))
    } else {
      "Not installed"
    }
  })

  version_df <- data.frame(
    package = names(versions),
    version = versions,
    stringsAsFactors = FALSE
  )

  # Add R version
  r_info <- data.frame(
    package = "R",
    version = paste(R.version$major, R.version$minor, sep = "."),
    stringsAsFactors = FALSE
  )

  version_df <- rbind(r_info, version_df)

  if (!is.null(output_file)) {
    write.csv(version_df, output_file, row.names = FALSE)
    cat(sprintf("Session info saved to: %s\n", output_file))
  }

  return(version_df)
}
