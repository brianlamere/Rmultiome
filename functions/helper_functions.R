#helper functions used by project, but not involved in the scientific workflow
update_provenance <- function(seurat_obj, milestone, params = NULL) {
  setwd(Rmultiome_path)
  git_hash <- tryCatch(system("git rev-parse HEAD", intern = TRUE),
                       error=function(e) NA)
  git_branch <- tryCatch(system("git rev-parse --abbrev-ref HEAD",
                                intern = TRUE), error=function(e) NA)
  git_tag <- tryCatch(system("git describe --tags --always",
                             intern = TRUE), error=function(e) NA)
  git_dirty <- tryCatch({
    if (length(system("git status --porcelain", intern = TRUE)) > 0) "YES" else "NO"
  }, error=function(e) NA)
  info <- list(
    step = milestone,
    timestamp = Sys.time(),
    git_commit = git_hash,
    git_branch = git_branch,
    git_tag = git_tag,
    git_dirty = git_dirty,
    params = params
  )
  if (is.null(seurat_obj@misc$provenance)) {
    seurat_obj@misc$provenance <- list()
  }
  seurat_obj@misc$provenance[[length(seurat_obj@misc$provenance) + 1]] <- info
  seurat_obj
}

#this will print the provenance log.
print_provenance <- function(seurat_obj) {
  prov <- seurat_obj@misc$provenance
  if (is.null(prov) || length(prov) == 0) {
    cat("No provenance information found.\n")
    return(invisible(NULL))
  }
  for (i in seq_along(prov)) {
    cat(sprintf("Step %d: %s\n", i, prov[[i]]$step))
    cat(sprintf("  Timestamp: %s\n", prov[[i]]$timestamp))
    cat(sprintf("  Git commit: %s\n", prov[[i]]$git_commit))
    cat(sprintf("  Branch: %s | Tag: %s | Dirty: %s\n",
                prov[[i]]$git_branch, prov[[i]]$git_tag, prov[[i]]$git_dirty))
    if (!is.null(prov[[i]]$params)) {
      cat("  Params:\n")
      print(prov[[i]]$params)
    }
    cat("\n")
  }
}

get_rds_path <- function(sample, milestone) {
  file.path(rdsdir, paste0(sample, "_", milestone, ".rds"))
}

#this just computes them mean percent of mitochondrial genes in an entire sample
meanMT <- function(samplename) {
  return(mean(samplename@meta.data[["percent.mt"]], na.rm = TRUE))
}

findElbow <- function(seurat_obj) {
  pc_sd <- seurat_obj[["pca"]]@stdev
  pc_table <- data.frame(PC = 1:length(pc_sd), StandardDeviation = pc_sd)
  
  # Custom ElbowPlot with grid
  ggplot(pc_table, aes(x = PC, y = StandardDeviation)) +
    geom_point() +
    geom_line() +
    theme_minimal() +
    scale_x_continuous(breaks = 1:length(pc_sd)) +
    scale_y_continuous(breaks = pretty(pc_sd)) +
    labs(title = "Elbow Plot of PCs", x = "Principal Component", y = "Standard Deviation") +
    theme(panel.grid.major = element_line(color = "gray", size = 0.5),
          panel.grid.minor = element_line(color = "lightgray", size = 0.2))
}

# Helper function to summarize results
summarize_results <- function(results) {
  # For each result, compute summary metrics
  metrics <- lapply(results, function(res) {
    # Example: Use variance of first 20 PCs
    variance_total <- sum(res$variance[1:20])
    # Example: Mean jackstraw p-value of first 20 PCs
    jackstraw_mean <- mean(res$jackstraw[1:20], na.rm = TRUE)
    # Example: Number of significant PCs (p < 0.05)
    jackstraw_sig <- sum(res$jackstraw[1:20] < 0.05, na.rm = TRUE)
    data.frame(
      variance_total = variance_total,
      jackstraw_mean = jackstraw_mean,
      jackstraw_sig = jackstraw_sig,
      n_clusters = res$n_clusters,
      n_singletons = res$n_singletons,
      # Combine parameter columns
      res$params
    )
  })
  do.call(rbind, metrics)
}

#' Initialize project structure and settings
#'
#' Creates directories and project_settings.R, then loads settings into global env
#' Call once at start of run_qc.R (stepping through interactively)
#'
#' @param random_seed Global random seed (default 42)
#' @param use_cellbender Use CellBender RNA (default TRUE)
#' @param use_scdblfinder Detect/remove doublets with scDblFinder (default TRUE)
#' @param doublet_rate_per_1000 Expected doublets per 1000 cells (default 8.0, 10X formula)
#' @param doublet_rate_sd Uncertainty in doublet rate (default 0.015; use 1.0 if very uncertain)
#' @param project_name Project name (default "multiome_project")
#' @param genome_build Genome build (default "hg38")
#' @return Invisible list of project settings (also loaded into .GlobalEnv)
init_project <- function(random_seed = 42,
                        use_cellbender = TRUE,
                        use_scdblfinder = TRUE,
                        doublet_rate_per_1000 = 8.0,
                        doublet_rate_sd = 0.015,
                        project_name = "multiome_project",
                        genome_build = "hg38") {

  cat("\n=== Initializing Project ===\n")

  # Require system_settings
  if (!exists("project_export")) {
    stop("Source system_settings.R first!")
  }

  # Create directories
  dirs <- c(project_export, rdsdir, tmpfiledir)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE)
      cat(sprintf("  Created: %s\n", d))
    }
  }

  # Create project_settings.R
  settings_path <- file.path(project_export, "project_settings.R")

  if (file.exists(settings_path)) {
    cat(sprintf("\nproject_settings.R exists, loading...\n"))
  } else {
    # Generate settings file
    settings_content <- sprintf(
'# ====================================================================
# PROJECT SETTINGS
# Created: %s
# Manually edit this file to change project-wide parameters
# After editing, re-source: source(file.path(project_export, "project_settings.R"))
# ====================================================================

# === Reproducibility ===
random_seed <- %d  # Global seed for all stochastic processes

# === Feature toggles ===
use_cellbender <- %s    # Use CellBender-corrected RNA counts
use_scdblfinder <- %s   # Detect and remove doublets with scDblFinder

# === Analysis scope ===
standard_chroms <- paste0("chr", c(1:22, "X", "Y"))  # Chromosomes to include

# === Doublet detection parameters ===
# Based on 10X Chromium documentation: ~1%% doublet rate per 1000 cells captured
doublet_rate_per_1000 <- %.1f  # Expected doublets per 1000 cells
doublet_rate_sd <- %.3f        # Uncertainty (0.015 = confident, 1.0 = use misclassification only)

# === Project metadata ===
project_name <- "%s"
genome_build <- "%s"
species <- "Homo sapiens"
tissue <- "brain"  # Edit as appropriate
analysis_version <- "1.0.0"
analysis_date <- "%s"

# ====================================================================
# END PROJECT SETTINGS
# ====================================================================
',
      Sys.Date(),
      random_seed,
      use_cellbender,
      use_scdblfinder,
      doublet_rate_pct,
      doublet_rate_sd,
      project_name,
      genome_build,
      Sys.Date()
    )

    writeLines(settings_content, settings_path)
    cat(sprintf("  Created: %s\n", settings_path))
  }

  # Load into global environment (like run_pipeline1.R does)
  source(settings_path, local = .GlobalEnv)

  cat("\nProject settings loaded:\n")
  cat(sprintf("  random_seed = %d\n", random_seed))
  cat(sprintf("  use_cellbender = %s\n", use_cellbender))
  cat(sprintf("  use_scdblfinder = %s\n", use_scdblfinder))
  cat(sprintf("  doublet_rate_per_1000 = %.1f\n", doublet_rate_per_1000))
  cat(sprintf("  project_name = %s\n", project_name))

  cat("\n=== Ready for QC ===\n\n")

  # Return settings invisibly
  invisible(list(
    random_seed = random_seed,
    use_cellbender = use_cellbender,
    use_scdblfinder = use_scdblfinder,
    doublet_rate_per_1000 = doublet_rate_per_1000,
    doublet_rate_sd = doublet_rate_sd,
    project_name = project_name,
    genome_build = genome_build,
    standard_chroms = standard_chroms
  ))
}

init_pipeline1_settings <- function(pipeline1_settings_file) {
  if (file.exists(pipeline1_settings_file)) {
    pipeline1_settings_init <- readRDS(pipeline1_settings_file)
  } else {
    # Use OLD column names to match existing code
    pipeline1_settings_init <- data.frame(
      sample = character(0),

      # === 1D Trimming parameters (OLD naming) ===
      min_nCount_ATAC = numeric(0),
      max_nCount_ATAC = numeric(0),
      min_nCount_RNA = numeric(0),
      max_nCount_RNA = numeric(0),
      min_nss = numeric(0),
      max_nss = numeric(0),
      max_percentMT = numeric(0),
      min_TSS = numeric(0),
      max_TSS = numeric(0),

      # === KDE parameters ===
      atac_percentile = numeric(0),
      rna_percentile = numeric(0),
      combine_method = character(0),
      kde_bandwidth = numeric(0),   # Future use

      # === Doublet detection ===
      n_cells_after_kde = numeric(0),
      expected_dbr = numeric(0),

      stringsAsFactors = FALSE
    )
  }

  return(pipeline1_settings_init)
}

#' Verify pipeline1 settings before updating
#'
#' Shows current settings vs proposed changes for a sample
#' Only compares columns present in my_settings (flexible for trim or KDE)
#'
#' @param pipeline1_settings The main settings dataframe
#' @param my_settings List of proposed settings (trim or KDE fields)
#' @param quiet If TRUE, only show changes (no full printout)
verify_pipeline1_settings <- function(pipeline1_settings, my_settings, quiet = FALSE) {

  sample_name <- my_settings$sample
  existing_row <- pipeline1_settings[pipeline1_settings$sample == sample_name, ]

  if (nrow(existing_row) == 0) {
    cat(sprintf("Sample '%s' is new. No previous settings to compare to.\n", sample_name))
  } else {

    # Only look at columns present in my_settings (ignoring "sample")
    fields_to_compare <- setdiff(names(my_settings), "sample")

    # Filter to only columns that exist in pipeline1_settings
    fields_to_compare <- intersect(fields_to_compare, colnames(pipeline1_settings))

    if (quiet == FALSE) {
      cat(sprintf("Current settings for sample '%s':\n", sample_name))
      print(existing_row[, c("sample", fields_to_compare), drop = FALSE])
      cat("Proposed new settings:\n")
      print(as.data.frame(my_settings, stringsAsFactors = FALSE)[, c("sample", fields_to_compare), drop = FALSE])
    }

    # Show changes (only for fields being updated)
    changed <- sapply(fields_to_compare, function(field) {
      old <- existing_row[[field]]
      new <- my_settings[[field]]
      !identical(old, new)
    })

    if (any(changed)) {
      cat("Fields that will change:\n")
      print(fields_to_compare[changed])
    } else {
      cat("No changes detected for this sample.\n")
    }
  }
}

#' Update pipeline1 settings for a sample
#'
#' Updates only the fields provided in my_settings, leaves others untouched
#' Works with any subset of settings (trim, KDE, or mix)
#'
#' @param pipeline1_settings The main settings dataframe
#' @param my_settings List of settings for one sample (any fields)
#' @return Updated pipeline1_settings dataframe
update_pipeline1_settings <- function(pipeline1_settings, my_settings) {

  sample_name <- my_settings$sample

  # Find existing row for this sample
  idx <- which(pipeline1_settings$sample == sample_name)

  if (length(idx) == 0) {
    # Sample doesn't exist - add new row with NAs, fill in what we have
    new_row <- pipeline1_settings[0, ]  # Empty row with correct structure
    new_row[1, ] <- NA  # Fill with NAs
    new_row[1, "sample"] <- sample_name

    # Fill in only the provided fields
    for (field in names(my_settings)) {
      if (field != "sample" && field %in% names(pipeline1_settings)) {
        new_row[1, field] <- my_settings[[field]]
      }
    }

    pipeline1_settings <- rbind(pipeline1_settings, new_row)
    cat(sprintf("Added new settings for sample '%s'\n", sample_name))

  } else {
    # Sample exists - update only the provided fields
    for (field in names(my_settings)) {
      if (field != "sample" && field %in% names(pipeline1_settings)) {
        pipeline1_settings[idx, field] <- my_settings[[field]]
      }
    }
    cat(sprintf("Updated settings for sample '%s'\n", sample_name))
  }

  rownames(pipeline1_settings) <- NULL
  return(pipeline1_settings)
}

read_cluster_settings <- function(cluster_settings_file = cluster_settings_file){
  if (!file.exists(cluster_settings_file)) {
    stop(sprintf("Cluster settings file not found: %s", cluster_settings_file))
  }
  cluster_settings_read <- readRDS(cluster_settings_file)
  return(cluster_settings_read)
}

read_celltype_settings <- function(celltype_settings_file = celltype_settings_file){
  if (!file.exists(celltype_settings_file)) {
    stop(sprintf("Celltype settings file not found: %s", celltype_settings_file))
  }
  celltype_settings_read <- readRDS(celltype_settings_file)
  return(celltype_settings_read)
}

read_harmony_settings <- function(harmony_settings_file = harmony_settings_file){
  if (!file.exists(harmony_settings_file)) {
    stop(sprintf("harmony settings file not found: %s", harmony_settings_file))
  }
  harmony_settings_read <- readRDS(harmony_settings_file)
  return(harmony_settings_read)
}

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

#' Open new graphics device if not in RStudio
#' @param width Window width in inches (default 10)
#' @param height Window height in inches (default 8)
maybe_new_device <- function(width = 10, height = 8) {
  if (!requireNamespace("rstudioapi", quietly = TRUE) ||
      !rstudioapi::isAvailable()) {
    dev.new(width = width, height = height)
  }
}

#' Open new graphics device with optional Hyprland workspace routing
#' @param width Window width in inches
#' @param height Window height in inches
#' @param workspace Hyprland workspace number (NULL = current workspace)
#' @param title Window title for Hyprland matching
maybe_new_device_workspace <- function(width = 10, height = 8,
                                      workspace = NULL,
                                      title = "R Graphics") {
  # Check if in RStudio
  in_rstudio <- requireNamespace("rstudioapi", quietly = TRUE) &&
                rstudioapi::isAvailable()

  if (in_rstudio) {
    # RStudio: plots go to Plots pane, workspace parameter ignored
    return(invisible(NULL))
  }

  # Not in RStudio, open new X11 device
  dev.new(width = width, height = height, title = title)

  # Handle Hyprland workspace routing if requested
  if (!is.null(workspace)) {
    # Check if running under Hyprland
    hypr_sig <- Sys.getenv("HYPRLAND_INSTANCE_SIGNATURE")

    if (hypr_sig == "") {
      warning(sprintf(
        "Workspace routing requested (workspace=%d) but not running under Hyprland.\nPlot window opened on current workspace.",
        workspace
      ))
      return(invisible(NULL))
    }

    # Running under Hyprland, attempt to move window
    Sys.sleep(0.2)  # Let window spawn

    # Use new Hyprland syntax: hyprctl dispatch movetoworkspacesilent
    cmd <- sprintf('hyprctl dispatch movetoworkspacesilent %d,title:%s',
                  workspace, title)
    result <- system(cmd, intern = TRUE, ignore.stderr = TRUE)

    # Optional: Check if window rule exists (helps with first-time setup)
    if (workspace == 9 && title == "R_ParamSweep") {
      # Check hyprland.conf for the expected rule
      hypr_conf <- path.expand("~/.config/hypr/hyprland.conf")
      if (file.exists(hypr_conf)) {
        conf_text <- readLines(hypr_conf, warn = FALSE)
        has_rule <- any(grepl("R_ParamSweep.*workspace.*9", conf_text))

        if (!has_rule) {
          message(
            "Note: For automatic workspace routing, add to ~/.config/hypr/hyprland.conf:\n",
            "  windowrule = match:title R_ParamSweep, workspace 9 silent\n",
            "Then reload Hyprland config with: hyprctl reload"
          )
        }
      }
    }
  }

  invisible(NULL)
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
