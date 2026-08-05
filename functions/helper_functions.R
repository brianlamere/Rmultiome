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
                        analysis_version = "1.1.0a",
                        species = "Homo sapiens",
                        tissue_type = "cortex",
                        genome_build = "hg38") {

  cat("\n=== Initializing Project ===\n")

  #TODO:  do something with genome_build or similar; have a list of supported
  # genome references, and use this setting to then select which is used.
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
species <- "%s"
tissue_type <- "%s"
analysis_version <- "%s"
analysis_date <- "%s"

# ====================================================================
# END PROJECT SETTINGS
# ====================================================================
',
      Sys.Date(),
      random_seed,
      use_cellbender,
      use_scdblfinder,
      doublet_rate_per_1000,
      doublet_rate_sd,
      project_name,
      genome_build,
      species,
      tissue_type,
      analysis_version,
      Sys.Date()
    )

    writeLines(settings_content, settings_path)
    cat(sprintf("  Created: %s\n", settings_path))
  }

  # Load into global environment (like run_pipeline1.R does)
  source(settings_path, local = .GlobalEnv)

  cat("\nProject settings loaded:\n")
  cat(sprintf("  random_seed = %d\n", .GlobalEnv$random_seed))
  cat(sprintf("  use_cellbender = %s\n", .GlobalEnv$use_cellbender))
  cat(sprintf("  use_scdblfinder = %s\n", .GlobalEnv$use_scdblfinder))
  cat(sprintf("  doublet_rate_per_1000 = %.3f\n", .GlobalEnv$doublet_rate_per_1000))
  cat(sprintf("  doublet_rate_sd = %.3f\n", .GlobalEnv$doublet_rate_sd))
  cat(sprintf("  tissue_type = %s\n", .GlobalEnv$tissue_type))
  cat(sprintf("  project_name = %s\n", .GlobalEnv$project_name))
  cat(sprintf("  genome_build = %s\n", .GlobalEnv$genome_build))

  cat("\n=== Ready for QC ===\n\n")

  # Return settings invisibly
  invisible(list(
    random_seed = random_seed,
    use_cellbender = use_cellbender,
    use_scdblfinder = use_scdblfinder,
    doublet_rate_per_1000 = doublet_rate_per_1000,
    doublet_rate_sd = doublet_rate_sd,
    tissue_type = tissue_type,
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

#' Read parameter settings
#'
#' Pulls in parameter settings from stored settings file
#'
#' @param sample sample name
#' @param seurat_obj Base Seurat Object
#why is this necessary?  Who should be at the point of calling this, and not already
#have this info? Unless it is purely to make sure you have the validated, stored info
#and not emphemeral info, before doing a Data modification.
get_params <- function(
	sample_name,
	seurat_obj
	)
{
  sample_name <- seurat_obj@project.name
  # Fetch pipeline1_settings if not provided
  if (is.null(pipeline1_settings)) {
    if (!exists("pipeline1_settings", envir = .GlobalEnv)) {
      stop("pipeline1_settings not found in global environment, and not provided as argument.")
    }
    pipeline1_settings <- get("pipeline1_settings", envir = .GlobalEnv)
  }

  # Look up parameters for this sample
  params <- pipeline1_settings[pipeline1_settings$sample == sample_name, ]
  if (nrow(params) == 0) {
    stop(sprintf("No KDE settings found for sample '%s'.", sample_name))
  }
}


#' Load the consolidated marker set for the current project's tissue type.
#'
#' Reads `tissue_type` from the global environment (set by init_project via
#' project_settings.R) and dispatches to the appropriate function in
#' functions/Consolidated_Markers.R.  Returns a list with two elements:
#'   $reference_table : data.frame, one row per marker gene
#'   $marker_lists    : named list of character vectors for identify_all_celltypes()
#'
#' Usage in run_qc_merged.R:
#'   tissue_markers <- load_tissue_markers()
#'   results <- identify_all_celltypes(all_markers, tissue_markers$marker_lists, ...)
#'
#' @return list(reference_table, marker_lists)
load_tissue_markers <- function() {
  tissue <- tolower(trimws(get("tissue_type", envir = .GlobalEnv)))

  if (!tissue %in% names(TISSUE_MARKER_FUNCTIONS)) {
    stop(sprintf(
      "Unsupported tissue_type: '%s'. Supported: %s",
      tissue,
      paste(names(TISSUE_MARKER_FUNCTIONS), collapse = ", ")
    ))
  }

  markers_call <- TISSUE_MARKER_FUNCTIONS[[tissue]]()

  cat(sprintf(
    "✓ Loaded %s markers: %d genes across %d cell types\n",
    tissue,
    nrow(markers_call$reference_table),
    length(markers_call$marker_lists)
  ))

  markers_call
}
