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

init_project <- function(){
  if (!dir.exists(project_outdir)) dir.create(project_outdir, recursive = TRUE)
  if (!dir.exists(rdsdir)) dir.create(rdsdir, recursive = TRUE)

  # Scratch space for interactive QC + temporary outputs
  if (!dir.exists(tmpfiledir)) dir.create(tmpfiledir, recursive = TRUE)

  # Output directory for per-sample CellBender merge reports (written by pipeline1)
  if (!dir.exists(cellbender_report_dir)) dir.create(cellbender_report_dir, recursive = TRUE)
}

init_trimming_settings <- function(trimming_settings_file){
  if (file.exists(trimming_settings_file)) {
    trimming_settings_init <- readRDS(trimming_settings_file)
  } else {
    trimming_settings_init <- data.frame(
      sample = character(0),
      min_nCount_ATAC = numeric(0),
      max_nCount_ATAC = numeric(0),
      min_nCount_RNA = numeric(0),
      max_nCount_RNA = numeric(0),
      min_nss = numeric(0),
      max_nss = numeric(0),
      max_percentMT = numeric(0),
      min_TSS = numeric(0),
      max_TSS = numeric(0),
      stringsAsFactors = FALSE
    )
  }
  return(trimming_settings_init)
}

init_kde_settings <- function(kde_settings_file){
  if (file.exists(kde_settings_file)) {
    kde_settings_init <- readRDS(kde_settings_file)
  } else {
    kde_settings_init <- data.frame(
      sample = character(0),
      atac_percentile = numeric(0),
      rna_percentile = numeric(0),
      combine_method = character(0)
    )
  }
  return(kde_settings_init)
}

read_trimming_settings <- function(trimming_settings_file){
  if (!file.exists(trimming_settings_file)) {
    stop(sprintf("Trimming settings file not found: %s", trimming_settings_file))
  }
  trimming_settings_read <- readRDS(trimming_settings_file)
  return(trimming_settings_read)
}

read_kde_settings <- function(kde_settings_file){
  if (!file.exists(kde_settings_file)) {
    stop(sprintf("KDE settings file not found: %s", kde_settings_file))
  }
  kde_settings_read <- readRDS(kde_settings_file)
  return(kde_settings_read)
}

verify_trimming_settings <- function(trimming_settings, my_trimming_settings,
                                     quiet = FALSE) {
  sample_name <- my_trimming_settings$sample
  existing_row <- trimming_settings[trimming_settings$sample == sample_name, ]
  if (nrow(existing_row) == 0) {
    cat(sprintf("Sample '%s' is new. No previous trimming settings to compare to.\n", sample_name))
  } else {
    if (quiet == FALSE) {
      cat(sprintf("Current settings for sample '%s':\n", sample_name))
      print(existing_row)
      cat("Proposed new settings:\n")
      # Reorder proposed settings to match the data frame's column order
      print(as.data.frame(my_trimming_settings, stringsAsFactors = FALSE)[, colnames(trimming_settings), drop = FALSE])
    }
    # Show changes
    changed <- sapply(colnames(trimming_settings), function(field) {
      old <- existing_row[[field]]
      new <- my_trimming_settings[[field]]
      !identical(old, new)
    })
    
    if (any(changed)) {
      cat("Fields that will change:\n")
      print(colnames(trimming_settings)[changed])
    } else {
      cat("No changes detected for this sample.\n")
    }
  }
}

# Updates or adds a sample's settings in the in-memory data.frame
update_trimming_settings <- function(trimming_settings, my_trimming_settings) {
  sample_name <- my_trimming_settings$sample
  if (is.null(trimming_settings)) {
    trimming_settings <- as.data.frame(my_trimming_settings, stringsAsFactors = FALSE)
  } else {
    idx <- which(trimming_settings$sample == sample_name)
    if (length(idx) == 0) {
      # Add new row
      trimming_settings <- rbind(trimming_settings, as.data.frame(my_trimming_settings, stringsAsFactors = FALSE))
      cat(sprintf("Added new settings for sample '%s'.\n", sample_name))
    } else {
      # Update existing row
      for (field in names(my_trimming_settings)) {
        trimming_settings[idx, field] <- my_trimming_settings[[field]]
      }
      cat(sprintf("Updated settings for sample '%s'.\n", sample_name))
    }
  }
  rownames(trimming_settings) <- NULL
  return(trimming_settings)
}

verify_kde_settings <- function(kde_settings, my_kde_settings) {
  sample_name <- my_kde_settings$sample
  existing_row <- kde_settings[kde_settings$sample == sample_name, ]
  if (nrow(existing_row) == 0) {
    cat(sprintf("Sample '%s' is new. No existing KDE settings.\n", sample_name))
  } else {
    cat(sprintf("Current KDE settings for sample '%s':\n", sample_name))
    print(existing_row)
    cat("Proposed new KDE settings:\n")
    print(as.data.frame(my_kde_settings, stringsAsFactors = FALSE)[, colnames(kde_settings), drop = FALSE])
    # Show changes
    changed <- sapply(colnames(kde_settings), function(field) {
      old <- existing_row[[field]]
      new <- my_kde_settings[[field]]
      !identical(old, new)
    })
    if (any(changed)) {
      cat("Fields that will change:\n")
      print(colnames(kde_settings)[changed])
    } else {
      cat("No changes detected for this sample.\n")
    }
  }
}

update_kde_settings <- function(kde_settings, my_kde_settings) {
  sample_name <- my_kde_settings$sample
  if (is.null(kde_settings)) {
    kde_settings <- as.data.frame(my_kde_settings, stringsAsFactors = FALSE)
  } else {
    idx <- which(kde_settings$sample == sample_name)
    if (length(idx) == 0) {
      kde_settings <- rbind(kde_settings, as.data.frame(my_kde_settings,
                      stringsAsFactors = FALSE)[, colnames(kde_settings), drop = FALSE])
      cat(sprintf("Added new KDE settings for sample '%s'.\n", sample_name))
    } else {
      for (field in names(my_kde_settings)) {
        kde_settings[idx, field] <- my_kde_settings[[field]]
      }
      cat(sprintf("Updated KDE settings for sample '%s'.\n", sample_name))
    }
  }
  rownames(kde_settings) <- NULL
  return(kde_settings)
}

generate_sample_report <- function(
    sample,
    base_obj,
    trim_obj,
    kde_obj,
    pipeline1_path,
    trimming_settings,
    kde_settings,
    report_dir = project_outdir
) {
  report_path <- file.path(report_dir, paste0(sample, "_pipeline1_report.txt"))
  
  # Gather trimming settings for this sample
  this_trim <- trimming_settings[trimming_settings$sample == sample, , drop = FALSE]
  
  # Compose trimming settings as readable lines
  trim_lines <- paste(sprintf("  %s = %s", names(this_trim), as.character(this_trim[1,])))
  
  # KDE trim settings (allow variable arguments)
  params_row <- kde_settings[kde_settings$sample == sample, , drop = FALSE]
  
  atac_percentile <- as.numeric(params_row$atac_percentile[1])
  rna_percentile  <- as.numeric(params_row$rna_percentile[1])
  combine_method  <- as.character(params_row$combine_method[1])
  
  kde_lines <- c(
    sprintf("  atac_percentile = %s", atac_percentile),
    sprintf("  rna_percentile = %s", rna_percentile),
    sprintf("  combine_method = %s", combine_method)
  )
  
  base_rna_counts <- if ("RNA" %in% names(base_obj)) GetAssayData(base_obj[["RNA"]], layer = "counts") else NULL
  base_atac_counts <- if ("ATAC" %in% names(base_obj)) GetAssayData(base_obj[["ATAC"]], layer = "counts") else NULL
  base_rna_features <- if (!is.null(base_rna_counts)) nrow(base_rna_counts) else NA
  base_atac_features <- if (!is.null(base_atac_counts)) nrow(base_atac_counts) else NA
  
  trim_rna_counts <- if ("RNA" %in% names(trim_obj)) GetAssayData(trim_obj[["RNA"]], layer = "counts") else NULL
  trim_atac_counts <- if ("ATAC" %in% names(trim_obj)) GetAssayData(trim_obj[["ATAC"]], layer = "counts") else NULL
  trim_rna_features <- if (!is.null(trim_rna_counts)) nrow(trim_rna_counts) else NA
  trim_atac_features <- if (!is.null(trim_atac_counts)) nrow(trim_atac_counts) else NA
  
  # RNA/ATAC cell and feature counts after KDE trim
  rna_cells <- ncol(GetAssayData(kde_obj[["RNA"]], layer = "counts"))
  rna_features <- nrow(GetAssayData(kde_obj[["RNA"]], layer = "counts"))
  atac_cells <- ncol(GetAssayData(kde_obj[["ATAC"]], layer = "counts"))
  atac_features <- nrow(GetAssayData(kde_obj[["ATAC"]], layer = "counts"))
  
  writeLines(c(
    sprintf("Pipeline1 Sample Report for: %s", sample),
    sprintf("Date: %s", Sys.time()),
    "",
    sprintf("Step 1: Cells at base (import): %d", nrow(base_obj@meta.data)),
    sprintf("Step 1: Features at base (import):"),
    sprintf("  RNA features (base): %s", ifelse(is.na(base_rna_features), "NA", as.character(base_rna_features))),
    sprintf("  ATAC features (base): %s", ifelse(is.na(base_atac_features), "NA", as.character(base_atac_features))),
    sprintf("Step 2: Cells after 1D trim: %d", nrow(trim_obj@meta.data)),
    sprintf("Step 2: Features after 1D trim:"),
    sprintf("  RNA features (after 1D trim): %s", ifelse(is.na(trim_rna_features), "NA", as.character(trim_rna_features))),
    sprintf("  ATAC features (after 1D trim): %s", ifelse(is.na(trim_atac_features), "NA", as.character(trim_atac_features))),
    sprintf("Step 3: Cells after KDE trim: %d", nrow(kde_obj@meta.data)),
    sprintf("Step 3: Features remaining after KDE trim (final QC):"),
    sprintf("  RNA features (after KDE trim): %d", rna_features),
    sprintf("  ATAC features (after KDE trim): %d", atac_features),
    "",
    "The following initial 1D trimming settings were used:",
    trim_lines,
    "",
    "KDE trim settings:",
    kde_lines,
    "",
    sprintf("Final object saved at: %s", pipeline1_path)
  ), con = report_path)
}

#' Compare CellBender merge reports across samples
#'
#' @param mode "qc" or "pipeline" (determines report directory)
#' @param samplelist Sample names to compare (if NULL, reads from 
#'   trimming_settings_file)
compare_cellbender_reports <- function(
  mode = c("qc", "pipeline"), 
  samplelist = NULL
) {
  mode <- match.arg(mode)
  
  # If samplelist not provided, read from trimming settings
  if (is.null(samplelist)) {
    trimming_settings <- read_trimming_settings(trimming_settings_file)
    samplelist <- trimming_settings$sample
    cat(sprintf(
      "Using samples from trimming_settings_file: %s\n", 
      trimming_settings_file
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
  cat(sprintf(
    "\n=== CellBender Merge Report Comparison (%s mode) ===\n", 
    mode
  ))
  cat(sprintf("Samples: %d\n", nrow(combined)))
  cat(sprintf("Report directory: %s\n\n", report_dir))
  
  print(combined, row.names = FALSE)
  
  # Write aggregated report to tmp for easy viewing
  agg_path <- file.path(report_dir, "aggregated_comparison.csv")
  write.csv(combined, agg_path, row.names = FALSE)
  cat(sprintf("\nAggregated report saved to: %s\n", agg_path))
  cat("View with: less -S <path> or import into your preferred tool\n\n")

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
