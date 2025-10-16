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
    atac_percentile = 0.95,
    rna_percentile = 0.95,
    combine_method = "intersection",
    report_dir = "/projects/opioid/project_export"
) {
  report_path <- file.path(report_dir, paste0(sample, "_pipeline1_report.txt"))
  
  # Gather trimming settings for this sample
  this_trim <- trimming_settings[trimming_settings$sample == sample, , drop = FALSE]
  
  # Compose trimming settings as readable lines
  trim_lines <- paste(sprintf("  %s = %s", names(this_trim), as.character(this_trim[1,])))
  
  # KDE trim settings (allow variable arguments)
  kde_lines <- c(
    sprintf("  atac_percentile = %s", atac_percentile),
    sprintf("  rna_percentile = %s", rna_percentile),
    sprintf("  combine_method = %s", combine_method)
  )
  
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
    sprintf("Step 2: Cells after 1D trim: %d", nrow(trim_obj@meta.data)),
    sprintf("Step 3: Cells after KDE trim: %d", nrow(kde_obj@meta.data)),
    "",
    sprintf("Features remaining after trimming:"),
    sprintf("  RNA features: %d", rna_features),
    sprintf("  ATAC features: %d", atac_features),
    sprintf("Cells remaining after trimming:"),
    sprintf("  RNA cells: %d", rna_cells),
    sprintf("  ATAC cells: %d", atac_cells),
    "",
    "The following initial trimming settings were used:",
    trim_lines,
    "",
    "KDE trim settings:",
    kde_lines,
    "",
    sprintf("Final object saved at: %s", pipeline1_path)
  ), con = report_path)
}