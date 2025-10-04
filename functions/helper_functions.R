#helper functions used by project, but not involved in the scientific workflow
update_provenance <- function(seurat_obj, milestone, params = NULL) {
  setwd(Rmultiome_files)
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

init_trimming_settings <- function(trimming_settings_file){
  if (file.exists(trimming_settings_file)) {
    trimming_settings_init <- readRDS(trimming_settings_file)
  } else {
    trimming_settings_init <- NULL
  }
  return(trimming_settings_init)
}

read_trimming_settings <- function(trimming_settings_file){
  if (!file.exists(trimming_settings_file)) {
    stop(sprintf("Trimming settings file not found: %s", trimming_settings_file))
  }
  trimming_settings_read <- readRDS(trimming_settings_file)
  return(trimming_settings_read)
}

verify_trimming_settings <- function(trimming_settings, my_trimming_settings) {
  sample_name <- my_trimming_settings$sample
  existing_row <- trimming_settings[trimming_settings$sample == sample_name, ]
  if (nrow(existing_row) == 0) {
    cat(sprintf("Sample '%s' is new. No existing trimming settings.\n", sample_name))
  } else {
    cat(sprintf("Current settings for sample '%s':\n", sample_name))
    print(existing_row)
    cat("Proposed new settings:\n")
    # Reorder proposed settings to match the data frame's column order
    print(as.data.frame(my_trimming_settings, stringsAsFactors = FALSE)[, colnames(trimming_settings), drop = FALSE])
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