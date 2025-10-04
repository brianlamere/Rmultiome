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

update_trimming_settings <- function(settings_df, sample_name, ...) {
  # Check sample column exists
  if (!"sample" %in% colnames(settings_df)) stop("settings_df must have a 'sample' column.")
  args <- list(...)
  # Only allow updating existing columns (except sample)
  valid_cols <- setdiff(colnames(settings_df), "sample")
  wrong_cols <- setdiff(names(args), valid_cols)
  if (length(wrong_cols)) stop("Unknown columns in update: ", paste(wrong_cols, collapse=", "))
  # Find row index for this sample
  idx <- which(settings_df$sample == sample_name)
  # If not present, add a new row
  if (length(idx) == 0) {
    new_row <- as.list(rep(NA, length(colnames(settings_df))))
    names(new_row) <- colnames(settings_df)
    new_row$sample <- sample_name
    for (nm in names(args)) new_row[[nm]] <- args[[nm]]
    settings_df <- rbind(settings_df, as.data.frame(new_row, stringsAsFactors = FALSE))
    message(sprintf("Added new entry for sample '%s': %s", sample_name, 
                    paste(sprintf("%s=%s", names(args), args), collapse=", ")))
  } else {
    # Update in place
    old_vals <- as.list(settings_df[idx, names(args), drop = FALSE])
    for (nm in names(args)) settings_df[idx, nm] <- args[[nm]]
    message(sprintf("Updated sample '%s': %s", sample_name, 
                    paste(sprintf("%s: %s -> %s", names(args), old_vals, args), collapse=", ")))
  }
  rownames(settings_df) <- NULL
  return(settings_df)
}

update_trimming_settings_in_file <- function(settings_file, my_trimming_settings) {
  # Source the settings file to get trimming_settings
  source(settings_file, local = TRUE)
  ts <- trimming_settings
  
  sample_name <- my_trimming_settings$sample
  
  # Build a data.frame row with the correct structure/order
  new_row <- as.data.frame(my_trimming_settings, stringsAsFactors = FALSE)
  # Ensure all columns present
  for (col in setdiff(colnames(ts), names(my_trimming_settings))) {
    new_row[[col]] <- NA
  }
  new_row <- new_row[, colnames(ts)]
  
  # Check if sample exists, update or add
  if (sample_name %in% ts$sample) {
    ts[ts$sample == sample_name, ] <- new_row
  } else {
    ts <- rbind(ts, new_row)
  }
  
  # Sorting helper for mixed alphanumeric (e.g., LG2, LG10, LG100)
  extract_num <- function(x) as.numeric(sub(".*?(\\d+)$", "\\1", x))
  extract_prefix <- function(x) sub("(.*?)(\\d+)$", "\\1", x)
  ts <- ts[order(extract_prefix(ts$sample), extract_num(ts$sample)), ]
  
  # Write out new trimming_settings block to settings_file
  # Read the whole file
  lines <- readLines(settings_file)
  start_idx <- grep("trimming_settings <- data.frame\\(", lines)
  end_idx <- grep("^\\s*\\)", lines[(start_idx+1):length(lines)]) + start_idx
  # Compose new data.frame block
  df_text <- capture.output(dput(ts))
  df_text[1] <- "trimming_settings <- "  # replace with assignment
  # Replace old block with new block
  new_lines <- c(lines[1:(start_idx-1)], df_text, lines[(end_idx+1):length(lines)])
  writeLines(new_lines, settings_file)
  cat(sprintf("Updated trimming_settings for sample %s in %s\n", sample_name, settings_file))
}

verify_trimming_settings_file_changes <- function(settings_file, my_trimming_settings) {
  # Source the settings file to get trimming_settings
  source(settings_file, local = TRUE)
  ts <- trimming_settings
  sample_name <- my_trimming_settings$sample
  
  # Check if sample exists
  if (!(sample_name %in% ts$sample)) {
    cat(sprintf("Sample '%s' is a new entry. Settings to be added:\n", sample_name))
    print(my_trimming_settings)
    return(invisible(NULL))
  }
  
  # Sample exists, compare fields
  existing_row <- ts[ts$sample == sample_name, ]
  changed_fields <- list()
  for (field in names(my_trimming_settings)) {
    if (field %in% colnames(existing_row)) {
      old_val <- existing_row[[field]]
      new_val <- my_trimming_settings[[field]]
      # Use identical for NA-safe comparison
      if (!identical(old_val, new_val)) {
        changed_fields[[field]] <- list(
          old = old_val,
          new = new_val
        )
      }
    }
  }
  if (length(changed_fields) == 0) {
    cat(sprintf("No changes for sample '%s'.\n", sample_name))
  } else {
    cat(sprintf("Changes for sample '%s':\n", sample_name))
    for (field in names(changed_fields)) {
      cat(sprintf("  %s: %s -> %s\n", field,
                  as.character(changed_fields[[field]]$old),
                  as.character(changed_fields[[field]]$new)))
    }
  }
  invisible(NULL)
}

