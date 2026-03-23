#trimming functions

#' trimSample: Trims cells and features from a Seurat object using sample-specific thresholds.
#'
#' This version expects the Seurat object to have its sample identity in the @project.name slot.
#' It looks up trimming parameters in the global `pipeline1_settings` data.frame.
#'
#' @param seurat_obj A Seurat object (with @project.name set to sample name)
#' @param pipeline1_settings (optional) Data frame with sample-specific trimming parameters.
#'        If not supplied, will look for `pipeline1_settings` in the global environment.
#' @return A filtered/truncated Seurat object
#' @export
trimSample <- function(seurat_obj, pipeline1_settings = NULL) {
  # Get sample name from Seurat object
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
    stop(sprintf("No trimming settings found for sample '%s'.", sample_name))
  }
  # Extract thresholds, with names matching your pipeline1_settings
  max_nCount_ATAC     <- params$max_nCount_ATAC
  max_nCount_RNA     <- params$max_nCount_RNA
  min_nCount_ATAC     <- params$min_nCount_ATAC
  min_nCount_RNA     <- params$min_nCount_RNA
  max_nss      <- params$max_nss
  min_nss      <- params$min_nss
  max_TSS   <- params$max_TSS
  min_TSS   <- params$min_TSS
  max_percentMT <- params$max_percentMT
  
  # Ensure required metadata columns exist
  if (!all(c("nCount_RNA", "nFeature_RNA", "percent.mt") %in% colnames(seurat_obj@meta.data))) {
    stop("Required metadata columns not found in Seurat object.")
  }
  
  # Cell Filtering (QC)
  trimmed <- subset(
    seurat_obj,
    subset = nCount_RNA < max_nCount_RNA &
      nCount_RNA > min_nCount_RNA &
      nCount_ATAC < max_nCount_ATAC &
      nCount_ATAC > min_nCount_ATAC &
      percent.mt < max_percentMT &
      TSS.enrichment < max_TSS &
      TSS.enrichment > min_TSS &
      nucleosome_signal > min_nss &
      nucleosome_signal < max_nss
    # Add further criteria for nslt, TSSesgt, etc., as needed
  )
  
  # filter genes with zero counts across all cells
  DefaultAssay(trimmed) <- "RNA"
  counts <- GetAssayData(trimmed, layer = "counts")
  nonzero_genes <- rowSums(counts) > 0
  nonzero_gene_names <- rownames(counts)[nonzero_genes]
  trimmed[["RNA"]] <- subset(trimmed[["RNA"]], features = nonzero_gene_names)
  
  #filter features with zero counts across all cells
  DefaultAssay(trimmed) <- "ATAC"
  counts <- GetAssayData(trimmed, layer = "counts")
  nonzero_peaks <- rowSums(counts) > 0
  nonzero_peak_names <- rownames(counts)[nonzero_peaks]
  trimmed[["ATAC"]] <- subset(trimmed[["ATAC"]], features = nonzero_peak_names)
  
  return(trimmed)
}

get_perc_level <- function(kde, kdepercent) {
  dz <- sort(as.vector(kde$z), decreasing = TRUE)
  cumprob <- cumsum(dz) / sum(dz)
  dz[which(cumprob >= kdepercent)[1]]
}

get_density_values <- function(x, y, kde) {
  ix <- findInterval(x, kde$x)
  iy <- findInterval(y, kde$y)
  ix <- pmax(pmin(ix, length(kde$x)-1), 1)
  iy <- pmax(pmin(iy, length(kde$y)-1), 1)
  sapply(seq_along(x), function(i) kde$z[ix[i], iy[i]])
}

kdeTrimSample <- function(seurat_obj, pipeline1_settings = NULL, qc_report = FALSE) {
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
  
  # Required KDE columns
  required_cols <- c("atac_percentile", "rna_percentile", "combine_method")
  missing_cols <- setdiff(required_cols, names(params))
  if (length(missing_cols) > 0) {
    stop(sprintf("KDE settings for sample '%s' are missing required columns: %s",
                 sample_name, paste(missing_cols, collapse = ", ")))
  }
  
  atac_percentile <- as.numeric(params$atac_percentile)
  rna_percentile <- as.numeric(params$rna_percentile)
  combine_method <- as.character(params$combine_method)
  
  # Check for missing or NA values
  if (is.na(atac_percentile) || is.na(rna_percentile) || is.na(combine_method) || combine_method == "") {
    stop(sprintf("Missing or invalid KDE settings for sample '%s': atac_percentile = %s, rna_percentile = %s, combine_method = %s",
                 sample_name, atac_percentile, rna_percentile, combine_method))
  }
  
  df <- seurat_obj@meta.data
  
  # For ATAC
  x_atac <- df$nCount_ATAC
  y_atac <- df$TSS.enrichment
  
  # For RNA
  x_rna <- df$nCount_RNA
  y_rna <- df$percent.mt
  
  #correlation of metrics
  corCount <- cor(x_atac, x_rna)
  corQual <- cor(y_atac, y_rna, use="complete.obs")
  cat(sprintf(
    "\nInformational Message:\nData correlation: %f correlation of atac to rna counts,\n %f correlation of TSS.enrichment to percent.mt\n",
    corCount, corQual
  ))
  
  # KDE objects for calculations
  kde_atac <- kde2d(x_atac, y_atac, n = 100)
  kde_rna  <- kde2d(x_rna, y_rna, n = 100)
  
  # Setting levels
  level_atac <- get_perc_level(kde_atac, kdepercent = atac_percentile)
  level_rna  <- get_perc_level(kde_rna, kdepercent = rna_percentile)
  
  # Assign density for each cell
  dens_atac <- get_density_values(x_atac, y_atac, kde_atac)
  dens_rna  <- get_density_values(x_rna, y_rna, kde_rna)
  
  # Which pass the threshold? (≥ level)
  pass_atac <- dens_atac >= level_atac
  pass_rna  <- dens_rna  >= level_rna
  top_cells <- pass_atac & pass_rna
  
  cor_top <- cor(x_atac[top_cells], x_rna[top_cells])
  cor_top_quality <- cor(y_atac[top_cells], y_rna[top_cells], use="complete.obs")
  cat(sprintf("\nInformational message:\nCorrelation in KDE-filtered subset: %f (counts), %f (quality)\n", cor_top, cor_top_quality))
  
  PMTmean_before <- mean(df$percent.mt, na.rm=TRUE)
  PMTmean_after <- mean(df$percent.mt[pass_atac & pass_rna], na.rm=TRUE)
  cat(sprintf("Average percent.mt before: %.3f, after KDE filter: %.3f\n", PMTmean_before, PMTmean_after))
  
  TSSmean_before <- mean(df$TSS.enrichment, na.rm=TRUE)
  TSSmean_after <- mean(df$TSS.enrichment[pass_atac & pass_rna], na.rm=TRUE)
  cat(sprintf("Average TSS.enrichment before: %.3f, after KDE filter: %.3f\n", TSSmean_before, TSSmean_after))
  
  # Combine according to method
  if (combine_method == "intersection") {
    keep_cells <- rownames(df)[pass_atac & pass_rna]
  } else if (combine_method == "union") {
    keep_cells <- rownames(df)[pass_atac | pass_rna]
  } else {
    stop(sprintf(
      "Invalid combine_method '%s'. Must be either 'intersection' or 'union'.",
      combine_method
    ))
  }
  
  if (qc_report == TRUE) {
    n_before <- nrow(df)
    n_pass_atac <- sum(pass_atac)
    n_pass_rna <- sum(pass_rna)
    n_after <- sum(top_cells)
    percent_retained <- 100 * n_after / n_before
  
    cat(sprintf(
      "\nCells before KDE filtering: %d\nCells passing ATAC KDE: %d\nCells passing RNA KDE: %d\nCells after combined KDE filtering: %d (%.2f%% retained)\n",
      n_before, n_pass_atac, n_pass_rna, n_after, percent_retained))
    cat(sprintf("\nKDE settings for sample '%s': atac_percentile = %s, rna_percentile = %s, combine_method = %s",
           sample_name, atac_percentile, rna_percentile, combine_method))
  }
  # Subset Seurat object
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  
  # filter genes with zero counts across all cells (RNA assay)
  DefaultAssay(seurat_obj) <- "RNA"
  counts <- GetAssayData(seurat_obj, layer = "counts")
  nonzero_genes <- rowSums(counts) > 0
  nonzero_gene_names <- rownames(counts)[nonzero_genes]
  seurat_obj[["RNA"]] <- subset(seurat_obj[["RNA"]], features = nonzero_gene_names)
  
  # filter peaks with zero counts across all cells (ATAC assay)
  DefaultAssay(seurat_obj) <- "ATAC"
  counts <- GetAssayData(seurat_obj, layer = "counts")
  nonzero_peaks <- rowSums(counts) > 0
  nonzero_peak_names <- rownames(counts)[nonzero_peaks]
  seurat_obj[["ATAC"]] <- subset(seurat_obj[["ATAC"]], features = nonzero_peak_names)
  
  return(seurat_obj)
}

#' Remove doublets using scDblFinder
#'
#' Detects and removes doublets from a Seurat object using scDblFinder.
#' Looks up expected doublet rate from pipeline1_settings dataframe.
#'
#' @param seurat_obj A Seurat object (with @project.name set to sample name)
#' @param pipeline1_settings (optional) Data frame with sample-specific settings
#'        If not supplied, will look for `pipeline1_settings` in global environment
#' @param qc_report If TRUE, prints detailed output and plots (for run_qc.R)
#'        If FALSE, minimal output (for run_pipeline1.R)
#' @return List with:
#'   - obj: Filtered Seurat object containing only singlets
#'   - stats: List of doublet detection statistics for reporting
#' @export
doubletRemoveSample <- function(seurat_obj, pipeline1_settings = NULL, qc_report = FALSE) {

  # Get sample name from Seurat object
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
    stop(sprintf("No doublet settings found for sample '%s'.", sample_name))
  }

  # Extract expected doublet rate
  expected_dbr <- params$expected_dbr
  if (is.na(expected_dbr)) {
    stop(sprintf("expected_dbr not set for sample '%s'. Run QC first to calculate this value.",
                sample_name))
  }

  # Get global settings
  if (!exists("random_seed", envir = .GlobalEnv)) {
    stop("random_seed not found. Source project_settings.R first.")
  }
  if (!exists("doublet_rate_sd", envir = .GlobalEnv)) {
    stop("doublet_rate_sd not found. Source project_settings.R first.")
  }

  random_seed <- get("random_seed", envir = .GlobalEnv)
  doublet_rate_sd <- get("doublet_rate_sd", envir = .GlobalEnv)

  # Calculate parameters
  n_cells <- ncol(seurat_obj)
  dbr_rate <- expected_dbr / n_cells

  if (qc_report) {
    cat(sprintf("\n=== scDblFinder: %s ===\n", sample_name))
    cat(sprintf("Cells before doublet removal: %d\n", n_cells))
    cat(sprintf("Expected doublets: %.1f (%.2f%%)\n", expected_dbr, 100*dbr_rate))
    cat(sprintf("Running scDblFinder with dbr=%.4f, dbr.sd=%.3f\n",
        dbr_rate, doublet_rate_sd))
  }

  # Convert to SingleCellExperiment (using Seurat's conversion method)
  sce <- as.SingleCellExperiment(seurat_obj, assay = "RNA")

  # Set seed for reproducibility
  set.seed(random_seed)

  # Run scDblFinder
  sce <- scDblFinder::scDblFinder(
    sce,
    clusters = FALSE,
    dbr = dbr_rate,
    dbr.sd = doublet_rate_sd,
    verbose = qc_report
  )

  # Extract results
  doublet_class <- sce$scDblFinder.class
  doublet_score <- sce$scDblFinder.score

  n_doublets <- sum(doublet_class == "doublet")
  n_singlets <- sum(doublet_class == "singlet")
  pct_doublets <- 100 * n_doublets / n_cells

  # Calculate statistics
  threshold <- min(doublet_score[doublet_class == "doublet"])
  singlet_score_median <- median(doublet_score[doublet_class == "singlet"])
  doublet_score_median <- median(doublet_score[doublet_class == "doublet"])

  # Report results
  if (qc_report) {
    cat(sprintf("\n=== Results ===\n"))
    cat(sprintf("Singlets: %d (%.1f%%)\n", n_singlets, 100*n_singlets/n_cells))
    cat(sprintf("Doublets: %d (%.1f%%)\n", n_doublets, pct_doublets))
    cat(sprintf("Expected: %.1f (%.2f%%)\n", expected_dbr, 100*dbr_rate))
    cat(sprintf("Difference: %+.1f doublets (%+.2f%%)\n\n",
               n_doublets - expected_dbr, pct_doublets - 100*dbr_rate))

    cat("Doublet score summary:\n")
    print(summary(doublet_score))
    cat("\nSinglet scores:\n")
    print(summary(doublet_score[doublet_class == "singlet"]))
    cat("Doublet scores:\n")
    print(summary(doublet_score[doublet_class == "doublet"]))

    # Two-panel visualization
    par(mfrow = c(1, 2))

    # Panel 1: Full distribution
    hist(doublet_score, breaks = 50,
         main = paste("All Scores -", sample_name),
         xlab = "scDblFinder Score",
         ylab = "Frequency",
         col = "lightblue")
    abline(v = threshold, col = "red", lwd = 2, lty = 2)
    legend("topright",
           legend = sprintf("Threshold: %.3f", threshold),
           col = "red", lty = 2, lwd = 2, cex = 0.8)

    # Panel 2: Doublet region zoomed (0.7-1.0)
    doublet_region <- doublet_score[doublet_score > 0.7]
    if (length(doublet_region) > 0) {
      hist(doublet_region, breaks = 30,
           main = "Doublet Region (0.7-1.0)",
           xlab = "scDblFinder Score",
           ylab = "Frequency",
           col = "salmon",
           xlim = c(0.7, 1.0))
      abline(v = threshold, col = "red", lwd = 2, lty = 2)
      abline(v = doublet_score_median, col = "darkred", lwd = 2, lty = 3)
      legend("topleft",
             legend = c(sprintf("Threshold: %.3f", threshold),
                       sprintf("Median: %.3f", doublet_score_median)),
             col = c("red", "darkred"),
             lty = c(2, 3),
             lwd = 2,
             cex = 0.7)
      text(0.85, max(hist(doublet_region, breaks = 30, plot = FALSE)$counts) * 0.9,
           sprintf("%d doublets\n(%.1f%%)", n_doublets, pct_doublets),
           cex = 0.9)
    }

    par(mfrow = c(1, 1))

    cat("\nFiltering to singlets only...\n")
  }

  # Filter Seurat object to singlets only
  singlets_idx <- which(doublet_class == "singlet")
  filtered_obj <- seurat_obj[, singlets_idx]

  if (qc_report) {
    cat(sprintf("Cells after doublet removal: %d\n", ncol(filtered_obj)))
  } else {
    cat(sprintf("Removed %d doublets from %s, %d cells remaining\n",
               n_doublets, sample_name, ncol(filtered_obj)))
  }

  # Compile statistics for reporting
  doublet_stats <- list(
    n_cells_before = n_cells,
    expected_dbr = expected_dbr,
    n_doublets = n_doublets,
    n_singlets = n_singlets,
    pct_doublets = pct_doublets,
    threshold = threshold,
    singlet_score_median = singlet_score_median,
    doublet_score_median = doublet_score_median
  )

  # Return both the filtered object and statistics
  return(list(
    obj = filtered_obj,
    stats = doublet_stats
  ))
}
