#trimming functions

#' trimSample: Trims cells and features from a Seurat object using sample-specific thresholds.
#'
#' This version expects the Seurat object to have its sample identity in the @project.name slot.
#' It looks up trimming parameters in the global `trimming_settings` data.frame.
#'
#' @param seurat_obj A Seurat object (with @project.name set to sample name)
#' @param trimming_settings (optional) Data frame with sample-specific trimming parameters.
#'        If not supplied, will look for `trimming_settings` in the global environment.
#' @return A filtered/truncated Seurat object
#' @export
trimSample <- function(seurat_obj, trimming_settings = NULL) {
  # Get sample name from Seurat object
  sample_name <- seurat_obj@project.name
  # Fetch trimming_settings if not provided
  if (is.null(trimming_settings)) {
    if (!exists("trimming_settings", envir = .GlobalEnv)) {
      stop("trimming_settings not found in global environment, and not provided as argument.")
    }
    trimming_settings <- get("trimming_settings", envir = .GlobalEnv)
  }
  # Look up parameters for this sample
  params <- trimming_settings[trimming_settings$sample == sample_name, ]
  if (nrow(params) == 0) {
    stop(sprintf("No trimming settings found for sample '%s'.", sample_name))
  }
  # Extract thresholds, with names matching your trimming_settings
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
  counts <- GetAssayData(trimmed, slot = "counts")
  nonzero_genes <- rowSums(counts) > 0
  nonzero_gene_names <- rownames(counts)[nonzero_genes]
  trimmed[["RNA"]] <- subset(trimmed[["RNA"]], features = nonzero_gene_names)
  
  #filter features with zero counts across all cells
  DefaultAssay(trimmed) <- "ATAC"
  counts <- GetAssayData(trimmed, slot = "counts")
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

kdeTrimSample <- function(seurat_obj,
                          atac_percentile = 0.75,
                          rna_percentile = 0.75,
                          combine_method = c("intersection", "union"),
                          ...) {
  combine_method <- match.arg(combine_method)
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
  }
  
  # Subset Seurat object
  seurat_obj <- subset(seurat_obj, cells = keep_cells)
  
  #message about the before and after counts, for sanity.
  #cat(sprintf(
  #  "KDE trimming: %d cells before, %d cells after (%.2f%% removed)\n",
  #  n_before, n_after, 100 * (n_before - n_after) / n_before
  #))
  # Optionally, attach pass/fail vectors for record-keeping
  #attr(seurat_obj, "kde_pass_atac") <- pass_atac
  #attr(seurat_obj, "kde_pass_rna")  <- pass_rna
  
  return(seurat_obj)
}