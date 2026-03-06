#Functions that are part of premerge processing.  Nothing here should be steps for a merged object

#we're using v86, though v114 is available it is from this month.
loadannotations <- function(ensdb = EnsDb.Hsapiens.v86) {
  annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  return(annotation)
}

#' Read CellBender-corrected RNA counts from H5 file
#'
#' @param sample Sample identifier
#' @return dgCMatrix of CellBender RNA counts (genes × barcodes)
read_cellbender_rna_counts <- function(sample) {
  cb_h5_path <- file.path(cb_datadir, sample, cellbender_rna_h5filename)
  
  if (!file.exists(cb_h5_path)) {
    stop(sprintf("CellBender H5 file not found for sample '%s': %s", sample, cb_h5_path))
  }
  
  # Check if scCustomize is available
  if (!requireNamespace("scCustomize", quietly = TRUE)) {
    stop("Package 'scCustomize' is required for reading CellBender output. Install with: install.packages('scCustomize')")
  }
  
  cat(sprintf("Reading CellBender RNA counts from: %s\n", cb_h5_path))
  cb_rna <- scCustomize::Read_CellBender_h5_Mat(file_name = cb_h5_path)
  
  cat(sprintf("Loaded CellBender RNA: %d genes, %d cells\n", nrow(cb_rna), ncol(cb_rna)))
  
  return(cb_rna)
}

#' Splice CellBender RNA into multiome data by intersecting barcodes
#'
#' @param rna_orig CellRanger RNA counts matrix
#' @param atac_orig CellRanger ATAC counts matrix
#' @param cb_rna CellBender RNA counts matrix
#' @param sample Sample name (for error messages)
#' @return Named list with rna_counts, atac_counts, common_cells, in_mat, cb_mat2
splice_cellbender_rna_into_multiome <- function(rna_orig, atac_orig, cb_rna, sample) {
  # Get barcode sets
  atac_cells <- colnames(atac_orig)
  cb_cells <- colnames(cb_rna)
  
  # Intersect barcodes
  common_cells <- intersect(atac_cells, cb_cells)
  
  if (length(common_cells) == 0) {
    stop(sprintf("Sample '%s': No common barcodes between ATAC (%d cells) and CellBender RNA (%d cells).", 
                 sample, length(atac_cells), length(cb_cells)))
  }
  
  # Report coverage
  prop_atac_covered <- length(common_cells) / length(atac_cells)
  cat(sprintf("Barcode intersection: %d common cells (%.1f%% of ATAC barcodes)\n", 
              length(common_cells), prop_atac_covered * 100))
  
  if (prop_atac_covered < 0.80) {
    warning(sprintf("Sample '%s': Only %.1f%% of ATAC barcodes found in CellBender output (< 80%%). Consider re-running CellBender.", 
                    sample, prop_atac_covered * 100))
  }
  
  # Subset both modalities to common cells
  rna_counts <- cb_rna[, common_cells, drop = FALSE]
  atac_counts <- atac_orig[, common_cells, drop = FALSE]
  
  # Also subset the original CellRanger RNA to common cells (for metrics comparison)
  in_mat <- rna_orig[, common_cells, drop = FALSE]
  cb_mat2 <- cb_rna[, common_cells, drop = FALSE]
  
  return(list(
    rna_counts = rna_counts,
    atac_counts = atac_counts,
    common_cells = common_cells,
    in_mat = in_mat,
    cb_mat2 = cb_mat2
  ))
}

#' Compute CellBender merge metrics
#'
#' @param sample Sample name
#' @param in_mat CellRanger RNA (subsetted to common cells)
#' @param cb_mat2 CellBender RNA (subsetted to common cells)
#' @param atac_orig Original ATAC counts (full)
#' @param cb_rna_full CellBender RNA counts (full)
#' @param common_cells Vector of intersected barcodes
#' @return 1-row data.frame with merge quality metrics
compute_cb_merge_metrics <- function(sample, in_mat, cb_mat2, atac_orig, cb_rna_full, common_cells) {
  # Barcode overlap metrics
  n_atac <- ncol(atac_orig)
  n_cb <- ncol(cb_rna_full)
  n_intersection <- length(common_cells)
  prop_atac_covered <- n_intersection / n_atac
  
  # Per-cell UMI totals (for correlation)
  cellranger_umis <- Matrix::colSums(in_mat)
  cellbender_umis <- Matrix::colSums(cb_mat2)
  
  # Correlation metrics
  pearson <- cor(cellranger_umis, cellbender_umis, method = "pearson")
  spearman <- cor(cellranger_umis, cellbender_umis, method = "spearman")
  
  # Linear regression (CellBender ~ CellRanger)
  lm_fit <- lm(cellbender_umis ~ cellranger_umis)
  r2 <- summary(lm_fit)$r.squared
  slope <- coef(lm_fit)[2]
  intercept <- coef(lm_fit)[1]
  
  # Removal metrics (weighted by original counts)
  removed_counts <- cellranger_umis - cellbender_umis
  weighted_removed <- sum(removed_counts * cellranger_umis) / sum(cellranger_umis^2)
  
  # Ratio distribution (CB / CellRanger per cell)
  ratios <- cellbender_umis / cellranger_umis
  ratio_median <- median(ratios)
  ratio_q05 <- quantile(ratios, 0.05)
  ratio_q95 <- quantile(ratios, 0.95)
  
  # Non-zero ratio (sparsity comparison)
  nnz_cellranger <- Matrix::nnzero(in_mat)
  nnz_cellbender <- Matrix::nnzero(cb_mat2)
  nnz_ratio <- nnz_cellbender / nnz_cellranger
  
  # Assemble metrics into 1-row data.frame
  metrics <- data.frame(
    sample = sample,
    n_atac = n_atac,
    n_cb = n_cb,
    n_intersection = n_intersection,
    prop_atac_covered = prop_atac_covered,
    pearson = pearson,
    spearman = spearman,
    r2 = r2,
    slope = slope,
    intercept = intercept,
    weighted_removed = weighted_removed,
    ratio_median = ratio_median,
    ratio_q05 = ratio_q05,
    ratio_q95 = ratio_q95,
    nnz_ratio = nnz_ratio,
    stringsAsFactors = FALSE
  )
  
  return(metrics)
}

#' Emit CellBender merge report (display and/or write)
#'
#' @param metrics 1-row data.frame from compute_cb_merge_metrics()
#' @param mode "write", "display", or "none"
#' @param sample Sample name (for filename)
emit_cb_report <- function(metrics, mode = c("write", "display", "none"), sample) {
  mode <- match.arg(mode)
  
  if (mode == "none") {
    return(invisible(NULL))
  }
  
  # Determine output directory based on mode
  if (mode == "display") {
    # QC mode: write to tmp and display
    out_dir <- file.path(tmpfiledir, "cellbender_merge_reports")
    
    # Ensure directory exists
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    # Write CSV
    out_path <- file.path(out_dir, paste0(sample, ".csv"))
    write.csv(metrics, out_path, row.names = FALSE)
    
    # Display to console
    cat("\n=== CellBender Merge Report ===\n")
    cat(sprintf("Sample: %s\n", sample))
    cat(sprintf("ATAC cells: %d | CellBender cells: %d | Intersection: %d (%.1f%%)\n",
                metrics$n_atac, metrics$n_cb, metrics$n_intersection, 
                metrics$prop_atac_covered * 100))
    cat(sprintf("Correlation - Pearson: %.3f | Spearman: %.3f | R²: %.3f\n",
                metrics$pearson, metrics$spearman, metrics$r2))
    cat(sprintf("Regression - Slope: %.3f | Intercept: %.1f\n",
                metrics$slope, metrics$intercept))
    cat(sprintf("Removal - Weighted: %.3f | Ratio median: %.3f [Q05: %.3f, Q95: %.3f]\n",
                metrics$weighted_removed, metrics$ratio_median, 
                metrics$ratio_q05, metrics$ratio_q95))
    cat(sprintf("Sparsity - NNZ ratio: %.3f\n", metrics$nnz_ratio))
    cat(sprintf("Report saved to: %s\n", out_path))
    cat("================================\n\n")
    
    ## View in RStudio (safe in interactive sessions)
    #if (interactive()) {
    #  utils::View(metrics, title = paste("CellBender Merge:", sample))
    #}
    
  } else if (mode == "write") {
    # Pipeline mode: write to export directory (no display)
    out_dir <- cellbender_report_dir
    
    # Ensure directory exists
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    
    out_path <- file.path(out_dir, paste0(sample, ".csv"))
    write.csv(metrics, out_path, row.names = FALSE)
    cat(sprintf("CellBender merge report written to: %s\n", out_path))
  }
  
  return(invisible(NULL))
}

#' Create the initial base seurat object
#'
#' @param samplename Sample identifier
#' @param cb_report CellBender merge report mode: "write", "display", or "none"
#' @return base Seurat5 object with percent.mt added
base_object <- function(samplename, cb_report = c("write", "display", "none")) {
  cb_report <- match.arg(cb_report)
  
  # Always read the multiome H5 (needed for ATAC, and for RNA if not using CellBender)
  fullrna <- paste(rawdatadir, samplename, h5filename, sep = "/")
  fullatac <- paste(rawdatadir, samplename, atacfilename, sep = "/")
  
  counts <- Read10X_h5(filename = fullrna)
  rna_counts_orig <- counts$`Gene Expression`
  atac_counts_orig <- counts$Peaks
  
  # Determine final RNA and ATAC counts based on use_cellbender setting
  if (use_cellbender) {
    cat("CellBender mode enabled: loading corrected RNA counts...\n")
    
    # Read CellBender RNA
    cb_rna <- read_cellbender_rna_counts(samplename)
    
    # Splice CellBender RNA into multiome (intersect barcodes, subset both modalities)
    spliced <- splice_cellbender_rna_into_multiome(
      rna_orig = rna_counts_orig,
      atac_orig = atac_counts_orig,
      cb_rna = cb_rna,
      sample = samplename
    )
    
    # Compute merge metrics
    metrics <- compute_cb_merge_metrics(
      sample = samplename,
      in_mat = spliced$in_mat,
      cb_mat2 = spliced$cb_mat2,
      atac_orig = atac_counts_orig,
      cb_rna_full = cb_rna,
      common_cells = spliced$common_cells
    )
    
    # Emit report (display/write/none)
    emit_cb_report(metrics, mode = cb_report, sample = samplename)
    
    # Use spliced (intersected) counts
    rna_counts <- spliced$rna_counts
    atac_counts <- spliced$atac_counts
    
    cat(sprintf("Using CellBender RNA counts: %d cells, %d genes\n", 
                ncol(rna_counts), nrow(rna_counts)))
    cat(sprintf("ATAC counts subsetted to common cells: %d cells, %d peaks\n", 
                ncol(atac_counts), nrow(atac_counts)))
    
  } else {
    # Use original CellRanger counts
    rna_counts <- rna_counts_orig
    atac_counts <- atac_counts_orig
  }
  
  # Create Seurat object (same logic whether using CellBender or not)
  print("Creating the RNA assay for the Seurat object...")
  baseSeuratObj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = samplename
  )
  
  print("Adding the ATAC assay for the Seurat object...")
  baseSeuratObj[["ATAC"]] <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fullatac,
    annotation = EnsDbAnnos
  )
  
  print("Calculating a slot for percent.mt for downstream QC")
  DefaultAssay(baseSeuratObj) <- "RNA"
  baseSeuratObj[["percent.mt"]] <- PercentageFeatureSet(baseSeuratObj, pattern = "^MT-")
  
  return(baseSeuratObj)
}

base_qc_object <- function(sample, EnsDbAnnos, save = FALSE, cb_report = "display") {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  base_path <- get_rds_path(sample, "base")
  
  # Create base object with CellBender reporting
  base_obj <- base_object(sample, cb_report = cb_report)
  
  cat("Adding chromosome mapping information to ATAC assay.\n")
  base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
  
  DefaultAssay(base_obj) <- "ATAC"
  
  cat("Calculating Nucleosome Signal...\n")
  base_obj <- NucleosomeSignal(base_obj)
  
  cat("Calculating TSS Enrichment...\n")
  base_obj <- TSSEnrichment(base_obj)
  
  if (save) {
    saveRDS(base_obj, base_path)
    cat(sprintf("Saved base object to: %s\n", base_path))
  }
  
  return(base_obj)
}

#' Annotate RNA features by EnsDb annotation and ATAC peaks by parsing peak names (hyphen format).
#'
#' @param seurat_obj Seurat object.
#' @param rna_annos GRanges object of gene annotation (from loadannotations/EnsDb).
#' @param warn_threshold Warn if fewer features mapped than this (RNA).
#' @return Seurat object with chromosome info in @misc$feature.info for both RNA and ATAC.
chromosome_mapping <- function(seurat_obj, rna_annos, warn_threshold = 16000) {
  # RNA: Map gene symbols to chromosomes using annotation
  if (!is.null(seurat_obj[["RNA"]])) {
    feature_names <- rownames(seurat_obj[["RNA"]])
    anno_symbols <- mcols(rna_annos)$gene_name
    #anno_ensembl <- mcols(rna_annos)$gene_id # or $gene_id depending on EnsDb version
    anno_chroms <- as.character(seqnames(rna_annos))
    
    chrom_vec <- anno_chroms[match(feature_names, anno_symbols)]
    names(chrom_vec) <- feature_names
    
    n_mapped <- sum(!is.na(chrom_vec))
    n_total <- length(chrom_vec)
    
    if (n_mapped < warn_threshold) {
      warning(sprintf("chromosome_mapping (RNA): Only %d/%d features mapped to a chromosome (expected 18,000-24,000).", n_mapped, n_total))
    } else {
      message(sprintf("chromosome_mapping (RNA): %d/%d features mapped to a chromosome.", n_mapped, n_total))
    }
    
    # Ensure feature.info is a data.frame with rownames = feature names
    feature_info <- seurat_obj[["RNA"]]@misc$feature.info
    if (is.null(feature_info) || !is.data.frame(feature_info) || nrow(feature_info) == 0) {
      feature_info <- data.frame(row.names = feature_names)
    }
    feature_info <- feature_info[feature_names, , drop=FALSE]
    feature_info$chromosome <- chrom_vec
    seurat_obj[["RNA"]]@misc$feature.info <- feature_info
  }
  
  # ATAC: Parse chromosome from peak names (chrN-...)
  if (!is.null(seurat_obj[["ATAC"]])) {
    peak_names <- rownames(seurat_obj[["ATAC"]])
    chrom_vec <- sub("-.*", "", peak_names)
    names(chrom_vec) <- peak_names
    
    feature_info <- seurat_obj[["ATAC"]]@misc$feature.info
    if (is.null(feature_info) || !is.data.frame(feature_info) || nrow(feature_info) == 0) {
      feature_info <- data.frame(row.names = peak_names)
    }
    feature_info <- feature_info[peak_names, , drop=FALSE]
    feature_info$chromosome <- chrom_vec
    seurat_obj[["ATAC"]]@misc$feature.info <- feature_info
    message(sprintf("chromosome_mapping (ATAC): Chromosome parsed for %d peaks.", length(peak_names)))
  }
  
  return(seurat_obj)
}

merge_sample_objects <- function(samplelist, suffix = "pipeline1", project_name = "opioid", path_fun = get_rds_path) {
  # Get file paths for all samples
  file_paths <- sapply(samplelist, function(sample) path_fun(sample, suffix))
  
  # Check file existence before proceeding
  missing_files <- file_paths[!file.exists(file_paths)]
  if (length(missing_files) > 0) {
    stop(sprintf("Missing RDS files for samples: %s", paste(missing_files, collapse=", ")))
  }
  
  # Read sample objects
  sample_objs <- lapply(file_paths, readRDS)
  
  # Merge sample objects
  merged_seurat <- merge(x = sample_objs[[1]], y = sample_objs[-1])
  
  # Optionally store input provenance for traceability
  merged_seurat@misc$input_provenance <- lapply(sample_objs, function(obj) obj@misc$provenance)
  
  # Update merged object provenance
  merged_seurat <- update_provenance(merged_seurat, "merged_object")
  
  # Save merged object
  #why are we doing this inside?  everything else did it in pipeline
  #saveRDS(merged_seurat, path_fun(project_name, "merge"))
  
  invisible(merged_seurat)
}
                                                



  
