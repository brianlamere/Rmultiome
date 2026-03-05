#Functions that are part of premerge processing.  Nothing here should be steps for a merged object

#we're using v86, though v114 is available it is from this month.
loadannotations <- function(ensdb = EnsDb.Hsapiens.v86) {
  annotation <- GetGRangesFromEnsDb(ensdb = ensdb)
  seqlevels(annotation) <- paste0('chr', seqlevels(annotation))
  return(annotation)
}

read_cellbender_rna_counts <- function(sample) {
  # Inputs: sample name (string)
  #  Uses settings: cb_datadir, cellbender_rna_h5filename
  #  Returns: dgCMatrix (genes × barcodes) from CellBender H5
  #  Behavior: Hard-fails if file missing or can't load scCustomize
}

splice_cellbender_rna_into_multiome <- function(rna_orig, atac_orig, cb_rna, sample) {
  # Inputs:
    # rna_orig: CellRanger RNA counts (matrix)
    # atac_orig: CellRanger ATAC counts (matrix)
    # cb_rna: CellBender RNA counts (matrix)
    # sample: sample name (for error messages)
  # Returns: named list:
    # list(
    # rna_counts = <subsetted cb_rna>,
    # atac_counts = <subsetted atac_orig>,
    # common_cells = <character vector of intersected barcodes>,
    # in_mat = <CellRanger RNA subsetted to common_cells>,  # for metrics
    # cb_mat2 = <CellBender RNA subsetted to common_cells>)  # for metrics
  # Behavior: Intersects barcodes, subsets both modalities; stops if no common cells
}

compute_cb_merge_metrics <- function(sample, in_mat, cb_mat2, atac_orig, cb_rna_full, common_cells) {
  #Inputs:
  #      sample: sample name
  #      in_mat: CellRanger RNA (subsetted to common cells)
  #      cb_mat2: CellBender RNA (subsetted to common cells)
  #      atac_orig: original ATAC (full, for barcode comparison)
  #      cb_rna_full: CellBender RNA (full, for barcode comparison)
  #      common_cells: intersected barcode vector
  #  Returns: 1-row data.frame with columns:
  #      sample, n_atac, n_cb, n_intersection, prop_atac_covered
  #      pearson, spearman, r2, slope, intercept
  #      weighted_removed, ratio_median, ratio_q05, ratio_q95, nnz_ratio

}

emit_cb_report <- function(metrics, mode = c("write", "display", "none"), sample) {
  #Inputs:
  #      metrics: 1-row data.frame from compute_cb_merge_metrics()
  #      mode: "write", "display", or "none"
  #      sample: sample name (for filename)
  #  Uses settings: tmpfiledir, cellbender_report_dir (determines output dir based on calling context)
  #  Behavior:
  #      "display": print to console + View() + write to tmpfiledir/cellbender_merge_reports/<sample>.csv
  #      "write": write to cellbender_report_dir/<sample>.csv
  #      "none": do nothing
}

#' Create the initial base seurat object
#'
#' @param samplename raw file from Read10X_h5
#' @return base Seurat5 object with percent.mt added
base_object <- function(samplename) {
  fullrna <- paste(rawdatadir, samplename, h5filename, sep = "/")
  fullatac <- paste(rawdatadir, samplename, atacfilename, sep = "/")
  counts <- Read10X_h5(filename = fullrna)
  rna_counts <- counts$`Gene Expression`
  atac_counts <- counts$Peaks
  
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

base_qc_object <- function(sample, EnsDbAnnos, save = FALSE) {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  base_path <- get_rds_path(sample, "base")
  
  base_obj <- base_object(sample)
  
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
                                                



  
