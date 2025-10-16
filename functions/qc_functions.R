#functions for QC of objects

#QC density plots against ATAC assay
QCDensity_ATAC <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  nCA_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nCount_ATAC',
                                    y = 'TSS.enrichment',
                                    log_x = TRUE,
                                    quantiles = TRUE
  ) +
    labs(
      title = paste("ATAC counts + TSS density plot for ",
                    seurat_obj@project.name,
                    " for debris and doublets")
    ) +
    stat_density_2d(
      color = "red",
      linewidth = 1,   # updated from size
      bins = 6,
      contour = TRUE
    )
  nFA_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nFeature_ATAC',
                                    y = 'TSS.enrichment',
                                    log_x = TRUE,
                                    quantiles = TRUE
  ) +
    labs(
    title = paste("ATAC features + TSS density plot for ",
                  seurat_obj@project.name,
                  " for debris and doublets")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 6, # Increase bins for smoother contours
      contour = TRUE
    )
  nss_TSS_density <- DensityScatter(seurat_obj,
                                    x = 'nucleosome_signal',
                                    y = 'TSS.enrichment',
                                    quantiles = TRUE
  ) +
    labs(
    title = paste("nucleosome signal vs TSS enrichment density plot for ",
                  seurat_obj@project.name,
                  " to find nucleosome-rich cells")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 6, # Increase bins for smoother contours
      contour = TRUE
    )
  print(nCA_TSS_density)
  print(nFA_TSS_density)
  print(nss_TSS_density)
}

QCDensity_RNA <- function(seurat_obj) {
  #DefaultAssay(seurat_obj) <- "RNA"
  pMT_nCR_density <- DensityScatter(seurat_obj,
                                    x = 'percent.mt',
                                    y = 'nCount_RNA',
                                    log_y = TRUE,
                                    quantiles = c(5, 10, 50, 90, 95)
  ) +
    labs(
      title = paste("percent.mt + RNA counts density plot for ",
                  seurat_obj@project.name,
                  " for cells with high mitochondrial content")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 10, # Increase bins for smoother contours
      contour = TRUE
    )
  pMT_nFR_density <- DensityScatter(seurat_obj,
                                    x = 'percent.mt',
                                    y = 'nFeature_RNA',
                                    log_y = TRUE,
                                    quantiles = c(5, 10, 50, 90, 95)
  ) +
    labs(
    title = paste("percent.mt + RNA features density plot for ",
                  seurat_obj@project.name,
                  " for cells with high mitochondrial content")) + 
    stat_density_2d(
      color = "red",
      linewidth = 1,
      bins = 10, # Increase bins for smoother contours
      contour = TRUE
    )
  nCR_nFR_density <- DensityScatter(seurat_obj,
                                    x = 'nCount_RNA',
                                    y = 'nFeature_RNA',
                                    log_x = TRUE,
                                    log_y = TRUE,
                                    quantiles = TRUE
  ) +
    ggtitle(paste("RNA count vs RNA feature density plot for ",
                  seurat_obj@project.name,
                  " for doublets and general QC"))
  print(nCR_nFR_density)
  print(pMT_nCR_density)
  print(pMT_nFR_density)
}

#QC Vln plots against ATAC assay
QCVlnA <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "ATAC"
  VPsn2 <- VlnPlot(
    object = seurat_obj,
    features = c("nCount_ATAC", "nCount_RNA", "nFeature_ATAC", "nFeature_RNA",
                 "TSS.enrichment", "nucleosome_signal", "percent.mt"),
    ncol = 4,
    pt.size = 0
  ) 
  print(VPsn2)
}

#QC Violin plots against RNA assay
QCVlnR <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  VP2sn2 <- VlnPlot(seurat_obj, features = "percent.mt") +
    ggtitle(paste("MT Density plot for", seurat_obj@project.name))
  print(VP2sn2)
}

QCpmt <- function(seurat_obj) {
  percent_mt_values <- seurat_obj@meta.data$percent.mt
  quantiles <- quantile(percent_mt_values, probs = c(0.05, 0.10, 0.90, 0.95), na.rm = TRUE)
  print(paste("The mean value of the MT percentages in this sample is",
              mean(percent_mt_values, na.rm = TRUE)))
  print(paste("The 5%, 10%, 90%, and 95% percentiles of percent.mt value are",
              paste(quantiles, collapse = ", ")))
}

plot_nCountRNA_percentMT <- function(seurat_obj) {
  df <- seurat_obj@meta.data
  ggplot(df, aes(x = nCount_RNA, y = percent.mt)) +
    geom_point(alpha = 0.3, color = "blue") +
    geom_density2d(color = "red") +
    labs(title = paste("nCount_RNA vs percent.mt for", seurat_obj@project.name),
         x = "nCount_RNA", y = "percent.mt") +
    theme_bw()
}

plotKDETrim <- function(seurat_obj, x_col, y_col, kde_percentile = 0.95,
                        retained_cells = NULL, main = NULL, ...) {
  library(MASS)
  df <- seurat_obj@meta.data
  
  x <- df[[x_col]]
  y <- df[[y_col]]
  
  kde <- kde2d(x, y, n = 100)
  
  # Find contour level for percentile
  dz <- sort(as.vector(kde$z), decreasing = TRUE)
  cumprob <- cumsum(dz) / sum(dz)
  level <- dz[which(cumprob >= kde_percentile)[1]]
  
  if (is.null(retained_cells)) {
    retained_cells <- rownames(df)
  }
  
  keep <- rownames(df) %in% retained_cells
  
  plot(x, y, col = ifelse(keep, "blue", "grey"), pch = 20, xlab = x_col,
       ylab = y_col, main = main, ...)
  contour(kde, levels = level, add = TRUE, col = "red", lwd = 2)
  legend("topright", legend = c("Retained", "Trimmed"),
         col = c("blue", "grey"), pch = 20)
}

plot_kde_filter_contours <- function(seurat_obj, kde_settings = NULL, sample_name = NULL, ...) {
  # Use project name if not specified
  if (is.null(sample_name)) sample_name <- seurat_obj@project.name
  
  # Load kde_settings if not provided
  if (is.null(kde_settings)) {
    if (!exists("kde_settings", envir = .GlobalEnv)) {
      stop("kde_settings not found in global environment, and not provided as argument.")
    }
    kde_settings <- get("kde_settings", envir = .GlobalEnv)
  }
  
  # Extract params for this sample
  params <- kde_settings[kde_settings$sample == sample_name, ]
  if (nrow(params) == 0) stop(sprintf("No KDE settings found for sample '%s'.", sample_name))
  atac_percentile <- as.numeric(params$atac_percentile)
  rna_percentile <- as.numeric(params$rna_percentile)
  
  df <- seurat_obj@meta.data
  
  # ATAC
  x_atac <- df$nCount_ATAC
  y_atac <- df$TSS.enrichment
  kde_atac <- kde2d(x_atac, y_atac, n = 100)
  level_atac <- get_perc_level(kde_atac, kdepercent = atac_percentile)
  
  plot(x_atac, y_atac, pch=20, cex=0.5,
       main=sprintf("ATAC KDE filtering at %.0f%%", 100 * atac_percentile),
       xlab="nCount_ATAC", ylab="TSS.enrichment", ...)
  contour(kde_atac, levels=level_atac, add=TRUE, col="red", lwd=2, drawlabels=FALSE)
  legend("topright", legend=sprintf("Contour at %.2f percentile", atac_percentile), col="red", lwd=2)
  
  # RNA
  x_rna <- df$percent.mt
  y_rna <- df$nCount_RNA

  kde_rna <- kde2d(x_rna, y_rna, n = 100)
  level_rna <- get_perc_level(kde_rna, kdepercent = rna_percentile)
  
  plot(x_rna, y_rna, pch=20, cex=0.5,
       main=sprintf("RNA KDE filtering at %.0f%%", 100 * rna_percentile),
       xlab="percent.MT", ylab="nCount_RNA", ...)
  contour(kde_rna, levels=level_rna, add=TRUE, col="blue", lwd=2, drawlabels=FALSE)
  legend("topright", legend=sprintf("Contour at %.2f percentile", rna_percentile), col="blue", lwd=2)
}

plot_kde_filter_combine_compare_atac <- function(seurat_obj, kde_settings = NULL, sample_name = NULL, ...) {
  if (is.null(sample_name)) sample_name <- seurat_obj@project.name
  if (is.null(kde_settings)) {
    if (!exists("kde_settings", envir = .GlobalEnv)) {
      stop("kde_settings not found in global environment, and not provided as argument.")
    }
    kde_settings <- get("kde_settings", envir = .GlobalEnv)
  }
  params <- kde_settings[kde_settings$sample == sample_name, ]
  if (nrow(params) == 0) stop(sprintf("No KDE settings found for sample '%s'.", sample_name))
  atac_percentile <- as.numeric(params$atac_percentile)
  rna_percentile <- as.numeric(params$rna_percentile)
  
  df <- seurat_obj@meta.data
  x_atac <- df$nCount_ATAC
  y_atac <- df$TSS.enrichment
  kde_atac <- kde2d(x_atac, y_atac, n = 100)
  level_atac <- get_perc_level(kde_atac, kdepercent = atac_percentile)
  dens_atac <- get_density_values(x_atac, y_atac, kde_atac)
  pass_atac <- dens_atac >= level_atac
  
  x_rna <- df$nCount_RNA
  y_rna <- df$percent.mt
  kde_rna <- kde2d(x_rna, y_rna, n = 100)
  level_rna <- get_perc_level(kde_rna, kdepercent = rna_percentile)
  dens_rna <- get_density_values(x_rna, y_rna, kde_rna)
  pass_rna <- dens_rna >= level_rna
  
  union_idx <- pass_atac | pass_rna
  inter_idx <- pass_atac & pass_rna
  
  plot(x_atac, y_atac, pch=20, cex=0.5, 
       main=sprintf("ATAC KDE filtering at %.0f%%\nUnion (red) vs Intersection (blue)", 100 * atac_percentile),
       xlab="nCount_ATAC", ylab="TSS.enrichment", col="grey80", ...)
  if (sum(union_idx) > 10) {
    kde_union <- kde2d(x_atac[union_idx], y_atac[union_idx], n = 100)
    contour(kde_union, levels=level_atac, add=TRUE, col="red", lwd=2, drawlabels=FALSE)
  }
  if (sum(inter_idx) > 10) {
    kde_inter <- kde2d(x_atac[inter_idx], y_atac[inter_idx], n = 100)
    contour(kde_inter, levels=level_atac, add=TRUE, col="blue", lwd=2, drawlabels=FALSE)
  }
  legend("topright", legend=c("Union", "Intersection"),
         col=c("red", "blue"), lwd=2, bty="n")
}

plot_kde_filter_combine_compare_rna <- function(seurat_obj, kde_settings = NULL, sample_name = NULL, ...) {
  if (is.null(sample_name)) sample_name <- seurat_obj@project.name
  if (is.null(kde_settings)) {
    if (!exists("kde_settings", envir = .GlobalEnv)) {
      stop("kde_settings not found in global environment, and not provided as argument.")
    }
    kde_settings <- get("kde_settings", envir = .GlobalEnv)
  }
  params <- kde_settings[kde_settings$sample == sample_name, ]
  if (nrow(params) == 0) stop(sprintf("No KDE settings found for sample '%s'.", sample_name))
  atac_percentile <- as.numeric(params$atac_percentile)
  rna_percentile <- as.numeric(params$rna_percentile)
  
  df <- seurat_obj@meta.data
  x_atac <- df$nCount_ATAC
  y_atac <- df$TSS.enrichment
  kde_atac <- kde2d(x_atac, y_atac, n = 100)
  level_atac <- get_perc_level(kde_atac, kdepercent = atac_percentile)
  dens_atac <- get_density_values(x_atac, y_atac, kde_atac)
  pass_atac <- dens_atac >= level_atac
  
  x_rna <- df$nCount_RNA
  y_rna <- df$percent.mt
  kde_rna <- kde2d(x_rna, y_rna, n = 100)
  level_rna <- get_perc_level(kde_rna, kdepercent = rna_percentile)
  dens_rna <- get_density_values(x_rna, y_rna, kde_rna)
  pass_rna <- dens_rna >= level_rna
  
  union_idx <- pass_atac | pass_rna
  inter_idx <- pass_atac & pass_rna
  
  plot(x_rna, y_rna, pch=20, cex=0.5, 
       main=sprintf("RNA KDE filtering at %.0f%%\nUnion (red) vs Intersection (blue)", 100 * rna_percentile),
       xlab="nCount_RNA", ylab="percent.mt", col="grey80", ...)
  if (sum(union_idx) > 10) {
    kde_union <- kde2d(x_rna[union_idx], y_rna[union_idx], n = 100)
    contour(kde_union, levels=level_rna, add=TRUE, col="red", lwd=2, drawlabels=FALSE)
  }
  if (sum(inter_idx) > 10) {
    kde_inter <- kde2d(x_rna[inter_idx], y_rna[inter_idx], n = 100)
    contour(kde_inter, levels=level_rna, add=TRUE, col="blue", lwd=2, drawlabels=FALSE)
  }
  legend("topright", legend=c("Union", "Intersection"),
         col=c("red", "blue"), lwd=2, bty="n")
}