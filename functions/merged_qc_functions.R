
#' Check PC correlations with technical covariates
#' @param seurat_obj Merged Seurat object with PCA reduction
#' @param n_pcs Number of PCs to test (default 10)
#' @param reduction Name of reduction to use (default "pca")
#' @return ggplot object showing correlation heatmap
check_pc_technical_bias <- function(seurat_obj, n_pcs = 10, reduction = "pca") {
  require(ggplot2)
  require(reshape2)
  
  # Technical covariates to test
  covariates <- c("percent.mt", "nCount_RNA", "nFeature_RNA", 
                  "TSS.enrichment", "nCount_ATAC", "nFeature_ATAC")
  
  # Check which covariates exist
  available_covariates <- covariates[covariates %in% colnames(seurat_obj@meta.data)]
  
  # Get PC embeddings
  pcs <- Embeddings(seurat_obj, reduction = reduction)[, 1:n_pcs]
  
  # Compute correlations
  cor_matrix <- matrix(NA, nrow = n_pcs, ncol = length(available_covariates))
  rownames(cor_matrix) <- paste0("PC", 1:n_pcs)
  colnames(cor_matrix) <- available_covariates
  
  for (i in 1:n_pcs) {
    for (j in seq_along(available_covariates)) {
      cov_name <- available_covariates[j]
      cor_matrix[i, j] <- cor(pcs[, i], seurat_obj@meta.data[[cov_name]], 
                              use = "complete.obs")
    }
  }
  
  # Melt for ggplot
  cor_df <- melt(cor_matrix)
  colnames(cor_df) <- c("PC", "Covariate", "Correlation")
  
  # Create heatmap
  p1 <- ggplot(cor_df, aes(x = Covariate, y = PC, fill = Correlation)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                         midpoint = 0, limits = c(-1, 1)) +
    geom_text(aes(label = sprintf("%.2f", Correlation)), size = 3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("PC Correlation with Technical Covariates -", 
                       seurat_obj@project.name),
         x = "Technical Covariate", y = "Principal Component")
  
  # Create line plot (absolute correlations)
  cor_df$AbsCorrelation <- abs(cor_df$Correlation)
  p2 <- ggplot(cor_df, aes(x = PC, y = AbsCorrelation, 
                           color = Covariate, group = Covariate)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    theme_minimal() +
    labs(title = "Absolute Correlation by PC (Technical Bias Decay)",
         x = "Principal Component", y = "Absolute Correlation",
         color = "Covariate") +
    geom_hline(yintercept = 0.1, linetype = "dashed", color = "gray50") +
    annotate("text", x = n_pcs * 0.8, y = 0.12, 
             label = "Threshold (|r| = 0.1)", color = "gray50")
  
  # Print summary stats
  cat("\n=== PC1 Correlations (Technical Bias Check) ===\n")
  pc1_cors <- cor_matrix[1, ]
  print(sort(abs(pc1_cors), decreasing = TRUE))
  
  cat("\n=== PC2 Correlations (Should be lower) ===\n")
  pc2_cors <- cor_matrix[2, ]
  print(sort(abs(pc2_cors), decreasing = TRUE))
  
  cat("\n=== Maximum Absolute Correlation by PC ===\n")
  max_cors <- apply(abs(cor_matrix), 1, max)
  print(max_cors)
  
  # Return both plots
  list(heatmap = p1, lineplot = p2, cor_matrix = cor_matrix)
}
