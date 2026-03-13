#manual mapping
FeaturePlot(merged_data4, features = c("ST18", "MOBP", "TMEM144",
                                       "CTNNA3", "KCNH8")) +
  ggtitle("Feature Algorith = SLM, resolution = 0.4")
FeaturePlot(merged_data, features = c("CBLN2", "MEG3", "CUX2", "HS6ST3", "CA10"))
VlnPlot(merged_data, features = c("ST18", "MOBP", "TMEM144", "CTNNA3", "KCNH8"), group.by = "seurat_clusters")
VlnPlot(merged_data, features = c("ST18"), group.by = "seurat_clusters")
VlnPlot(merged_data, features = "ST18", group.by = "seurat_clusters", pt.size = 0)



# Extract log-normalized expression matrix from Seurat
expr_mat <- GetAssayData(merged_data, assay = "RNA", slot = "data")

# Convert to data frame
expr_mat_df <- as.data.frame(as.matrix(expr_mat))

# Write to CSV (genes = rows, cells = columns)
write.csv(expr_mat_df, file = "expression_matrix_for_sctype.csv", row.names = TRUE)

DefaultAssay(merged_data) <- "RNA"
merged_data <- FindVariableFeatures(merged_data)
top_genes <- VariableFeatures(merged_data)
expr_mat <- GetAssayData(merged_data, assay = "RNA", slot = "data")[top_genes, ]
expr_mat_df <- as.data.frame(as.matrix(expr_mat))
write.csv(expr_mat_df, file = "sctype_topgenes.csv", row.names = TRUE)