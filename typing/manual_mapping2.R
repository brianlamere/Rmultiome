#manual mapping
# find these markers associated with HIV stress: NfL, MAP2, Synaptophysin, PSD-95, 8-OHdG

## Oligodendrocytes
DefaultAssay(labeled_obj) <- "RNA"
DimPlot(labeled_obj,label=T, raster=FALSE)
FeaturePlot(labeled_obj, features = c("MBP","MOBP","PLP1","MOG"), raster=FALSE)
DotPlot(object = labeled_obj, features = c("MBP","MOBP","PLP1","MOG"))

labeled_obj$celltypes[labeled_obj$celltypes == "0"] <- "Oligodendrocyte"

## Astrocytes
# consider adding:  YKL-40 (CHI3L1), S100B
DefaultAssay(labeled_obj) <- "RNA"
FeaturePlot(labeled_obj, features = c("GFAP","ALDH1L1","GLUL","AQP4","SLC1A2","SLC4A4"), raster=FALSE)
DotPlot(object = labeled_obj, features = c("GFAP","ALDH1L1","GLUL","AQP4","SLC1A2","SLC4A4"))

labeled_obj$celltypes[labeled_obj$celltypes=="2"] <- "Astrocyte"

## Microglia
# consider adding: CD163, sCD163, Iba-1, CD68, TREM2, GBP2, ACSL1, NAMPT, DPYD
DefaultAssay(labeled_obj) <- "RNA"
FeaturePlot(labeled_obj, features = c("CSF1R","APBB1IP","P2RY12"), raster=FALSE)
DotPlot(object = labeled_obj, features = c("CSF1R","APBB1IP","P2RY12"))

labeled_obj$celltypes[labeled_obj$celltypes=="3"] <- "Microglia"

## Macrophages
# consider adding:  CD163, sCD163, CD14, sCD14 
dp <- DotPlot(object = chosen_obj, features = c("CD163","MRC1","LYVE1"))

# Lymphocytes
# "CD96","CD3D","IL2RG"
DotPlot(object = labeled_obj, features = c("CD96","CD3D","IL2RG"))

##########################################
#focus on above cell types
#################################
## OPC
FeaturePlot(labeled_obj, features = c("VCAN","PDGFRA","PCDH15"), raster=FALSE)
DotPlot(object = labeled_obj, features = c("VCAN","PDGFRA","PCDH15"))
labeled_obj$celltypes <- as.character(labeled_obj$seurat_clusters)
labeled_obj$celltypes[labeled_obj$celltypes=="6"] <- "OPC"

## Interneurons
DefaultAssay(labeled_obj) <- "RNA"
FeaturePlot(labeled_obj,
            features = c("GAD1","GAD2","SLC2A1","LAMP5","PAX6","VIP","SST",
                         "PVALB","ADARB2"), raster=FALSE)
DotPlot(object = labeled_obj,
        features = c("GAD1","GAD2","SLC2A1","LAMP5","PAX6","VIP","SST","PVALB",
                     "ADARB2","CALB2"))

## Endothelial cells 
DefaultAssay(labeled_obj) <- "RNA"
FeaturePlot(labeled_obj, features = c("FLT1","CLDN5","KDR"), raster=FALSE)
DotPlot(object = labeled_obj, features = c("FLT1","CLDN5","KDR"))

labeled_obj$celltypes[labeled_obj$celltypes=="7"] <- "Endothelial"


###########
#looking for identity of cluster 9

markers_9 <- FindMarkers(
  labeled_obj,
  ident.1 = 9,                # Cluster number
  only.pos = TRUE,            # Only positive markers
  min.pct = 0.25,             # Expressed in at least 25% of cells in cluster
  logfc.threshold = 0.25      # Minimum log-fold change
)

# Write out full marker tables for reference
write.csv(markers_9[order(markers_9$p_val_adj), ], "cluster9_markers_pval.csv")
write.csv(markers_9[order(-markers_9$avg_log2FC), ], "cluster9_markers_log2FC.csv")

# Top 10 markers by significance
top_genes_pval <- rownames(head(markers_9[order(markers_9$p_val_adj), ], 10))

# Top 10 markers by log2FC
top_genes_log2FC <- rownames(head(markers_9[order(-markers_9$avg_log2FC), ], 10))

# Plots
FeaturePlot(labeled_obj, features = top_genes_pval, raster = FALSE)
FeaturePlot(labeled_obj, features = top_genes_log2FC, raster = FALSE)
DotPlot(labeled_obj, features = top_genes_pval)
DotPlot(labeled_obj, features = top_genes_log2FC)

FeaturePlot(tagged_obj, features = c("JUNB", "ARC", "FOS", "NR4A1", "EGR1",
                                     "NPAS4", "DUSP1", "AURKAIP1", "MRPL20",
                                     "SSU72", "GABRD", "TPRG1L"))
