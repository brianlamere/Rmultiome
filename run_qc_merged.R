source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#load object created at end of run_pipeline1.R
merged_obj <- readRDS(file.path(rdsdir,"merged_preharmony.Rds"))

# ============================================================================
# PHASE 1: Clustering Parameter Selection
# ============================================================================

# === STEP 1: Check for technical bias in PCs ===
pc_check <- check_pc_technical_bias(merged_obj, n_pcs = 30)

maybe_new_device(width = 10, height = 8)
print(pc_check$heatmap)

maybe_new_device(width = 12, height = 6)
print(pc_check$lineplot)

# Decision: Exclude PC1 based on high technical correlation

# === STEP 2: Decide harmony settings ===
# Based on PC bias check, decide which PCs to use

harmony_settings <- list(
  dims_use = 2:50,
  random_seed = 1984,
  max_iter = 50,
  project_dim = FALSE
)

saveRDS(harmony_settings, harmony_settings_file)

# === STEP 3: Run Harmony for parameter sweep ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = harmony_settings$random_seed,
  harmony_max_iter = harmony_settings$max_iter,
  harmony_project.dim = harmony_settings$project_dim,
  harmony_dims = harmony_settings$dims_use
)

saveRDS(harmony_obj, file.path(rdsdir, "harmonized.rds"))
#next line is only present due to debugging the code and restarting here
#harmony_obj <- readRDS(file.path(rdsdir, "harmonized.rds"))

# === STEP 4: Find elbow (informative but not prescriptive) ===
# Note: Elbow point is informative but not necessarily optimal for clustering.
maybe_new_device(width = 10, height = 6)
print(findElbow(harmony_obj))

#BELOW IS VERY INTERESTING, DETERMINE WHAT TOP CELLS IN THE TRANSITIONING GROUPS ARE
# === STEP 5: Run parameter sweep with metrics ===
sweep_results <- run_parameter_sweep_plots(
  seurat_obj = harmony_obj,
  dims_range = list(c(1:43),c(1:44),c(1:45)),
  knn_values = c(40),
  res_values = c(0.16,0.17,0.18),
  alg = 3,                   # SLM algorithm (required parameter)
  cluster_seed = 1984       # Reproducibility (required parameter)
)

# Save full sweep results table
write.csv(sweep_results, file.path(tmpfiledir, "param_sweep_results.csv"),
         row.names = FALSE)

# Print results sorted by cluster count
cat("Results sorted by cluster count:\n")
print(sweep_results[order(sweep_results$n_clusters), ])

# === STEP 6: USER INPUT - Set chosen parameters ===
# After reviewing plots, set these to your chosen values:
chosen_dims_min <- 1      # CHANGE THIS
chosen_dims_max <- 44     # CHANGE THIS
chosen_knn <- 40          # CHANGE THIS
chosen_resolution <- 0.16 # CHANGE THIS

# === STEP 7: Save cluster settings ===
cluster_settings <- data.frame(
  dims_min = chosen_dims_min,
  dims_max = chosen_dims_max,
  knn = chosen_knn,
  resolution = chosen_resolution,
  algorithm = 3,
  random_seed = 1984
)

saveRDS(cluster_settings, cluster_settings_file)

# === STEP 8: Create final clustered object from settings ===
cat("\n=== STEP 8: Applying cluster settings ===\n")

# Read the settings you just saved
dims <- cluster_settings$dims_min:cluster_settings$dims_max

# Apply FMMN
chosen_obj <- FMMN_task(harmony_obj,
                       knn = cluster_settings$knn,
                       dims = dims)

# Apply clustering
chosen_obj <- cluster_data(
  chosen_obj,
  alg = cluster_settings$algorithm,
  res = cluster_settings$resolution,
  cluster_seed = cluster_settings$random_seed,
  singleton_handling = "keep",
  run_umap = TRUE
)

# Quick summary
n_clusters <- length(unique(Idents(chosen_obj))) -
              (if("singleton" %in% levels(Idents(chosen_obj))) 1 else 0)
n_singletons <- sum(Idents(chosen_obj) == "singleton")

cat(sprintf("Result: %d clusters, %d singletons (%.2f%%)\n",
           n_clusters, n_singletons,
           100 * n_singletons / ncol(chosen_obj)))

# Quick viz
maybe_new_device(width = 10, height = 8)
print(DimPlot(chosen_obj, reduction = "wnn.umap", label = TRUE, raster = FALSE))

# Save
saveRDS(chosen_obj, file.path(rdsdir, "chosen_clustered_obj.rds"))
#adding as a checkpoint during code development
# chosen_obj <- readRDS(file.path(rdsdir, "chosen_clustered_obj.rds"))

# ============================================================================
# PHASE 2: Cell Type Annotation
# ============================================================================

# Note: refer to lymphocyte_investigation.R for exclusion reason

# === STEP 9: Define and save marker panels ===
# markers from https://pubmed.ncbi.nlm.nih.gov/40204703/
# Zillich, Cell type-specific multiomics analysis of cocaine use disorder
# in the human caudate nucleus
cocaine_markers_list <- split(cocaine_paper_markers$gene,
                             cocaine_paper_markers$celltype)

# === STEP 10: Finding cluster markers ===

DefaultAssay(chosen_obj) <- "RNA"

# CRITICAL: Join layers before FindAllMarkers (Seurat v5)
chosen_obj <- JoinLayers(chosen_obj)

# Exclude singletons from marker finding
clusters_to_analyze <- setdiff(levels(Idents(chosen_obj)), "singleton")

cat("Running FindAllMarkers...\n")
all_markers <- FindAllMarkers(
  chosen_obj,
  idents = clusters_to_analyze,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

# Save all markers
write.csv(all_markers, 
         file.path(tmpfiledir, "all_cluster_markers.csv"),
         row.names = FALSE)

cat(sprintf("Found %d marker genes across %d clusters\n",
           nrow(all_markers), length(unique(all_markers$cluster))))

# checkpoint if restarting
# all_markers <- read.csv(file.path(tmpfiledir, "all_cluster_markers.csv"))

# === STEP 11: save top markers per cluster ===
# useful for cell typing

top_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top_markers,
         file.path(tmpfiledir, "top10_markers_per_cluster.csv"),
         row.names = FALSE)

# === STEP 12:  programmatic cell type assessment ===

results <- identify_all_celltypes(
  all_markers,
  cocaine_markers_list,
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# Save results
write.csv(results$assignments,
         file.path(tmpfiledir, "celltype_assignments.csv"),
         row.names = FALSE)

# Show summary
cat("\n=== Assigned clusters ===\n")
print(results$assignments %>%
      filter(!is.na(cluster)) %>%
      select(cluster, celltype, confidence, n_markers, score))

# Check for unassigned
all_clusters <- setdiff(levels(Idents(chosen_obj)), "singleton")
assigned <- results$assignments$cluster[!is.na(results$assignments$cluster)]
unassigned <- setdiff(all_clusters, assigned)

if (length(unassigned) > 0) {
  cat(sprintf("\n%d unassigned clusters: %s\n",
             length(unassigned),
             paste(unassigned, collapse = ", ")))

# === STEP 12: Manual cell type assignment ===

# USER INPUT: Edit this mapping based on marker genes
# This is an EXAMPLE - adjust based on your actual markers
celltype_mapping <- data.frame(
  cluster = 0:18,
  celltype = c(
    "Oligodendrocytes",           # 0 nature assignment
    "Excitatory_Neurons",         # 1
    "Astrocytes",                 # 2 nature assignment
    "Microglia",                  # 3 nature assignment
    "GABAergic_SST",              # 4
    "OPCs",                       # 5 nature assignment
    "Stressed_Dying_Cells",       # 6 REMOVE
    "GABAergic_VIP",              # 7
    "Excitatory_Neurons",         # 8
    "Excitatory_Neurons",         # 9
    "Excitatory_Neurons",         # 10
    "Excitatory_Neurons",         # 11
    "Pericytes",                  # 12
    "Endothelial",                # 13 nature assignment
    "GABAergic_LAMP5",            # 14
    "Excitatory_Neurons",         # 15
    "Excitatory_Neurons",         # 16
    "GABAergic_PVALB",            # 17 nature assignment
    "Singleton"			  # singletons - remove
  ),
  markers_used = c(
    "MBP, MOBP, PLP1",
    "SLC17A7, CAMK2A, SATB2, TBR1, NRGN",
    "GFAP, ALDH1L1, GLUL, AQP4, SLC1A2, SLC4A4",
    "CSF1R, APBB1IP, P2RY12",
    "GAD1, GAD2, SST, TAC1, ARX",
    "VCAN, PDGFRA, PCDH15",
    "VIM, GFAP, NEFM, NEFL, NRGN",
    "GAD1, GAD2, VIP, CALB2, TAC3",
    "SLC17A7, CAMK2A, SATB2, TBR1, FOXP2",
    "SLC17A7, CAMK2A, SATB2, TBR1, FOXP2",
    "SLC17A7, CAMK2A, SATB2, TBR1, FOXP2",
    "SLC17A7, CAMK2A, SATB2, BCL11B, FEZF2",
    "PDGFRB, RGS5, NOTCH3, CARMN",
    "FLT1, CLDN5",
    "GAD1, GAD2, LAMP5, KIT",
    "SLC17A7, CAMK2A, SATB2, TBR1",
    "SLC17A7, CAMK2A, SATB2, BCL11B, FEZF2, PCP4",
    "GAD1, GAD2, PVALB, PTHLH",
    "QC_flag"
  ),
  confidence = c(
    "high", "medium", "high", "high",
    "high", "high", "low", "high",
    "medium", "medium", "medium", "medium",
    "high", "high", "high",
    "medium", "medium", "high","high"
  ),
  stringsAsFactors = FALSE
)

#celltype_mapping$celltype[celltype_mapping$cluster == 3] <- "Microglia"
#
#celltype_mapping$markers_used[celltype_mapping$cluster == 3] <-
#  "CSF1R, APBB1IP, P2RY12; elevated inflammatory markers (ITGAX, TLR2, SLC11A1, C3)"
#
#celltype_mapping$notes[celltype_mapping$cluster == 3] <-
#  "CNS-resident microglia showing gradient of activation states. Elevated expression of inflammatory markers (TLR2, ITGAX, SLC11A1, C3) in subset of cells indicates neuroinflammatory response to HIV/opioid exposure. Continuous activation spectrum rather than discrete subtypes."
#
#celltype_mapping$confidence[celltype_mapping$cluster == 3] <- "high"

# a reminder of the DimPlot:
# DimPlot(chosen_obj, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, raster = FALSE)

# === STEP 13: Save cell type mapping ===
saveRDS(celltype_mapping, celltype_settings_file)

# === STEP 14: Preview cell type assignment ===

# Apply cell types to preview
#cluster_ids <- as.numeric(as.character(Idents(chosen_obj)))
#chosen_obj$celltypes <- celltype_mapping$celltype[match(cluster_ids, 
#                                                        celltype_mapping$cluster)]
# Add celltype to metadata BEFORE filtering
chosen_obj$celltype <- celltype_mapping$celltype[match(chosen_obj$seurat_clusters,
                                                        celltype_mapping$cluster)]

# Report before filtering
cat("\n=== BEFORE FILTERING ===\n")
cat(sprintf("Total cells: %d\n", ncol(chosen_obj)))
table(chosen_obj$celltype) %>% print()

# Filter out cluster 6 and singleton
cells_to_keep <- !chosen_obj$celltype %in% c("Stressed_Dying_Cells", "Singleton")

chosen_obj_filtered <- subset(chosen_obj, cells = colnames(chosen_obj)[cells_to_keep])

# Report after filtering
cat("\n=== AFTER FILTERING ===\n")
cat(sprintf("Total cells: %d\n", ncol(chosen_obj_filtered)))
cat(sprintf("Cells removed: %d (%.1f%%)\n",
           sum(!cells_to_keep),
           100 * sum(!cells_to_keep) / ncol(chosen_obj)))

table(chosen_obj_filtered$celltype) %>% print()

# Update the main object name
chosen_obj <- chosen_obj_filtered
rm(chosen_obj_filtered)

# ============================================================================
# VISUALIZATION: CELL TYPE DIMPLOT
# ============================================================================

# Define nice colors for cell types
celltype_colors <- c(
  "Oligodendrocytes" = "#E69F00",        # Orange
  "OPCs" = "#F0E442",                    # Yellow
  "Astrocytes" = "#56B4E9",              # Light blue
  "Microglia" = "#CC79A7",               # Pink
  "Excitatory_Neurons" = "#009E73",      # Green
  "GABAergic_SST" = "#D55E00",           # Red-orange
  "GABAergic_VIP" = "#0072B2",           # Dark blue
  "GABAergic_LAMP5" = "#CC0000",         # Red
  "GABAergic_PVALB" = "#9900CC",         # Purple
  "Endothelial" = "#999999",             # Gray
  "Pericytes" = "#666666"                # Dark gray
)

# Create the main DimPlot
p1 <- DimPlot(chosen_obj,
             reduction = "wnn.umap",
             group.by = "celltype",
             cols = celltype_colors,
             label = TRUE,
             label.size = 5,
             repel = TRUE,
             raster = FALSE) +
  ggtitle("Cell Type Annotation - Prefrontal Cortex") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))

# Create a version without labels for cleaner look
p2 <- DimPlot(chosen_obj,
             reduction = "wnn.umap",
             group.by = "celltype",
             cols = celltype_colors,
             label = FALSE,
             raster = FALSE) +
  ggtitle("Cell Type Annotation - Prefrontal Cortex") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
       legend.text = element_text(size = 10))

# Save both versions
pdf(file.path(project_outdir, "UMAP_celltype_labeled_font5.pdf"), width = 12, height = 10)
print(p1)
dev.off()

pdf(file.path(project_outdir, "UMAP_celltype_clean.pdf"), width = 12, height = 10)
print(p2)
dev.off()

# Display in R
print(p1)

cat("\n✓ DimPlots saved to:\n")
cat(sprintf("  - %s\n", file.path(project_outdir, "UMAP_celltype_labeled.pdf")))
cat(sprintf("  - %s\n", file.path(project_outdir, "UMAP_celltype_clean.pdf")))

# ============================================================================
# CELL TYPE SUMMARY STATISTICS
# ============================================================================

# Calculate summary by cell type
celltype_summary <- chosen_obj@meta.data %>%
  group_by(celltype) %>%
  summarize(
    n_cells = n(),
    percent = 100 * n() / ncol(chosen_obj),
    median_nFeature = median(nFeature_RNA),
    median_nCount = median(nCount_RNA),
    median_pct_mt = median(percent.mt)
  ) %>%
  arrange(desc(n_cells))

# Add cell class grouping
celltype_summary <- celltype_summary %>%
  mutate(
    major_class = case_when(
      grepl("Excitatory", celltype) ~ "Neurons_Excitatory",
      grepl("GABAergic", celltype) ~ "Neurons_GABAergic",
      celltype %in% c("Oligodendrocytes", "OPCs") ~ "Oligodendrocyte_Lineage",
      celltype == "Astrocytes" ~ "Astrocytes",
      celltype == "Microglia" ~ "Microglia",
      celltype %in% c("Endothelial", "Pericytes") ~ "Vascular",
      TRUE ~ "Other"
    )
  )

# Display
cat("\n=== CELL TYPE COMPOSITION ===\n")
print(celltype_summary, n = 20)

# Summary by major class
major_class_summary <- celltype_summary %>%
  group_by(major_class) %>%
  summarize(
    n_cells = sum(n_cells),
    percent = sum(percent),
    n_subtypes = n()
  ) %>%
  arrange(desc(percent))

cat("\n=== MAJOR CELL CLASSES ===\n")
print(major_class_summary)

# Save
write.csv(celltype_summary,
         file.path(project_outdir, "celltype_summary_final.csv"),
         row.names = FALSE)

write.csv(major_class_summary,
         file.path(project_outdir, "major_class_summary_final.csv"),
         row.names = FALSE)

cat("\n✓ Summary statistics saved\n")

# ============================================================================
# SAVE FILTERED OBJECT
# ============================================================================

# Save the filtered object with cell type annotations
saveRDS(chosen_obj, file = file.path(rdsdir, "chosen_obj_filtered_annotated.rds"))

