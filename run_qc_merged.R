source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

init_project()
#load object created at end of run_pipeline1.R
merged_obj <- readRDS(file.path(rdsdir,"merged_preharmony.Rds"))

# here in case you're restarting and already have these settings
# set to TRUE then skip to where you were
if (FALSE) {
	pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
	cluster_settings <- read_cluster_settings(cluster_settings_file)
	celltype_settings <- read_celltype_settings(celltype_settings_file)
	harmony_settings <- read_harmony_settings(harmony_settings_file)
}

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
  random_seed = random_seed,
  max_iter = 50,
  project_dim = FALSE
)

saveRDS(harmony_settings, harmony_settings_file)

# === STEP 3: Run Harmony for parameter sweep ===
harmony_obj <- harmonize_both(
  merged_obj,
  random_seed = random_seed,
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
# I like 1:30knn40res0.16, 1:40k40r.160 
sweep_results <- run_parameter_sweep_plots(
  seurat_obj = harmony_obj,
  dims_range = list(c(1:30),c(1:40)),
  knn_values = c(30, 40),
  res_values = c(0.16,0.18),
  alg = 3,                     # SLM algorithm (required parameter)
  cluster_seed = random_seed   # Reproducibility (required parameter)
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
chosen_dims_max <- 40     # CHANGE THIS
chosen_knn <- 44          # CHANGE THIS
chosen_resolution <- 0.18 # CHANGE THIS

# === STEP 7: Save cluster settings ===
cluster_settings <- data.frame(
  dims_min = chosen_dims_min,
  dims_max = chosen_dims_max,
  knn = chosen_knn,
  resolution = chosen_resolution,
  algorithm = 3,
  random_seed = random_seed
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
  cluster_seed = random_seed,
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
#maybe_new_device(width = 10, height = 8)
#print(DimPlot(chosen_obj, reduction = "wnn.umap", label = TRUE, raster = FALSE))

# Save
saveRDS(chosen_obj, file.path(rdsdir, "clustered_obj.rds"))
# ISKIPPEDTOTHISLINE
# chosen_obj <- readRDS(file.path(rdsdir, "clustered_obj.rds"))

# ============================================================================
# PHASE 2: Cell Type Annotation
# ============================================================================

# Note: refer to lymphocyte_investigation.R for exclusion reason

# Load consolidated markers
cortex_markers <- readRDS(file.path(Rmultiome_path, "Cortex_Consolidated_Markers.rds"))

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
  min.pct = 0.10,
  logfc.threshold = 0.25,
  test.use = "wilcox"
)

#adding just to speed up reproducing the error.  remove after.
# saveRDS(all_markers, file.path(tmpfiledir, "saved_all_markers.rds"))
# all_markers <- readRDS(file.path(tmpfiledir, "saved_all_markers.rds"))

# Use for automated typing
results <- identify_all_celltypes(
  all_markers,
  cortex_markers$marker_lists,  # Uses the simplified lists
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)


# Quick microglia vs macrophage check
microglia_genes <- cortex_markers$marker_lists$Microglia
macrophage_genes <- cortex_markers$marker_lists$Macrophages

DotPlot(chosen_obj, 
        features = c(microglia_genes, macrophage_genes),
        idents = suspected_myeloid_cluster) +
  ggtitle("Microglia vs Macrophage: Sankowski 2019 markers")

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
  cluster = c(0:29, "singleton"),
  celltype = c(
    "Oligodendrocytes",           # 0 - high confidence (8.60)
    "Excitatory_Neurons",         # 1 - high confidence (9.31)
    "Astrocytes",                 # 2 - high confidence (32.03)
    "Microglia_Macrophages",      # 3 - MIXED (Micro: 35.27, Macro: 38.89)
    "OPCs",                       # 4 - high confidence (89.60)
    "Astrocytes",                 # 5 - potential match (GFAP, AQP4)
    "GABAergic_VIP",              # 6 - high confidence (35.57)
    "Excitatory_Neurons",         # 7 - unassigned, likely neuronal
    "GABAergic_SST",              # 8 - high confidence (36.59)
    "GABAergic_PVALB",            # 9 - high confidence (38.50)
    "Excitatory_Neurons",         # 10 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 11 - potential match (7.63)
    "Excitatory_Neurons",         # 12 - unassigned, likely neuronal
    "Pericytes",                  # 13 - HIGH confidence (74.00) ✅ CLEAN!
    "GABAergic_LAMP5",            # 14 - high confidence (32.64)
    "Endothelial",                # 15 - HIGH confidence (210.29) ✅ CLEAN!
    "Excitatory_Neurons",         # 16 - potential match (8.04)
    "Excitatory_Neurons",         # 17 - unassigned, likely neuronal
    "GABAergic_PVALB",            # 18 - potential match (PVALB, GAD1)
    "Oligodendrocytes",           # 19 - potential match (MBP, PLP1, MOBP)
    "Excitatory_Neurons",         # 20 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 21 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 22 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 23 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 24 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 25 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 26 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 27 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 28 - unassigned, likely neuronal
    "Excitatory_Neurons",         # 29 - unassigned, likely neuronal
    "QC_Remove"                   # singleton
  ),
  action = c(
    rep("keep", 30),              # Keep clusters 0-29
    "remove"                      # Remove singleton
  ),
  markers_used = c(
    "MBP, MOBP, PLP1",                                    # 0
    "SLC17A7, CAMK2A, SATB2, TBR1",                      # 1
    "GFAP, ALDH1L1, GLUL, AQP4, SLC1A2, SLC4A4",        # 2
    "P2RY12, CX3CR1 (microglia); CD163, MRC1 (macrophages)", # 3
    "VCAN, PDGFRA, PCDH15",                              # 4
    "GFAP, AQP4 (weak)",                                 # 5
    "GAD1, GAD2, VIP",                                   # 6
    "Excitatory markers (check manually)",               # 7
    "GAD1, GAD2, SST",                                   # 8
    "GAD1, GAD2, PVALB",                                 # 9
    "Excitatory markers (check manually)",               # 10
    "SLC17A7, CAMK2A, SATB2, TBR1",                     # 11
    "Excitatory markers (check manually)",               # 12
    "PDGFRB, RGS5, NOTCH3",                             # 13 PERICYTES!
    "GAD1, GAD2, LAMP5",                                 # 14
    "FLT1, CLDN5",                                       # 15 ENDOTHELIAL!
    "SATB2, TBR1, SLC17A7, CAMK2A",                     # 16
    "Excitatory markers (check manually)",               # 17
    "PVALB, GAD1",                                       # 18
    "MBP, PLP1, MOBP (weak)",                           # 19
    "Excitatory markers (check manually)",               # 20
    "Excitatory markers (check manually)",               # 21
    "Excitatory markers (check manually)",               # 22
    "Excitatory markers (check manually)",               # 23
    "Excitatory markers (check manually)",               # 24
    "Excitatory markers (check manually)",               # 25
    "Excitatory markers (check manually)",               # 26
    "Excitatory markers (check manually)",               # 27
    "Excitatory markers (check manually)",               # 28
    "Excitatory markers (check manually)",               # 29
    "QC_flag"                                            # singleton
  ),
  confidence = c(
    "high",    # 0
    "high",    # 1
    "high",    # 2
    "high",    # 3 (biological mixture)
    "high",    # 4
    "low",     # 5 (weak markers)
    "high",    # 6
    "medium",  # 7
    "high",    # 8
    "high",    # 9
    "medium",  # 10
    "high",    # 11
    "medium",  # 12
    "high",    # 13 PERICYTES
    "high",    # 14
    "high",    # 15 ENDOTHELIAL
    "high",    # 16
    "medium",  # 17
    "medium",  # 18
    "low",     # 19 (weak oligo markers)
    "medium",  # 20-29 (unassigned excitatory)
    "medium",
    "medium",
    "medium",
    "medium",
    "medium",
    "medium",
    "medium",
    "medium",
    "medium",
    "high"     # singleton (remove)
  ),
  notes = c(
    "Mature myelinating oligodendrocytes",
    "Excitatory neurons, cortical layer unspecified",
    "Protoplasmic astrocytes",
    "MIXED: CNS-resident microglia (P2RY12+ 26%, CX3CR1+ 14%) + infiltrating macrophages (CD163+ 13%, MRC1+ 13%). Biologically expected in HIV+ tissue. Macrophage score slightly higher (38.89 vs 35.27).",
    "Oligodendrocyte precursor cells",
    "Weak astrocyte markers; possible reactive astrocytes or mixed population",
    "VIP+ interneurons",
    "Likely excitatory neurons; verify with top markers",
    "SST+ Martinotti interneurons",
    "PVALB+ fast-spiking basket interneurons",
    "Likely excitatory neurons; verify with top markers",
    "Excitatory neurons, likely upper cortical layers",
    "Likely excitatory neurons; verify with top markers",
    "Pericytes - mural cells surrounding blood vessels. CLEANLY SEPARATED from endothelial!",
    "LAMP5+ neurogliaform interneurons",
    "Vascular endothelial cells with tight junction markers. CLEANLY SEPARATED from pericytes!",
    "Excitatory neurons, likely deep cortical layers",
    "Likely excitatory neurons; verify with top markers",
    "Additional PVALB+ interneuron population",
    "Weak oligodendrocyte markers; possible stressed/transitioning cells",
    rep("Likely excitatory neurons; verify with top markers", 10),  # 20-29
    "Singletons - QC flagged for removal"
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

# After creating celltype_mapping, test it:
chosen_obj <- apply_celltype_labels(
  chosen_obj,
  celltype_settings = celltype_mapping,
  remove_flagged = TRUE,     # Remove cluster 6 and singleton
  verbose = TRUE
)

# Visualize
DimPlot(chosen_obj, group.by = "celltype", label = TRUE, raster = FALSE)

# Update the main object name
# NO.  This creates a non-idempotent process.
#chosen_obj <- chosen_obj_filtered
#rm(chosen_obj_filtered)
saveRDS(chosen_obj, file.path(rdsdir, "final_qc_obj.rds"))

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
p1 <- DimPlot(chosen_obj_filtered,
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
p2 <- DimPlot(chosen_obj_filtered,
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
celltype_summary <- chosen_obj_filtered@meta.data %>%
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
