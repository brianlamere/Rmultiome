source("/projects1/opioid/Rmultiome/system_settings.R")
source("/projects1/opioid/Rmultio:me/Rmultiome-main.R")

all_markers <- read.csv(file.path(tmpfiledir, "all_cluster_markers.csv"))

# TBD clusters
tbd_clusters <- c(1, 4, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16)

# === Step 1: Check for excitatory neurons (HIGHEST PRIORITY) ===
cat("=== Checking for excitatory neurons ===\n")

excitatory_markers <- c("SLC17A7", "CAMK2A", "SATB2", "TBR1", "NEUROD6", "NRGN")

excitatory_results <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = excitatory_markers,
  celltype_name = "Excitatory_Neurons",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# === Step 2: Check cluster 6 (related to oligodendrocytes) ===
cat("\n=== Cluster 6: Oligodendrocyte-related ===\n")

oligo_extended_markers <- c("OLIG1", "OLIG2", "SOX10", "MBP", "PLP1", "CNP")

clust6_markers <- all_markers %>%
  filter(cluster == 6) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)

cat("Top 10 markers for cluster 6:\n")
print(clust6_markers %>% select(gene, avg_log2FC, pct.1, pct.2))

# Check if oligodendrocyte markers present
oligo_in_6 <- clust6_markers %>% filter(gene %in% oligo_extended_markers)
if (nrow(oligo_in_6) > 0) {
  cat("\nOligodendrocyte markers in cluster 6:\n")
  print(oligo_in_6)
}

# === Step 3: Check cluster 12 (related to endothelial) ===
cat("\n=== Cluster 12: Endothelial-related (pericytes?) ===\n")

pericyte_markers <- c("PDGFRB", "RGS5", "ACTA2", "CSPG4", "ANPEP")

pericyte_result <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = pericyte_markers,
  celltype_name = "Pericytes",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# === Step 4: Check for GABAergic subtypes ===
cat("\n=== Checking for GABAergic neuron subtypes ===\n")

# D2 MSN
d2_result <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = c("DRD2", "ADORA2A", "GPR6"),
  celltype_name = "GABAergic_D2_MSN",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# PVALB+
pvalb_result <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = c("PVALB", "GAD1", "GAD2"),
  celltype_name = "GABAergic_PVALB",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# === Step 5: Show top markers for all remaining TBD ===
cat("\n=== Top markers for all TBD clusters ===\n")

for (clust in tbd_clusters) {
  cat(sprintf("\n--- Cluster %d ---\n", clust))
  
  top5 <- all_markers %>%
    filter(cluster == clust) %>%
    arrange(desc(avg_log2FC)) %>%
    head(5)
  
  print(top5)
}


# Check if it's doublets
# Check doublet scores if you have them

# Or check if it's immature astrocytes
immature_astro_markers <- c("VIM", "GFAP", "GLAST", "FABP7", "HES1")

# Or radial glia
radial_glia_markers <- c("VIM", "GFAP", "NES", "SOX2", "PAX6")


imm_ast <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = immature_astro_markers,
  celltype_name = "Immature Astrocytes",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

radial_glia <- identify_celltype(
  all_markers %>% filter(cluster %in% tbd_clusters),
  marker_genes = radial_glia_markers,
  celltype_name = "Radial Glia",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)


microglia_specific <- c(
  "P2RY12",   # THE definitive microglia marker (you used this!)
  "TMEM119",  # Microglia-specific transmembrane protein
  "CX3CR1",   # High in microglia, lower in macrophages
  "SALL1",    # Transcription factor, microglia identity
  "FCRLS",    # Fc receptor-like, microglia-specific
  "SIGLECH"   # Sialic acid binding, microglia
)

macrophage_specific <- c(
  "CD163",    # Scavenger receptor, macrophages
  "MRC1",     # Mannose receptor (CD206), macrophages
  "LYVE1",    # Lymphatic vessel marker, perivascular
  "F13A1",    # Factor XIII, macrophages
  "APOE",     # High in macrophages, lower in microglia
  "CD163L1",  # Macrophage marker
  "MS4A7"     # Monocyte/macrophage marker
)

shared_myeloid <- c(
  "CD68",     # Pan-myeloid
  "AIF1",     # IBA1 - both express
  "CSF1R",    # Both depend on CSF1 (you used this!)
  "CD14",     # Monocyte/myeloid
  "ITGAM"     # CD11b - both express
)

clust3_markers <- all_markers %>%
  filter(cluster == 3)

# Microglia markers in cluster 3
cat("=== Microglia-specific markers in cluster 3 ===\n")
microglia_check <- clust3_markers %>%
  filter(gene %in% microglia_specific) %>%
  arrange(desc(avg_log2FC))

if (nrow(microglia_check) > 0) {
  print(microglia_check %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2))
} else {
  cat("No specific microglia markers found in top DE genes\n")
}

# Macrophage markers in cluster 3
cat("\n=== Macrophage-specific markers in cluster 3 ===\n")
macrophage_check <- clust3_markers %>%
  filter(gene %in% macrophage_specific) %>%
  arrange(desc(avg_log2FC))

if (nrow(macrophage_check) > 0) {
  print(macrophage_check %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2))
  cat("\n⚠️ MACROPHAGE MARKERS PRESENT - May have infiltrating macrophages!\n")
} else {
  cat("No macrophage markers found - likely pure microglia\n")
}

# Visualize with DotPlot
cat("\n=== Generating DotPlot for microglia vs macrophage markers ===\n")
DotPlot(chosen_obj, 
        features = c(microglia_specific, macrophage_specific),
        idents = 3,
        cols = "RdYlBu") + 
  RotatedAxis() +
  ggtitle("Cluster 3: Microglia vs Macrophage markers")



clust13_markers <- all_markers %>%
  filter(cluster == 13) %>%
  arrange(desc(avg_log2FC)) %>%
  head(20)

print(clust13_markers %>% dplyr::select(gene, avg_log2FC, pct.1, pct.2))
