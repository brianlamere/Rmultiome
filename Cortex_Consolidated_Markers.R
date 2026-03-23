# ============================================================================
# Consolidated Cortical Cell Type Markers - High-Impact References
# ============================================================================
# 
# Primary Sources:
# 1. Mathys et al. 2019 Nature (DOI: 10.1038/s41586-019-1195-2)
#    - Prefrontal cortex snRNA-seq, Alzheimer's disease
#    - Comprehensive cell type markers for all major cortical types
#
# 2. Sankowski et al. 2019 Nat Neurosci (DOI: 10.1038/s41593-019-0532-y)
#    - Microglia vs macrophage distinction in human brain
#
# 3. Zillich et al. 2025 Nat Commun (DOI: 10.1038/s41467-024-55467-0)
#    - Caudate nucleus multiome, validation reference
# ============================================================================

source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
#init_project()

# === OLIGODENDROCYTES ===
oligodendrocyte_markers <- data.frame(
  gene = c("MBP", "MOBP", "PLP1", "MOG", "MAG"),
  source = "Mathys_2019",
  confidence = c("high", "high", "high", "medium", "medium"),
  notes = "Mature myelinating oligodendrocytes"
)

# === OLIGODENDROCYTE PRECURSOR CELLS (OPCs) ===
opc_markers <- data.frame(
  gene = c("VCAN", "PDGFRA", "PCDH15", "CSPG4", "SOX10"),
  source = "Mathys_2019",
  confidence = c("high", "high", "high", "medium", "medium"),
  notes = "Oligodendrocyte precursors, proliferative"
)

# === ASTROCYTES ===
astrocyte_markers <- data.frame(
  gene = c("GFAP", "ALDH1L1", "GLUL", "AQP4", "SLC1A2", "SLC4A4"),
  source = c(rep("Mathys_2019", 4), "Zillich_2025", "Zillich_2025"),
  confidence = rep("high", 6),
  notes = "Mature protoplasmic astrocytes; GFAP may be low in homeostatic state"
)

# === MICROGLIA (CNS-resident) ===
microglia_markers <- data.frame(
  gene = c("P2RY12", "TMEM119", "CX3CR1", "GPR34", "OLFML3", 
           "APBB1IP", "CSF1R", "SALL1"),
  source = c(rep("Sankowski_2019", 5), "Zillich_2025", 
             "Mathys_2019", "Sankowski_2019"),
  confidence = c("high", "high", "high", "high", "high", 
                "medium", "low_specificity", "high"),
  notes = c(
    "THE definitive microglia marker (homeostatic)",
    "Microglia-specific transmembrane protein",
    "Fractalkine receptor, high in microglia",
    "GPCR, microglia-specific",
    "Olfactomedin-like, microglia-specific",
    "Microglia marker from Zillich",
    "Pan-myeloid, NOT specific (both microglia and macrophages)",
    "Transcription factor, microglia identity"
  )
)

# === MACROPHAGES (CNS-infiltrating / border-associated) ===
macrophage_markers <- data.frame(
  gene = c("CD163", "MRC1", "LYVE1", "CD14", "FCGR3A", "MS4A7", "F13A1"),
  source = rep("Sankowski_2019", 7),
  confidence = c("high", "high", "high", "high", "medium", "medium", "low"),
  notes = c(
    "Scavenger receptor, perivascular macrophages (CRITICAL for distinguishing from microglia)",
    "Mannose receptor (CD206), macrophages",
    "Lymphatic marker, border-associated macrophages",
    "Monocyte marker, blood-derived macrophages",
    "CD16, non-classical monocytes",
    "Monocyte/macrophage marker",
    "Factor XIII, tissue macrophages"
  )
)

# === EXCITATORY NEURONS ===
excitatory_markers <- data.frame(
  gene = c("SLC17A7", "CAMK2A", "SATB2", "TBR1", "NRGN"),
  source = c(rep("Mathys_2019", 4), "General"),
  confidence = rep("high", 5),
  notes = "Glutamatergic neurons; SATB2/TBR1 for layer identity"
)

# Layer-specific excitatory markers (optional, for subtyping)
excitatory_layer_markers <- data.frame(
  gene = c("FOXP2", "RORB", "PCP4", "BCL11B", "FEZF2"),
  layer = c("L6", "L4", "L5", "L5", "L5_PT"),
  source = rep("Mathys_2019", 5),
  confidence = c("high", "medium", "medium", "high", "high"),
  notes = c(
    "Layer 6 corticothalamic neurons",
    "Layer 4 granular neurons",
    "Layer 5 marker",
    "Layer 5 intratelencephalic neurons",
    "Layer 5 pyramidal tract (PT) neurons"
  )
)

# === GABAergic INTERNEURONS ===
# Core GABAergic
gabaergic_core <- data.frame(
  gene = c("GAD1", "GAD2", "SLC32A1"),
  source = rep("Mathys_2019", 3),
  confidence = rep("high", 3),
  notes = "Pan-GABAergic markers (all interneurons)"
)

# GABAergic subtypes
gabaergic_subtypes <- data.frame(
  gene = c("PVALB", "SST", "VIP", "LAMP5", "CALB2"),
  subtype = c("PVALB", "SST", "VIP", "LAMP5", "VIP/LAMP5"),
  source = rep("Mathys_2019", 5),
  confidence = rep("high", 5),
  notes = c(
    "Parvalbumin+ fast-spiking interneurons",
    "Somatostatin+ Martinotti cells",
    "Vasoactive intestinal peptide+ interneurons",
    "LAMP5+ neurogliaform cells",
    "CALB2 (calretinin) in VIP/LAMP5 subsets"
  )
)

# Additional subtype markers (from Zillich, striatal-specific)
gabaergic_striatal <- data.frame(
  gene = c("DRD1", "DRD2", "ADORA2A", "TAC1", "ADARB2", "PTHLH"),
  subtype = c("D1_MSN", "D2_MSN", "D2_MSN", "D1_MSN", "D1_ADARB2", "PTHLH_PVALB"),
  source = rep("Zillich_2025", 6),
  confidence = rep("medium", 6),
  notes = "Striatal medium spiny neuron markers; may not be present in frontal cortex"
)

# === ENDOTHELIAL CELLS ===
endothelial_markers <- data.frame(
  gene = c("FLT1", "CLDN5", "KDR", "VWF"),
  source = c(rep("Mathys_2019", 3), "Mathys_2019"),
  confidence = c("high", "high", "high", "medium"),
  notes = "Vascular endothelial cells; CLDN5 = tight junctions (BBB)"
)

# === PERICYTES ===
pericyte_markers <- data.frame(
  gene = c("PDGFRB", "RGS5", "NOTCH3", "CARMN", "ACTA2"),
  source = rep("Mathys_2019", 5),
  confidence = c("high", "high", "high", "medium", "low"),
  notes = c(
    "PDGF receptor beta, definitive pericyte marker",
    "Regulator of G-protein signaling 5",
    "NOTCH3 mutations cause CADASIL (pericyte dysfunction)",
    "Cardiac mesoderm marker",
    "Smooth muscle actin (shared with smooth muscle)"
  )
)

# === LYMPHOCYTES (T-cells, if present) ===
lymphocyte_markers <- data.frame(
  gene = c("CD96", "CD3D", "IL2RG"),
  source = rep("Zillich_2025", 3),
  confidence = c("medium", "high", "high"),
  notes = c(
    "NK cell and T-cell marker",
    "T-cell receptor component (definitive)",
    "Common gamma chain, lymphocyte marker"
  )
)

# === STRESSED / LOW-QUALITY CELLS (for QC flagging) ===
stress_qc_markers <- data.frame(
  gene = c("VIM", "GFAP", "NEFM", "NEFL", "NRGN"),
  source = "Your_pipeline",
  confidence = "QC_flag",
  notes = "Combination suggests stressed/dying cells with mixed identity (remove from analysis)"
)

# ============================================================================
# COMBINE INTO MASTER REFERENCE
# ============================================================================

cortex_marker_reference <- bind_rows(
  oligodendrocyte_markers %>% mutate(celltype = "Oligodendrocytes"),
  opc_markers %>% mutate(celltype = "OPCs"),
  astrocyte_markers %>% mutate(celltype = "Astrocytes"),
  microglia_markers %>% mutate(celltype = "Microglia"),
  macrophage_markers %>% mutate(celltype = "Macrophages"),
  excitatory_markers %>% mutate(celltype = "Excitatory_Neurons"),
  excitatory_layer_markers %>% mutate(celltype = "Excitatory_Layer_Specific"),
  gabaergic_core %>% mutate(celltype = "GABAergic_Core"),
  gabaergic_subtypes %>% mutate(celltype = paste0("GABAergic_", subtype)),
  endothelial_markers %>% mutate(celltype = "Endothelial"),
  pericyte_markers %>% mutate(celltype = "Pericytes"),
  lymphocyte_markers %>% mutate(celltype = "Lymphocytes"),
  stress_qc_markers %>% mutate(celltype = "QC_Flag")
)

# Create simplified lists for quick plotting/typing
marker_lists <- list(
  Oligodendrocytes = c("MBP", "MOBP", "PLP1"),
  OPCs = c("VCAN", "PDGFRA", "PCDH15"),
  Astrocytes = c("GFAP", "ALDH1L1", "AQP4", "SLC1A2"),
  Microglia = c("P2RY12", "TMEM119", "CX3CR1"),
  Macrophages = c("CD163", "MRC1", "LYVE1"),  # ADD THIS for biologist
  Excitatory_Neurons = c("SLC17A7", "CAMK2A", "SATB2", "TBR1"),
  GABAergic_PVALB = c("GAD1", "GAD2", "PVALB"),
  GABAergic_SST = c("GAD1", "GAD2", "SST"),
  GABAergic_VIP = c("GAD1", "GAD2", "VIP"),
  GABAergic_LAMP5 = c("GAD1", "GAD2", "LAMP5"),
  Endothelial = c("FLT1", "CLDN5"),
  Pericytes = c("PDGFRB", "RGS5", "NOTCH3")
)

# ============================================================================
# SAVE
# ============================================================================

cortex_marker_file <- file.path(Rmultiome_path, "Cortex_Consolidated_Markers.rds")
saveRDS(list(
  reference_table = cortex_marker_reference,
  marker_lists = marker_lists
), cortex_marker_file)

write.csv(cortex_marker_reference, 
          file.path(project_export, "Cortex_Consolidated_Markers.csv"),
          row.names = FALSE)

cat("✓ Cortex marker reference created\n")
cat(sprintf("  - %d markers across %d cell types\n", 
           nrow(cortex_marker_reference),
           length(unique(cortex_marker_reference$celltype))))
cat(sprintf("  - Primary sources: Mathys 2019, Sankowski 2019\n"))
cat(sprintf("  - Saved to: %s\n", cortex_marker_file))
