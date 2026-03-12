# Markers from Zillich et al. 2025 Nature Communications
# "Multiome analysis of cocaine use disorder"
# Caudate nucleus (similar to your frontal cortex)

cocaine_paper_markers <- data.frame(
  celltype = c(
    rep("Oligodendrocytes", 3),
    rep("OPCs", 3),
    rep("Astrocytes", 6),
    rep("Microglia", 3),
    rep("GABAergic_Neurons", 10),
    rep("Endothelial", 3),
    rep("Lymphocytes", 3)
  ),
  gene = c(
    # Oligodendrocytes
    "MBP", "MOBP", "PLP1",
    # OPCs
    "VCAN", "PDGFRA", "PCDH15",
    # Astrocytes
    "GFAP", "ALDH1L1", "GLUL", "AQP4", "SLC1A2", "SLC4A4",
    # Microglia
    "CSF1R", "APBB1IP", "P2RY12",
    # GABAergic (they had subtypes, but core markers)
    "GAD1", "GAD2", "DRD1", "DRD2", "ADORA2A", "GPR6", 
    "ADARB2", "PTHLH", "PVALB", "CALB2",
    # Endothelial
    "FLT1", "CLDN5", "KDR",
    # Lymphocytes
    "CD96", "CD3D", "IL2RG"
  ),
  source = "Zillich_2025_NatComm",
  notes = "Caudate nucleus multiome, cocaine study",
  status_in_your_data = NA,  # To fill in during validation
  stringsAsFactors = FALSE
)
