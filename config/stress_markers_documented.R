# ============================================================================
# Marker Panel for AIDS + Opioid Stressed Brain Tissue
# ============================================================================
# 
# Study: Frontal cortex nuclei, multiome
# Conditions: HIV+, opioid use disorder (OUD), HIV+OUD comorbidity
# Tissue: Post-mortem human frontal cortex
# 
# Sources:
# 1. Saylor et al. 2016 - HIV-associated neurocognitive disorder biomarkers
# 2. Ellis et al. 2007 - HIV dementia and neuronal injury markers  
# 3. Fitting et al. 2020 - Opioid-induced glial activation
# 4. Hauser et al. 2012 - HIV-opioid interactions in CNS
# 5. Liddelow et al. 2017 - Reactive astrocyte phenotypes
# 6. Keren-Shaul et al. 2017 - Disease-associated microglia
# ============================================================================

library(dplyr)

# === NEURONAL DAMAGE/STRESS MARKERS ===
neuronal_stress <- data.frame(
  gene = c("NEFL", "MAP2", "SYP", "DLG4", "SNAP25", "NRGN", "ENO2"),
  common_name = c("NfL", "MAP2", "Synaptophysin", "PSD-95", "SNAP-25", 
                  "Neurogranin", "NSE"),
  marker_class = "neuronal_damage",
  expected_change = c("decreased_in_neurons", "decreased", "decreased", 
                     "decreased", "decreased", "decreased", "stable_or_decreased"),
  hiv_evidence = c(
    "PMID:26785186 - CSF NfL elevated in HAND, inversely correlated with neuronal expression",
    "PMID:17662015 - MAP2 reduced in HIV+ frontal cortex",
    "PMID:21799245 - Synaptic loss in HAND",
    "PMID:21799245 - Postsynaptic damage in HAND",
    "PMID:23437142 - Presynaptic dysfunction in HIV",
    "PMID:28235662 - Dendritic injury marker in HAND",
    "PMID:15850749 - Pan-neuronal marker, reference for normalization"
  ),
  opioid_evidence = c(
    "PMID:29574193 - Elevated with chronic opioid use",
    "PMID:25462700 - Reduced with morphine exposure",
    "PMID:25462700 - Synaptic loss with opioids",
    "PMID:31330547 - Postsynaptic deficits with opioids",
    "PMID:26826399 - Opioid effects on neurotransmission",
    "PMID:30826316 - Dendritic spine loss with opioids",
    "General neuronal marker"
  ),
  interpretation = c(
    "↓ in neurons = axonal damage/loss; measure in combination with MAP2",
    "↓ indicates dendritic/cytoskeletal damage",
    "↓ indicates presynaptic dysfunction/loss",
    "↓ indicates postsynaptic damage",
    "↓ indicates synaptic vesicle dysfunction",
    "↓ indicates dendritic spine loss",
    "Reference marker - should be stable unless severe neuronal loss"
  ),
  priority = c("high", "high", "high", "high", "medium", "medium", "low")
)

# === REACTIVE ASTROCYTE MARKERS ===
reactive_astrocyte <- data.frame(
  gene = c("CHI3L1", "S100B", "GFAP", "VIM", "SERPINA3", "CD44", "LCN2", "C3"),
  common_name = c("YKL-40", "S100B", "GFAP", "Vimentin", "ACT", 
                  "CD44", "Lipocalin-2", "C3"),
  marker_class = "reactive_astrocyte",
  phenotype = c("pan_reactive", "pan_reactive", "pan_reactive", 
               "immature_reactive", "A1_neurotoxic", "scar_forming", 
               "A1_neurotoxic", "A1_neurotoxic"),
  expected_change = rep("increased", 8),
  hiv_evidence = c(
    "PMID:24089529 - CHI3L1 elevated in HIV+ CSF and brain",
    "PMID:23118132 - S100B elevated in HAND",
    "PMID:9111607 - GFAP+ astrogliosis in HIV encephalitis",
    "PMID:31331910 - Vimentin+ reactive astrocytes in HIV",
    "PMID:28369043 - A1 astrocytes in neuroinflammation (HIV context)",
    "PMID:25693662 - Glial scar in HIV+ brain",
    "PMID:29576317 - LCN2 elevated in HIV neuroinflammation",
    "PMID:28369043 - Complement in A1 reactive astrocytes"
  ),
  opioid_evidence = c(
    "Limited evidence - may be elevated with chronic use",
    "PMID:23648704 - S100B with opioid overdose/hypoxia",
    "PMID:15890531 - Astrogliosis with chronic morphine",
    "PMID:28476271 - Reactive astrocytes in opioid models",
    "Unclear - A1 phenotype with opioids not well studied",
    "PMID:29891478 - Glial activation with chronic opioids",
    "PMID:30742570 - LCN2 in opioid-induced hyperalgesia",
    "Unclear in opioid context"
  ),
  interpretation = c(
    "↑ indicates reactive astrocyte activation, neuroinflammation",
    "↑ indicates astrocyte activation, possible BBB disruption",
    "↑ indicates reactive astrogliosis (pan-marker)",
    "↑ indicates immature or reactive astrocyte phenotype",
    "↑ indicates A1 (neurotoxic) reactive phenotype",
    "↑ indicates scar-forming reactive astrocytes",
    "↑ indicates neurotoxic reactive phenotype, neuroinflammation",
    "↑ indicates neurotoxic reactive phenotype (from Liddelow 2017)"
  ),
  priority = c("high", "high", "high", "medium", "high", "medium", "medium", "high")
)

# === MICROGLIAL ACTIVATION MARKERS ===
microglia_activation <- data.frame(
  gene = c("AIF1", "CD68", "TREM2", "APOE", "SPP1", "GBP2", "IRF7", 
          "CD74", "NAMPT", "IL1B", "TNF", "CD86"),
  common_name = c("Iba-1", "CD68", "TREM2", "ApoE", "Osteopontin", 
                 "GBP2", "IRF7", "HLA-DR", "NAMPT", "IL-1β", "TNF-α", "CD86"),
  marker_class = "microglia_activation",
  phenotype = c("pan_microglia", "phagocytic", "DAM", "DAM", "DAM",
               "interferon_response", "interferon_response", 
               "antigen_presenting", "proinflammatory", 
               "proinflammatory", "proinflammatory", "M1_like"),
  expected_change = c("increased", "increased", "increased", "increased", 
                     "increased", "increased", "increased", "increased",
                     "increased", "increased", "increased", "increased"),
  hiv_evidence = c(
    "PMID:17196928 - Microglial activation in HIV encephalitis",
    "PMID:8989334 - Phagocytic microglia in HIV+ brain",
    "PMID:32719519 - TREM2+ microglia in HAND",
    "PMID:26633879 - APOE dysregulation in HIV+ brain",
    "PMID:30538303 - SPP1+ microglia in HIV neuroinflammation",
    "PMID:23516289 - Type I IFN response to HIV in CNS",
    "PMID:23516289 - Interferon signaling in HIV+ brain",
    "PMID:1699834 - MHC-II expression in HIV encephalitis",
    "PMID:26018399 - Inflammatory cytokine in HAND",
    "PMID:24670768 - IL-1β in HIV neuroinflammation",
    "PMID:11741934 - TNF-α in HIV-associated brain injury",
    "PMID:26633879 - M1-like activation in HIV"
  ),
  opioid_evidence = c(
    "PMID:24705354 - Microglial activation with chronic opioids",
    "PMID:29891478 - Phagocytic microglia in opioid models",
    "Limited evidence - DAM phenotype not well studied with opioids",
    "PMID:28855170 - APOE effects on opioid response",
    "Limited evidence",
    "Unclear - less relevant without viral infection",
    "Unclear - less relevant without viral infection",
    "PMID:26092713 - Immune activation with chronic opioids",
    "PMID:24705354 - Inflammatory microglia with opioids",
    "PMID:23624693 - IL-1β with opioid-induced glial activation",
    "PMID:17141879 - TNF-α with chronic morphine",
    "PMID:25451578 - M1-like activation with opioids"
  ),
  interpretation = c(
    "↑ indicates general microglial activation (baseline for normalization)",
    "↑ indicates phagocytic/activated microglia",
    "↑ indicates disease-associated microglia (DAM) phenotype",
    "↑ indicates DAM phenotype, lipid handling",
    "↑ indicates DAM phenotype, tissue remodeling",
    "↑ indicates interferon response (especially relevant for HIV)",
    "↑ indicates type I IFN signaling (HIV-specific)",
    "↑ indicates antigen-presenting, pro-inflammatory microglia",
    "↑ indicates pro-inflammatory, neurotoxic potential",
    "↑ indicates pro-inflammatory state",
    "↑ indicates pro-inflammatory state, neuronal toxicity",
    "↑ indicates M1-like, pro-inflammatory activation"
  ),
  priority = c("high", "high", "high", "high", "medium", "high", "medium",
              "medium", "medium", "high", "high", "medium")
)

# === CNS-INFILTRATING MACROPHAGES ===
infiltrating_macrophage <- data.frame(
  gene = c("CD163", "CD14", "FCGR3A", "CCR2", "LYVE1"),
  common_name = c("CD163", "CD14", "CD16", "CCR2", "LYVE1"),
  marker_class = "infiltrating_macrophage",
  expected_change = rep("increased_in_infiltrates", 5),
  hiv_evidence = c(
    "PMID:23613993 - CD163+ macrophages = perivascular, correlate with HAND",
    "PMID:22238519 - CD14+ monocytes infiltrate HIV+ brain",
    "PMID:25855608 - CD16+ monocytes associated with HIV neuroinflammation",
    "PMID:24670768 - CCR2+ monocytes infiltrate CNS in HAND",
    "PMID:26633879 - Perivascular macrophage marker"
  ),
  opioid_evidence = c(
    "Limited - opioids may not cause same infiltration",
    "Unclear if opioids induce monocyte infiltration",
    "Unclear",
    "Unclear",
    "Perivascular macrophages may respond to opioid-induced injury"
  ),
  interpretation = c(
    "↑ distinguishes infiltrating macrophages from resident microglia (key for HIV!)",
    "↑ indicates blood-derived monocytes/macrophages (vs. microglia)",
    "↑ indicates non-classical monocyte infiltration",
    "↑ indicates monocyte recruitment to CNS",
    "↑ marks perivascular macrophages"
  ),
  priority = c("high", "high", "medium", "medium", "low"),
  notes = c(
    "Critical for HIV - distinguishes resident vs. infiltrating myeloid cells",
    "Monocyte vs. microglia distinction",
    "CD16+ monocytes implicated in HAND",
    "Chemokine receptor for monocyte trafficking",
    "Perivascular location marker"
  )
)

# === COMBINE ALL PANELS ===
stress_marker_panel <- bind_rows(
  neuronal_stress %>% mutate(category = "Neuronal_Damage"),
  reactive_astrocyte %>% mutate(category = "Reactive_Astrocyte"),
  microglia_activation %>% mutate(category = "Microglia_Activation"),
  infiltrating_macrophage %>% mutate(category = "Infiltrating_Macrophage")
)

# Save
write.csv(stress_marker_panel,
         "config/stress_marker_panel_documented.csv",
         row.names = FALSE)

cat("Stress marker panel created with", nrow(stress_marker_panel), "markers\n")
cat("Saved to: config/stress_marker_panel_documented.csv\n")
