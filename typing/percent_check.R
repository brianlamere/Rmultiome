# Current cluster assignments (without cluster 6)
celltype_mapping <- data.frame(
  cluster = 0:17,
  celltype = c(
    "Oligodendrocytes",           # 0
    "Excitatory_Neurons",         # 1
    "Astrocytes",                 # 2
    "Microglia",                  # 3
    "GABAergic_SST",              # 4
    "OPCs",                       # 5
    "Stressed_Dying_Cells",       # 6 - EXCLUDE
    "GABAergic_VIP",              # 7
    "Excitatory_Neurons",         # 8
    "Excitatory_Neurons",         # 9
    "Excitatory_Neurons",         # 10
    "Excitatory_Neurons",         # 11
    "Pericytes",                  # 12
    "Endothelial",                # 13
    "GABAergic_LAMP5",            # 14
    "Excitatory_Neurons",         # 15
    "Excitatory_Neurons",         # 16
    "GABAergic_PVALB"             # 17
  ),
  action = c(
    rep("keep", 6), "remove", rep("keep", 11)
  ),
  stringsAsFactors = FALSE
)

# Add celltype to metadata
chosen_obj$celltype <- celltype_mapping$celltype[match(chosen_obj$seurat_clusters, 
                                                        celltype_mapping$cluster)]

# Calculate proportions (EXCLUDING cluster 6)
cells_to_keep <- chosen_obj$seurat_clusters != "6"

celltype_summary <- chosen_obj@meta.data[cells_to_keep, ] %>%
  group_by(celltype) %>%
  summarize(n_cells = n()) %>%
  mutate(
    percent = 100 * n_cells / sum(n_cells),
    percent_formatted = sprintf("%.1f%%", percent)
  ) %>%
  arrange(desc(n_cells))

# Add expected ranges for frontal cortex
expected_ranges <- data.frame(
  celltype = c(
    "Excitatory_Neurons",
    "GABAergic_SST",
    "GABAergic_VIP", 
    "GABAergic_LAMP5",
    "GABAergic_PVALB",
    "Oligodendrocytes",
    "Astrocytes",
    "Microglia",
    "OPCs",
    "Endothelial",
    "Pericytes"
  ),
  expected_min = c(15, 2, 1, 0.5, 2, 15, 15, 5, 2, 1, 0.5),
  expected_max = c(30, 5, 3, 2, 5, 25, 20, 10, 5, 3, 2),
  source = c(
    rep("Literature (cortex)", 5),
    rep("Literature (cortex)", 6)
  ),
  stringsAsFactors = FALSE
)

# Combine GABAergic for overall check
gabaergic_total <- celltype_summary %>%
  filter(grepl("GABAergic", celltype)) %>%
  summarize(
    celltype = "GABAergic_TOTAL",
    n_cells = sum(n_cells),
    percent = sum(percent)
  )


# Merge with expected
celltype_comparison <- celltype_summary %>%
  left_join(expected_ranges, by = "celltype") %>%
  mutate(
    expected_range = ifelse(!is.na(expected_min),
                           sprintf("%.1f-%.1f%%", expected_min, expected_max),
                           "N/A"),
    status = case_when(
      is.na(expected_min) ~ "No reference",
      percent >= expected_min & percent <= expected_max ~ "✓ Within range",
      percent < expected_min ~ "⚠ Lower than expected",
      percent > expected_max ~ "⚠ Higher than expected"
    )
  ) %>%
  dplyr::select(celltype, n_cells, percent, percent_formatted, expected_range, status)

cat("\n=== Cell Type Composition (excluding cluster 6) ===\n\n")
print(celltype_comparison, n = 20)

# Summary for GABAergic
cat("\n=== GABAergic Neuron Summary ===\n")
cat(sprintf("Total GABAergic: %d cells (%.1f%% of dataset)\n",
           gabaergic_total$n_cells, gabaergic_total$percent))
cat("Expected range: 15-25% for frontal cortex\n")
if (gabaergic_total$percent >= 15 && gabaergic_total$percent <= 25) {
  cat("✓ Within expected range\n")
} else if (gabaergic_total$percent < 15) {
  cat("⚠ Lower than expected (nuclei prep may under-represent interneurons)\n")
} else {
  cat("⚠ Higher than expected\n")
}

# Overall summary
cat("\n=== Summary Statistics ===\n")
cat(sprintf("Total cells analyzed: %d (%.1f%% of original dataset)\n",
           sum(celltype_summary$n_cells),
           100 * sum(celltype_summary$n_cells) / ncol(chosen_obj)))
cat(sprintf("Cells excluded (cluster 6): %d (%.1f%% of original dataset)\n",
           sum(!cells_to_keep),
           100 * sum(!cells_to_keep) / ncol(chosen_obj)))

# Major cell classes
major_classes <- celltype_summary %>%
  mutate(
    major_class = case_when(
      grepl("Excitatory", celltype) ~ "Excitatory_Neurons",
      grepl("GABAergic", celltype) ~ "GABAergic_Neurons",
      celltype %in% c("Oligodendrocytes", "OPCs", "Astrocytes") ~ "Glia",
      celltype == "Microglia" ~ "Microglia",
      celltype %in% c("Endothelial", "Pericytes") ~ "Vascular",
      TRUE ~ "Other"
    )
  ) %>%
  group_by(major_class) %>%
  summarize(
    n_cells = sum(n_cells),
    percent = 100 * sum(n_cells) / sum(celltype_summary$n_cells)
  ) %>%
  arrange(desc(percent))

cat("\n=== Major Cell Classes ===\n")
print(major_classes)

# Save to file
write.csv(celltype_comparison,
         file.path(project_outdir, "celltype_composition_summary.csv"),
         row.names = FALSE)

cat(sprintf("\nSummary saved to: %s\n",
           file.path(project_outdir, "celltype_composition_summary.csv")))
