# === Generate lymphocyte investigation report ===
cat("\n=== Generating lymphocyte investigation report ===\n")

# Create report content
report_content <- sprintf("
# Lymphocyte Investigation Report
**Date:** %s  
**Dataset:** Frontal cortex nuclei, HIV+ and opioid-exposed subjects  
**Analysis:** Cell type identification via marker-based typing  

---

## Objective
Investigate the presence of lymphocytes in the dataset using established markers from the literature (Zillich et al. 2025, Nature Communications).

---

## Methods

### Markers tested:
- **CD96**: TACTILE, NK cell and T-cell marker
- **CD3D**: T-cell receptor component (definitive T-cell marker)
- **IL2RG**: Common gamma chain, lymphocyte marker

### Analysis approach:
1. **Programmatic cell typing**: Used `identify_all_celltypes()` function with Zillich et al. markers
2. **Visual inspection**: Generated DotPlot to visualize marker expression across clusters
3. **Quantitative assessment**: Extracted numerical values for percent expression and average expression per cluster

---

## Results

### 1. Programmatic typing result:
No markers (CD96, CD3D, IL2RG) were present in the top differentially expressed genes (FindAllMarkers output with thresholds: min.pct = 0.25, logfc.threshold = 0.25).

### 2. Visual inspection (DotPlot):
- **CD96**: Faint scattered expression across multiple clusters
- **CD3D**: Essentially absent across all clusters
- **IL2RG**: Essentially absent across all clusters

### 3. Quantitative assessment:

**Expression levels by cluster:**

| Cluster | CD96 pct | CD96 avg | CD3D pct | CD3D avg | IL2RG pct | IL2RG avg |
|---------|----------|----------|----------|----------|-----------|-----------|
| 0       | 0.16%%    | -0.76    | 0.01%%    | NA       | 0.01%%     | NA        |
| 1       | 3.43%%    | 0.48     | 0.01%%    | NA       | 0.01%%     | NA        |
| 2       | 0.29%%    | -0.74    | 0.00%%    | NA       | 0.04%%     | NA        |
| 3       | 1.17%%    | 2.50     | 0.11%%    | NA       | 0.34%%     | NA        |
| 4       | 0.98%%    | -0.33    | 0.01%%    | NA       | 0.00%%     | NA        |
| 5       | 2.27%%    | 2.03     | 0.00%%    | NA       | 0.00%%     | NA        |
| 6       | 0.22%%    | -0.62    | 0.00%%    | NA       | 0.08%%     | NA        |
| 7       | 0.91%%    | -0.39    | 0.00%%    | NA       | 0.02%%     | NA        |
| 8       | 2.17%%    | -0.30    | 0.05%%    | NA       | 0.00%%     | NA        |
| 9       | 2.36%%    | 0.09     | 0.04%%    | NA       | 0.04%%     | NA        |
| 10      | 2.85%%    | 0.15     | 0.04%%    | NA       | 0.00%%     | NA        |
| 11      | 0.65%%    | -0.71    | 0.00%%    | NA       | 0.00%%     | NA        |
| 12      | 0.83%%    | 0.19     | 0.00%%    | NA       | 0.14%%     | NA        |
| 13      | 0.95%%    | -0.12    | 0.05%%    | NA       | 0.10%%     | NA        |
| 14      | 3.12%%    | 0.79     | 0.00%%    | NA       | 0.00%%     | NA        |
| 15      | 0.96%%    | -0.48    | 0.11%%    | NA       | 0.00%%     | NA        |
| 16      | 1.15%%    | -0.44    | 0.00%%    | NA       | 0.00%%     | NA        |
| 17      | 0.33%%    | -0.81    | 0.00%%    | NA       | 0.00%%     | NA        |

**Key observations:**
- **CD96**: Maximum 3.43%% of cells in any cluster, mostly with negative average expression (below dataset mean = background)
- **CD3D**: Maximum 0.11%% of cells, average expression = NA (too few cells expressing for calculation)
- **IL2RG**: Maximum 0.34%% of cells, average expression = NA

**Expected for true lymphocyte population:**
- CD3D/IL2RG: >50%% of cells expressing
- CD96: >30%% of cells expressing
- Average expression: Strongly positive

---

## Interpretation

### Why lymphocytes are absent:

#### 1. **Biological reasons:**
- **Tissue type**: Frontal cortex parenchyma is relatively immune-privileged
- **Lymphocyte localization**: Brain lymphocytes primarily found in meninges, choroid plexus, and perivascular spaces (not sampled in cortical parenchyma)
- **Disease state**: Chronic HIV and opioid exposure without acute encephalitis
- **Expected frequency**: <0.1%% of cells in healthy cortex parenchyma

#### 2. **Technical reasons:**
- **Protocol**: Nuclei extraction (not whole cells) may preferentially lose small, fragile immune cells
- **Sample type**: Post-mortem tissue; mobile immune cells may migrate out during PMI
- **Enrichment**: No immune cell enrichment performed

#### 3. **Comparison to reference:**
- **Zillich et al. (cocaine paper)**: Caudate nucleus (more vascular, different immune environment) with possible whole-cell protocol
- **Our data**: Frontal cortex with nuclei extraction

---

## Conclusion

**No evidence of lymphocytes in this dataset.**

This finding is:
- ✓ Biologically plausible (cortical parenchyma in non-inflamed brain)
- ✓ Technically expected (nuclei prep of post-mortem tissue)
- ✓ Consistent across multiple analysis methods

**Recommendation:** Discontinue searching for lymphocytes. Focus cell typing efforts on the 6 other major brain cell types that were successfully identified:
1. Oligodendrocytes
2. Oligodendrocyte precursor cells (OPCs)
3. Astrocytes
4. Microglia
5. GABAergic neurons
6. Endothelial cells

---

## Methods Note
For publication: \"Lymphocyte markers (CD96, CD3D, IL2RG) from Zillich et al. 2025 were evaluated but not detected (CD3D <0.11%%, IL2RG <0.35%% of cells across all clusters), consistent with the relatively immune-privileged nature of cortical parenchyma in the absence of acute neuroinflammation.\"

---

**Analysis conducted:** %s
**Report generated:** %s
",
Sys.Date(),
Sys.time(),
Sys.time()
)

# Write report
report_file <- file.path(project_outdir, "Lymphocyte_Investigation_Report.md")
writeLines(report_content, report_file)

cat(sprintf("Report saved: %s\n", report_file))

# Also save the numerical data
write.csv(lymph_table,
         file.path(project_outdir, "lymphocyte_marker_expression_data.csv"),
         row.names = FALSE)

cat("Supporting data saved: lymphocyte_marker_expression_data.csv\n")
