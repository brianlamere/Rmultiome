source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
init_project()

pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
cluster_settings <- read_cluster_settings(cluster_settings_file)
celltype_settings <- read_celltype_settings(celltype_settings_file)
harmony_settings <- read_harmony_settings(harmony_settings_file)


EnsDbAnnos <- loadannotations()

obj_assigned <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))

# =============================================================================
# export_for_sarah.R
# Run AFTER typing_validation.R has loaded obj_assigned (annotated_obj_final.rds)
# Exports NATIVE celltype labels (11 types: Excitatory_Neurons, GABAergic_LAMP5/
# PVALB/SST/VIP, Astrocytes, Oligodendrocytes, OPCs, Microglia_Macrophages,
# Endothelial, Pericytes) - no collapsing to Sarah's Exc/Inh/Astro/... panel.
# She collapses subtypes on her end as needed.
# Produces:
#   1. RNA_avg_expr_pct_by_celltype.csv          (gene x celltype, + pct expressing)
#   2. ATAC_peak_accessibility_by_celltype.tsv   (peak_id,chr,start,end,<celltypes>)
#   3. ATAC_peak_pct_accessible_by_celltype.tsv  (same peaks, % nuclei accessible)
# =============================================================================

library(Seurat)
library(Signac)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# 0. No collapsing - export native celltype labels as-is.
#    Sarah's requested panel (Exc/Inh/Astro/Oligo/OPC/Micro/Endo) doesn't match
#    this object's actual typing 1:1 - GABAergic splits into 4 stable subtypes,
#    and Pericytes/Endothelial are distinct populations with no combined "Endo"
#    bucket here. Decision (per Sarah, 2026-06-28): send native labels, let her
#    collapse however she wants on her end rather than deciding it for her here.
# -----------------------------------------------------------------------------
obj_assigned$celltype_major <- as.character(obj_assigned$celltype)

# Known Seurat bug (satijalab/seurat#7893): AverageExpression silently replaces
# "_" with "-" in group.by labels when building output column names. Rather
# than chase that sanitization through every downstream call, confirm no
# label already contains a literal "-" (so the substitution is unambiguously
# reversible), then reverse it once, immediately after each AverageExpression
# call, so everything downstream operates on the real metadata values.
stopifnot(!any(grepl("-", unique(obj_assigned$celltype_major))))

# =============================================================================
# 1. RNA: average normalized expression + percent expressing, by major cell type
# =============================================================================
DefaultAssay(obj_assigned) <- "RNA"

# Layers are unjoined at this checkpoint (still split per-sample: counts.LG05, etc)
# pipeline2 joins them at Step 5 before DE - do the same here.
obj_assigned <- JoinLayers(obj_assigned)

avg_expr <- AverageExpression(
  obj_assigned,
  assays    = "RNA",
  group.by  = "celltype_major",
  layer     = "data"   # log-normalized data slot, standard "normalized expression"
)$RNA

# Reverse Seurat's "_" -> "-" sanitization (see guard above re: safety of this).
colnames(avg_expr) <- gsub("-", "_", colnames(avg_expr))
col_order <- sort(colnames(avg_expr))

avg_expr_df <- as.data.frame(avg_expr) %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, all_of(col_order))

# Percent of cells expressing each gene per cell type (count > 0 in raw counts)
pct_expr_by_group <- function(obj, group_col, assay = "RNA") {
  groups <- sort(unique(obj[[group_col, drop = TRUE]]))
  counts_mat <- GetAssayData(obj, assay = assay, layer = "counts")
  pct_list <- lapply(groups, function(g) {
    cells_in_group <- which(obj[[group_col, drop = TRUE]] == g)
    Matrix::rowMeans(counts_mat[, cells_in_group, drop = FALSE] > 0) * 100
  })
  names(pct_list) <- groups
  as.data.frame(pct_list)
}

pct_expr_df <- pct_expr_by_group(obj_assigned, "celltype_major") %>%
  tibble::rownames_to_column("gene") %>%
  dplyr::select(gene, all_of(col_order))

names(pct_expr_df)[-1] <- paste0(names(pct_expr_df)[-1], "_pct")

rna_export <- avg_expr_df %>%
  left_join(pct_expr_df, by = "gene")

write.csv(rna_export, file.path(project_export, "RNA_avg_expr_pct_by_celltype.csv"),
          row.names = FALSE)
cat(sprintf("Wrote RNA export: %d genes x %d celltypes (+pct columns)\n",
            nrow(rna_export), length(col_order)))

# =============================================================================
# 2. ATAC: peak-level average accessibility by major cell type (required file)
# =============================================================================
DefaultAssay(obj_assigned) <- "ATAC"

avg_acc <- AverageExpression(
  obj_assigned,
  assays    = "ATAC",
  group.by  = "celltype_major",
  layer     = "data"   # Signac normalized accessibility; will fall back to "counts" if no "data" layer exists - check warning on run
)$ATAC

# Reverse Seurat's "_" -> "-" sanitization (same as RNA block above).
colnames(avg_acc) <- gsub("-", "_", colnames(avg_acc))

# Sanity check: confirm ATAC's AverageExpression produced the same celltype
# set as RNA's, now that both are back in underscore form.
if (!setequal(colnames(avg_acc), col_order)) {
  stop(sprintf(
    "ATAC celltype columns don't match RNA's:\n  RNA: %s\n  ATAC: %s",
    paste(col_order, collapse = ", "),
    paste(sort(colnames(avg_acc)), collapse = ", ")
  ))
}

peak_ids <- rownames(avg_acc)

# chr-start-end peak names -> split columns
peak_coords <- data.frame(
  peak_id = peak_ids,
  chr     = sub("^(chr[^-]+)-.*$", "\\1", peak_ids),
  start   = as.integer(sub("^chr[^-]+-([0-9]+)-[0-9]+$", "\\1", peak_ids)),
  end     = as.integer(sub("^chr[^-]+-[0-9]+-([0-9]+)$", "\\1", peak_ids))
)

atac_export <- peak_coords %>%
  bind_cols(as.data.frame(avg_acc) %>% dplyr::select(all_of(col_order)))

write.table(atac_export,
            file.path(project_export, "ATAC_peak_accessibility_by_celltype.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote ATAC accessibility export: %d peaks x %d celltypes\n",
            nrow(atac_export), length(col_order)))

# =============================================================================
# 3. ATAC: percent of nuclei accessible per peak per cell type (optional file)
# =============================================================================
pct_acc_by_group <- function(obj, group_col, assay = "ATAC") {
  groups <- sort(unique(obj[[group_col, drop = TRUE]]))
  counts_mat <- GetAssayData(obj, assay = assay, layer = "counts")
  pct_list <- lapply(groups, function(g) {
    cells_in_group <- which(obj[[group_col, drop = TRUE]] == g)
    Matrix::rowMeans(counts_mat[, cells_in_group, drop = FALSE] > 0) * 100
  })
  names(pct_list) <- groups
  as.data.frame(pct_list)
}

pct_acc_df <- pct_acc_by_group(obj_assigned, "celltype_major") %>%
  dplyr::select(all_of(col_order))
names(pct_acc_df) <- paste0(col_order, "_pct")

atac_pct_export <- peak_coords %>%
  bind_cols(pct_acc_df)

write.table(atac_pct_export,
            file.path(project_export, "ATAC_peak_pct_accessible_by_celltype.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Wrote ATAC %% accessible export: %d peaks x %d celltypes\n",
            nrow(atac_pct_export), length(col_order)))

cat("\nDone. Files in: ", project_export, "\n")
