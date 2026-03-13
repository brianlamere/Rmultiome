# Edit this file to match your local environment.
# The settings in this file are very specific to installation and project!

# Working directory that has everything in it.
project_base_dir <- "/projects1/opioid"

# If TRUE, base_object() will replace the CellRanger RNA counts with CellBender-corrected
# RNA.  This requires that the CellBender output exists for each sample; otherwise
# the pipeline will stop.
use_cellbender <- TRUE

# if you downloaded a cellmarker csv from here:
# http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarkerSearch.jsp?index_species=Human&index_tissue=Brain
# define location below
CellMarker_file <- file.path(project_base_dir, "references/Cellmarker_Human_Brain.csv")

##################################################################################
# If following guide precisely, nothing below to end of file would need to change.
##################################################################################
# Main 10X multiome h5 file (RNA + ATAC peaks) from CellRanger ARC
# This file is still used even when CellBender is enabled, because it contains Peaks.
# expected at:  paste(rawdatadir, samplename, h5filename, sep = "/")
h5filename <- "filtered_feature_bc_matrix.h5"

# CellBender-corrected RNA file (per-sample)
# expected at:  paste(cb_datadir, samplename, cellbender_rna_h5filename, sep = "/")
cellbender_rna_h5filename <- "cellbender_gex_filtered.h5"

# ATAC fragments file name
atacfilename <- "atac_fragments.tsv.gz"

# Project source code directory
Rmultiome_path <- file.path(project_base_dir, "Rmultiome")

# Project raw data directory
rawdatadir <- file.path(project_base_dir, "rawdata")

# Project cellbender data directory
cb_datadir <- file.path(project_base_dir, "cb_data")

# Vault of rds files saved at different milestones
rdsdir <- file.path(project_base_dir, "vault")

# Project output directory for reports, settings, etc.
project_outdir <- file.path(project_base_dir, "project_export")

# Temporary directory (ok to delete; used for scratch outputs during QC, etc.)
tmpfiledir <- file.path(project_base_dir, "tmp")

# Reference files directory
referencedir <- file.path(project_base_dir, "references")

# Directory where pipeline1 writes per-sample CellBender merge reports (one CSV per sample)
cellbender_report_dir <- file.path(project_outdir, "cellbender_merge_reports")

# Project 1D settings file
trimming_settings_file <- file.path(project_outdir, "trimming_settings.Rds")

# Project KDE settings file
kde_settings_file <- file.path(project_outdir, "kde_settings.Rds")

#Project Cluster settings file: dims, KNN, resolution
cluster_settings_file <- file.path(project_outdir, "cluster_settings.Rds")

#Project Celltype Mapping settings file
celltype_settings_file <- file.path(project_outdir, "celltype_settings.Rds")

#Project Celltype marker panel file.  CONTENTS ARE VERY SPECIFIC TO EACH PROJECT
marker_panel_file <- file.path(referencedir, "stress_marker_panel_documented.csv")

#Project Harmony Settings file
harmony_settings_file <- file.path(project_outdir, "harmony_settings.Rds")

#directory to store temporary files used during parameter sweep
sweep_dir <- file.path(tmpfiledir, "param_sweep")

# Multiomic, needing to have two different reductions available
# Defining here only so that they aren't defined in a script
reduction.save.RNA = "harmony_RNA"
reduction.save.ATAC = "harmony_ATAC"
