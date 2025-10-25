# Edit this file to match your local environment.
# The settings in this file are very specific to installation and project!

# Working directory that has everything in it.
project_base_dir <- "/projects/opioid"

# If following guide precisely, nothing below to end of file would need to change.

# Project source code directory
Rmultiome_path <- file.path(project_base_dir, "Rmultiome")

# Project raw data directory
rawdatadir <- file.path(project_base_dir, "rawdata")

# Vault of rds files saved at different milestones
rdsdir <- file.path(project_base_dir, "vault")

# Project output directory for reports, settings, etc.
project_outdir <- file.path(project_base_dir, "project_export")

# Project 1D settings file
trimming_settings_file <- file.path(project_outdir, "trimming_settings.Rds")

# Project KDE settings file
kde_settings_file <- file.path(project_outdir, "kde_settings.Rds")

# Reference files directory
referencedir <- file.path(project_base_dir, "references")

# Main 10X h5 file (RNA)
h5filename <- "filtered_feature_bc_matrix.h5"

# ATAC fragments file name
atacfilename <- "atac_fragments.tsv.gz"



