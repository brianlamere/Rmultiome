source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#Step 1-1: set up your space and list the options for sample names
init_project(
  random_seed = 42,
  use_cellbender = TRUE,
  use_scdblfinder = TRUE,
  doublet_rate_per_1000 = 8.0,
  doublet_rate_sd = 0.015,
  project_name = "opioid_hiv_multiome",
  genome_build = "hg38"
)
# init_project() # if resuming

# the contents of the below directories will entirely dictate what samples can be
# run through this pipeline.
list.files(path = rawdatadir)
list.files(path = cb_datadir)
# Define samples, likely as contents of rawdatadir but maybe as cb_datadir instead
# if you are using cellbender.  We will not autopopulate this because you need to know
# what samples you want to include, but it needs to be a subset or full set of rawdatadir.
# note if cb_datadir isn't a subset of exact names from rawdata dir, things will break later
samplelist <- c("LG05", "LG08", "LG22", "LG23", "LG25", "LG26",
                "LG300", "LG301", "LG31", "LG33", "LG38")

# Initialize pipeline1 settings
pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)

EnsDbAnnos <- loadannotations()

#step 1-2: pick your sample name, from the listing of files in rawdatadir
mysample <- "LG05"

# Step 1-3: Create base QC object
qc_obj <- base_qc_object(mysample, EnsDbAnnos, cb_report="display")

# Step 1-4: Generate QC plots (before trimming)
QCVlnA(qc_obj)
QCVlnR(qc_obj)
QCDensity_ATAC(qc_obj)
QCDensity_RNA(qc_obj)

# Step 1-5: Define local trimming settings for session
#my_trimming_settings <- trimming_settings[trimming_settings$sample == mysample, ]
my_trimming_settings <- list(
  sample = mysample,
  # ATAC counts
  min_nCount_ATAC = 1200,
  max_nCount_ATAC = 40000,
  # RNA counts
  min_nCount_RNA = 200,
  max_nCount_RNA = 30000,
  # Nucleosome signal (nss)
  min_nss = 0.2,
  max_nss = 1.5,
  # % mitochondrial
  max_percentMT = 8,
  # TSS enrichment
  min_TSS = 2.5, # was 2.5
  max_TSS = 9  # was 9
)

# Step 1-6: Apply trimming.  trimSample reads from trimming_settings and to ensure
# congruence between QC and the pipeline, we need to update the dataframe
# TODO: modify verify step to warn when a change is beyond a particular threshold
verify_pipeline1_settings(pipeline1_settings, my_trimming_settings)
pipeline1_settings <- update_pipeline1_settings(pipeline1_settings, my_trimming_settings)

trimmed_obj <- trimSample(qc_obj)

# Step 1-7: Generate QC plots after trimming
QCVlnA(trimmed_obj)
QCVlnR(trimmed_obj)
QCDensity_ATAC(trimmed_obj)
QCDensity_RNA(trimmed_obj)

# step 1-8: save trimming settings
#if the new plots made in step6 still show change is needed, start back at step4
#Once you are done, write to disk to save the settings
saveRDS(pipeline1_settings, pipeline1_settings_file)

###############KDE settings################################

#step 2-1: define local KDE settings for sample
#can re-print the plots, but your "before" is the same as the last set from 1D-trim above
# KDE filtering combine methods:
#   - "intersection": Cell must pass BOTH modality thresholds
#                     (more stringent, more trimming)
#   - "union":        Cell must pass EITHER modality threshold
#                     (more permissive, less trimming)

my_kde_settings <- list(
  sample = mysample,
  atac_percentile = 0.98,
  rna_percentile = 0.98,
  combine_method = "intersection"
)

#set 2-2: save to the settings dataframe
verify_pipeline1_settings(pipeline1_settings, my_kde_settings)
pipeline1_settings <- update_pipeline1_settings(pipeline1_settings, my_kde_settings)

#step 2-3: Visualize via contours
plot_kde_filter_contours(trimmed_obj, pipeline1_settings)

#optional: Step 2-4: Visualize the difference between union and intersection
plot_kde_filter_combine_compare_atac(trimmed_obj, pipeline1_settings)
plot_kde_filter_combine_compare_rna(trimmed_obj, pipeline1_settings)

#Repeat steps 2-1 to 2-4 as desired until you find the percentile and combine
# method you want to use.

#step 2-5: save the KDE trimming setting
saveRDS(pipeline1_settings, pipeline1_settings_file)

cat("\n=== Applying KDE Trimming ===\n")
kde_obj <- kdeTrimSample(trimmed_obj, qc_report = TRUE)

#Step 3: scDblFinder
# Record cell count after KDE (for doublet rate calculation)
n_cells <- ncol(kde_obj)
cat(sprintf("Cells after KDE trim: %d\n", n_cells))

# Calculate expected doublet rate for this sample
expected_dbr <- (n_cells / 1000) * (doublet_rate_per_1000 / 10) * (n_cells / 100)
cat(sprintf("Expected doublets: %.1f (%.2f%% of %d cells)\n",
           expected_dbr, doublet_rate_pct, n_cells))

################REMOVE ONCE RUN_QC IS REPEATED############
# CHEAT: start from here and loop over samples, since you already have QC settings
samplelist <- c("LG05", "LG08", "LG22", "LG23", "LG25", "LG26",
                "LG300", "LG301", "LG31", "LG33", "LG38")
pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)
EnsDbAnnos <- loadannotations()
mysample <- "LG05"
qc_obj <- base_qc_object(mysample, EnsDbAnnos, cb_report="display")
trimmed_obj <- trimSample(qc_obj)
kde_obj <- kdeTrimSample(trimmed_obj, qc_report = TRUE)
##########################################################

# Store in pipeline1_settings
my_doublet_settings <- list(
  sample = mysample,
  n_cells_after_kde = n_cells,
  expected_dbr = expected_dbr
)

verify_pipeline1_settings(pipeline1_settings, my_doublet_settings)
pipeline1_settings <- update_pipeline1_settings(pipeline1_settings, my_doublet_settings)

# === STEP 4: Test scDblFinder ===
if (use_scdblfinder) {
  doublet_obj <- doubletRemoveSample(kde_obj, qc_report = TRUE)

  cat("\n=== Evaluation ===\n")
  cat("Review the results:\n")
  cat("  - Is the detected rate close to expected?\n")
  cat("  - Is the score distribution bimodal?\n")
  cat("  - Do the doublets make biological sense?\n\n")
  cat("If results look reasonable, settings will be used in run_pipeline1.R\n\n")
}

#step 4-1: start almost everything over
#protect yourself from stepping on yourself
rm(mysample,qc_obj,my_trimming_settings,trimmed_obj,my_kde_settings)
#Stop at this point, then repeat steps 1-2 to 3-1 for each sample, starting by
# changing the "mysample" setting and looping back to here

#step 4-1: compare cellbender reports
#compare all the cellbender reports as an aggregated list.  Call allows for
# "samplelist=" but defaults to samplelist=trimming_settings$sample.
# Example use: samplelist=c("LG05","LG08") as argument in call
compare_cellbender_reports("qc")

#if everything looks good to this point, you're ready to use run_pipeline1.R
