source("/projects1/opioid2/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#Step 1-1: set up your space and list the options for sample names
init_project(
  random_seed = 42,
  use_cellbender = TRUE,
  use_scdblfinder = TRUE,
  doublet_rate_per_1000 = 0.8,
  doublet_rate_sd = 0.015,
  project_name = "opioid_hiv_multiome",
  genome_build = "hg38"
)
# init_project() # if resuming

hgd(port=8777, token=FALSE)

# the contents of the below directories will entirely dictate what samples can be
# run through this pipeline.
list.files(path = cra_outdir)
list.files(path = cb_datadir)
# Define samples, likely as contents of cra_outdir but maybe as cb_datadir instead
# if you are using cellbender.  We will not autopopulate this because you need to know
# what samples you want to include, but it needs to be a subset or full set of cra_outdir.
# note if cb_datadir isn't a subset of exact names from cr-arc out dir, things will break later
samplelist <- c("LG22", "LG26", "LG300", "LG301", "LG31", "LG38")

# Initialize pipeline1 settings
pipeline1_settings <- init_pipeline1_settings(pipeline1_settings_file)

EnsDbAnnos <- loadannotations()

#step 1-2: pick your sample name, from the listing of files in cra_outdir
mysample <- "LG22"

# Step 1-3: Create base QC object
qc_obj <- base_qc_object(mysample, EnsDbAnnos, cb_report="display")

ncol(qc_obj[["ATAC"]])

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
  min_nCount_ATAC = 200,
  max_nCount_ATAC = 10000,
  # RNA counts
  min_nCount_RNA = 125,
  max_nCount_RNA = 3000,
  # Nucleosome signal (nss)
  min_nss = 0.1,
  max_nss = 1.4,
  # % mitochondrial
  max_percentMT = 30,
  # TSS enrichment
  min_TSS = .75, # 
  max_TSS = 7.5  # 
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
  atac_percentile = 0.965,
  rna_percentile = 0.965,
  combine_method = "intersection"
)

######set 2-2: save to the settings dataframe
verify_pipeline1_settings(pipeline1_settings, my_kde_settings)
pipeline1_settings <- update_pipeline1_settings(pipeline1_settings, my_kde_settings)

######step 2-3: Visualize via contours
plot_kde_filter_contours(trimmed_obj, pipeline1_settings)

######optional: Step 2-4: Visualize the difference between union and intersection
plot_kde_filter_combine_compare_atac(trimmed_obj, pipeline1_settings)
plot_kde_filter_combine_compare_rna(trimmed_obj, pipeline1_settings)

######Repeat steps 2-1 to 2-4 as desired until you find the percentile and combine
# method you want to use.

######step 2-5: save the KDE trimming setting
saveRDS(pipeline1_settings, pipeline1_settings_file)

cat("\n=== Applying KDE Trimming ===\n")
kde_obj <- kdeTrimSample(trimmed_obj, qc_report = TRUE)

######Step 3: scDblFinder
n_cells <- ncol(kde_obj)
cat(sprintf("Cells after KDE trim: %d\n", n_cells))

# Calculate expected doublet rate for this sample
# doublet rate is based on physical cells captured, so use pre-QC
n_true_cells <- ncol(qc_obj)
expected_dbr <- ((n_true_cells / 1000) * doublet_rate_per_1000)
#example: if you have need to alter n_cells due to cellbender, etc, do so like such
#expected_dbr <- ((56664 / 1000) * (doublet_rate_per_1000))
cat(sprintf("Expected doublets: %.1f (%.2f%% of %d cells)\n",
           (n_cells * expected_dbr / 100), expected_dbr, n_cells))

# Store in pipeline1_settings
my_doublet_settings <- list(
  sample = mysample,
  n_cells_scdbl = n_cells,
  expected_dbr = expected_dbr
)

verify_pipeline1_settings(pipeline1_settings, my_doublet_settings)
pipeline1_settings <- update_pipeline1_settings(pipeline1_settings, my_doublet_settings)
saveRDS(pipeline1_settings, pipeline1_settings_file)

###### === STEP 4: Test scDblFinder ===
if (use_scdblfinder) {
  result <- doubletRemoveSample(kde_obj, qc_report = TRUE)
  doublet_obj <- result$obj
  doublet_stats <- result$stat

  cat("\n=== Evaluation ===\n")
  cat("Review the results:\n")
  cat("  - Is the detected rate close to expected?\n")
  cat("  - Is the score distribution bimodal?\n")
  cat("  - Do the doublets make biological sense?\n\n")
  cat("If results look reasonable, settings will be used in run_pipeline1.R\n\n")
}

######step 5: start almost everything over
#protect yourself from stepping on yourself
rm(mysample,qc_obj,trimmed_obj,my_kde_settings,my_doublet_settings,
 n_cells, expected_dbr, doublet_obj)
#Stop at this point, then repeat steps 1-2 to 5 for each sample, starting by
# changing the "mysample" setting and looping back to here

######step 4-1: compare cellbender reports
#compare all the cellbender reports as an aggregated list.  Call allows for
# "samplelist=" but defaults to samplelist=trimming_settings$sample.
# Example use: samplelist=c("LG05","LG08") as argument in call
compare_cellbender_reports(samplelist=samplelist)

#if everything looks good to this point, you're ready to use run_pipeline1.R
