source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#Step 1-1: set up your space and list the options for sample names
init_project()

trimming_settings <- init_trimming_settings(trimming_settings_file)
kde_settings <- init_kde_settings(kde_settings_file)

EnsDbAnnos <- loadannotations()

list.files(path = rawdatadir)

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
  min_nCount_ATAC = 500,
  max_nCount_ATAC = 50000,
  # RNA counts
  min_nCount_RNA = 200,
  max_nCount_RNA = 30000,
  # Nucleosome signal (nss)
  min_nss = 0.4,
  max_nss = 4,
  # % mitochondrial
  max_percentMT = 8,
  # TSS enrichment
  min_TSS = 2,
  max_TSS = 10
)

# Step 1-6: Apply trimming.  trimSample reads from trimming_settings and to ensure
# congruence between QC and the pipeline, we need to update the dataframe
verify_trimming_settings(trimming_settings, my_trimming_settings)
trimming_settings <- update_trimming_settings(trimming_settings, my_trimming_settings)
trimmed_obj <- trimSample(qc_obj)

# Step 1-7: Generate QC plots after trimming
QCVlnA(trimmed_obj)
QCVlnR(trimmed_obj)
QCDensity_ATAC(trimmed_obj)
QCDensity_RNA(trimmed_obj)

# step 1-8: save trimming settings
#if the new plots made in step6 still show change is needed, start back at step4
#Once you are done, write to disk to save the settings
saveRDS(trimming_settings, trimming_settings_file)

###############KDE settings################################

#step 2-1: define local KDE settings for sample
#we can re-print the plots, but your "before" is the same as the last set of
# 4 from 1D trimming, above
my_kde_settings <- list(
  sample = mysample,
  atac_percentile = 0.98,
  rna_percentile = 0.98,
  combine_method = "intersection"
)

#set 2-2: save to the settings dataframe
verify_kde_settings(kde_settings, my_kde_settings)
kde_settings <- update_kde_settings(kde_settings, my_kde_settings)

#step 2-3: Visualize via contours
plot_kde_filter_contours(trimmed_obj, kde_settings)

#optional: Step 2-4: Visualize the difference between union and intersection
plot_kde_filter_combine_compare_atac(trimmed_obj, kde_settings)
plot_kde_filter_combine_compare_rna(trimmed_obj, kde_settings)

#Repeat steps 2-1 to 2-4 as desired until you find the percentile and combine
# method you want to use.

#step 2-6: save the KDE trimming setting
saveRDS(kde_settings, kde_settings_file)

#step 3-1: start almost everything over
#protect yourself from stepping on yourself
rm(mysample,qc_obj,my_trimming_settings,trimmed_obj,kde_settings,my_kde_settings)
#Stop at this point, then repeat steps 1-2 to 3-1 for each sample, starting by
# changing the "mysample" setting and looping back to here

#step 4-1: compare cellbender reports
#compare all the cellbender reports as an aggregated list.  Call allows for
# "samplelist=" but defaults to samplelist=trimming_settings$sample.
# Example use: samplelist=c("LG05","LG08") as argument in call
compare_cellbender_reports("qc")

#if everything looks good to this point, you're ready to use run_pipeline1.R
