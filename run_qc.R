source("/projects/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

trimming_settings <- init_trimming_settings(trimming_settings_file)

EnsDbAnnos <- loadannotations()

mysample <- "LG05"

# Step 1-2: Create base QC object
qc_obj <- base_qc_object(mysample, EnsDbAnnos)

# Step 1-3: Generate QC plots (before trimming)
QCVlnA(qc_obj)
QCVlnR(qc_obj)
QCDensity_ATAC(qc_obj)
QCDensity_RNA(qc_obj)

# Step 1-4: Define local trimming settings for session
my_trimming_settings <- list(
  sample = mysample,
  # ATAC counts
  min_nCount_ATAC = 1000,
  max_nCount_ATAC = 40000,
  # RNA counts
  min_nCount_RNA = 200,
  max_nCount_RNA = 25000,
  # Nucleosome signal (nss)
  min_nss = 0.4,
  max_nss = 1.5,
  # % mitochondrial
  max_percentMT = 8,
  # TSS enrichment
  min_TSS = 2.5,
  max_TSS = 10
)

# Step 1-5: Apply trimming.  trimSample reads from trimming_settings and to ensure
# congruence between QC and the pipeline, we need to update the dataframe
verify_trimming_settings(trimming_settings, my_trimming_settings)
trimming_settings <- update_trimming_settings(trimming_settings, my_trimming_settings)
trimmed_obj <- trimSample(qc_obj)

# Step 1-6: Generate QC plots after trimming
QCDensity_ATAC(trimmed_obj)
QCDensity_RNA(trimmed_obj)
QCVlnA(trimmed_obj)
QCVlnR(trimmed_obj)

#if the new plots made in step6 still show change is needed, start back at step4
#Once you are done, write to disk to save the settings
saveRDS(trimming_settings, trimming_settings_file)

###############KDE settings################################

#step 2-1: read any existing/current settings
kde_settings <- init_kde_settings(kde_settings_file)

#we can re-print the plots, but your "before" is the same as the last set of
# 4 from 1D trimming, above

#step 2-2: define local KDE settings for sample
my_kde_settings <- list(
  sample = mysample,
  atac_percentile = 0.96,
  rna_percentile = 0.96,
  combine_method = "intersection"
)

#set 2-3: save to the settings dataframe
verify_kde_settings(kde_settings, my_kde_settings)
kde_settings <- update_kde_settings(kde_settings, my_kde_settings)

#step 2-3: Visualize 
#there should be happy visuals here that have contour lines at the percentages
#for now, a basic eyeballing from previous is necessary

kde_obj <- kdeTrimSample(trimmed_obj)
