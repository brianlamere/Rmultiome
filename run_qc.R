source("/projects/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))
source(project_settings_file)

EnsDbAnnos <- loadannotations()

mysample <- "LG05"

# Step 2: Create base QC object
qc_obj <- base_qc_object(mysample, EnsDbAnnos)

# Step 3: Generate QC plots (before trimming)
QCVlnA(qc_obj)
QCVlnR(qc_obj)
QCDensity_ATAC(qc_obj)
QCDensity_RNA(qc_obj)

# Step 4: Define local trimming settings for session
my_trimming_settings <- list(
  sample = "LG05",
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

# Step 5: Apply trimming
trimmed_obj <- trimSample(qc_obj, trimming_settings = my_trimming_settings)

# Step 6: Generate QC plots after trimming
QCDensity_ATAC(trimmed_obj)
QCDensity_RNA(trimmed_obj)
QCVlnA(trimmed_obj)
QCVlnR(trimmed_obj)

verify_trimming_settings_file_changes("/projects/opioid/Rmultiome/settings.R", my_trimming_settings)
update_trimming_settings_in_file("/projects/opioid/Rmultiome/settings.R", my_trimming_settings)