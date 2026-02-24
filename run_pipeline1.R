source("/projects1/opioid/Rmultiome/system_settings.R")
source(file.path(Rmultiome_path, "Rmultiome-main.R"))

#The intent of this script is to be able to be run after you've selected all
#settings, done all QC, and this script will re-create your DE/DA/etc results
#from the raw output of cellranger-atac.  Step1 is not this file!  Step1 is to
#use run_base_QC.R and establish your pipeline1 settings, as stored in settings.R
#If you use run_all.R multiple times with settings tailored for your data, and
#the end results change, then you either have a data stability problem you haven't
#resolved yet, your 1D/KDE trimming is not correct, or you need to select different
#dims/knn/resolution during FMMN and clustering.

# If your biological goal is data that is "good enough" to inform another process,
# and the observed changes between runs fall within a threshold acceptable for your
# research context, then some small degree of variability may be tolerable.
# However, if your intent is a fully deterministic pipeline, any change in results
# between runs is a strong indicator of unresolved issues in data stability,
# preprocessing, or parameter selection.

standard_chroms <- paste0("chr", c(1:22, "X", "Y"))

trimming_settings <- read_trimming_settings(trimming_settings_file)
kde_settings <- init_kde_settings(kde_settings_file)

samplelist <- trimming_settings$sample

EnsDbAnnos <- loadannotations()

#########copy FROM here, to run in an IDE as a full block
for (sample in samplelist) {
  cat(sprintf("\nProcessing sample: %s\n", sample))
  
  # STEP 1: RAW (import)
  base_path <- get_rds_path(sample, "base")
  if (!file.exists(base_path)) {
    base_obj <- base_object(sample)
    print("Adding chromosome mapping information to ATAC assay.")
    base_obj <- chromosome_mapping(base_obj, rna_annos = EnsDbAnnos)
    #base_obj <- update_provenance(base_obj, "raw_import")
    saveRDS(base_obj, base_path)
  } else {
    print("Reading previous base file from vault\n")
    base_obj <- readRDS(base_path)
  }
  
  # STEP 2: 1D TRIM
  trimmed_path <- get_rds_path(sample, "trimmed")
  if (!file.exists(trimmed_path)) {
    trim_obj <- base_obj
    print("Calculating NucleosomeSignal and TSSEnrichment for ATAC data.")
    DefaultAssay(trim_obj) <- "ATAC"
    trim_obj <- NucleosomeSignal(trim_obj)
    trim_obj <- TSSEnrichment(trim_obj)
    print("Trimming based on QC")
    trim_obj <- trimSample(trim_obj)
    #trim_obj <- update_provenance(trim_obj, "trim_data")
    saveRDS(trim_obj, trimmed_path)
  } else {
    print("Reading previous trim file from vault\n")
    trim_obj <- readRDS(trimmed_path)
  }
  
  # STEP 3: 2D TRIM
  kde_path <- get_rds_path(sample, "kdetrim")
  if (!file.exists(kde_path)) {
    kde_obj <- trim_obj
    print("Doing n-Dimensional KDE trimming.")
    cat("Cells before trimming:", nrow(kde_obj@meta.data), "\n")
    kde_obj <- kdeTrimSample(kde_obj, qc_report = FALSE)
    cat("Cells after KDE trimming:", nrow(kde_obj@meta.data), "\n")
    saveRDS(kde_obj, kde_path)
  } else {
    print("Reading previous kdetrim file from vault\n")
    kde_obj <- readRDS(kde_path)
  }
    
  # STEP 4: preRNA
  preRNA_path <- get_rds_path(sample, "preRNA")
  if (!file.exists(preRNA_path)) {
    obj <- kde_obj
    DefaultAssay(obj) <- "RNA"
    obj <- NormalizeData(obj)
    obj <- FindVariableFeatures(obj, selection.method = "vst",
                                nfeatures = 2500)
    #obj <- update_provenance(obj, "pre-merge_rna")
    saveRDS(obj, preRNA_path)
  } else {
    print("Reading previous preRNA file from vault\n")
    obj <- readRDS(preRNA_path)
  }
  
  # STEP 5: preATAC/pipeline1
  pipeline1_path <- get_rds_path(sample, "pipeline1")
  if (!file.exists(pipeline1_path)) {
    DefaultAssay(obj) <- "ATAC"
    obj <- RunTFIDF(obj)
    obj <- FindTopFeatures(obj)
    #obj <- update_provenance(obj, "pre-merge_atac")
    saveRDS(obj, pipeline1_path)
  } else {
    print("Files already present for this sample for pipeline1\n")
    obj <- readRDS(pipeline1_path)
  }
  
  # Step 6: Report and cleanup
  generate_sample_report(
    sample = sample,
    base_obj = base_obj,
    trim_obj = trim_obj,
    kde_obj = kde_obj,
    pipeline1_path = pipeline1_path,
    trimming_settings = trimming_settings,
    kde_settings = kde_settings
  )
  
  gc() #R is obnoxious
  cat(sprintf("Sample %s completed successfully.\n", sample))
}
#############copy TO here, to run in an IDE as a full block

merged_data <- merge_sample_objects(samplelist)

saveRDS(merged_data, file.path(rdsdir,"merged_samples.rds"))
#merged_data <- readRDS("/projects/opioid/vault/merged_samples95.rds")

#post-merge RNA modality
merged_data <- post_merge_rna(merged_data)

#post-merge ATAC modality
merged_data <- post_merge_atac(merged_data)

DefaultAssay(merged_data) <- "RNA"

#time for harmony and FindMultiModalNeighbors
harmony_obj <- harmonize_both(merged_data, harmony_max_iter = 50)
saveRDS(harmony_obj, file.path(rdsdir, "harmonized.rds"))

