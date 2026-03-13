# Quick test script for POC
source("/projects1/opioid/Rmultiome/system_settings.R")
source("/projects/opioid/Rmultiome/Rmultiome-main.R")

#adding as a checkpoint during code development
# chosen_obj <- readRDS(file.path(rdsdir, "chosen_clustered_obj.rds"))

# Load your all_markers
all_markers <- read.csv(file.path(tmpfiledir, "all_cluster_markers.csv"))

# Test with just oligodendrocytes
oligo_result <- identify_celltype(
  all_markers,
  marker_genes = c("MBP", "MOBP", "PLP1"),
  celltype_name = "Oligodendrocytes",
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# Check result
print(oligo_result$assignment)
