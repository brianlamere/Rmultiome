source("/projects1/opioid/Rmultiome/system_settings.R")
source("/projects/opioid/Rmultiome/Rmultiome-main.R")


cocaine_markers_list <- split(cocaine_paper_markers$gene,
                             cocaine_paper_markers$celltype)

# Load your data
all_markers <- read.csv(file.path(tmpfiledir, "all_cluster_markers.csv"))

results <- identify_all_celltypes(
  all_markers,
  cocaine_markers_list,
  min_markers = 2,
  min_score = 5,
  verbose = TRUE
)

# Save results
write.csv(results$assignments,
         file.path(tmpfiledir, "celltype_assignments.csv"),
         row.names = FALSE)

# Show summary
cat("\n=== Assigned clusters ===\n")
print(results$assignments %>%
      filter(!is.na(cluster)) %>%
      select(cluster, celltype, confidence, n_markers, score))

# Check for unassigned
all_clusters <- setdiff(levels(Idents(chosen_obj)), "singleton")
assigned <- results$assignments$cluster[!is.na(results$assignments$cluster)]
unassigned <- setdiff(all_clusters, assigned)

if (length(unassigned) > 0) {
  cat(sprintf("\n%d unassigned clusters: %s\n",
             length(unassigned),
             paste(unassigned, collapse = ", ")))
}
