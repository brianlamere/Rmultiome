# Let DotPlot do the work
dp <- DotPlot(chosen_obj, features = c("CUX2","SATB2","LAMP5"))

# Extract the data
cluster_table <- dp$data %>%
  dplyr::select(
    cluster = id,
    gene = features.plot,
    pct_expressed = pct.exp,
    avg_expression = avg.exp.scaled
  )

# View
print(cluster_table)

# Show only clusters with >1% expression
print(cluster_table %>% filter(pct_expressed > 1))

# Or pivot to wide format
cluster_wide <- cluster_table %>%
  tidyr::pivot_wider(
    names_from = gene,
    values_from = c(pct_expressed, avg_expression)
  )

print(cluster_wide)
