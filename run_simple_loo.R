loo_jobs <- data.frame(
  comparison   = c("Low_vs_No_HIV", "Low_vs_No_HIV", "Low_vs_No_HIV",
                   "Acute_vs_Low",  "Acute_vs_Low",  "Acute_vs_Low",
                   "Acute_vs_Low",  "Acute_vs_Low",  "Acute_vs_Low",
                   "Chronic_vs_Low","Chronic_vs_Low","Chronic_vs_Low",
                   "Chronic_vs_Low","Chronic_vs_Low","Chronic_vs_Low"),
  ident.1      = c("Low","Low","Low",
                   "Acute","Acute","Acute",
                   "Acute","Acute","Acute",
                   "Chronic","Chronic","Chronic",
                   "Chronic","Chronic","Chronic"),
  ident.2      = c("No_HIV","No_HIV","No_HIV",
                   "Low","Low","Low",
                   "Low","Low","Low",
                   "Low","Low","Low",
                   "Low","Low","Low"),
  excluded     = c("LG22","LG25","LG38",
                   "LG05","LG26","LG31",
                   "LG22","LG25","LG38",
                   "LG08","LG23","LG33",
                   "LG22","LG25","LG38"),
  excluded_group = c("Low","Low","Low",
                     "Acute","Acute","Acute",
                     "Low","Low","Low",
                     "Chronic","Chronic","Chronic",
                     "Low","Low","Low"),
  stringsAsFactors = FALSE
)

celltypes_list <- c("Oligodendrocytes", "Microglia_Macrophages", "Astrocytes")

loo_jobs$filename_de <- sprintf("LOO_DE_%s_excluding_%s_%s.csv",
                                 loo_jobs$comparison,
                                 loo_jobs$excluded_group,
                                 loo_jobs$excluded)
loo_jobs$filename_da <- sprintf("LOO_DA_%s_excluding_%s_%s.csv",
                                 loo_jobs$comparison,
                                 loo_jobs$excluded_group,
                                 loo_jobs$excluded)
loo_jobs_full <- merge(loo_jobs,
                        data.frame(celltype = celltypes_list,
                                   stringsAsFactors = FALSE))

loo_jobs_full$filename_de <- sprintf("LOO_DE_%s_%s_excluding_%s_%s.csv",
                                      loo_jobs_full$celltype,
                                      loo_jobs_full$comparison,
                                      loo_jobs_full$excluded_group,
                                      loo_jobs_full$excluded)

loo_jobs_full$filename_da <- sprintf("LOO_DA_%s_%s_excluding_%s_%s.csv",
                                      loo_jobs_full$celltype,
                                      loo_jobs_full$comparison,
                                      loo_jobs_full$excluded_group,
                                      loo_jobs_full$excluded)

full_obj <- readRDS(file.path(rdsdir, "annotated_obj_final.rds"))
full_obj <- JoinLayers(full_obj)

for (i in 1:nrow(loo_jobs_full)) {
  job <- loo_jobs_full[i, ]
  cat(sprintf("[%d/%d] %s | %s | excluding %s (%s)\n",
              i, nrow(loo_jobs_full),
              job$celltype, job$comparison,
              job$excluded, job$excluded_group))

  loo_obj <- subset(full_obj, subset = orig.ident != job$excluded)

  de <- run_pseudobulk_DE(
    seurat_obj = loo_obj,
    celltype_col = "celltype",
    celltype_value = job$celltype,
    sample_col = "orig.ident",
    group_col = "group",
    ident.1 = job$ident.1,
    ident.2 = job$ident.2,
    min_cells_per_sample = 50
  )
  if (nrow(de) > 0)
    write.csv(de, file.path(tmpfiledir, job$filename_de), row.names = FALSE)

  da <- run_pseudobulk_DA(
    seurat_obj = loo_obj,
    celltype_col = "celltype",
    celltype_value = job$celltype,
    sample_col = "orig.ident",
    group_col = "group",
    ident.1 = job$ident.1,
    ident.2 = job$ident.2,
    min_cells_per_sample = 50
  )
  if (nrow(da) > 0)
    write.csv(da, file.path(tmpfiledir, job$filename_da), row.names = FALSE)
}
