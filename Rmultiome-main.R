library(SeuratObject)
library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
#library(SeuratDisk)
library(dplyr)
library(qlcMatrix)
library(BSgenome.Hsapiens.UCSC.hg38)
library(data.table)
library(future)
library(MASS)
source(file.path(Rmultiome_path, "functions", "premerge_processing_functions.R"))
source(file.path(Rmultiome_path, "functions", "postmerge_processing_functions.R"))
source(file.path(Rmultiome_path, "functions", "postmerge_parameter_functions.R"))
source(file.path(Rmultiome_path, "functions", "premerge_parameter_functions.R"))
source(file.path(Rmultiome_path, "functions", "helper_functions.R"))
source(file.path(Rmultiome_path, "functions", "DE_functions.R"))
source(file.path(Rmultiome_path, "functions", "qc_functions.R"))

if (FALSE) {
  #dependencies that should be outside the main libdir, but I don't use venv.
  dir.create("/projects/Seurat4", showWarnings = FALSE)
  install.packages("remotes") # if not already installed
  remotes::install_github("satijalab/seurat@v4.3.0",
                          lib = "/projects/Seurat4")
  remotes::install_github("satijalab/seurat-object@v4.1.3",
                          lib = "/projects/Seurat4")
  remotes::install_github("mojaveazure/seurat-disk",
                          lib = "/projects/Seurat4")
  install.packages("/projects/scratch/signac-1.11.0.tar.gz",
                   repos = NULL, type = "source",
                   lib = "/projects/Seurat4")
  install.packages("/projects/scratch/harmony_1.1.0.tar.gz",
                   repos = NULL, type = "source",
                   lib = "/projects/Seurat4")
  BiocManager::install("HGNChelper")
}
