#library(SeuratObject, lib.loc = "/projects/Seurat4")
#library(Seurat, lib.loc = "/projects/Seurat4")
#library(Signac, lib.loc = "/projects/Seurat4")
#library(ggplot2, lib.loc = "/projects/Seurat4")
library(SeuratObject)
library(Seurat)
library(Signac)
library(ggplot2)
library(EnsDb.Hsapiens.v86)
library(AnnotationFilter)
library(SeuratDisk)
library(dplyr)
library(qlcMatrix)
library(BSgenome.Hsapiens.UCSC.hg38)
#library(ggplot2)
library(data.table)
library(harmony, lib.loc = "/projects/Seurat4")
library(future)
library(MASS)
source("/projects/opioid/Rmultiome/functions/premerge_processing_functions.R")
source("/projects/opioid/Rmultiome/functions/postmerge_processing_functions.R")
source("/projects/opioid/Rmultiome/functions/postmerge_parameter_functions.R")
source("/projects/opioid/Rmultiome/functions/premerge_parameter_functions.R")
source("/projects/opioid/Rmultiome/functions/helper_functions.R")
source("/projects/opioid/Rmultiome/functions/DE_functions.R")

#TODO:  Seurat5 was painful because all documentation available for 10X multiome
#   uses Seurat4 and RunHarmony, and we have multiple samples.  Seurat5 doesn't
#   "merge" in the same way, as it keeps each sample separate (ie...not merged)
#   Sans-vignette, https://satijalab.org/seurat/reference/harmonyintegration is
#   part of the new process for this, and is also a great example of poor docs.
#   Intention, however, is to move this tool to using Seurat5, to future-proof.

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
