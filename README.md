# Rmultiome
A small tool set that is being designed for 10X multiome data from brain nuclei

Rebased after learning many of the vignettes that seem to apply, don't really
for my particular data set (very high depth, brain nuclei, thus very strong 
batch effects but no need to do things like redoing peaks, etc).

## Configuration

Project settings are stored in `settings.R`.  
Edit this file to match your local paths and data files before running any
scripts.

**Example:**
```r
sourcedir <- "/path/to/data"
h5filename <- "filtered_feature_bc_matrix.h5"
atacfilename <- "atac_fragments.tsv.gz"
referencedir <- "/path/to/data"
referencefile <- "referencefile.rds"
scratchdir1 <- "/path/to/scratch/dir"
```

### RDS Milestone Files

Milestone files are saved as `rds/{sample}_{milestone}.rds`.
This allows you to revert or restart analysis at any stage.

Example:
- `sampleA_qc.rds`
- `sampleA_filtered.rds`
- `sampleA_norm.rds`

### External Dependencies
Usually I like to be able to just have lines that say:
library("Seurat") and then move on, without there being version requirements, and assume people 
use a somewhat recent version.  However, Seurat5 added layers to the SeuratObject,
and many of the core Seurat tools I need (RunPCA, RunHarmony, etc) don't work
with those layers.  I honestly don't understand why Seurat5 was even released,
given it won't let you collapse the layers for ATAC (JoinLayers works for RNA,
but the tool doesn't work for ChromatinAssay objects), nor will it let you 
prevent the layers during merge(collapse=TRUE).  You're simply shut out of using 
numerous core tools of Seurat with 10X multiome data that includes ATAC.  With
that in mind, here is the environment I ended up using:

> packageVersion("Seurat")
[1] ‘4.3.0’

> packageVersion("SeuratObject")
[1] ‘4.1.3’

> packageVersion("Signac")
[1] ‘1.11.0’

Newer versions of Signac required Seurat5.
