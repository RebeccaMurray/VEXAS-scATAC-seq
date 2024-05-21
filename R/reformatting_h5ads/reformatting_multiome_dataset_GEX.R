# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# library(remotes)
# remotes::install_github("satijalab/seurat", "feat/dictionary")
# remotes::install_github("mojaveazure/seurat-disk")
# BiocManager::install("biovizBase")
library(Seurat)
library(SeuratDisk)
library(Signac)
# library(EnsDb.Hsapiens.v86)
# library(dplyr)
# library(ggplot2)
# library(GenomeInfoDb)
# library(GenomicRanges)
# library(BSgenome.Hsapiens.UCSC.hg38)
library(future)

## Memory/tech specs
set.seed(1234)
plan("multicore", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2) # for 100 Gb RAM


## Convert from h5ad format to h5seurat format
multiome.data.file <- "/gpfs/commons/home/rmurray/rscripts/VEXAS_GoTChA_signac/data/multiome_geo/GEX_only.h5ad"

print(paste0("File name: ", multiome.data.file))
print("Beginning conversion...")

Convert(multiome.data.file, dest = "h5seurat", assay = "RNA", overwrite = FALSE, verbose = TRUE)

print("Complete!")
