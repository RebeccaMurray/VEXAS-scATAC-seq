if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
library(remotes)
remotes::install_github("satijalab/seurat", "feat/dictionary")
remotes::install_github("mojaveazure/seurat-disk")
BiocManager::install("biovizBase")
library(Seurat)
library(SeuratDisk)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(GenomeInfoDb)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(future)

## Memory/tech specs
set.seed(1234)
plan("multicore", workers = 4)
options(future.globals.maxSize = 300000 * 1024^2) # for 300 Gb RAM


#### OPTION 1: Attempting to convert entire file at once

## Convert from h5ad format to h5seurat format
multiome.data.file <- "/gpfs/commons/home/rmurray/rscripts/VEXAS_GoTChA_signac/data/multiome_geo/GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad"

print("Beginning conversion...")

Convert(multiome.data.file, dest = "h5seurat", overwrite = FALSE, verbose = TRUE) #  This one succeeded without error for entire file

print("Complete!")
