library(Signac)
library(Seurat)
library(JASPAR2020)
library(TFBSTools)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(parallel)
set.seed(1234)

theme_set(theme_bw())

args = commandArgs(trailingOnly=TRUE)
current.cell.type <- do.call(paste, as.list(args))

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")

## Get a list of motif position frequency matrices from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)

DefaultAssay(atac.obj) <- "ATAC"

print(paste0("Running for cluster ", current.cell.type))

#################################### Calling peaks + running chromvar separately per cell type ######################

## run this using gotcha_signac_env_3_17_23
print("Running CallPeaks...")

  se <- subset(atac.obj, cluster_celltype == current.cell.type)
  DefaultAssay(se) <- "ATAC"
  print("Running for object:")
  print(se)
  rm(atac.obj)
  gc()
  temp.dir.location <- paste0("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/macs2_temp/", gsub(" ", "_", current.cell.type), "/")
  dir.create(temp.dir.location)
  macs2.peaks <- CallPeaks(se,
                           outdir = temp.dir.location, 
                           fragment.tempdir = temp.dir.location, 
                           macs2.path = "/nfs/sw/macs2/macs2-2.2.7.1/python/bin/macs2")
  macs2.counts <- FeatureMatrix(fragments = Fragments(se), features = macs2.peaks, cells = colnames(se))
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
  genome(annotations) <- "hg38"
  
  print("Peak calling complete")
  print(se)
  se[["macs2_peaks_celltype"]] <- CreateChromatinAssay(counts = macs2.counts, annotation = annotations)
  
  
  print("Adding motif information...")
  # add motif information
  DefaultAssay(se) <- "macs2_peaks_celltype"
  se <- AddMotifs(object = se, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "macs2_peaks_celltype", pfm = pfm)
  
  print("Running chromvar...")
  # Compute a per-cell motif activity score
  DefaultAssay(se) <- "macs2_peaks_celltype"
  se@assays$chromvar <- NULL
  se <- RunChromVAR(
    object = se,
    new.assay.name = "chromvar_celltype",
    assay = "macs2_peaks_celltype",
    genome = BSgenome.Hsapiens.UCSC.hg38
  )
  
  ## Name to ID matrix
  id.to.name <- data.frame()
  for (item.name in names(pfm)) {
    print(item.name)
    id.to.name <- rbind(id.to.name, data.frame(ID = pfm[[item.name]]@ID, name = pfm[[item.name]]@name))
  }
  
  print("Updating ID names...")
  
  ## Replace chromvar assay ID rownames with gene names
  rownames(se@assays$chromvar_celltype@data) <- plyr::mapvalues(rownames(se@assays$chromvar_celltype@data), from = id.to.name$ID, to = id.to.name$name)

  print("Done!")
  print(se)
    
  ## Return the chromvar matrix

## Save the chromvar data
print("Saving chromvar matrix...")

saveRDS(se@assays$chromvar_celltype@data, paste0("data/seurat_objects/", current.cell.type, "_chromvar.rds"))
