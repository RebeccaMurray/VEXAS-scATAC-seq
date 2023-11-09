library(Signac)
library(Seurat)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
library(ggrepel)
set.seed(1234)

theme_set(theme_bw())

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_20231103.RDS")

## Get a list of motif position frequency matrices from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)


################# Adding motif information (only need to do once) ##################

# add motif information
DefaultAssay(atac.obj) <- "ATAC"
atac.obj <- AddMotifs(object = atac.obj, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "ATAC", pfm = pfm)

# Compute a per-cell motif activity score
DefaultAssay(atac.obj) <- "ATAC"
atac.obj <- RunChromVAR(
  object = atac.obj,
  assay = "ATAC",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

################## Create an ID to name matrix for motifs #############

## Name to ID matrix
id.to.name <- data.frame()
for (item.name in names(pfm)) {
  print(item.name)
  id.to.name <- rbind(id.to.name, data.frame(ID = pfm[[item.name]]@ID, name = pfm[[item.name]]@name))
}

## Replace chromvar assay ID rownames with gene names
rownames(atac.obj@assays$chromvar@data) <- plyr::mapvalues(rownames(atac.obj@assays$chromvar@data), from = id.to.name$ID, to = id.to.name$name)

#################### Save updated ATAC obj ####################

saveRDS(atac.obj, "~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_20231103.RDS")
