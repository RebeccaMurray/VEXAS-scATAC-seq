# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(version = "3.14")
# # library(remotes)
# remotes::install_github("satijalab/seurat", "feat/dictionary")
#remotes::install_github("mojaveazure/seurat-disk")
# BiocManager::install("biovizBase")
library(Seurat)
#library(SeuratDisk)
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
# plan("multicore", workers = 4)
# options(future.globals.maxSize = 80000 * 1024^2) # for 80 Gb RAM

data.path <- paste0("~/rscripts/VEXAS_GoTChA_signac/data/bridge_integration_metadata/")
dir.create(data.path)
# 
# ### PART 1: Load and set up the multiome data set, to be used as a bridge
# obj.multi <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/multiome_geo/multi_ref_2000.rds")
# obj.multi <- UpdateSeuratObject(obj.multi)
# 
# 
# ## PART 2: Load and set up the ATAC seq query
# obj.atac <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/all_samples_50_svd.RDS")
# 
# # Requantify query ATAC to have same features as multiome ATAC dataset
# requant_multiome_ATAC <- FeatureMatrix(
#   fragments = Fragments(obj.atac@assays$ATAC),
#   features = granges(obj.multi[['ATAC']]),
#   cells = Cells(obj.atac@assays$ATAC)
# )
# 
# # Create assay with requantified ATAC data
# ATAC_assay <- CreateChromatinAssay(
#   counts = requant_multiome_ATAC,
#   annotation = obj.atac@assays$ATAC@annotation
# )
# 
# # Create Seurat object
# obj.atac  <- CreateSeuratObject(counts = ATAC_assay, assay = 'ATAC')
# obj.atac[['peak.orig']] <- obj.atac@assays$ATAC
# 
# 
# ### PART 3: Load and set up the RNA reference
# 
# ## bone marrow reference
# obj.rna = readRDS("~/rscripts/VEXAS_GoTChA_signac/data/rna_reference/ref.Rds") # entire azimuth reference
# obj.rna <- UpdateSeuratObject(obj.rna)
# #obj.rna = readRDS("~/rscripts/VEXAS_GoTChA_signac/data/rna_reference/bmcite_demo.rds") # azimuth demo reference
# # obj.rna = SCTransform(object = obj.rna) %>% RunPCA() %>% RunUMAP(dims = 1:50, return.model = TRUE)
# 
# 
# 
# ### Preprocessing/normalization for all datasets
# 
# # normalize multiome RNA
# DefaultAssay(obj.multi) <- "RNA"
# obj.multi <- SCTransform(obj.multi, verbose = FALSE)
# 
# # normalize multiome ATAC
# DefaultAssay(obj.multi) <- "ATAC"
# obj.multi <- RunTFIDF(obj.multi)
# obj.multi <- FindTopFeatures(obj.multi, min.cutoff = "q0")
# 
# # normalize query
# obj.atac <- RunTFIDF(obj.atac)
# 
# 
# ## Save and load
# 
# saveRDS(obj.multi, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_multi.rds")
# saveRDS(obj.rna, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_rna.rds")
# saveRDS(obj.atac, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_atac.rds")

obj.multi <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_multi.rds")
obj.rna <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_rna.rds")
obj.atac <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_atac.rds")


## Perform the mapping!

# Drop first dimension for ATAC reduction
dims.atac <- 2:50
dims.rna <- 1:50

# print('Creating extended reference...')
# # Prepare "extended" reference
# DefaultAssay(obj.multi) <-  "RNA"
# obj.rna.ext <- PrepareBridgeReference(reference = obj.rna,
#                                       bridge = obj.multi,
#                                       reference.reduction = "refDR", # 'spca' for pbmc ref, 'pca' for small BM ref
#                                       reference.dims = dims.rna,
#                                       normalization.method = "SCT"
# )
# 
# ## Save and load
# saveRDS(obj.rna.ext, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_rna_ext.rds")
obj.rna.ext <- readRDS("~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/bridge_integration_07_25_23/obj_rna_ext.rds")

print('Transferring bridge anchors...')
DefaultAssay(obj.atac) <- "ATAC"
DefaultAssay(obj.rna.ext) <- "ATAC"
# Find transfer anchors between the extended reference and our query object
bridge.anchor <- FindBridgeTransferAnchors(extended.reference = obj.rna.ext,
                                           query = obj.atac,
                                           reduction = "lsiproject",
                                           dims = dims.atac
)


print('Mapping the query onto the reference...')
# Map the query dataset onto the reference
DefaultAssay(obj.atac) <- "ATAC"
obj.atac <- MapQuery(anchorset = bridge.anchor,
                     reference = obj.rna.ext,
                     query = obj.atac,
                     refdata = list(
                       l1 = "celltype.l1",
                       l2 = "celltype.l2"),
                     reduction.model = "refUMAP"
)

print('Mapping complete!')

md <- obj.atac@meta.data
md$bridge_integration_barcodes <- rownames(md)
write.csv(md, file = "data/bridge_integration_metadata/all_samples_through_BM18_07_24_23.csv")

print('Mapped query saved')
