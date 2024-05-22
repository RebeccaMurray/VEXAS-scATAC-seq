library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
library(tidyverse)
library(Matrix)

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")

## Load gene scores that were computed in ArchR
archr.gene.scores <- readRDS(file = "/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_ArchR/data/gene_score_matrices/minTSS_2_minFrags_500_arrow_files.rds")
colnames(archr.gene.scores) <- gsub("#", "_", colnames(archr.gene.scores))
archr.gene.scores[1:5, 1:5]
dim(archr.gene.scores)

## Only keep barcodes that are in ATAC object
atac.bcs <- atac.obj@meta.data %>% rownames()
archr.bcs <- colnames(archr.gene.scores)
bcs.keep <- intersect(atac.bcs, archr.bcs)
archr.gene.scores <- archr.gene.scores[, bcs.keep] ## Any missing barcodes will be added as NA rows

## Add as a new assay
archr.gene.score.assay <- CreateAssayObject(data = archr.gene.scores)
atac.obj[["RNA_ArchR"]] <- archr.gene.score.assay

# ## Plot?
# DefaultAssay(atac.obj) <- "RNA_ArchR"
# FeaturePlot(atac.obj, features = "AVP", min.cutoff = "q10", max.cutoff = "q90")
# FeaturePlot(atac.obj, features = "MPO", min.cutoff = "q10", max.cutoff = "q90")
# FeaturePlot(atac.obj, features = "FLT3", min.cutoff = "q10", max.cutoff = "q90")

saveRDS(atac.obj, "~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")


############ Extra #############

# ## Option 1 for missing barcodes: just remove them from ATAC object
# atac.obj <- subset(atac.obj, cells = bcs.keep) ## Removes 1 cell

# ## Option 2 for missing barcodes: add NA rows to ArchR gene score matrix for them
# bcs.missing <- setdiff(atac.bcs, archr.bcs)
# missing_matrix <- Matrix(0, nrow=nrow(archr.gene.scores), ncol=length(bcs.missing), dimnames=list(rownames(archr.gene.scores), bcs.missing))
# archr.gene.scores <- cbind(archr.gene.scores, missing_matrix)