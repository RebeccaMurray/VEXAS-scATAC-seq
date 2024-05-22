library(ggh4x)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
library(scales)
theme_set(theme_classic())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "CD14 Mono", "B", "Plasma", "pre-mDC", "pDC", "MAIT", "CD4 T", "CD8 T", "NK") ## Ignore the Stromal cluster
atac.obj$CellType <- factor(atac.obj$cluster_celltype, levels = celltype.order)


###################################### Set up the ADTs - pull average expression per cluster ################################################

current.assay <- "ADT_DSB"

## Create an seurat object with only cells with non-zero ADT expression
adt.counts.per.bc <- colSums(atac.obj@assays[[current.assay]]@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
atac.obj.adt <- subset(atac.obj, cells = bcs.keep)

## Scale the DSB assay
atac.obj.adt <- ScaleData(atac.obj.adt, assay = current.assay)
atac.obj.adt@assays$ADT_DSB@data[1:5, 1:5]

## Compile a data frame with the average expression per cluster
avg.exp <- lapply(levels(atac.obj.adt$CellType), function(x) {
  print(x)
  bc.list <- atac.obj.adt@meta.data %>% filter(CellType == x) %>% rownames()
  if (length(bc.list) < 50) {
    return(NA_real_)
  }
  df <- data.frame(expression_mean = rowMeans(atac.obj.adt[[current.assay]]@scale.data[, bc.list], na.rm = T))
  colnames(df) <- x
  return(df)
})
avg.exp <- avg.exp[!is.na(avg.exp)]
avg.exp.matrix <- do.call(cbind, avg.exp)
avg.exp.df <- avg.exp.matrix %>% rownames_to_column("ADT") %>% pivot_longer(cols = colnames(avg.exp.matrix), names_to = "cluster", values_to = "average_expression")


###################################### ADT annotation plot - manually selected tags ################################################

## Manually selecting tags)
tags.keep <- c("Hu.CD90", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD36", "Hu.CD41", "Hu.CD123", "Hu.CD11c", "Hu.CD45RA",  "Hu.CD14-M5E2",  "Hu.CD20-2H7", "Hu.CD19", "Hu.CD161", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD3-UCHT1", "Hu.CD56", "Hu.CD16", "Hu.CD7")
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "MkP", "BaEoMa", "CLP", "GMP", "CD14 Mono", "B", "Plasma", "pre-mDC", "pDC", "CD4 T", "CD8 T", "NK") ## Ignore the Stromal cluster
dev.off()
pdf("figures/current_figure_drafts/ATAC_adt_heatmap.pdf", width = 4, height = 4)
pheatmap::pheatmap(t(avg.exp.matrix[tags.keep, celltype.order]), 
                   # scale = "column",
                   main = "ADT expression",
                   cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-2, 2, length.out = 50))
dev.off()
