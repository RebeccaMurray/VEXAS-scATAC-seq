library(Seurat)
library(Signac)
library(tidyverse)


## Load the bridge integration annotations, the atac object
df <- read_csv("data/bridge_integration_metadata/all_samples_through_BM18_07_24_23.csv") %>% 
  select(predicted.l1, predicted.l1.score, predicted.l2, predicted.l2.score, bridge_integration_barcodes) %>% 
  column_to_rownames("bridge_integration_barcodes")

atac.obj <-  readRDS("data/seurat_objects/all_samples_50_svd.RDS")

## Add the bridge annotations
atac.obj <- AddMetaData(atac.obj, df)

## Make some plots
DimPlot(atac.obj, group.by = "predicted.l2", label = T, label.box = T, repel = T) + NoLegend()
atac.obj@meta.data %>% filter(predicted.l2 == "HSC"& predicted.l2.score > 0.5) %>% count(orig.ident) 
atac.obj@meta.data %>% filter(orig.ident == "VEX_BM11_LIN_NEG_ATAC") %>% count(predicted.l2)
