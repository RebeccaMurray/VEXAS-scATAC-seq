library(Seurat)
library(Signac)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(ggsci)
library(viridis)
library(msigdbr)
theme_set(theme_bw())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits.RDS")
md <- atac.obj@meta.data 

########################################### UMAPs with cell type labels ########################################
# 
# ## Plot Azimuth labels
# p1 <- DimPlot(atac.obj, group.by = "seurat_clusters", label = T, label.box = T, repel = T) + NoLegend() + coord_fixed()
# p2 <- DimPlot(atac.obj, group.by = "predicted.l2", label = T, label.box = T, repel = T) + NoLegend() + coord_fixed()
# p3 <- DimPlot(atac.obj, group.by = "orig.ident", label = T, label.box = T, repel = T, shuffle = T) + NoLegend() + coord_fixed()
# p1 | p2
# p1 | p2 | p3

########################################### Adjust clustering if necessary ##################################

DefaultAssay(atac.obj) <- "ATAC"
atac.obj <- FindClusters(object = atac.obj, algorithm = 3, resolution = 2, random.seed = 1)
md <- atac.obj@meta.data 

################################## Automatic assignment ######################################################

cell.type.mapping <- atac.obj@meta.data %>% 
  group_by(seurat_clusters) %>% 
  # filter(predicted.l2.score > 0.5) %>%
  dplyr::count(predicted.l2) %>% 
  slice_max(n) %>% 
  select(seurat_clusters, predicted.l2) %>% 
  mutate(predicted.l2 = gsub(".*CD8.*", "CD8 T", predicted.l2)) %>% 
  mutate(predicted.l2 = gsub(".*CD4.*", "CD4 T", predicted.l2)) %>% 
  mutate(predicted.l2 = gsub("Memory ", "", predicted.l2)) %>% 
  mutate(predicted.l2 = gsub("Prog Mk", "MkP", predicted.l2)) %>% 
  mutate(predicted.l2 = gsub("cDC.*", "cDC", predicted.l2))
 

# ## Plot some QC metrics, to see if there are any clusters we need to remove
# atac.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize(mean_ncount = mean(nCount_ATAC), mean_pct_reads_peaks = mean(pct_reads_in_peaks), mean_tss = mean(TSS.enrichment)) %>% 
#   ggplot(aes(x = mean_tss, y = mean_pct_reads_peaks, label = seurat_clusters)) +
#   geom_point() +
#   geom_label_repel()
# DimPlot(atac.obj, split.by = "seurat_clusters", group.by = "seurat_clusters", ncol = 6) + NoLegend()
# DimPlot(atac.obj, cells.highlight = atac.obj@meta.data %>% filter(seurat_clusters == 30) %>% rownames())
# atac.obj$mito_fragment_percent = atac.obj$mitochondrial / atac.obj$total
# VlnPlot(atac.obj, features = "nCount_ATAC", pt.size = 0) + NoLegend()
# VlnPlot(atac.obj, features = "nCount_ATAC", pt.size = 0) + NoLegend() + coord_cartesian(ylim=c(500, 5000))
# atac.obj@meta.data %>% 
#   group_by(seurat_clusters) %>% 
#   summarize(mean_mito_frag = mean(mito_fragment_percent)) %>% 
#   ggplot(aes(x = mean_mito_frag)) +
#   geom_histogram()
# atac.obj@meta.data %>% 
#     group_by(seurat_clusters) %>% 
#     summarize(mean_mito_frag = mean(mito_fragment_percent), mean_ncount = mean(nCount_ATAC), mean_pct_reads_peaks = mean(pct_reads_in_peaks), mean_tss = mean(TSS.enrichment)) %>% 
#     ggplot(aes(x = mean_mito_frag, y = mean_ncount, label = seurat_clusters)) +
#     geom_point() +
#     geom_label_repel()
# atac.obj@meta.data %>% filter(predicted.l2 == "Prog Mk") %>% dplyr::count(seurat_clusters) %>% arrange(-n) %>% head

## Adjust clustering slightly
cell.type.mapping$predicted.l2[cell.type.mapping$seurat_clusters == 17] <- "EMP"

## Add the cluster annotations
atac.obj$cluster_celltype <- plyr::mapvalues(atac.obj$seurat_clusters, from = cell.type.mapping$seurat_clusters, to = cell.type.mapping$predicted.l2)
atac.obj$CellType <- atac.obj$cluster_celltype

# p2 <- DimPlot(atac.obj, group.by = "predicted.l2", label = T, label.box = T, repel = T) + NoLegend()
# p3 <- DimPlot(atac.obj, group.by = "cluster_celltype", label = T, label.box = T, repel = T, cols = cell.type.palette) + NoLegend()
# 
# p2 | p3

saveRDS(atac.obj, "~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")



########################## Extra #########################
# 
# adt.counts.per.bc <- colSums(atac.obj@assays$ADT_DSB@data)
# bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
# atac.obj.adt <- subset(atac.obj, cells = bcs.keep)
# hspc.clusters <- c("6", "4", "11", "1", "30", "9", "10", "23", "25", "27", "12", "20")
# p3 <- VlnPlot(atac.obj.adt, features = "Hu.CD90", group.by = "seurat_clusters", pt.size = 0, idents = hspc.clusters) + NoLegend() + coord_cartesian(ylim=c(0, 25))
# p4 <- VlnPlot(atac.obj.adt, features = "Hu.CD71", group.by = "seurat_clusters", pt.size = 0, idents = hspc.clusters) + NoLegend() + coord_cartesian(ylim=c(0, 25))
# p3 / p4
# 
# current.assay <- "ADT_DSB"
# 
# ## Scale the DSB assay
# atac.obj.adt <- ScaleData(atac.obj.adt, assay = current.assay)
# atac.obj.adt@assays$ADT_DSB@data[1:5, 1:5]
# atac.obj.adt@assays$ADT_DSB@scale.data[1:5, 1:5]
# 
# ## Compile a data frame with the average expression per cluster
# avg.exp <- lapply(levels(atac.obj.adt$seurat_clusters), function(x) {
#   print(x)
#   bc.list <- atac.obj.adt@meta.data %>% filter(seurat_clusters == x) %>% rownames()
#   df <- data.frame(expression_mean = rowMeans(atac.obj.adt[[current.assay]]@scale.data[, bc.list], na.rm = T))
#   colnames(df) <- x
#   return(df)
# })
# avg.exp.matrix <- do.call(cbind, avg.exp)
# avg.exp.df <- avg.exp.matrix %>% rownames_to_column("ADT") %>% pivot_longer(cols = colnames(avg.exp.matrix), names_to = "cluster", values_to = "average_expression")
# 
# 
# ## Manually selecting tags)
# tags.keep = c("Hu.CD90", "Hu.CD38-HIT2", "Hu.CD71", "Hu.CD36", "Hu.CD41", "Hu.CD45RA", "Hu.CD7", "Hu.CD123")
# pheatmap::pheatmap(avg.exp.matrix[tags.keep, hspc.clusters], 
#                    # scale = "row",
#                    cluster_rows = F, cluster_cols = F, color=colorRampPalette(c("navy", "white", "red"))(50), breaks = seq(-1, 1, length.out = 50))
# 
# 
# ################################## Extra analysis ######################################################
# 
# DimPlot(atac.obj, group.by = "seurat_clusters", split.by = "seurat_clusters", ncol = 6)
# DimPlot(atac.obj, group.by = "seurat_clusters", split.by = "orig.ident", ncol = 3)
# DimPlot(atac.obj, cells.highlight = atac.obj@meta.data %>% filter(seurat_clusters == 23) %>% rownames())
# atac.obj@meta.data %>% 
#   filter(seurat_clusters == 23) %>% 
#   count(orig.ident) %>% arrange(-n) %>% head()
# 
# ## Plot some QC metrics
# FeaturePlot(atac.obj, features = c("nCount_ATAC")) + NoLegend()
# VlnPlot(atac.obj, features = c("nCount_ATAC"), group.by = "seurat_clusters") + NoLegend()
# VlnPlot(atac.obj, features = c("nFeature_ATAC"), group.by = "seurat_clusters") + NoLegend()
# VlnPlot(atac.obj, features = c("pct_reads_in_peaks"), group.by = "seurat_clusters") + NoLegend()
# 
# ## Make a bar plot of the number of cell types per cluster
# md %>%
#   filter(predicted.l2.score > 0.5) %>%
#   group_by(predicted.l2, seurat_clusters) %>%
#   count() %>% 
#   group_by(seurat_clusters) %>% 
#   # top_n(n, 5) %>%
#   ggplot(aes(x = reorder(predicted.l2, -n), y = n, fill = predicted.l2)) +
#   geom_col() +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
#   facet_wrap(~ seurat_clusters, scales = "free") +
#   scale_fill_manual(values = cell.type.palette) +
#   theme(legend.position = "none") +
#   labs(x = "Azimuth prediction", y = "Cell count") +
#   ggtitle("Prediction score > 0.5")
# 
# ## Make a histogram of the celltype prediction scores
# md %>% 
#   ggplot(aes(x = predicted.l2.score, y = ..count..)) +
#   geom_histogram() +
#   ggtitle("Predicted celltype score, all cells")+ 
#   geom_vline(xintercept = 0.6)
# 
# md %>% 
#   filter(predicted.l2 %in% c("NK", "CD8 Effector_1", "CD8 Effector_2", "HSC")) %>% 
#   ggplot(aes(x = predicted.l2.score, y = ..count.., fill = predicted.l2)) +
#   geom_histogram() +
#   facet_grid(rows = vars(predicted.l2), scales = "free") +
#   scale_fill_manual(values = cell.type.palette) +
#   ggtitle("Predicted celltype score") + 
#   geom_vline(xintercept = 0.6)
# 
# ## List cell types in each cluster
# seurat_cluster = 30
# md %>% filter(seurat_clusters == seurat_cluster & predicted.l2.score > 0.1) %>% count(predicted.l2) %>% arrange(-n) %>% head(n=10)
# md %>% filter(seurat_clusters == seurat_cluster) %>% count(predicted.l2) %>% arrange(-n) %>% head(n=10)
# md %>% filter(predicted.l2 == "LMPP") %>% count(seurat_clusters) %>% arrange(-n) %>% head
# DimPlot(atac.obj, cells.highlight = md %>% filter(seurat_clusters == seurat_cluster) %>% rownames())
# 
# ## Plot cell type predictions
# predicted_celltype = 'HSC'
# DimPlot(atac.obj, cells.highlight = md %>% filter(predicted.l2 == predicted_celltype & predicted.l2.score > 0.6) %>% rownames())
# 
# 
# ## Resolution 1
# cluster.to.celltype.mapping <- list(
#   `0` = "CD8 T",
#   `1` = "LMPP",
#   `2` = "CD14 Mono",
#   `3` = "CD14 Mono",
#   `4` = "CD8 T",
#   `5` = "HSC",
#   `6` = "HSC",
#   `7` = "Early Eryth",
#   `8` = "CD14 Mono",
#   `9` = "GMP",
#   `10` = "EMP",
#   `11` = "MkP",
#   `12` = "pDC",
#   `13` = "HSC",
#   `14` = "CD4 T",
#   `15` = "pre-mDC",
#   `16` = "Early Eryth",
#   `17` = "Late Eryth",
#   `18` = "GMP",
#   `19` = "GMP",
#   `20` = "B",
#   `21` = "NK",
#   `22` = "CLP",
#   `23` = "BaEoMa",
#   `24` = "Plasma",
#   `25` = "CD14 Mono",
#   `26` = "Early Eryth",
#   `27` = "CD8 T",
#   `28` = "LMPP"
# )
# 
# 
# atac.obj$cluster_celltype <- plyr::mapvalues(atac.obj$seurat_clusters, from = names(cluster.to.celltype.mapping), to = unlist(cluster.to.celltype.mapping))
# atac.obj$CellType <- atac.obj$cluster_celltype
# 
# p1 <- DimPlot(atac.obj, group.by = "cluster_celltype", label = T, label.box = T, repel = T, cols = cell.type.palette) + NoLegend() + theme(aspect.ratio = 1)
# p2 <- DimPlot(atac.obj, group.by = "predicted.l2", label = T, label.box = T, repel = T) + NoLegend() + theme(aspect.ratio = 1)
# p1 | p2
# 
# saveRDS(atac.obj, "~/rscripts/VEXAS_GoTChA_signac/data/seurat_objects/atac_obj_updated.rds")
# 
# 
# 
# 
# 
# 
# 
# ######################################### Extra ##########################################
# 
# ## Identifying neutrophil progenitors (not included in azimuth, but in GoT + CH papers)
# ## CH paper
# ch.paper.modules <- readxl::read_excel("data/CH_paper_references/41588_2022_1179_MOESM3_ESM.xlsx", skip = 1, sheet = 2) %>% 
#   filter(!grepl("et al", `Cell Cycle module`))
# 
# DefaultAssay(atac.obj) <- "RNA"
# atac.obj <- NormalizeData(atac.obj)
# 
# ## Adding CH paper modules
# for (module.name in colnames(ch.paper.modules)) {
#   module <- ch.paper.modules[[module.name]]
#   module <- module[!is.na(module)]
#   atac.obj <- AddModuleScore(atac.obj, features = list(module), name = gsub(" ", "_", module.name), assay = "RNA", search = TRUE)
# }
# 
# module.names <- colnames(atac.obj@meta.data)[44:51]
# module.names <- colnames(atac.obj@meta.data)[54:69]
# 
# FeaturePlot(atac.obj, features = module.names, cols = c("blue", "yellow"), ncol = 4)
# 
# 
# ######################## Using ADT for help ###############################
# 
# ## Subset for samples we've run ADT for
# se.subset <- subset(atac.obj, orig.ident %in% c("BM1_CD34_plus", "BM2_CD34_plus", "VX_BM12_A", "VX_BM13_B", "VX_BM14_Fpos", "VX_BM15a", "VX_BM15b"))
# DefaultAssay(se.subset) <- "CLR_CITE"
# tags <- c("Hu.CD14-M5E2", "Hu.CD19", "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8")
# FeaturePlot(se.subset, features = tags, cols = c("blue", "red"), ncol = 3, min.cutoff = "q5", max.cutoff = "q95")
# 
# DefaultAssay(se.subset) <- "CLR_CITE"
# avg.exp <- list()
# for (celltype in unique(se.subset$cluster_celltype)) {
#   cluster.barcodes <- se.subset@meta.data %>% filter(cluster_celltype == celltype) %>% rownames()
#   avg.exp.cluster <- se.subset@assays$CLR_CITE@data[, cluster.barcodes] %>% rowMeans(na.rm = T)
#   avg.exp[[celltype]] <- avg.exp.cluster
# }
# avg.exp <- as.data.frame(avg.exp)
# colnames(avg.exp) <- gsub("X", "", colnames(avg.exp))
# 
# ## Markers to plot
# marker.list <- c(
#   "Hu.CD14-M5E2", 
#   "Hu.CD19", 
#   "Hu.CD3-UCHT1", 
#   "Hu.CD4-RPA.T4", 
#   "Hu.CD8",
#   "Hu.CD38-HIT2",
#   "Hu.CD71",
#   "Hu.CD41",
#   "Hu.CD123",
#   "Hu.CD11c",
#   "Hu.CD45RA"
# )
# 
# pheatmap::pheatmap(t(avg.exp[marker.list, ]),
#                    breaks = seq(0, 30, length.out = 100))
# 
