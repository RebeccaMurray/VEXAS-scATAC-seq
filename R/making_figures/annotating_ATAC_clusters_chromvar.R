library(ggh4x)
library(Signac)
library(Seurat)
library(tidyverse)
library(gridExtra)
library(ggrepel)
library(ggsci)
library(ggpubr)
library(viridis)
library(pheatmap)
library(patchwork)
library(scales)
theme_set(theme_classic())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "CD14 Mono", "B", "Plasma", "pDC", "pre-mDC", "MAIT", "CD4 T", "CD8 T", "NK") ## Ignore the Stromal cluster
atac.obj$CellType <- factor(atac.obj$cluster_celltype, levels = celltype.order)
table(atac.obj$CellType, useNA = "always")

####################################### Load per-cluster chromvar results #################################

chromvar.dfs <- lapply(list.files("~/rscripts/VEXAS_ATAC_signac/data/cluster_chromvar_results", full.names = T), readRDS)

chromvar.df <- do.call(cbind, chromvar.dfs)
dim(chromvar.df)


####################################### Plot specified genes - chromvar run separately #################################

# tf.motifs <- c("HOXB7", , "SPI1", "CEBPA", "PAX5", "EOMES", "TCF7L2", "HNF4G")
tf.motifs <- c("NFIC", "STAT3", "REL", "GATA1", "SPI1", "CEBPA", "EOMES", "TCF7L2")

# m <- merge(atac.obj@meta.data, t(chromvar.df), by = "row.names")
m <- merge(atac.obj@meta.data, atac.obj@assays$chromvar_celltype@data %>% t(), by = "row.names")
m[1:5, 1:5]

plist <- lapply(tf.motifs, function(x) {
  m %>% 
    dplyr::rename(current_tf = x) %>% 
    mutate(CellType = factor(CellType, levels = celltype.order)) %>% 
    ggplot(aes(x = CellType, y = current_tf, fill = CellType)) +
    geom_boxplot(outlier.size = 0.1) +
    ylab("chromVAR score") +
    xlab(element_blank()) +
    ggtitle(x) +
    scale_fill_manual(values = cell.type.palette) + theme(legend.position = "none") +
    geom_hline(yintercept = 0, linetype = "longdash") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
})

names(plist) <- tf.motifs
p.tfs <- wrap_plots(plist, ncol = 2)
p.tfs
ggsave("figures/current_figure_drafts/ATAC_cluster_markers_chromvar.pdf", plot = p.tfs, dpi = 300, width = 8, height = 10)



 ####################################### Plot specified genes - chromvar run together #################################

tf.motifs <- c("PBX1", "GATA1", "SPI1", "CEBPA", "PAX5", "EOMES", "TCF7L2", "HNF4G")

m <- atac.obj@assays$chromvar_celltype@data %>% as.matrix() %>% t()
m <- merge(m, atac.obj@meta.data, by = "row.names")
m[1:5, 1:5]

plist <- lapply(tf.motifs, function(x) {
  m %>% 
    dplyr::rename(current_tf = x) %>% 
    mutate(CellType = factor(CellType, levels = celltype.order)) %>% 
    ggplot(aes(x = CellType, y = current_tf, fill = CellType)) +
    geom_violin(outer())+
    ylab(x) +
    xlab(element_blank()) +
    ggtitle(x) +
    scale_fill_manual(values = cell.type.palette) + theme(legend.position = "none") +
    geom_hline(yintercept = 0, linetype = "longdash") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
})
names(plist) <- tf.motifs
plist$GATA1
p.tfs <- wrap_plots(plist, ncol = 2)
p.tfs
ggsave("figures/current_figure_drafts/ATAC_cluster_markers_chromvar.pdf", plot = p.tfs, dpi = 300, width = 8, height = 17)




####################################### Extra ###########################################
da.tf.motifs <- read_csv("data/chromVAR/JASPAR2020/cluster_markers/landau_and_NIH_09_13_23_resolution_1_2_cell_types_chromvar.csv")

## Pull the top N unique tfs for each cluster
num.markers.to.pull = 5
top.genes <- da.tf.motifs %>%
  arrange(p_val) %>% 
  distinct(gene, .keep_all = T) %>%
  mutate(cluster = factor(cluster, levels = celltype.order)) %>%
  group_by(cluster) %>% 
  arrange(p_val) %>% 
  slice_min(n = num.markers.to.pull, order_by = p_val) %>% 
  slice_max(n = num.markers.to.pull, order_by = avg_diff) %>% View
  pull(gene)

## Downsample to 50 BCs per cluster
downsampled.bc.list <- atac.obj@meta.data %>% 
  filter(cluster_celltype != "cycling") %>% 
  filter(!is.na(cluster_celltype)) %>% 
  rownames_to_column("Barcode") %>% 
  group_by(cluster_celltype) %>% 
  slice_sample(n = 50, replace = T) %>% 
  pull(Barcode)

## Pull metadata for annotation row
annotation.df <- atac.obj@meta.data[downsampled.bc.list, ] %>% 
  select(cluster_celltype)

## Heatmap
cluster.colors <- cell.type.palette[levels(atac.obj$CellType)]
pheatmap::pheatmap(atac.obj@assays$chromvar@data[top.genes, downsampled.bc.list], 
                     annotation_col = annotation.df,
                     annotation_colors = list(`cluster_celltype` = cluster.colors),
                     show_colnames = F, cluster_rows = F, cluster_cols = F,
                     color=colorRampPalette(c("navy", "white", "red"))(50),  breaks = seq(-2, 2, length.out = 50))
