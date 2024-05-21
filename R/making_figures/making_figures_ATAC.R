library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(ggh4x)
library(ggsci)
library(gridExtra)
library(tidyverse)
library(Matrix)
library(ggpubr)
library(ggrepel)
library(patchwork)
theme_set(theme_classic())

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")
source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")
source("~/rscripts/general_helpers/custom_hypergeometric_helpers.R")

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")


############################# Making a metadata data frame with coordinates #############################

## Set an order for the cluster annotations
celltype.order <- c("HSC", "LMPP", "EMP", "Early Eryth", "Late Eryth", "MkP", "BaEoMa", "CLP", "GMP", "CD14 Mono", "B", "Plasma", "pre-mDC", "pDC", "MAIT", "CD4 T", "CD8 T", "NK") ## Ignore the Stromal cluster
atac.obj$CellType <- factor(atac.obj$cluster_celltype, levels = celltype.order)

## Formatting the sample and donor columns
meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_2024-03-29.csv")
atac.obj$Sample <- plyr::mapvalues(atac.obj$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$updated_name)
atac.obj$Sample <- factor(atac.obj$Sample, levels = meta.data.mapping$updated_name %>% sort())
atac.obj@meta.data %>% count(orig.ident, Sample)
atac.obj$Donor <- gsub("_BMMC|_CD34.*", "", atac.obj$Sample) %>% gsub("_[A|B]", "", .)
atac.obj$Donor <- factor(atac.obj$Donor, levels = unique(atac.obj$Donor) %>% sort())
atac.obj@meta.data %>% count(orig.ident, Sample, Donor)

## Get metadata plot with coordinates
md <- atac.obj[[]]
coords <- Embeddings(atac.obj[["umap"]]) 
md <- cbind(md, coords) %>% dplyr::rename(`UMAP1` = umap_1) %>% dplyr::rename(`UMAP2` = umap_2)
md <- md[sample(1:nrow(md)), ] ## Shuffle the rows in the dataframe, so plotting order is random

## To make small arrows in UMAP corner only
axis <- ggh4x::guide_axis_truncated(
  trunc_lower = unit(0, "npc"),
  trunc_upper = unit(1, "in")
)

## Fix the het-NA problem in the Genotyping column
atac.obj$Genotype <- atac.obj$genotype_pred_manual
atac.obj$Genotype[!(atac.obj$Genotype %in% c("MUT", "WT")) ] <- "NA"
atac.obj$Genotype <- factor(atac.obj$Genotype, levels = c("WT", "MUT", "NA"))
atac.obj@meta.data %>% 
  count(Genotype, genotype_pred_manual)


##################################### Saving metadata for tables ####################################################

## Donor, genotype
atac.obj@meta.data %>% 
  group_by(Sample, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped) %>% 
  write_csv("data/metadata/cell_counts_donor_genotype.csv")

## Cell type, genotype
atac.obj@meta.data %>% 
  group_by(CellType, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  mutate(total = sum(MUT, WT, `NA`, na.rm = T)) %>% 
  mutate(frac_genotyped = sum(MUT, WT, na.rm = T) / total) %>% 
  mutate(mut_cell_frequency = MUT / sum(MUT, WT)) %>% 
  arrange(-frac_genotyped) %>% 
  write_csv("data/metadata/cell_counts_celltype_genotype.csv")


## Cell type, genotype, donor
atac.obj@meta.data %>% 
  filter(CellType == "HSC") %>% 
  group_by(Donor, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  mutate(mut_cell_frequency = MUT / sum(MUT, WT)) %>% 
  ungroup() %>% 
  arrange(-total) %>% 
  write_csv("data/metadata/cell_counts_celltype_genotype_HSCs.csv")

#################################### Making plots ########################################

# ## Plot orig.ident
# p.sample.origin <- ggplot(md, aes(x = UMAP1, y = UMAP2, color = Sample )) +
#   geom_point(size = 1) + theme_classic() + coord_fixed() +
#   guides(x = axis, y = axis) +
#   scale_x_continuous(breaks = NULL) +
#   scale_y_continuous(breaks = NULL) +
#   theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))
#   # scale_color_manual(values = ggsci::pal_aaas()(length(unique(md$Sample))))
# p.sample.origin
# ggsave("figures/current_figure_drafts/UMAP_ATAC_donor.pdf", plot = p.sample.origin, device = "pdf", dpi = 300, width = 7, height = 5, unit = "in")

## Sample of origin, split by donor
p.sample.origin.split <- DimPlot(atac.obj, split.by = "Sample", group.by = "Sample", ncol = 5, pt.size = 0.1) + 
  coord_fixed()+
  ggtitle(element_blank()) + theme(legend.position = "none")
  # scale_color_manual(values = ggsci::pal_aaas()(length(unique(md$Sample))))
p.sample.origin.split
ggsave("figures/current_figure_drafts/UMAP_ATAC_sample_split_20240328.pdf", plot = p.sample.origin.split, device = "pdf", dpi = 300, width = 10, height = 8, unit = "in")


## Plot QC metrics - number of fragments
p.fragments <- md %>% 
  ggplot(aes(x = reorder(Donor, -log10(nCount_ATAC)), y = log10(nCount_ATAC), fill = Donor)) +
  geom_boxplot() + 
  theme(legend.position = "none") +
  # geom_hline(yintercept = log10(1500), linetype = "longdash") + 
  ylab("ATAC fragments (log10)") +
  xlab("Patient") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p.fragments
ggsave("figures/current_figure_drafts/ATAC_num_fragments_sorted_20240420.pdf", plot = p.fragments, device = "pdf", dpi = 300, width = 4, height = 3, unit = "in")

## Plot QC metrics - nucleosome signal
p.nucleosome <- md %>% 
  ggplot(aes(x = reorder(Donor, -nucleosome_signal), y = nucleosome_signal, fill = Donor)) +
  geom_boxplot() + 
  theme(legend.position = "none") +
  ylab("Nucleosome signal") +
  xlab("Patient") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
p.nucleosome
ggsave("figures/current_figure_drafts/ATAC_nucleosome_signal_20240328.pdf", plot = p.nucleosome, device = "pdf", dpi = 300, width = 4, height = 3, unit = "in")


## Plot QC metrics - fragment length
p.frag.length <- FragmentHistogram(atac.obj, group.by = "Donor", assay = "ATAC") + NoLegend()
p.frag.length
ggsave("figures/current_figure_drafts/ATAC_fragment_length_20240328.pdf", plot = p.frag.length, device = "pdf", dpi = 300, width = 5, height = 5, unit = "in")

## Plot the azimuth labels
to.plot.tags <- atac.obj@meta.data %>% group_by(predicted.l2) %>% count() %>% filter(n > 50) %>% pull(predicted.l2)
atac.obj$predicted.celltype.l2.plot <- atac.obj$predicted.l2
atac.obj$predicted.celltype.l2.plot[!(atac.obj$predicted.celltype.l2.plot %in% to.plot.tags)] <- "Other"
p.azimuth <- DimPlot(atac.obj, group.by = "predicted.celltype.l2.plot", label = T, label.box = T, repel = T, pt.size = 0.1) + NoLegend() + 
  coord_fixed() +
  ggtitle(element_blank()) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) 
ggsave("figures/current_figure_drafts/UMAP_ATAC_azimuth.pdf", plot = p.azimuth, device = "pdf", dpi = 300, width = 7, height = 5, unit = "in")


## Plot the azimuth label/predicted celltype label overlap
predicted.l2.order <- c("HSC" ,  "LMPP"   ,   "EMP" , "Early Eryth", "Late Eryth",  "MkP", "BaEoMa",   "CLP", "GMP", "CD14 Mono", "B", "Plasma",  "pre-mDC", "pre-pDC", "pDC" , "cDC", "CD4 T", "CD8 T", "NK", "CD16 Mono", "Macrophage", "ASDC", "MAIT", "NK CD56+", "Stromal")
p.azimuth.comparison <- md %>% 
  filter(!is.na(predicted.l2)) %>% ## Excludes 12 cells
  mutate(predicted.l2 = gsub(".*CD8.*", "CD8 T", predicted.l2)) %>%
  mutate(predicted.l2 = gsub(".*CD4.*", "CD4 T", predicted.l2)) %>%
  mutate(predicted.l2 = gsub(".* B", "B", predicted.l2)) %>% 
  mutate(predicted.l2 = gsub("Memory ", "", predicted.l2)) %>%
  mutate(predicted.l2 = gsub("Prog Mk", "MkP", predicted.l2)) %>%
  mutate(predicted.l2 = gsub("cDC.*", "cDC", predicted.l2)) %>%
  group_by(CellType) %>% 
  count(predicted.l2) %>% 
  mutate(percent = (n / sum(n))*100) %>% select(-n) %>% 
  pivot_wider(names_from = CellType, values_from = percent, values_fill = 0) %>% 
  column_to_rownames("predicted.l2")
pdf("figures/current_figure_drafts/ATAC_cell_type_labels_heatmap.pdf", width = 5, height = 5)
pheatmap::pheatmap(p.azimuth.comparison[predicted.l2.order, ], cluster_rows = F, cluster_cols = F, main = "Assigned cell type labels", color = inferno(50))
dev.off()

## Highlighting Phospho-seq samples
cols.highlight = list(
  `PT17_all` =  "#A6C9E0",
  `PT17_HSPC` = "#3A71B1",
  `PT16_all` =  "#7F7EB7",
  `PT16_HSPC` = "#66539E"
)
hspc.clusters <- c("HSC", "EMP", "LMPP", "MkP")
md.phospho.hspcs <- md %>% filter(Sample %in% c("PT16", "PT17") & (cluster_celltype %in% hspc.clusters)) %>% mutate(category = paste0(Sample, "_", "HSPC"))
md.phospho.all <- md %>% filter(Sample %in% c("PT16", "PT17") & !(cluster_celltype %in% hspc.clusters)) %>% mutate(category = paste0(Sample, "_", "Other"))
md.na <- md %>% filter(!(Sample %in% c("PT16", "PT17")))
p.phospho <- ggplot(md.na, aes(x = UMAP1, y = UMAP2)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.phospho.all, aes(fill = category), pch = 21, size = 1, stroke = 0.1, color = "black") +
  geom_point(data = md.phospho.hspcs, aes(fill = category), pch = 21, size = 1, stroke = 0.1, color = "black") +
  scale_fill_manual(values = c(PT16_Other = "#CEADD1", PT17_Other = "#A5BFC5", PT16_HSPC = "#7C3B86", PT17_HSPC = "#2C858C")) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.phospho
ggsave("figures/current_figure_drafts/UMAP_ATAC_phospho_sample_highlight_20240328.tiff", device = "tiff", width = 5, height = 4)

#################################### Making QC figures ########################################

ordered.sample.list <- c("BM5_LLL_CD34_plus", "BM6_LLL_CD34_plus", "BM7_LLL_CD34_plus", "VEX_BM8_POS_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM11_LIN_NEG_ATAC", "SG_VX16_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC")
## Print QC scatter plot
plist <- lapply(ordered.sample.list, function(x) {
  if (x == "SG_VX17_ATAC") {
    md %>% 
      filter(orig.ident == x) %>%
      filter(TSS.enrichment < 20) %>% 
      ggplot(aes(x = log10(nCount_ATAC), y = TSS.enrichment)) +
      geom_bin2d(bins = 70) +
      scale_fill_continuous(type = "viridis") +
      ggtitle(x)
  } else {
    md %>% 
      filter(orig.ident == x) %>% 
      ggplot(aes(x = log10(nCount_ATAC), y = TSS.enrichment)) +
      geom_bin2d(bins = 70) +
      scale_fill_continuous(type = "viridis") +
      ggtitle(x)
  }

  
})
p.combined <- wrap_plots(plist, ncol = 3)
ggsave(filename = "figures/current_figure_drafts/QC_heatmaps_final_object.pdf", plot = p.combined, width = 10, height = 7, dpi = 300)

## Print QC violin plots
plist <- lapply(ordered.sample.list, function(x) {
  p1 <- md %>% 
    filter(orig.ident == x) %>% 
    ggplot(aes(x = "Sample", y = pct_reads_in_peaks)) +
    geom_violin(fill = "dodgerblue3") + xlab(element_blank()) + ggtitle("Percent reads\nin peaks") + coord_cartesian(ylim = c(0, 100))
  p2 <- md %>% 
    filter(orig.ident == x) %>% 
    ggplot(aes(x = "Sample", y = nucleosome_signal)) +
    geom_violin(fill = "dodgerblue3") + xlab(element_blank()) + ggtitle("Nucleosome\nsignal")
  p3 <- md %>% 
    filter(orig.ident == x) %>% 
    ggplot(aes(x = "Sample", y = blacklist_ratio)) +
    geom_violin(fill = "dodgerblue3") + xlab(element_blank()) + ggtitle("Blacklist\nratio") 
  p.combined <- (p1 | p2 | p3) + plot_annotation(title = x)
  return(patchwork::wrap_elements(p.combined))
})
p.combined <- wrap_plots(plist, ncol = 3)
ggsave(filename = paste0("figures/current_figure_drafts/QC_violins_final_object.pdf"), plot = p.combined, width = 15.5, height = 8)



#################################### Making plots - genotyping ########################################

## Print the genotyping fraction (overall)
atac.obj@meta.data %>% 
  count(Genotype) %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped)

## Print the genotyping fraction (by sample)
atac.obj@meta.data %>% 
  group_by(Donor, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped)

## Print the genotyping fraction (by sample, stat summaries)
atac.obj@meta.data %>% 
  group_by(Donor, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  ungroup() %>% 
  summarize(mean_frac_genotyped = mean(frac_genotyped), sd = sd(frac_genotyped))

## Print the genotyping fraction (HSCs)
atac.obj@meta.data %>%
  filter(CellType == "HSC") %>% 
  group_by(Donor, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped) %>% write_csv("data/genotyping_metrics_ATAC_HSCs.csv")

## Print the genotyping fraction (HSPCs)
atac.obj@meta.data %>%
  filter(CellType %in% c("HSC", "EMP", "LMPP", "MkP", "GMP")) %>% 
  group_by(Donor, Genotype) %>% 
  count() %>% 
  pivot_wider(names_from = Genotype, values_from = n, values_fill = 0) %>% 
  mutate(total = sum(MUT, WT, `NA`)) %>% 
  mutate(frac_genotyped = sum(MUT, WT) / total) %>% 
  arrange(-frac_genotyped) %>% write_csv("data/genotyping_metrics_ATAC_HSPCs.csv")



## Genotyping UMAP (points)
mut.count <- md %>% filter(Genotype == "MUT") %>% nrow()
wt.count <- md %>% filter(Genotype == "WT") %>% nrow()
na.or.het.count <- md %>% filter(Genotype == "NA") %>% nrow()
mut.title <- paste0("MUT (n = ", format(mut.count, big.mark = ","), ")")
wt.title <- paste0("WT (n = ", format(wt.count, big.mark = ","), ")")
na.title <- paste0("NA (n = ", format(na.or.het.count, big.mark = ","), ")")
md.wt <- md[md$Genotype %in% c("WT"), ]
md.mut <- md[md$Genotype %in% c("MUT"), ]
md.genotyped <- md[md$Genotype %in% c("MUT", "WT"), ]
md.na <- md[md$Genotype %in% c("NA", "Control"), ]

p.genotyping <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, fill = Genotype)) +
  geom_point(shape = 19, size = 1.5, color = "grey80") + 
  geom_point(data = md.mut, shape = 21, size = 1, fill = vexas.genotyping.palette[["MUT"]], stroke = 0.1, color = "black") +
  geom_point(data = md.wt, shape = 21, size = 1, fill = vexas.genotyping.palette[["WT"]], stroke = 0.1, color = "black") +
  # geom_point(data = md.genotyped, shape = 21, size = 1) +
  # scale_color_manual(values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.genotyping
ggsave("figures/current_figure_drafts/UMAP_ATAC_genotyping.tiff", plot = p.genotyping, device = "tiff", dpi = 300, width = 5, height = 4, unit = "in")

p.genotyping <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, fill = Genotype)) +
  geom_point(shape = 19, size = 0.25, color = "grey80") + 
  geom_point(data = md.mut, shape = 21, size = 0.7, fill = vexas.genotyping.palette[["MUT"]], size = 0.5) +
  geom_point(data = md.wt, shape = 21, size = 0.7, fill = vexas.genotyping.palette[["WT"]], size = 0.5) +
  # geom_point(data = md.genotyped, shape = 21, size = 1) +
  scale_color_manual(values = unlist(vexas.genotyping.palette)) +
  theme_classic() + coord_fixed() +
  theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.genotyping
ggsave("figures/current_figure_drafts/UMAP_ATAC_genotyping.pdf", plot = p.genotyping, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Genotyping UMAP, split by patient
per.patient.plots <- lapply(levels(atac.obj$Donor), function(x) {
  print(x)
  md.mut.donor <- md.mut %>% filter(Donor == x)
  md.wt.donor <- md.wt %>% filter(Donor == x)
  fraction.gentotyped <- md %>% filter(Donor == x) %>% filter(Genotype %in% c("MUT", "WT")) %>% nrow() /  md %>% filter(Donor == x) %>% nrow()
  percent.genotyped <- round(fraction.gentotyped, digits = 4) *100
  p.donor.genotyped <- ggplot(md.na, aes(x = UMAP1, y = UMAP2, color = Genotype)) +
    geom_point(shape = 19, size = 1.5, color = "black") + 
    geom_point(shape = 19, size = 1, color = "grey80")+
    geom_point(data = md.mut.donor, shape = 19, size = 0.5, color = vexas.genotyping.palette[["MUT"]]  ) +
    geom_point(data = md.wt.donor, shape = 19, size = 0.5, color = vexas.genotyping.palette[["WT"]]) +
    scale_color_manual(values = unlist(vexas.genotyping.palette)) +
    theme_classic() + coord_fixed() +
    theme(legend.position = c(1, 0.9), legend.title = element_blank(), legend.text = element_text(size = 12)) +
    guides(x = axis, y = axis) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    ggtitle(x, subtitle = paste0("Cells genotyped = ", percent.genotyped, "%")) +
    theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) + coord_fixed() 
  return(p.donor.genotyped)
}) 
p.combined <- wrap_plots(per.patient.plots, nrow = 2)
ggsave("figures/current_figure_drafts/ATAC_UMAP_genotyping_split_20240328.tiff", plot = p.combined, device = "tiff", dpi = 300, width = 12, height = 6)


## Print the overall genotyping fraction
atac.obj@meta.data %>% 
  count(Genotype) %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(genotyping_fraction = (MUT + WT) / (MUT + WT + `NA`)) %>% 
  arrange(-genotyping_fraction)

atac.obj@meta.data %>% 
  count(Donor, Genotype) %>% 
  pivot_wider(names_from = Genotype, values_from = n) %>% 
  mutate(genotyping_fraction = (MUT + WT) / (MUT + WT + `NA`)) %>% 
  arrange(-genotyping_fraction)

## Plot genotyping - heatmaps
wt<-ggplot(md[which(md$Genotype =="WT"),], aes(x=UMAP1, y=UMAP2))+
  stat_density2d(aes(alpha=..level.., fill=..level..), size=1,
                 bins=10, geom="polygon") +
  scale_fill_gradient(low = "white", high = "#253494") +
  scale_alpha(range = c(0.00, 1), guide = FALSE) +
  xlim(-11,11)+
  ylim(-11,11)+
  guides(alpha=FALSE) +
  theme_classic(base_size = 17) +
  theme(legend.position="none",axis.line=element_blank(), axis.text=element_blank())
mut<-ggplot(md[which(md$Genotype =="MUT"),], aes(x=UMAP_1, y=UMAP_2))+
  stat_density2d(aes(alpha=..level.., fill=..level..), size=1,
                 bins=10, geom="polygon") +
  scale_fill_gradient(low = "white", high = "red") +
  scale_alpha(range = c(0.00, 1), guide = FALSE) +
  xlim(-11,11)+
  ylim(-11,11)+
  guides(alpha=FALSE) +
  theme_classic(base_size = 17) +
  theme(legend.position="none",axis.line=element_blank(), axis.text=element_blank())
grid.arrange(wt,mut,ncol=2)

## Plot celltypes, with labels
p.celltypes <- DimPlot(atac.obj, group.by = "CellType", cols = cell.type.palette, label = T, label.box = T)  + NoLegend() +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0)) + coord_fixed()
p.celltypes
ggsave("figures/current_figure_drafts/UMAP_ATAC_celltypes.pdf", plot = p.celltypes, device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Plot celltypes (with outlines)
p.celltypes <- ggplot(md, aes(x = UMAP1, y = UMAP2, fill = CellType)) +
  geom_point(pch = 21, size = 1, stroke = 0.1, color = "black") + ## size controls width of point, stroke controls width of border, color is border color
  scale_fill_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0), legend.position = "none")  + coord_fixed()
p.celltypes
ggsave("figures/current_figure_drafts/UMAP_ATAC_celltypes_no_labels_colors_v2.pdf", device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")

## Plot celltypes (no labels, outlines)
p.celltypes <- ggplot(md, aes(x = UMAP1, y = UMAP2, color = CellType)) +
  geom_point(pch = 19, size = 1) + ## size controls width of point, stroke controls width of border, color is border color
  scale_color_manual(values = cell.type.palette) +
  theme_classic() + coord_fixed() +
  theme(legend.position = "right", legend.title = element_blank(), legend.text = element_text(size = 12)) +
  guides(x = axis, y = axis) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(axis.line = element_line(arrow = arrow()), axis.title = element_text(hjust = 0))  + coord_fixed()
p.celltypes
ggsave("figures/current_figure_drafts/UMAP_ATAC_celltypes_no_labels.pdf", device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")


########################################### Making coverage plots #############################################

## Plotting coverage track of a few important genes
atac.obj.subset <- subset(atac.obj, CellType == "HSC" & Genotype %in% c("MUT", "WT"))
atac.obj.subset$Genotype_format <- gsub("WT", "_WT", atac.obj.subset$Genotype)
Idents(atac.obj.subset) <- "Genotype_format"
# vexas.genotyping.palette.custom <- c(MUT = vexas.genotyping.palette$MUT, `_WT` = vexas.genotyping.palette$

CoveragePlot(atac.obj.subset, region = "chrX-47198691-47199456", group.by = "Genotype", assay = "ATAC", extend.upstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_UBA1_mutation_site.pdf", device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")

## Plot DDIT4
CoveragePlot(atac.obj.subset, region = "chr10-72273000-72274500", group.by = "Genotype", assay = "ATAC", ymax = 24, peaks = F) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_DDIT4.pdf", device = "pdf", dpi = 300, width = 4, height = 5, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr10-72273000-72274500", group.by = "Genotype", assay = "ATAC", idents = c("MUT"), peaks = F, ymax = 22) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_DDIT4_MUT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr10-72273000-72274500", group.by = "Genotype", assay = "ATAC", idents = c("WT"), peaks = F, ymax = 22) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_DDIT4_WT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")

## Plot CXCL3
CoveragePlot(atac.obj.subset, region = "chr4-74038000-74039500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 10) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3.pdf", device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr4-74038000-74039500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 10, idents = c("MUT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3_MUT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr4-74038000-74039500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 10, idents = c("WT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3_WT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")

## Plot IFITM3
CoveragePlot(atac.obj.subset, region = "chr11-326000-328000", group.by = "Genotype", assay = "ATAC", peaks = F) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3.pdf", device = "pdf", dpi = 300, width = 4, height = 4, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr4-74038000-74039500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 10, idents = c("MUT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3_MUT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr4-74038000-74039500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 10, idents = c("WT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_CXCL3_WT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in")


## Plot IL1B
CoveragePlot(atac.obj.subset, region = "chr2-112836000-112837500", group.by = "Genotype", assay = "ATAC", peaks = F, ymax = 8) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_IL1B.pdf", device = "pdf", dpi = 300, width = 4, height = 5, unit = "in")
CoveragePlot(atac.obj.subset, region = "chr2-112836000-112837500", group.by = "Genotype_format", assay = "ATAC", peaks = F, ymax = 8, idents = c("MUT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_IL1B_MUT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in", aspect.ratio = 1)
CoveragePlot(atac.obj.subset, region = "chr2-112836000-112837500", group.by = "Genotype_format", assay = "ATAC", peaks = F, ymax = 8, idents = c("WT")) & scale_fill_manual(values = vexas.genotyping.palette)
ggsave("figures/current_figure_drafts/coverage_plot_IL1B_WT.pdf", device = "pdf", dpi = 300, width = 4, height = 3.5, unit = "in", )


## Extra plots of other regions of interest
CoveragePlot(atac.obj.subset, region = "chr13-28098592-28100592", group.by = "Genotype", assay = "ATAC", extend.downstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette) # FLT3
CoveragePlot(atac.obj.subset, region = "NLRP1", group.by = "Genotype", assay = "ATAC", extend.downstream = 2000) ## NLRP1
CoveragePlot(atac.obj.subset, region = "chr17-5615424-5619424", group.by = "Genotype", assay = "ATAC", extend.downstream = 2000) ## NLRP1
CoveragePlot(atac.obj.subset, region = "NLRP3", group.by = "Genotype", assay = "ATAC", extend.upstream = 5000)
CoveragePlot(atac.obj.subset, region = "ATF4", group.by = "Genotype", assay = "ATAC", extend.upstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette)
CoveragePlot(atac.obj.subset, region = "chr22-39518000-39520000", group.by = "Genotype", assay = "ATAC", extend.upstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette)
CoveragePlot(atac.obj.subset, region = "chr22-39518000-39521000", group.by = "Genotype", assay = "ATAC") & scale_fill_manual(values = vexas.genotyping.palette)
CoveragePlot(atac.obj.subset, region = "MCL1", group.by = "Genotype", assay = "ATAC", extend.downstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette)
CoveragePlot(atac.obj.subset, region = "HLA-DRA", group.by = "Genotype", assay = "ATAC", extend.upstream = 2000) & scale_fill_manual(values = vexas.genotyping.palette)
CoveragePlot(atac.obj.subset, region = "chr5-52987000-52991000", group.by = "Genotype", assay = "ATAC") & scale_fill_manual(values = vexas.genotyping.palette)


###################################### Plots - cell type proportion bar plots ################################################

## Bar plot of cell type proportions - per smample
p.celltype.fractions.per.donor <- md %>% 
  group_by(Sample) %>% 
  mutate(total_per_sample = n()) %>% 
  group_by(Sample, total_per_sample) %>% 
  count(CellType) %>% 
  mutate(fraction = n / total_per_sample) %>% 
  ggplot(aes(x = Sample, y = fraction, fill = CellType)) +
  geom_col() +
  scale_fill_manual(values = cell.type.palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  xlab(element_blank()) +
  ylab("Percent of total cells") +
  scale_y_continuous(labels = scales::percent)
p.celltype.fractions.per.donor
ggsave("figures/current_figure_drafts/UMAP_ATAC_celltype_fractions_per_donor.pdf", device = "pdf", dpi = 300, width = 7, height = 7, unit = "in")


############################################### MUT cell frequency per cell type, comparing with RNA ##########################################


## Summarize the MUT cell frequency per donor
atac.mut.freqs <- md %>% 
  filter(CellType != "Other") %>%
  filter(Genotype %in% c("MUT", "WT")) %>%
  group_by(CellType) %>% 
  filter(n() > 50) %>% ## Only including cell types with > 20 genotyped cells
  group_by(CellType, Donor) %>% 
  filter(n() > 10) %>% ## Only including patient data points with > 10 genotyped cells
  group_by(CellType, Genotype, Donor) %>%
  count() %>%
  group_by(CellType, Donor) %>% 
  mutate(mut_fraction = n / sum(n)) %>% 
  filter(Genotype == "MUT") %>% 
  group_by(CellType) %>%
  filter(n_distinct(Donor) >= 2) %>% 
  summarize(mean_val = mean(mut_fraction), sd_val = sd(mut_fraction), n = n(), se_val = sd_val / sqrt(n)) %>% 
  mutate(Assay = "ATAC")

## Make the bar plot
p.mutant_cell_fraction.bar.plot <- atac.mut.freqs %>% 
  ggplot(aes(x = reorder(CellType, -mean_val), y = mean_val, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(x = CellType, ymin = mean_val-sd_val, ymax = mean_val+sd_val), width=0.4) +
  scale_fill_manual(values = cell.type.palette) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  labs(y = "Mutant cell fraction", x = element_blank()) + coord_fixed(9)
p.mutant_cell_fraction.bar.plot
ggsave("figures/current_figure_drafts/ATAC_MUT_cell_freq_bar_plot.pdf", plot = p.mutant_cell_fraction.bar.plot, device = "pdf", dpi = 300, width = 6, height = 4, unit = "in")

## Loading RNA data
rna.mut.freqs <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/metadata/mut_cell_fractions_RNA.csv") %>% 
  mutate(Assay = "RNA") %>% 
  select(-median_pseudotime)

mut.freqs <- rbind(rna.mut.freqs, atac.mut.freqs)
celltype.order <-  rna.mut.freqs %>% arrange(-mean_val) %>% pull(CellType)

## Make a scatter plot
rna.mut.freqs <- rna.mut.freqs %>% dplyr::rename(mean_val_RNA = mean_val) %>% dplyr::rename(sd_val_RNA = sd_val) %>% select(-Assay)
atac.mut.freqs <- atac.mut.freqs %>% dplyr::rename(mean_val_ATAC = mean_val) %>% dplyr::rename(sd_val_ATAC = sd_val) %>% select(-Assay)
mut.freqs.merged <- merge(rna.mut.freqs, atac.mut.freqs, by = "CellType")
adjusted.colors <- cell.type.palette
adjusted.colors[["BaEoMa"]] <- "gold3"
adjusted.colors[["MkP"]] <- "goldenrod3"
p.mutant_cell_fraction.scatter.plot <- mut.freqs.merged %>% 
  ggplot(aes(x = mean_val_RNA, y = mean_val_ATAC)) +
  geom_point(size = 0) +
  # geom_abline(linetype = "longdash", color = "grey70") +
  stat_cor(aes(x = mean_val_RNA, y = mean_val_ATAC), method = "pearson") +
  geom_smooth(method='lm', formula= y~x, color = "black") +
  geom_point(aes(color = CellType), size = 4) +
  scale_color_manual(aes(color = CellType, label = CellType), values = adjusted.colors) +
  geom_text_repel(aes(color = CellType, label = CellType), point.padding = 15, size = 4) +
  # geom_crossbar(mapping = aes(x = mean_val_RNA,
  #                      ymin = mean_val_ATAC - sd_val_ATAC,
  #                      ymax = mean_val_ATAC + sd_val_ATAC, color = CellType)) +
  # geom_crossbar(mapping = aes(y = mean_val_ATAC,
  #                       xmin = mean_val_RNA - sd_val_RNA,
  #                       xmax = mean_val_RNA + sd_val_RNA, color = CellType)) +
  theme(legend.position = "none") +
  ggtitle("Mutant cell frequencies", subtitle = "Genotyping from RNA (GoT) vs. ATAC (GoTChA)") +
  coord_fixed() +
  xlab("Mean MUT cell fraction,\nRNA genotyping") +
  ylab("Mean MUT cell fraction,\nATAC genotyping")
  # coord_cartesian(xlim=c(0, 1), ylim = c(0, 1))
p.mutant_cell_fraction.scatter.plot
ggsave("figures/current_figure_drafts/MUT_cell_freq_scatter_plot.pdf", plot = p.mutant_cell_fraction.scatter.plot, device = "pdf", dpi = 300, width = 4.5, unit = "in")


############################################### HSC - gene score volcano plot #############################

# de.results <- read_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_DiffLMM2_expressed_only_5p_RNA_obj.csv") %>% 
#   mutate(avg_log2FC = log2(fc))

de.results <- read_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_genotyping_updated_DiffLMM2_20240304.csv") %>%
  mutate(avg_log2FC = log2(fc))

plot_volcano(de.results, effect_line = 0.25, stat_line = 0.05, stat_column = "pval", max.overlaps = 20) + ggtitle("HSC gene scores, MUT vs. WT") | 
  plot_volcano(de.results, effect_line = 0.25, stat_line = 0.1, stat_column = "fdr") + ggtitle("HSC gene scores, MUT vs. WT")
ggsave("figures/current_figure_drafts/gene_score_gsea_20240125.pdf", width = 12, height = 10)

plot_volcano(de.results, effect_line = 0.25, stat_line = 0.05, stat_column = "fdr") + ggtitle("HSC gene scores, MUT vs. WT") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0, 5)) + theme(legend.position = "none")
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_gene_score_fdr_20240223.pdf", plot = last_plot(), device = "pdf", width = 6, height = 6)
plot_volcano(de.results, effect_line = 0.25, stat_line = 0.05, stat_column = "pval") + ggtitle("HSC gene scores, MUT vs. WT") +
  coord_cartesian(xlim = c(-1.5, 1.5), ylim = c(0, 5)) + theme(legend.position = "none")
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_gene_score_pval_20240223.pdf", plot = last_plot(), device = "pdf", width = 6, height = 6)

atac.obj.hsc <- subset(atac.obj, cluster_celltype == "HSC" & genotype_pred_manual %in% c("MUT", "WT") & orig.ident %in% c("VEX_BM8_POS_ATAC", "VEX_BM10_POS_ATAC", "SG_VX16_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC"))
atac.obj.hsc <- subset(atac.obj.hsc, nCount_ATAC > 1500 & pct_reads_in_peaks > 50)
atac.obj.hsc@meta.data %>% count(Sample, Genotype) %>% print()
atac.obj.hsc@meta.data %>% count(Sample) %>% print()
VlnPlot(atac.obj.hsc, features = c("nCount_ATAC", "TSS.enrichment", "pct_reads_in_peaks"), group.by = "orig.ident")
VlnPlot(atac.obj.hsc, assay = "RNA_ArchR", feature = "DDIT4", split.by = "genotype_pred_manual", group.by = "orig.ident")

fgsea.out <- run_fgsea_for_list(de.results)

## Try re-ranking
de_results.ranked <- de.results %>%
  filter(!is.na(pval)) %>% 
  # mutate(rank = -log10(pval)*sign(avg_log2FC)) %>%
  mutate(rank = avg_log2FC) %>%
  # mutate(rank = -log10(pval)) %>%
  arrange(-rank)
current.ranks.list <- de_results.ranked$rank
names(current.ranks.list) <- de_results.ranked$feature

fgsea.out <- run_fgsea_for_list(de.results, ranks.list = current.ranks.list)
## Save the gsea output
fgsea.out %>% 
     write_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/gsea_output/HSC_DiffLMM2_fold_change_rank_gsea_2024-02-23.csv")

plot_enrichment(ranks.list = current.ranks.list, pathway.name = "HALLMARK_INFLAMMATORY_RESPONSE", pval_col = "pval", pval_show = "padj")

## Adjust with fdr
fgsea.out$fdr <- p.adjust(fgsea.out$pval, method = "fdr")
fgsea.out %>% 
  arrange(pval) %>% head()
# fgsea.out %>% 
#   write_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/gsea_output/HSC_DiffLMM2_10p_20231206_gsea.csv")

## Use a custom set of msigdb gene sets, if not specified
m_t2g <- tibble()
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% filter(gs_name %in% fgsea.pathways))
pathways = split(x = m_t2g$gene_symbol, f = m_t2g$gs_name) # format for fgsea
## Add GoT paper modules
got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv") %>% as.list()
got.modules <- lapply(got.modules, function(x) return(x[!is.na(x)]))
pathways <- c(pathways, got.modules)

## Add Han et al modules
han.modules <- list(
  `ATF4 targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
  `CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% pull(Feature),
  `Shared ATF4/CHOP targets (Han et al.)` = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_and_ATF4.txt") %>% filter(Feature != "N/A") %>% pull(Feature)
)
pathways <- c(pathways, han.modules)

volcano.gene.list <- list(
  `MUT` = de.results %>% slice_min(pval, n = 20) %>% filter(avg_log2FC > 0) %>% pull(feature),
  `WT` = de.results %>% slice_min(pval, n = 20) %>% filter(avg_log2FC < 0) %>% pull(feature)
  # `ATF4` = c(intersect(pathways$`ATF4 targets (Han et al.)`, de.results %>% filter(pval < 0.05) %>% pull(feature))),
  # `Inflammation` = intersect(pathways$HALLMARK_INFLAMMATORY_RESPONSE, de.results %>% filter(pval < 0.05) %>% pull(feature)),
  # `Interferon` = c(union(
  #   intersect(pathways$HALLMARK_INTERFERON_ALPHA_RESPONSE, de.results %>% filter(pval < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
  #   intersect(pathways$HALLMARK_INTERFERON_GAMMA_RESPONSE, de.results %>% filter(pval < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature))
  # )),
  # `Antigen presentation` = c(de.results %>% filter(pval < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature) %>% grep("HLA", ., value = T), "CD74", "B2M", "CD99", "CIITA"),
  # `Translation` = intersect(pathways$REACTOME_TRANSLATION, de.results %>% filter(pval < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
  # `CHOP` = c(intersect(pathways$`CHOP targets (Han et al.)`, de.results %>% filter(pval < 0.05) %>% pull(feature)))
)
# volcano.gene.list <- list(
#   `P53` = intersect(fgsea.out %>% filter(pathway == "HALLMARK_P53_PATHWAY") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#   `Translation` = intersect(fgsea.out %>% filter(pathway == "REACTOME_TRANSLATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature)),
#   `Myeloid bias` = c("CD99", intersect(fgsea.out %>% filter(pathway == "GOBP_LEUKOCYTE_DIFFERENTIATION") %>% pull(leadingEdge) %>% unlist(), de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature))),
#   `Antigen presentation` = c(de.results %>% filter(fdr < 0.05 & abs(avg_log2FC) > 0.2) %>% pull(feature) %>% grep("HLA", ., value = T), "CD74", "B2M")
# )

## Pick colors + plot volcano
volcano.gene.list.colors <- c(
  `all_remaining` = "#5F4690",
  # `Interferon` = "#CC503E",
  `Inflammation` = "#94346E",
  `ATF4` = "#0F8554", 
  # `Translation` = "#5F4690",
  `CHOP` = "#1D6996"
)
p.HSC.volcano <- plot_volcano(de.results, stat_line = 0.05, effect_line = 0.2, stat_column = "pval", title = paste0("HSCs"), 
            only_genes = volcano.gene.list, only_genes_colors = vexas.genotyping.palette, max.overlaps = 20) + 
  theme(legend.position = "none") +
  coord_cartesian(ylim=c(0, 5), xlim = c(-1.5, 1.5)) + ggtitle("HSC differential gene scores")
p.HSC.volcano 
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_gene_scores_20240304.pdf", plot = p.HSC.volcano, device = "pdf", width = 6, height = 5)

## Plot volcano, without labels
p.HSC.volcano <- plot_volcano(de.results, metadata.df =atac.obj@meta.data %>% dplyr::filter(cluster_celltype == "HSC") %>% dplyr::filter(orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC")), 
                              stat_line = 0.05, effect_line = 0.2, stat_column = "pval", title = paste0("HSCs"), 
                              genotyping_column = "Genotype", 
                              max.overlaps = 20, label_genes = F) + 
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 5), xlim = c(-2, 2)) + ggtitle("HSC differential gene scores")
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_gene_score_no_labels.pdf", plot = p.HSC.volcano, device = "pdf", width = 7, height = 6)

## Enrichment plots
p1 <- plot_enrichment(de.results, ranks.list = current.ranks.list, pathway = "HALLMARK_INFLAMMATORY_RESPONSE", pval_show = "padj") + theme_classic() + geom_line(color = "#BA242A", size = 1)
p2 <- plot_enrichment(de.results, pathway = "CHOP_01", pval_show = "pval") + theme_classic() + geom_line(color = "#BA242A", size = 1) + ggtitle("CHOP target genes")
p3 <- plot_enrichment(de.results, pathway = "ATF4 targets (Han et al.)", pval_show = "pval") + geom_line(color = "#BA242A", size = 1) + ggtitle("ATF4 targets (Han et al.)")
p4 <- plot_enrichment(de.results, pathway = "HALLMARK_P53_PATHWAY", pval_show = "pval") + theme_classic()
p.combined <- (p1 / p3 / p2)
p.combined
ggsave("figures/current_figure_drafts/ATAC_gene_score_gsea_inflammatory_response.pdf", plot = p.combined, device = "pdf", width = 4, height = 7)
ggsave("figures/current_figure_drafts/ATAC_gene_score_gsea_inflammatory_response.pdf", plot = p1, device = "pdf", width = 4, height = 2.5)

## Run a hypergeometric test
gene.list <- de.results %>% filter(avg_log2FC < 0 & pval < 0.05) %>% pull(feature)
gene.list <- de.results %>% filter(avg_log2FC > 0) %>% slice_min(order_by = pval, n = 1000) %>% pull(feature)

m_t2g <- tibble()
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "H", subcategory = NULL) %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C5") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:GTRD") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- rbind(m_t2g, msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT:TFT_Legacy") %>% filter(gs_name %in% fgsea.pathways))
m_t2g <- m_t2g %>% select(gs_name, human_entrez_gene) %>% dplyr::rename(TermID = gs_name) %>% dplyr::rename(geneID = human_entrez_gene)

## Add GoT paper modules
got.modules <- read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/GoT_paper_references/modules.csv")
additional.modules <- got.modules %>% pivot_longer(cols = colnames(got.modules), names_to = "TermID", values_to = "geneID") %>% 
  filter(!is.na(geneID))
## Add Han et al modules
atf4.df <- data.frame(geneID = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/ATF4_only.txt") %>% filter(Feature != "N/A") %>% dplyr::rename(geneID = Feature))
atf4.df$TermID = "ATF4 only (Han et al.)"
atf4.df <- atf4.df %>% select(TermID, geneID)
chop.df <- data.frame(geneID = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% dplyr::rename(geneID = Feature))
chop.df$TermID = "CHOP only (Han et al.)"
chop.df <- chop.df %>% select(TermID, geneID)

additional.modules <- rbind(additional.modules, atf4.df)
additional.modules <- rbind(additional.modules, chop.df)

library(org.Hs.eg.db)
library(clusterProfiler)
hs <- org.Hs.eg.db
mapping <- AnnotationDbi::select(hs, 
                                 keys = additional.modules$geneID,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")
additional.modules <- merge(additional.modules, mapping, by.x = "geneID", by.y = "SYMBOL")
additional.modules <- additional.modules %>% dplyr::select(TermID, ENTREZID) %>% dplyr::rename(geneID = ENTREZID) %>% 
  filter(!is.na(geneID))
m_t2g <- rbind(m_t2g, additional.modules)
gene.list.entrez.ids <- AnnotationDbi::select(hs, 
                                 keys = gene.list,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")
res.custom <- enricher(gene.list.entrez.ids$ENTREZID, gson = NULL, TERM2GENE = m_t2g, pvalueCutoff = 1)
res.custom@result %>% head

res.custom@result %>% View

atac.obj.HSC <- subset(atac.obj, subset = cluster_celltype == "HSC")
VlnPlot(atac.obj.HSC, features = "DDIT4", assay = "RNA_ArchR", split.by = "Genotype", group.by = "Sample", cols = c("blue", "red", "grey80"))


enrichr.outs.1 <- run_enrichr_for_list(de.results %>% filter(avg_log2FC > 0 & pval < 0.05) %>% pull(feature))
enrichr.outs.2 <- run_enrichr_for_list(de.results %>% filter(avg_log2FC < 0 & pval < 0.05) %>% pull(feature))

print("Up in MUT:")
enrichr.outs.1 %>% select(GeneRatio, BgRatio, pvalue, p.adjust, Count) %>% filter(pvalue < 0.05) %>% print
print("Up in WT:")
enrichr.outs.2 %>% select(GeneRatio, BgRatio, pvalue, p.adjust, Count) %>% filter(pvalue < 0.05) %>% print


############################################### HSC - TF motif volcano plots ##########################################

## This is the file currently used in the manuscript: HSC_chromvar_analysis_with_lmm_20231205.csv
da.motifs <- read_csv("data/chromVAR/JASPAR2020/HSC_chromvar_analysis_with_lmm_20231205.csv")

## This is a line for "experimental" files
da.motifs.new <- read_csv("data/chromVAR/JASPAR2020/HSC_chromvar_analysis_with_lmm_chromvar_celltype_20240411.csv")

plot_volcano(da.motifs, effect_column = "delta", stat_column = "fdr", stat_line = 0.05, effect_line = 0.1)

## Subtitle is based on filtering for metadata column below
p.motifs <- plot_volcano(da.motifs, metadata.df = atac.obj@meta.data %>% dplyr::filter(cluster_celltype == "HSC") %>% dplyr::filter(orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC")),
             effect_column = "delta", stat_column = "fdr", effect_line = 0.1, stat_line = 0.05, 
             title = "",
             genotyping_column = "genotype_pred_manual") +
             # only_genes = volcano.TF.list, only_genes_colors = volcano.TF.list.colors) + 
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 6), xlim = c(-0.5, 0.5)) +
  ggtitle("HSC differential TF motifs")
p.motifs 
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_tf_motifs_volcano_20240126.pdf", plot = p.motifs, device = "pdf", width = 7, height = 7)

p.motifs <- plot_volcano(da.motifs, metadata.df = atac.obj@meta.data %>% dplyr::filter(cluster_celltype == "HSC") %>% dplyr::filter(orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC")),
                         effect_column = "delta", stat_column = "fdr", effect_line = 0.1, stat_line = 0.05, 
                         title = "",
                         genotyping_column = "genotype_pred_manual", label_genes = F) +
  # only_genes = volcano.TF.list, only_genes_colors = volcano.TF.list.colors) + 
  theme(legend.position = "bottom") +
  coord_cartesian(ylim=c(0, 13), xlim = c(-0.5, 0.5)) +
  ggtitle("HSC differential TF motifs")
p.motifs 
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_tf_motifs_volcano_no_labels.pdf", plot = p.motifs, device = "pdf", width = 7, height = 7)


############################################## Other progenitors - TF motif volcano plots ##############################

## Subtitle is based on filtering for metadata column below
plist <-  lapply(c("LMPP", "MkP", "EMP"), function(x) {
  da.motifs <- read_csv(paste0("data/chromVAR/JASPAR2020/", x, "_chromvar_analysis_with_lmm.csv"))
  p.motifs <- plot_volcano(da.motifs, metadata.df = atac.obj@meta.data %>% dplyr::filter(orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC") & cluster_celltype == x),
                           effect_column = "delta", stat_column = "fdr", effect_line = 0.2, stat_line = 0.2, 
                           title = "",
                           genotyping_column = "genotype_pred_manual") +
    # only_genes = volcano.TF.list, only_genes_colors = volcano.TF.list.colors) + 
    theme(legend.position = "bottom") +
    # coord_cartesian(ylim=c(0, 6), xlim = c(-0.7, 0.7)) +
    ggtitle(x, subtitle = "Differential TF motifs")
  return(p.motifs)
})
wrap_plots(plist, ncol = 1)
ggsave("figures/current_figure_drafts/ATAC_HSC_diff_tf_motifs_volcano.pdf", plot = p.motifs, device = "pdf", width = 7, height = 7)


############################################### HSC - TF motif heatmap ##########################################

atac.obj.subset <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_HSC_chromvar_20231205.rds")
atac.obj.subset.donor <- subset(atac.obj.subset, CellType == "HSC" & orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC"))
# atac.obj.subset.donor$Genotype <- atac.obj.subset.donor$genotype_pred_manual
# atac.obj.subset.donor$Genotype[is.na(atac.obj.subset.donor$Genotype)] <- "NA"
# atac.obj.subset.donor$Genotype <- factor(atac.obj.subset.donor$Genotype, levels = c("WT", "MUT", "NA"))

## Formatting the sample and donor columns
meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_2024-03-29.csv")
atac.obj.subset.donor$Sample <- plyr::mapvalues(atac.obj.subset.donor$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$updated_name)
atac.obj.subset.donor$Sample <- factor(atac.obj.subset.donor$Sample, levels = meta.data.mapping$updated_name %>% sort())
atac.obj.subset.donor@meta.data %>% count(orig.ident, Sample)
atac.obj.subset.donor$Donor <- gsub("_BMMC|_CD34.*", "", atac.obj.subset.donor$Sample) %>% gsub("_[A|B]", "", .)
atac.obj.subset.donor$Donor <- factor(atac.obj.subset.donor$Donor, levels = unique(atac.obj.subset.donor$Donor) %>% sort())
atac.obj.subset.donor@meta.data %>% count(orig.ident, Sample, Donor)

da.motifs <- read_csv("data/chromVAR/JASPAR2020/HSC_chromvar_analysis_with_lmm_20231205.csv")
tf.list <- c(da.motifs %>% dplyr::filter(delta > 0) %>% slice_min(pval, n = 10) %>% pull(feature), da.motifs %>% dplyr::filter(delta < 0) %>% slice_min(pval, n = 10) %>% pull(feature))

tf.data.frames <- lapply(unique(atac.obj.subset.donor$Sample), function(x) {
  bcs.mut <- atac.obj.subset.donor@meta.data %>% dplyr::filter(Sample == x) %>% dplyr::filter(Genotype == "MUT") %>% rownames()
  bcs.wt <- atac.obj.subset.donor@meta.data %>% dplyr::filter(Sample == x) %>% dplyr::filter(Genotype == "WT") %>% rownames()
  difference.matrix <- atac.obj.subset.donor@assays$chromvar_HSC@data[tf.list, bcs.mut] %>% rowMeans() - atac.obj.subset.donor@assays$chromvar_HSC@data[tf.list, bcs.wt] %>% rowMeans()
  
  df <- data.frame(mean_chromvar_score = difference.matrix, sample = x)
  df <- df %>% rownames_to_column("TF") 
  return(df)
})
tf.motifs.per.sample <- do.call(rbind, tf.data.frames)
tf.matrix <- tf.motifs.per.sample %>% pivot_wider(names_from = sample, values_from = mean_chromvar_score) %>% column_to_rownames("TF")

dev.off()
pdf(file = "figures/current_figure_drafts/TF_motif_heatmap_20240328.pdf", width = 4, height = 4)
pheatmap::pheatmap(tf.matrix, color=colorRampPalette(c("navy", "white", "red"))(50),  
                   cluster_cols = F,
                   breaks = seq(-1, 1, length.out = 50), main = "Delta z-score for\nmotif accessibility")
dev.off()

## Plot number of cells genotyped for each donor
atac.obj.subset.donor@meta.data%>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(CellType == "HSC") %>% 
  count(Donor) %>% ggplot(aes(x = Donor, y = log10(n))) +
  geom_col(color = "black", fill = "grey40") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Cells Genotyped\n(log10)") +
  xlab("Patient sample")
ggsave(paste0("figures/current_figure_drafts/HSC_TF_motif_donor_genotyping_bar_plot.pdf"), dpi = 300, device = "pdf", height = 1.5, width = 2.25)

######################################### TF motif chromvar score boxplots ###################################

atac.obj.subset <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_HSC_chromvar_20231205.rds")
atac.obj.subset.donor <- subset(atac.obj.subset, CellType %in% c("HSC") & orig.ident %in% c("SG_VX16_ATAC", "SG_VX17_ATAC"))

## rename stat1::stat2
rownames(atac.obj.subset.donor@assays$chromvar_celltype@data) <- gsub("::", "_", rownames(atac.obj.subset.donor@assays$chromvar_celltype@data))

tf.data.frame <- atac.obj.subset.donor@assays$chromvar_HSC@data %>% as.matrix() %>% t()
tf.data.frame <- merge(tf.data.frame, atac.obj.subset.donor@meta.data, by = "row.names")

plot_TF_motif_boxplot <- function(x) {
  min.quantile <- tf.data.frame %>%
    filter(Genotype %in% c("MUT", "WT")) %>%
    dplyr::rename(module_plotting = x) %>%
    group_by(Genotype) %>% 
    summarize(quantile = quantile(module_plotting, c(0.05))) %>% 
    slice_min(quantile) %>% 
    pull(quantile)
  max.quantile <- tf.data.frame %>%
    filter(Genotype %in% c("MUT", "WT")) %>%
    dplyr::rename(module_plotting = x) %>%
    group_by(Genotype) %>% 
    summarize(quantile = quantile(module_plotting, c(0.95))) %>% 
    slice_max(quantile) %>% 
    pull(quantile)
  
  p1 <- tf.data.frame %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    mutate(genotype_category = paste0(Sample, "_", Genotype)) %>% 
    ggplot(aes_string(x = "Genotype", y = x, fill = "Genotype")) +
    geom_boxplot(outlier.shape = NA) +
    theme(legend.position = "none") +
    ggtitle(x) +
    # scale_fill_manual(values = c("PT16_WT" = "#B7B6E3", "PT16_MUT" = "#66539E", "PT17_WT" = "#A6C9E0", "PT17_MUT" = "#3A71B1")) +
    scale_fill_manual(values = vexas.genotyping.palette) +
    # scale_y_continuous(expand = c(0.2, 0)) +
    # coord_cartesian(ylim=c(-5, 5)) +
    scale_y_continuous(expand = c(0.2, 0), limits = c(min.quantile, max.quantile)) +
    stat_compare_means(label = "p.format", label.y = Inf, vjust = 1.5) +
    ylab(paste0("Motif accessibility\n(z-score)")) +
    xlab("")
  return(p1)
}

## Plot all donors combined for a few tags of interest
plist <- lapply(c("SPIB", "ATF4", "CEBPD", "GATA1_TAL1"), plot_TF_motif_boxplot)
p.combined <- wrap_plots(plist, nrow = 1)
p.combined
ggsave("figures/current_figure_drafts/ATAC_TF_motif_boxplot.pdf", plot = p.combined, device = "pdf", width = 8, height = 2.5, dpi = 300)

## Plot BM16, BM17 for CEBPA, ATF4
p.atf4 <- plot_TF_motif_boxplot("ATF4") + facet_grid(cols = vars(Sample))
p.atf4
ggsave("figures/current_figure_drafts/ATAC_TF_motif_boxplot_BM17_BM18_highlight_formatted.pdf", plot = p.atf4, device = "pdf", width = 2, height = 3, dpi = 300)


########################################## ATAC- ADT ##############################################################

## Make a scatter plot separating NK cells, T cells via ADT expression
atac.obj.adt.NK <- subset(atac.obj, CellType %in% c("CD4 T cells", "CD8 T cells", "NK"))
table(atac.obj.adt.NK$cluster_celltype)
geno.df <- atac.obj.adt.NK@meta.data %>% dplyr::select(genotype_pred_manual, Sample)
m <- atac.obj.adt.NK@assays$ADT_DSB@data %>% t() %>% as.data.frame()
m.genotyped <- merge(m, geno.df, by = "row.names", all = F)
m.genotyped %>% 
  dplyr::filter(genotype_pred_manual %in% c("MUT", "WT")) %>%
  ggplot(aes(x = `Hu.CD3-UCHT1`, y = `Hu.CD7`, color = genotype_pred_manual)) +
  geom_point()
  # coord_fixed() +
  # facet_grid(cols = vars(Sample), scales = "free") +
  # scale_color_manual(values = vexas.genotyping.palette) 

##################################### Accessibility quantile vs. genotyping plot ##############################

## UBA1 mutation site peak: chrX-47198691-47199456
atac.obj@assays$ATAC@data[1:5, 1:5]
rownames(atac.obj@assays$ATAC@data) %>% grep("chrX-4719", ., value = T)

## Stratify cells based on accessilibity at that peak
mutation.site.reads <- data.frame(UBA1_mutation_site_reads = atac.obj@assays$ATAC@counts["chrX-47198691-47199456", ])
md.mut.site <- merge(md, mutation.site.reads, by = "row.names")
table(md.mut.site$UBA1_mutation_site_reads)
md.mut.site %>% 
  ggplot(aes(x = UBA1_mutation_site_reads, y = log10(..count..))) +
  geom_histogram(binwidth = 1) + 
  xlab("Read count") +
  ylab("Log10(cell count)") + 
  ggtitle("Number of reads covering\nUBA1 mutation site", subtitle = "chrX-47198691-47199456")
ggsave("figures/current_figure_drafts/ATAC_UBA1_site_coverage.pdf", device = "pdf", width = 2.5, height = 3, dpi = 300)

##################################### Saving HSC reads for mutation calling ###########################

lapply(c("SG_VX17_ATAC", "SG_VX18_ATAC", "VEX_BM10_POS_ATAC"), function(x) {
  atac.obj@meta.data %>% 
    filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP") & orig.ident == x) %>% 
    filter(Genotype == "MUT") %>% 
    rownames_to_column("barcode") %>% 
    mutate(barcode = gsub(paste0(x, "_"), "", barcode)) %>% 
    select(barcode) %>% 
    write_csv(file = paste0("/gpfs/commons/groups/landau_lab/VEXAS/checking_mutation_rates/barcode_lists/", x, "_MUT.txt"), col_names = F)
  
  atac.obj@meta.data %>% 
    filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP") & orig.ident == x) %>% 
    filter(Genotype == "WT") %>% 
    rownames_to_column("barcode") %>% 
    mutate(barcode = gsub(paste0(x, "_"), "", barcode)) %>% 
    select(barcode) %>% 
    write_csv(file = paste0("/gpfs/commons/groups/landau_lab/VEXAS/checking_mutation_rates/barcode_lists/", x, "_WT.txt"), col_names = F)
})

##################################### Comparing old vs. new genotyping protocol ########################

## Print averages
md %>% 
  group_by(Sample) %>% 
  summarize(genotyping_rate = sum(Genotype %in% c("MUT", "WT")) / n()) %>% 
  mutate(Protocol = if_else(Sample %in% c("PT16", "PT17", "PT18"), "Updated", "Original")) %>% 
  group_by(Protocol) %>% 
  summarize(mean = mean(genotyping_rate))

## Bar plot
md %>% 
  group_by(Sample) %>% 
  summarize(genotyping_rate = sum(Genotype %in% c("MUT", "WT")) / n()) %>% 
  mutate(Protocol = if_else(Sample %in% c("PT16", "PT17", "PT18"), "Updated", "Original")) %>% 
  group_by(Protocol) %>% 
  mutate(average_val = mean(genotyping_rate)) %>% 
  ggplot(aes(x = Sample, y = genotyping_rate, fill = Protocol, label = average_val)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c("grey80", "goldenrod2")) +
  theme(legend.position = "none") +
  xlab(element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Genotyping rate") +
  ggtitle("GoTChA Genotyping Protocol", subtitle = "Mean genotyping rates, original protocol: 8.33%,\nupdated protocol: 17.0%") +
  geom_hline(color = "black", aes(yintercept = average_val), linetype = "longdash") +
  # annotate("text", x = 0, y = 0.25) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  facet_grid(cols = vars(Protocol), drop = TRUE, scales = "free_x", space = "free")
  # geom_hline(color = "goldenrod2", yintercept = 0.170)

ggsave("figures/current_figure_drafts/GoTChA_genotyping_protocols_20240328.pdf", device = "pdf", width = 4.5, height = 3.5)

###################################### Phospho-seq plots ################################################

## Normalize within cells
sample.plot = "PT16"
total.tsb.counts <- atac.obj@assays$ADT_TSB@counts %>% colSums()
atac.obj <- AddMetaData(atac.obj, total.tsb.counts, col.name = "total_tsb_counts")
bcs.keep <- atac.obj@meta.data %>% 
  filter(!is.na(total_tsb_counts)) %>% 
  filter(Sample == sample.plot) %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>% rownames()
atac.obj.sub <- subset(atac.obj, cells = bcs.keep)

atac.obj.sub <- NormalizeData(atac.obj.sub, normalization.method = "CLR", margin = 2, assay = "ADT_TSB")
atac.obj.sub@assays$ADT_TSB@counts[1:5, 1:5]
atac.obj.sub@assays$ADT_TSB@data[1:5, 1:5]
atac.obj.sub <- AddMetaData(atac.obj.sub, atac.obj.sub@assays$ADT_TSB@data %>% t())

table(atac.obj.sub$orig.ident)

## Make a rank plot
tag.list <- rownames(atac.obj.sub@assays$ADT_TSB@data)
p.rank <- atac.obj.sub@meta.data %>% 
  filter(!is.na(IgG1)) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  filter(cluster_celltype == "HSC") %>%
  group_by(Genotype) %>% 
  summarize(across(tag.list, mean)) %>% 
  pivot_longer(cols = tag.list, names_to = "TSB", values_to = "expression") %>% 
  pivot_wider(names_from = "Genotype", values_from = "expression") %>% 
  mutate(diff = MUT - WT, abs_diff = abs(MUT - WT)) %>% 
  arrange(-diff) %>% 
  mutate(rank = row_number()) %>%
  mutate(direction = ifelse(diff > 0, "Up in MUT", "Up in WT")) %>% 
  mutate(direction = ifelse(rank > 10 & rank < 40, "NA", direction)) %>% 
  mutate(label = ifelse(direction != "NA" | TSB == "CHOP" | TSB == "ATF4", TSB, NA_character_)) %>%
  ggplot(aes(x = reorder(TSB, rank), y = diff, label = label, color = direction)) +
  geom_point(size = 2) +
  scale_color_manual(values = c(`Up in MUT` = vexas.genotyping.palette[["MUT"]], `Up in WT` = vexas.genotyping.palette[["WT"]])) +
  # theme(axis.text.x = element_blank()) +
  # scale_x_discrete(breaks = c(1, 10, 20, 30, 40, 50)) +
  xlab("Rank") + ylab("Average(MUT) - Average(WT)") + 
  ggtitle(paste0(sample.plot, ", TSB expression"), subtitle = "HSC clusters") + 
  geom_hline(yintercept = 0, linetype = "longdash")
ggsave(paste0("figures/current_figure_drafts/ATAC_HSC_", sample.plot, "_phospho_seq_rank_plot_no_labels_20240409.pdf"), plot = p.rank, dpi = 300, device = "pdf", width = 7, height = 5)
p.rank <- p.rank + geom_text_repel(max.overlaps = 25, point.padding = 0.5, box.padding = 0.5)
p.rank
ggsave(paste0("figures/current_figure_drafts/ATAC_HSC_", sample.plot, "_phospho_seq_rank_plot_20240409.pdf"), plot = p.rank, dpi = 300, device = "pdf", width = 7, height = 5)
# ggsave(paste0("figures/current_figure_drafts/ATAC_HSC_EMP_LMPP_MkP_", sample.plot, "_phospho_seq_rank_plot_20240326.pdf"), plot = p.rank, dpi = 300, device = "pdf", width = 7, height = 5)

## Save the rank plot metrics to a table for the supplementary
atac.obj.sub@meta.data %>% 
  filter(!is.na(IgG1)) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  # filter(cluster_celltype == "HSC") %>%
  group_by(Genotype) %>% 
  summarize(across(tag.list, mean)) %>% 
  pivot_longer(cols = tag.list, names_to = "TSB", values_to = "expression") %>% 
  pivot_wider(names_from = "Genotype", values_from = "expression") %>% 
  mutate(diff = MUT - WT, abs_diff = abs(MUT - WT)) %>% 
  arrange(-diff) %>% 
  mutate(rank = row_number()) %>% 
  write_csv(file = paste0("data/phospho_seq/ATAC_HSC_", sample.plot, "_phospho_seq_rank_plot_metrics_20240409.csv"))

## Making boxplots of all tags
making_boxplot_for_adt <-  function(x) {
  plabel.position <- atac.obj.sub@meta.data %>% 
    filter(!is.na(IgG1)) %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    pull(x) %>% max()
    
  p.boxplot <- atac.obj.sub@meta.data %>% 
    filter(!is.na(IgG1)) %>% 
    # filter(!! sym(x) != 0) %>% 
    filter(Genotype %in% c("MUT", "WT")) %>% 
    ggplot(aes(x = Genotype, y = !! sym(x), fill = Genotype)) +
    geom_boxplot() +
    # coord_cartesian(ylim = c(min(quantiles.current.mut[[1]], quantiles.current.wt[[1]]), max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]))) +
    theme_classic() +
    ylab("ADT expression") +
    theme(axis.title.x = element_blank()) +
    # stat_compare_means(label = "p.format", label.y =max(quantiles.current.mut[[2]], quantiles.current.wt[[2]]) * 1.1) +
    stat_compare_means(label = "p.format", label.y = plabel.position*1.2) +
    scale_y_continuous(expand = c(0.3, 0)) +
    theme(legend.position = "none") +
    scale_fill_manual(values = vexas.genotyping.palette) +
    ggtitle(x)
  return(p.boxplot)
}
tag.list <- rownames(atac.obj.sub@assays$ADT_TSB@data)
plist <- lapply(rownames(atac.obj.sub@assays$ADT_TSB@data), making_boxplot_for_adt)
names(plist) <- rownames(atac.obj.sub@assays$ADT_TSB@data)

p.select.tags <- (plist$CHOP | plist$ATF4 | plist$CEBPA | plist$MCL1) + plot_annotation("BM17, TSB quantification", subtitle = "HSC, EMP, LMPP, MkP clusters")
p.select.tags
ggsave("figures/current_figure_drafts/ATAC_boxplots_select_tags.pdf", plot = p.select.tags, dpi = 300, device = "pdf", width = 7, height = 4)

## Plot all tags together
p.combined <- wrap_plots(plist, nrow = 5)
p.combined
ggsave("figures/current_figure_drafts/ATAC_HSC_EMP_LMPP_MkP_BM17_phospho_seq_boxplots.pdf", plot = p.combined, dpi = 300, device = "pdf", width = 15, height = 10)

## Making heatmaps by tag
target.tags <- c("CEBPA", "ATF4", "CHOP", "MCL1", "SPI-B")
m.phospho <- atac.obj.sub@meta.data %>% 
  select(Donor, CellType, Genotype, all_of(tag.list)) %>% 
  rownames_to_column("Barcode") %>% 
  pivot_longer(cols = all_of(tag.list), names_to = "tag", values_to = "value") %>% 
  group_by(Donor, CellType, tag, Genotype) %>% 
  summarize(average_expression = mean(value)) %>% 
  pivot_wider(names_from = Genotype, values_from = average_expression) %>% 
  mutate(difference_average_expression = MUT - WT) %>%
  ungroup() %>% 
  select(CellType, tag, difference_average_expression) %>%
  filter(tag %in% target.tags) %>% 
  pivot_wider(names_from = tag, values_from = difference_average_expression) %>% 
  column_to_rownames("CellType") %>% 
  # column_to_rownames("Donor") %>% 
  t()

atac.obj.sub@meta.data %>% 
  filter(CellType == "HSC" & Genotype %in% c("MUT", "WT")) %>% 
  ggplot(aes(x = Genotype, y = ATF4)) +
  geom_boxplot() +
  geom_point() + stat_compare_means()

atac.obj.sub@meta.data %>% 
  filter(CellType == "HSC" & Genotype %in% c("MUT", "WT")) %>% 
  group_by(Genotype) %>% 
  summarize(mean = mean(ATF4))

pdf(paste0("figures/current_figure_drafts/per_patient_phosphoseq_tags_", sample.plot, "_no_scaling.pdf"), width = 4, height = 3)
pheatmap(m.phospho, cellheight=10, cellwidth = 10, 
         scale = "none", 
         cluster_rows = F, cluster_cols = F, 
         main = sample.plot,
         color=colorRampPalette(c("navy", "white", "red"))(50),  breaks = seq(-0.1, 0.1, length.out = 50))
dev.off()

m.phospho.all <- atac.obj.sub@meta.data %>% 
  select(Donor, CellType, Genotype, all_of(tag.list)) %>% 
  rownames_to_column("Barcode") %>% 
  pivot_longer(cols = all_of(tag.list), names_to = "tag", values_to = "value") %>% 
  group_by(Donor, CellType, tag, Genotype) %>% 
  summarize(average_expression = mean(value)) %>% 
  pivot_wider(names_from = Genotype, values_from = average_expression) %>% 
  mutate(difference_average_expression = MUT - WT) %>%
  ungroup() %>% 
  select(CellType, tag, difference_average_expression) %>%
  pivot_wider(names_from = tag, values_from = difference_average_expression) %>% 
  column_to_rownames("CellType") %>% 
  # column_to_rownames("Donor") %>% 
  t()

pdf(paste0("figures/current_figure_drafts/per_patient_phosphoseq_tags_", sample.plot, "_no_scaling_all_tags.pdf"), width = 3, height = 10)
pheatmap(m.phospho.all, cellheight=10, cellwidth = 10, 
         scale = "none", 
         cluster_rows = T, cluster_cols = F, 
         main = sample.plot,
         color=colorRampPalette(c("navy", "white", "red"))(50),  breaks = seq(-0.2, 0.2, length.out = 50))
dev.off()

atac.obj.sub@meta.data %>% 
  select(Donor, CellType, Genotype, all_of(tag.list)) %>% 
  rownames_to_column("Barcode") %>% 
  mutate(ATF4_CHOP_ratio = ATF4 / CHOP) %>% 
  ggplot(aes(x = Genotype, y = ATF4_CHOP_ratio)) +
  geom_boxplot() +
  facet_grid(cols = vars(CellType))
  
## Correlate chromVAR scores with protein levels
atac.obj.hsc.chromvar <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_HSC_chromvar_20231205.rds")
atf4.chromvar <- atac.obj.hsc.chromvar@assays$chromvar_HSC$data["ATF4", ]
atac.obj.sub <- AddMetaData(atac.obj.sub, atf4.chromvar, col.name = "ATF4_chromvar")

## Chromvar run with HSCs only
atac.obj.sub@meta.data %>% 
  filter(!is.na(IgG1)) %>% 
  filter(cluster_celltype == "HSC") %>%
  ggplot(aes(x = ATF4_chromvar, y = ATF4, color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  geom_point() +
  facet_grid(cols = vars(Genotype)) +
  ggtitle("HSCs", subtitle = "chromvar run with HSCs only") +
  stat_cor(label.y = 1.5)

## Chromvar run with all cells, HSPCs
atf4.chromvar <- atac.obj.sub@assays$chromvar_celltype$data["ATF4", ]
atac.obj.sub$ATF4_chromvar <- NULL
atac.obj.sub <- AddMetaData(atac.obj.sub, atf4.chromvar, col.name = "ATF4_chromvar")
atac.obj.sub@meta.data %>% 
  filter(!is.na(IgG1)) %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP")) %>%
  ggplot(aes(x = ATF4_chromvar, y = ATF4, color = Genotype)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  geom_point() +
  facet_grid(cols = vars(Genotype)) +
  ggtitle("HSPCs", subtitle = "chromvar run with all cells") +
  stat_cor(label.y = 1.5)

## Chromvar run with all cells, all cells
bcs.keep <- atac.obj@meta.data %>% 
  filter(!is.na(total_tsb_counts)) %>% 
  filter(Sample == sample.plot) %>% rownames()
atac.obj.sub <- subset(atac.obj, cells = bcs.keep)
atac.obj.sub <- NormalizeData(atac.obj.sub, normalization.method = "CLR", margin = 2, assay = "ADT_TSB")
atac.obj.sub <- AddMetaData(atac.obj.sub, atac.obj.sub@assays$ADT_TSB@data %>% t())
atf4.chromvar <- atac.obj.sub@assays$chromvar_celltype$data["ATF4", ]
atac.obj.sub$ATF4_chromvar <- NULL
atac.obj.sub <- AddMetaData(atac.obj.sub, atf4.chromvar, col.name = "ATF4_chromvar")
atac.obj.sub@meta.data %>% 
  filter(!is.na(IgG1)) %>% 
  ggplot(aes(x = ATF4_chromvar, y = ATF4)) +
  scale_color_manual(values = vexas.genotyping.palette) +
  geom_point(alpha = 0.1) +
  # facet_grid(cols = vars(Genotype)) +
  ggtitle("All cells", subtitle = "chromvar run with all cells") +
  stat_cor(label.y = 1.8)
  
  
atac.obj@meta.data$Donor %>% table()


atac.obj.sub <- subset(atac.obj, Donor %in% c("PT16","PT17", "PT18"))
Idents(atac.obj) <- "Donor"
CoveragePlot(atac.obj, region = "chrY-2786855-2787682", assay = "ATAC")

 ############# Checking cell type assignments with ADT ################

to.plot.tags <- c("Isotype-RTK2071", "Isotype-RTK2758", "Isotype-RTK4174",  "Isotype-G0114F7", "Isotype-HTK888", 
                  "Hu.CD19", "Hu.CD3-UCHT1", "Hu.CD4-RPA.T4", "Hu.CD8", "Hu.CD14-M5E2",  "Hu.CD16", "Hu.CD56",
                  "Hu.CD7", "Hu.CD41", "Hu.CD71", "Hu.HLA.DR", "Hu.CD36", "Hu.CD99", "Hu.CD90", "Hu.CD49b", "HuMs.CD49f", "Hu.CD150", "Hu.CD38-HIT2", "Hu.CD45RA", "Hu.CD123")

## CLR normalization
adt.counts.per.bc <- colSums(atac.obj@assays$ADT_CLR@counts)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
atac.obj.adt <- subset(atac.obj, cells = bcs.keep)
DotPlot(atac.obj.adt, assay = "ADT_CLR", features = to.plot.tags, group.by = "cluster_celltype", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
atac.obj.adt.hpscs <- subset(atac.obj.adt, cluster_celltype %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk"))
DotPlot(atac.obj.adt.hpscs, assay = "ADT_CLR", features = to.plot.tags, group.by = "cluster_celltype", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
atac.obj.adt.hpscs <- subset(atac.obj.adt, cluster_celltype %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk"))
DotPlot(atac.obj.adt.hpscs, assay = "ADT_CLR", features = to.plot.tags, group.by = "predicted.l2", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## DSB normalization
adt.counts.per.bc <- colSums(atac.obj@assays$ADT_DSB@data)
bcs.keep <- adt.counts.per.bc[!is.na(adt.counts.per.bc)] %>% names()
atac.obj.adt <- subset(atac.obj, cells = bcs.keep)
atac.obj.adt <- subset(atac.obj.adt, cluster_celltype != "CD14 Mono")
all.tags <- atac.obj.adt@assays$ADT_DSB %>% rownames()
DotPlot(atac.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[1:75], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(atac.obj.adt, assay = "ADT_DSB", group.by = "cluster_celltype", features = all.tags[75:164], col.min = -1, col.max = 2) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
DotPlot(atac.obj.adt, assay = "ADT_DSB", features = to.plot.tags, group.by = "cluster_celltype", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
atac.obj.adt <- subset(atac.obj.adt, predicted.l2.score > 0.5)
DotPlot(atac.obj.adt, assay = "ADT_DSB", features = to.plot.tags, group.by = "predicted.l2") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

atac.obj.adt.hpscs <- subset(atac.obj.adt, cluster_celltype %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk"))
DotPlot(atac.obj.adt.hpscs, assay = "ADT_DSB", features = to.plot.tags, group.by = "cluster_celltype", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
atac.obj.adt.hpscs <- subset(atac.obj.adt, predicted.l2 %in% c("HSC", "EMP", "LMPP", "GMP", "Prog Mk") & predicted.l2.score > 0.5)
DotPlot(atac.obj.adt.hpscs, assay = "ADT_DSB", features = to.plot.tags, group.by = "predicted.l2", ) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Idents(atac.obj.adt) <- "predicted.l2"
VlnPlot(atac.obj.adt, assay = "ADT_DSB", features = "Hu.CD14-M5E2", group.by = "predicted.l2", idents = c("HSC", "EMP", "LMPP", "GMP", "CD14 Mono", "CD8 Effector_2", "Memory B")) + NoLegend()

###################### Combine plots in grid format and save #######################

### Using patchwork
# p.combined <- ((p.celltypes + p.genotyping) / plot_spacer()) | (p.genotype.rates / p.uba1.expression) + plot_layout(widths = c(5, 1), ncol = 2, nrow = 2)
pane.1 <- (p.celltypes | p.genotyping | p.mut.freq.barplots) # + plot_layout(widths = c(2, 2, 1))
pane.2 <- ((p.genotype.rates / p.uba1.expression) | p.hsc | p.lmpp)
# p.combined <- ((pane.1 / pane.2) | pane.3) + plot_layout(widths = c(3, 1))
p.combined <- pane.1 / pane.2

ggsave("figures/vexas_genotyping_umaps_2_02_23.pdf", plot = p.combined, width = 40, height = 30, units = "cm", dpi = "retina")

### Using gridExtra
# grobs.list <- list(`1` = p.celltypes, `2` = p.cellcycle.module, `3` = p.hspc.module, `4` = p.genotyping)
# 
# lay <- rbind(c(1,1,2,2, 4, 4),
#              c(1,1,3,3, 4, 4))
# p.combined <- grid.arrange(grobs = grobs.list, layout_matrix = lay)

# lay <- rbind(c(0,0,0,0,1,1,1),
#              c(0,0,0,0,1,1,1),
#              c(4,4,2,2,5,6,7),
#              c(4,4,3,3,8,9,10))

# p.combined <- grid.arrange(grobs = grobs.list, layout_matrix = lay)

# 
# design <- "
#   11122244
#   11122244
#   33333355
# "
# p.combined <- p.celltypes + p.genotyping + plot_spacer() + p.genotype.rates + p.uba1.expression + plot_layout(design = design)
# 
# p.combined <- p.celltypes + p.genotyping + p.genotype.rates + plot_spacer() + plot_spacer() + p.uba1.expression + plot_layout(ncol = 3, nrow = 2, widths = c(1.5, 1.5, 1), heights = c(2, 1))

#################### Looking at TF motifs across progenitors #####################

atac.obj@meta.data %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "CLP")) %>% 
  mutate(CellType = factor(CellType, levels = c("HSC", "EMP", "LMPP", "MkP", "CLP"))) %>% 
  filter(Genotype %in% c("MUT", "WT")) %>% 
  count(CellType, Genotype) %>% 
  pivot_wider(names_from = Genotype, values_from = n)

atac.obj@assays$chromvar_celltype %>% rownames() %>% grep("ATF", ., value = T)

tf.data.frame <- atac.obj@assays$chromvar_celltype@data %>% as.matrix() %>% t()
tf.data.frame <- merge(tf.data.frame, atac.obj@meta.data, by = "row.names")

tf.data.frame %>% 
  filter(cluster_celltype %in% c("HSC", "EMP", "LMPP", "MkP", "CLP")) %>% 
  mutate(CellType = factor(CellType, levels = c("HSC", "EMP", "LMPP", "MkP", "CLP"))) %>% 
  filter(Genotype == "MUT") %>% 
  ggplot(aes(x = CellType, y = ATF4)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("HSC", "CLP"))) + ggtitle("MUT cells only")


###################### Updating sample labels ########################

df.2 <- read_csv("data/formatting_patient_ids/patient_id_mappings_2024-03-28.csv")

orig.ident.to.label.mapping <- df.2$`Study ID`
names(orig.ident.to.label.mapping) <- df.2$`Patient ID`
head(orig.ident.to.label.mapping)

## Do an initial replace of the sample IDs
df.1 <- data.frame(orig.ident = atac.obj@meta.data$orig.ident %>% unique %>% sort())
df.1$updated_name <- df.1$orig.ident
for (sample in names(orig.ident.to.label.mapping)) {
  df.1 <- df.1 %>% mutate(updated_name = gsub("VX", "BM", updated_name) %>% str_replace_all(., paste0(sample, "[a|b]*_"), paste0(orig.ident.to.label.mapping[[sample]], "_")))
}

## Save the list
df.1 %>% 
  mutate(updated_name = gsub("SG_|VEX_", "", updated_name)) %>% 
  mutate(updated_name = gsub("_[A-Za-z].*", "", updated_name)) %>% 
  write_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_2024-03-28.csv")

