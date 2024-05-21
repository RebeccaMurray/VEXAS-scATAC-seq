library(Seurat)
library(tidyverse)
library(patchwork)
library(ggrepel)
library(fgsea)
library(msigdbr)

source("~/rscripts/general_helpers/volcano_plots.R")
source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

atac.obj <-  readRDS("data/seurat_objects/all_samples_50_svd.RDS")
# atac.obj.subset <- subset(atac.obj, seurat_clusters == 1)

## Pull metadata 
meta.data.df <- atac.obj@meta.data

de.results <- list()

## DE results - HSCs
de.results[["HSC"]] <- read_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/10p_genotyped_cells_sample_filtered/HSC_DiffLMM2_expressed_only.csv") %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = "HSC")
p.hsc <- de.results[["HSC"]] %>% plot_volcano(., metadata.df = NULL, stat_line = 0.05, effect_line = 0.2, stat_column = "pval", title = paste0("HSCs"), 
                                           genotyping_column = NA, subtitle_addition = "Linear mixed model") + theme(legend.position = "none")
p.hsc + coord_cartesian(xlim = c(-2, 2))
p.hsc

### running GSEA
genes.tested <- de.results[["HSC"]]$feature
fgsea.out <- run_fgsea_for_list(de.list = de.results[["HSC"]], genes.tested = genes.tested)

p1 <- plot_enrichment(de.results[["HSC"]], pathway.name = "HALLMARK_INFLAMMATORY_RESPONSE")

p2 <- plot_enrichment(de.results[["HSC"]], pathway.name = "HALLMARK_P53_PATHWAY")
p3 <- plot_enrichment(de.results[["HSC"]], pathway.name = "KEGG_CHEMOKINE_SIGNALING_PATHWAY")

p1 / p3

fgsea.out %>% write_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/10p_genotyped_cells_sample_filtered/gsea_results.csv")

## running GSEA - manual ranking
de_results.ranked <- de.results[["HSC"]] %>%
  filter(!is.na(pval)) %>% 
  filter(is.finite(-log10(pval))) %>% 
  mutate(rank = sign(avg_log2FC) * -log10(pval)) %>%
  arrange(-rank)

ranks.list <- de_results.ranked$rank
names(ranks.list) <- de_results.ranked$feature
head(ranks.list)
tail(ranks.list)

genes.tested <- atac.obj.subset@assays$RNA_ArchR@data %>% rownames()
fgsea.out <- run_fgsea_for_list(ranks.list = ranks.list, genes.tested = genes.tested)
plot_enrichment(ranks.list = ranks.list, pathway.name = "HALLMARK_P53_PATHWAY")
plot_enrichment(ranks.list = ranks.list, pathway.name = "HALLMARK_INFLAMMATORY_RESPONSE")


################################ Save some GSEA results #######################

p53.genes <- fgsea.out %>% filter(pathway == "HALLMARK_P53_PATHWAY") %>% pull(leadingEdge) %>% unlist()
data.frame(p53_leading_edge_genes = p53.genes) %>% 
  write_csv(., file = "data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/gsea/HSCs_HALLMARK_P53_PATHWAY_leading_edge.csv")






## DE results - HSPCs
de.results[["HSPC"]] <- read_csv("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/celltype_labels_hspcs_donor_cutoff_all_genes_expressed_5p_hspcs_only.csv") %>% mutate(avg_log2FC = log2(fc)) %>% mutate(cluster = "HSC")
p.hsc <- de.results[["HSPC"]] %>% plot_volcano(., metadata.df = NULL, stat_cutoff = 0.01, effect_cutoff = 0.5, effect_line = 0.01, stat_column = "pval", title = paste0("HSCs"), 
                                              genotyping_column = NA, subtitle_addition = "Linear mixed model") + theme(legend.position = "none")
p.hsc + coord_cartesian(xlim = c(-10, 10))
p.hsc









##################################### Writing to a temp file for use with GoST ###########################

de.results[["HSPC"]] %>% 
  filter(avg_log2FC > 0 & pval < 0.05) %>% 
  select(feature) %>% 
  write_csv("tmp.txt")

de.results.combined <- do.call(rbind, de.results)
de.results.combined %>% 
  filter(avg_log2FC > 0 & pval < 0.05) %>% 
  select(feature, cluster) %>% 
  count(feature) %>% 
  arrange(-n) %>% 
  filter(n > 3) %>% 
  View()
  select(feature) %>% 
  write_csv("tmp.txt")

de.results.combined %>% 
  filter(avg_log2FC < 0 & pval < 0.05) %>% 
  select(feature, cluster) %>% 
  count(feature) %>% 
  arrange(-n) %>% 
  filter(n > 2) %>% 
  select(feature) %>% 
  write_csv("tmp.txt")

gene.list <- c("MYC", "TXNIP", "GADD45B", "DDIT4","DDIT3", "UCHL1", "CALR", "NPM1", "CAT", "DLK1", "CRYGD", "CDK6")

m <- de.results.combined %>% 
  # mutate(negative_log10_fdr = -log10(fdr)) %>% 
  filter(pval < 0.05) %>% 
  filter(feature %in% gene.list) %>% 
  select(cluster, feature, avg_log2FC) %>% 
  pivot_wider(names_from = feature, values_from = avg_log2FC) %>% 
  column_to_rownames("cluster")

m <- t(m)
m[is.na(m)] = 0

# pheatmap::pheatmap(m[, c("HSPC", "LMPP", "EMP", "MkP")], cluster_cols = F, breaks = seq(-1, 1, length.out = 25), color = colorRampPalette(brewer.pal(8, "Blues"))(25)))

################################## Making combined dotplot ###########################

fgsea.res.all <- lapply(names(de.results), function(x) { 
  fgsea.res <- run_fgsea_for_list(de.results[[x]], genes.tested = genes.tested)
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
  }
)

fgsea.res.combined <- do.call(rbind, fgsea.res.all)

## Plot the fgsea dotplot - all results
fgsea.res.combined %>% 
  mutate(pathway.name = gsub("_", " ", pathway)) %>% 
  mutate(pathway.name = str_to_title(pathway.name)) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  # mutate(cluster = factor(cluster, levels = unique(cluster))) %>% 
  mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSPC", "LMPP", "EMP", "MkP", "GMP")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_color_gradient2(low = "blue", high = "red", mid = "grey80", na.value = "grey") +
  geom_point(stat = "identity") +
  # facet_grid(rows = vars(category), scales = "free") +
  theme(panel.spacing = unit(2, "lines")) +
  xlab("") + ylab("") + theme_classic()


## Load the fgsea pathways list from bottom of file
fgsea.pathways.df <- data.frame(category = names(fgsea.pathways.specific), pathway = unlist(fgsea.pathways.specific))

## Plot the fgsea dotplot
fgsea.res.combined %>% 
  filter(pathway %in% fgsea.pathways.specific) %>% 
  merge(., fgsea.pathways.df, by.x = "pathway", by.y = "pathway", all.x = T, all.y = T) %>% 
  mutate(pathway.name = gsub("_", " ", pathway)) %>% 
  mutate(pathway.name = str_to_title(pathway.name)) %>% 
  mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  # mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSPC", "LMPP", "EMP", "MkP", "GMP")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_color_gradient2(low = "blue", high = "red", mid = "grey80", na.value = "grey") +
  geom_point(stat = "identity") +
  facet_grid(rows = vars(category), scales = "free") +
  theme(panel.spacing = unit(2, "lines")) +
  xlab("") + ylab("") + theme_classic()


fgsea.pathways.specific = list(
  Inflammation = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  Inflammation = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  Inflammation = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  Inflammation = "HALLMARK_TGF_BETA_SIGNALING",
  `Metabolism` = "HALLMARK_GLYCOLYSIS",
  `Metabolism` = "KEGG_OXIDATIVE_PHOSPHORYLATION",
  `Cellular regulation` = "REACTOME_TRANSLATION",
  `Cellular regulation` ="HALLMARK_MYC_TARGETS_V1",
  `Cellular regulation` ="HALLMARK_MYC_TARGETS_V2",
  `Cellular regulation` = "HALLMARK_G2M_CHECKPOINT",
  `Cellular regulation` = "HALLMARK_APOPTOSIS",
  `Cellular regulation` = "HALLMARK_P53_PATHWAY",
  # `Cellular regulation` = "KEGG_P53_SIGNALING_PATHWAY",
  Proteostasis = "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
  Proteostasis = "KEGG_REGULATION_OF_AUTOPHAGY",
  Proteostasis = "KEGG_RNA_DEGRADATION",
  Proteostasis = "PERK module",
  Proteostasis = "ATF6 module",
  Proteostasis = "XBP1 target module"
)


# fgsea.pathways.specific = list(
#   Inflammation = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
#   Inflammation = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
#   Inflammation = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
#   Inflammation = "HALLMARK_TGF_BETA_SIGNALING",
#   `Cellular regulation` = "REACTOME_TRANSLATION",
#   `Cellular regulation` ="HALLMARK_MYC_TARGETS_V1",
#   `Cellular regulation` ="HALLMARK_MYC_TARGETS_V2",
#   `Cellular regulation` = "HALLMARK_G2M_CHECKPOINT",
#   `Cellular regulation` = "KEGG_OXIDATIVE_PHOSPHORYLATION",
#   `Cellular regulation` = "HALLMARK_APOPTOSIS",
#   `Cellular regulation` = "HALLMARK_P53_PATHWAY",
#   `Cellular regulation` = "KEGG_P53_SIGNALING_PATHWAY",
#   Proteostasis = "KEGG_UBIQUITIN_MEDIATED_PROTEOLYSIS",
#   Proteostasis = "KEGG_REGULATION_OF_AUTOPHAGY",
#   Proteostasis = "KEGG_RNA_DEGRADATION",
#   Proteostasis = "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS",
#   Proteostasis = "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE",
#   Proteostasis = "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE"
# )


################################## Making combined dotplot - hallmark only ###########################

fgsea.res.all.h <- lapply(names(de.results), function(x) { 
  fgsea.res <- run_fgsea_for_list_hallmark_only(de.results[[x]])
  fgsea.res[["cluster"]] <- x
  return(fgsea.res)
}
)

fgsea.res.combined.h <- do.call(rbind, fgsea.res.all.h)

## Plot the fgsea dotplot - all results
fgsea.res.combined.h %>% 
  mutate(pathway.name = gsub("_", " ", pathway)) %>% 
  mutate(pathway.name = str_to_title(pathway.name)) %>% 
  # mutate(pathway.name = str_wrap(pathway.name, width = 30)) %>% 
  # mutate(cluster = factor(cluster, levels = unique(cluster))) %>% 
  mutate(NES = if_else(pval < 0.05, NES, NA_real_)) %>%
  ggplot(aes(x = factor(cluster, levels = c("HSPC", "LMPP", "EMP", "MkP", "GMP")), y = pathway.name, size = -log10(pval), color = NES)) +
  scale_color_gradient2(low = "blue", high = "red", mid = "grey80", na.value = "grey") +
  geom_point(stat = "identity") +
  # facet_grid(rows = vars(category), scales = "free") +
  theme(panel.spacing = unit(2, "lines")) +
  xlab("") + ylab("") + theme_classic()







# hsc.clusters <- c(9, 11, 2)
# samples.vexas.and.control@meta.data %>% 
#   filter(seurat_clusters %in% hsc.clusters) %>% 
#   count(orig.ident, seurat_clusters)
# Idents(samples.vexas.and.control) <- "seurat_clusters"
# DimPlot(samples.vexas.and.control, group.by = "seurat_clusters", split.by = "seurat_clusters", ncol = 4)
# VlnPlot(samples.vexas.and.control, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), idents = hsc.clusters)
# DimPlot(samples.vexas.and.control, cells.highlight = samples.vexas.and.control@meta.data %>% filter(seurat_clusters %in% hsc.clusters) %>% rownames())
# DimPlot(subset.rna.obj, group.by = "seurat_clusters", shuffle = T) | DimPlot(subset.rna.obj, group.by = "SampleType", shuffle =T)
# DimPlot(subset(samples.vexas.and.control, seurat_clusters %in% hsc.clusters), size = 1.5, 
#         cols = c(`VEXAS MUT` = "red", `VEXAS WT` = "blue", `VEXAS unknown` = "green", `Control` = "orange"), 
#         group.by = "Genotype", order = c("VEXAS MUT", "VEXAS WT"), na.value = "grey80", shuffle = T) +
#   ggtitle("") + coord_fixed() 
# samples.to.plot <- samples.vexas.and.control@meta.data %>% 
#   filter(seurat_clusters %in% hsc.clusters) %>% 
#   count(orig.ident) %>% filter(n > 50) %>% pull(orig.ident)
# samples.vexas.and.control@meta.data %>% 
#   filter(orig.ident %in% samples.to.plot) %>% 
#   filter(seurat_clusters %in% hsc.clusters) %>% 
#   count(seurat_clusters, orig.ident) %>% 
#   ggplot(aes(x = orig.ident, y = n, fill = orig.ident)) +
#   geom_col(position = "dodge") +
#   facet_grid(cols = vars(seurat_clusters)) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") 
# samples.vexas.and.control@meta.data %>% 
#   filter(seurat_clusters %in% hsc.clusters) %>% 
#   count(seurat_clusters, Sample, SampleType) %>% 
#   ggplot(aes(x = Sample, y = n, fill = SampleType)) +
#   geom_col(position = "dodge") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none") 
# 
# DefaultAssay(samples.vexas.and.control) <- "RNA"
# samples.vexas.and.control <- NormalizeData(samples.vexas.and.control)
# samples.vexas.and.control <- ScaleData(samples.vexas.and.control)
# 
# ## Optional: load gene list to highlight
# genes.to.highlight <- read_csv("data/gene_lists_to_highlight/genes_to_highlight_vexas.csv", col_names = "gene")
# 
# ################### Subset for HSCs ################
# 
# ## Subset HSCs
# subset.rna.obj <- subset(samples.vexas.and.control, subset = seurat_clusters %in% hsc.clusters & Sample != "N-02-T1-PB-L")
# # subset.rna.obj <- subset(samples.vexas.and.control, subset = seurat_clusters %in% hsc.clusters)
# Idents(subset.rna.obj) <- "predicted.celltype.l2"
# VlnPlot(subset.rna.obj, features = "HSPCModule1", idents = c("HSC", "LMPP", "EMP"))
# VlnPlot(subset.rna.obj, features = "EPModule1", idents = c("HSC", "LMPP", "EMP"))
# 
# ## Only test non mito, ribo genes
# features.to.test <- rownames(subset.rna.obj@assays$RNA@data)
# features.to.test <- features.to.test[!(grepl("^MT-|^RP|^MALAT1|^MTRNR", features.to.test))]
# 
# ## Remove genes with zero expression
# gene.count.totals <- rowSums(subset.rna.obj@assays$RNA@counts)
# zero.count.genes <- names(gene.count.totals[gene.count.totals == 0])
# features.to.test <- features.to.test[!(features.to.test %in% zero.count.genes)]
# 
# DefaultAssay(subset.rna.obj) <- "RNA"
# Idents(subset.rna.obj) <- "SampleType"
# de_genes <- FindMarkers(subset.rna.obj, features = features.to.test, group.by = "SampleType", ident.1 = "VEXAS", ident.2 = "Control", assay = "RNA", test.use = "wilcox", logfc.threshold = 0, max.cells.per.ident = 300)
# plot_volcano(de_genes, subset.rna.obj@meta.data, pval_field = "p_val", title = "N-02 removed, cluster 2 included, downsampled to 300 cells per group") + coord_cartesian(xlim = c(-4, 4))
# 
# write_de_gene_csv(de_genes, filename = "clusters_2_9_11_N-02_removed_downsampled_300.csv")
# 
# de_analysis.ranked <- de_genes %>%
#     mutate(rank = avg_log2FC * -log10(p_val)) %>%
#     filter(is.finite(rank)) %>% 
#     arrange(-rank)
#   
# ranks.list <- de_analysis.ranked$rank
# names(ranks.list) <- rownames(de_analysis.ranked)
# print(head(ranks.list))
# print(tail(ranks.list))
#   
# pathways <- gmtPathways("~/rscripts/reference/gsea/h.all.v7.5.1.symbols.gmt")
#   
# fgseaRes <- fgsea(pathways = pathways,
#                     stats = ranks.list)
# 
# plotEnrichment(pathways$HALLMARK_INFLAMMATORY_RESPONSE, ranks.list) + ggtitle("HALLMARK_INFLAMMATORY_RESPONSE")
# plotEnrichment(pathways$HALLMARK_TNFA_SIGNALING_VIA_NFKB, ranks.list) + ggtitle("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
# plotEnrichment(pathways$HALLMARK_UNFOLDED_PROTEIN_RESPONSE, ranks.list) + ggtitle("HALLMARK_UNFOLDED_PROTEIN_RESPONSE")
# 
# pathways <- gmtPathways("~/rscripts/reference/gsea/c5.all.v7.5.1.symbols.gmt")
# 
# fgseaRes <- fgsea(pathways = pathways,
#                   stats = ranks.list)
# plotEnrichment(pathways$GOBP_RESPONSE_TO_GROWTH_FACTOR, ranks.list) + ggtitle("GOBP_RESPONSE_TO_GROWTH_FACTOR")
# 
# VlnPlot(subset.rna.obj, features = c("CD81", "IFITM1", "GADD45B", "DNAJB1", "LYZ", "HLA-DRB5"), group.by = "Sample", assay = "RNA", ncol = 3)
# VlnPlot(subset.rna.obj, features = c("CD81", "IFITM1", "GADD45B", "XBP1", "ATF4", "DNAJB9"), group.by = "Sample", assay = "RNA", ncol = 3)
# VlnPlot(subset.rna.obj, features = c("ATF3", "JUN", "FOS", "IL1B", "MCL1", "KLF6"), group.by = "Sample", assay = "RNA", ncol = 3)
# 
# ################### Subset for HSCs ################
# 
# ## Subset Monocytes
# monocyte.clusters <- c(0)
# subset.rna.obj <- subset(samples.vexas.and.control, subset = seurat_clusters %in% monocyte.clusters)
# 
# ## Only test non mito, ribo genes
# features.to.test <- rownames(samples.vexas.and.control@assays$RNA@data)
# features.to.test <- features.to.test[!(grepl("^MT-|^RP|^MALAT1|^MTRNR", features.to.test))]
# 
# ## Remove genes with zero expression
# gene.count.totals <- rowSums(subset.rna.obj@assays$RNA@counts)
# zero.count.genes <- names(gene.count.totals[gene.count.totals == 0])
# features.to.test <- features.to.test[!(features.to.test %in% zero.count.genes)]
# 
# DefaultAssay(subset.rna.obj) <- "RNA"
# Idents(subset.rna.obj) <- "SampleType"
# de_genes.mono <- FindMarkers(subset.rna.obj, features = features.to.test, group.by = "SampleType", ident.1 = "VEXAS", ident.2 = "Control", assay = "RNA", test.use = "wilcox", logfc.threshold = 0, max.cells.per.ident = 1000)
# plot_volcano(de_genes.mono, subset.rna.obj@meta.data)
# 
# 
# de_analysis.ranked <- de_genes.mono %>%
#   mutate(rank = avg_log2FC * -log10(p_val)) %>%
#   filter(is.finite(rank)) %>% 
#   arrange(-rank)
# 
# ranks.list <- de_analysis.ranked$rank
# names(ranks.list) <- rownames(de_analysis.ranked)
# print(head(ranks.list))
# print(tail(ranks.list))
# 
# pathways <- gmtPathways("~/rscripts/reference/gsea/h.all.v7.5.1.symbols.gmt")
# 
# fgseaRes <- fgsea(pathways = pathways,
#                   stats = ranks.list)
# 
# plotEnrichment(pathways$HALLMARK_INFLAMMATORY_RESPONSE, ranks.list)
# 
# 

# DefaultAssay(samples.vexas.and.control) <- "RNA"
# Idents(samples.vexas.and.control) <- "cluster_celltype"
# VlnPlot(samples.vexas.and.control.sub, features = "IL1B", split.by = "SampleType", group.by = "cluster_celltype", idents = c("HSC", "EMP", "GMP", "LMPP", "CD14 Mono"), cols = c(VEXAS = vexas.color, Control = control.color)) + xlab("") +
#   stat_compare_means(method = "wilcox.test", aes(label = paste("Wilcoxon,\np =", ..p.format..)))
