#Import packages
message("Import packages")
library(Seurat)
library(Signac)
library(ggplot2)
library(msigdbr)
library(fgsea)
library(ggrepel)
library(stringr)
library(reshape2)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(patchwork)
theme_set(theme_classic(base_size = 20))

source("~/rscripts/general_helpers/VEXAS_color_palettes.R")

atac_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_genotyping_updated_DiffLMM2_20240221.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
rna_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
combined.df <- merge(atac_difflmm, rna_difflmm, by = "feature", keep = F, suffixes = c("_atac", "_rna"))

top.genes.mut <- combined.df %>% 
  filter(sign(avg_log2FC_rna) == sign(avg_log2FC_atac)) %>% 
  filter(avg_log2FC_rna > 0) %>% 
  mutate(combined_rank = avg_log2FC_atac * avg_log2FC_rna) %>% 
  arrange(-combined_rank) %>% 
  head(n=10) %>% pull(feature)

top.genes.wt <- combined.df %>% 
  filter(sign(avg_log2FC_rna) == sign(avg_log2FC_atac)) %>% 
  filter(avg_log2FC_rna < 0) %>% 
  mutate(combined_rank = avg_log2FC_atac * avg_log2FC_rna) %>% 
  arrange(-combined_rank) %>% 
  head(n=10) %>% pull(feature)

atac_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_genotyping_updated_DiffLMM2_20240221.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval)) %>% 
  mutate(assay = "ATAC") %>% 
  filter(feature %in% c(top.genes.mut, top.genes.wt))
rna_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval)) %>% 
  mutate(assay = "RNA") %>% 
  filter(feature %in% c(top.genes.mut, top.genes.wt))

combined.long.format <- rbind(atac_difflmm, rna_difflmm)
combined.long.format %>% 
  mutate(assay = factor(assay, levels = c("RNA", "ATAC"))) %>% 
  mutate(size_stat = abs(avg_log2FC)) %>%
  ggplot(aes(x = assay, y = feature, size = size_stat, color = stat*sign(avg_log2FC))) +
  geom_point() +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", limits=c(-20, 20), oob=scales::squish)

p.mut <- combined.long.format %>% 
  filter(feature %in% c(top.genes.mut)) %>% 
  mutate(feature = factor(feature, levels = rev(top.genes.mut))) %>% 
  mutate(size_stat = abs(avg_log2FC)) %>%
  ggplot(aes(x = assay, y = feature, fill = avg_log2FC)) +
  geom_tile() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1.5, 1.5), oob=scales::squish)

p.wt <- combined.long.format %>% 
  filter(feature %in% c(top.genes.wt)) %>% 
  mutate(feature = factor(feature, levels = rev(top.genes.wt))) %>% 
  mutate(size_stat = abs(avg_log2FC)) %>%
  ggplot(aes(x = assay, y = feature, fill = avg_log2FC)) +
  geom_tile() +
  xlab(element_blank()) +
  ylab(element_blank()) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", limits=c(-1.5, 1.5), oob=scales::squish)

p.mut | p.wt | plot_layout(guides = "collect")
ggsave("figures/current_figure_drafts/top_genes_rna_atac.pdf", device = "pdf", width = 7, height = 4, dpi = 300)

combined.df %>% 
  filter(feature %in% c(top.genes.mut, top.genes.wt)) %>% 
  select(feature, stat_rna, avg_log2FC_rna, stat_atac, avg_log2FC_atac) %>% 
  pivot_longer(cols = c(stat_rna, stat_atac), names_to = "assay", values_to = "stat")

combined.df %>% 
  filter(feature %in% c(top.genes.mut, top.genes.wt)) %>% 
  ggplot(aes(x = stat_rna, y = stat_atac)) +
  geom_point() +
  coord_cartesian(xlim = c(0, 20))


## Plot by p-val
combined.df %>% 
  filter(pval_rna < 0.05 & avg_log2FC_rna > 0) %>% 
  filter(pval_atac < 0.05 & avg_log2FC_atac > 0) %>% 
  filter(feature %in% genes.expressed) %>% 
  ggplot(aes(x = stat_atac, y = stat_rna)) + 
  geom_point(alpha = 0.5, color = vexas.genotyping.palette[["MUT"]]) +
  # coord_cartesian(ylim = c(0, 20), xlim = c(0, 20)) +
  geom_text_repel(aes(label = feature), max.overlaps = 5, color = "darkred") +
  stat_cor() +
  xlab("-log10(p-value), ATAC Gene Scores") +
  ylab("-log10(p-value), RNA Gene Expression")
ggsave("figures/current_figure_drafts/RNA_vs_ATAC_up_MUT_pval.pdf", device = "pdf", width = 3.5, height = 3.5)


## Plot by p-val
combined.df %>% 
  filter(pval_rna < 0.05 & avg_log2FC_rna > 0) %>% 
  filter(pval_atac < 0.05 & avg_log2FC_atac > 0) %>% 
  filter(feature %in% genes.expressed) %>% 
  mutate(total_rank = stat_rna*stat_atac) %>% 
  top_n(n=15, wt = total_rank)
  ggplot(aes(x = stat_atac, y = stat_rna)) + 
  geom_point(alpha = 0.5, color = vexas.genotyping.palette[["MUT"]]) +
  # coord_cartesian(ylim = c(0, 20), xlim = c(0, 20)) +
  geom_text_repel(aes(label = feature), max.overlaps = 5, color = "darkred") +
  stat_cor() +
  xlab("-log10(p-value), ATAC Gene Scores") +
  ylab("-log10(p-value), RNA Gene Expression")
ggsave("figures/current_figure_drafts/RNA_vs_ATAC_up_MUT_pval.pdf", device = "pdf", width = 3.5, height = 3.5)

## Plot by log2FC
combined.df %>% 
  filter(pval_rna < 0.05 & avg_log2FC_rna > 0) %>% 
  filter(pval_atac < 0.05 & avg_log2FC_atac > 0) %>% 
  filter(feature %in% genes.expressed) %>% 
  ggplot(aes(x = avg_log2FC_rna, y = avg_log2FC_atac)) + 
  geom_point(alpha = 0.25, color = vexas.genotyping.palette[["MUT"]]) +
  # coord_cartesian(ylim = c(0, 20), xlim = c(0, 20)) +
  geom_label_repel(aes(label = feature), max.overlaps = 5) +
  stat_cor()

## Print ranks
m <- combined.df %>% 
  filter(pval_rna < 0.05 & avg_log2FC_rna > 0) %>% 
  filter(pval_atac < 0.05 & avg_log2FC_atac > 0) %>% 
  filter(feature %in% genes.expressed) %>% 
  mutate(rank = stat_rna * stat_atac) %>% 
  arrange(-rank) %>% 
  head(n=10) %>% 
  select(feature, avg_log2FC_rna, avg_log2FC_atac) %>% 
  column_to_rownames("feature") %>% 
  as.matrix()

pheatmap::pheatmap(m, color=colorRampPalette(c("white", "red"))(50),  breaks = seq(0, 1, length.out = 50), cluster_rows = F)


############### look at overlap ##################
atac_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_DiffLMM2_10p_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
rna_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
rownames(atac_difflmm) <- atac_difflmm$feature
rownames(rna_difflmm) <- rna_difflmm$feature

atac.up.mut <- atac_difflmm %>% 
  filter(pval < 0.05) %>% 
  filter(avg_log2FC > 0) %>% 
  pull(feature)

rna.up.mut <- rna_difflmm %>% 
  filter(pval < 0.05) %>% 
  filter(avg_log2FC > 0) %>% 
  pull(feature)

genes.overlap <- intersect(atac.up.mut, rna.up.mut)
atac_difflmm %>% 
  filter(feature %in% genes.overlap) %>% 
  arrange(pval)


atac.up.mut <- atac_difflmm %>% 
  filter(pval < 0.05) %>% 
  filter(avg_log2FC > 0) %>% 
  pull(feature)

rna.up.mut <- rna_difflmm %>% 
  filter(pval < 0.05) %>% 
  filter(avg_log2FC > 0) %>% 
  pull(feature)

genes.overlap <- intersect(atac.up.mut, rna.up.mut)
atac_difflmm %>% 
  filter(feature %in% genes.overlap) %>% 
  arrange(pval)

############### pval correlation plots ##################
atac_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_DiffLMM2_10p_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
rna_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  mutate(stat = -log10(pval))
rownames(atac_difflmm) <- atac_difflmm$feature
rownames(rna_difflmm) <- rna_difflmm$feature

common.genes <- intersect(atac_difflmm$feature, rna_difflmm$feature)

combined.df <- merge(atac_difflmm, rna_difflmm, by = "feature", keep = F, suffixes = c("_atac", "_rna"))
rownames(combined.df)<-combined.df$feature

combined.df <- filter(combined.df, filter(pval_atac < 0.05 & pval_rna < 0.05))

## Add rank
combined.df <- combined.df %>% arrange(avg_log2FC_atac)
combined.df$ATAC_rank <- NA
combined.df$ATAC_rank <- seq(1, nrow(combined.df))
combined.df <- combined.df %>% arrange(avg_log2FC_rna)
combined.df$RNA_rank <- NA
combined.df$RNA_rank <- seq(1, nrow(combined.df))

combined.df %>% 
  filter(feature %in% genes.overlap) %>% 
  filter(feature %in% genes.expressed) %>% 

combined.df %>% 
  filter(feature %in% genes.overlap) %>% 
  filter(feature %in% genes.expressed) %>% 
  ggplot(aes(x = stat_atac, y = stat_rna)) + 
  geom_point(alpha = 0.25) +
  geom_label_repel(aes(label = feature))
  # coord_cartesian(ylim = c(0, 20), xlim = c(0, 10))
  

combined.df %>% 
  filter(feature %in% genes.overlap) %>% 
  ggplot(aes(x = avg_log2FC_atac, y = avg_log2FC_rna)) + 
  geom_point(alpha = 0.25) +
  geom_text(aes(label = feature))
  
combined.df %>% 
  filter(feature %in% genes.overlap) %>% 
    ggplot(aes(x = ATAC_rank, y = RNA_rank)) + 
    geom_point(alpha = 0.25) +
    geom_text(aes(label = feature))

############### LFC correlation plots ##################
atac_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/HSC_DiffLMM2_10p_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  filter(pval < 0.05) %>% filter(avg_log2FC > 0)
rna_difflmm <- read_csv("/gpfs/commons/home/rmurray/rscripts/VEXAS_RNA_seurat/data/differential_expression/VEXAS_MUT_vs_WT/5p_genotyped_cells_no_RP_MT_renormalized/HSC_5cells_genotyped_20231206.csv") %>% 
  mutate(avg_log2FC = log2(fc)) %>% 
  filter(pval < 0.05) %>% filter(avg_log2FC > 0)
rownames(atac_difflmm) <- atac_difflmm$feature
rownames(rna_difflmm) <- rna_difflmm$feature

common.genes <- intersect(genes.use,union(atac_difflmm$feature,rna_difflmm$feature))

xmodal <- data.frame(feature = common.genes,
                     ATAC = atac_difflmm[common.genes,]$avg_log2FC,
                     RNA = rna_difflmm[common.genes,]$avg_log2FC)
rownames(xmodal)<-xmodal$feature


## Add rank
xmodal <- xmodal %>% arrange(ATAC)
xmodal$ATAC_rank <- NA
xmodal$ATAC_rank <- seq(1, nrow(xmodal))
xmodal <- xmodal %>% arrange(RNA)
xmodal$RNA_rank <- NA
xmodal$RNA_rank <- seq(1, nrow(xmodal))

xmodal$label <- ifelse(abs(xmodal$ATAC) > 0.2 & abs(xmodal$RNA) > 0.2 & sign(xmodal$ATAC) == sign(xmodal$RNA),xmodal$feature,NA)

lm_eqn <- function(df){
  m <- lm(RNA ~ ATAC, df);
  eq <- substitute(italic(RNA) == a + b %.% italic(ATAC)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.expression(eq);                 
}


cortest <- cor.test(xmodal$ATAC,xmodal$RNA,method = "spearman")
cortest.pearson <- cor.test(xmodal$ATAC,xmodal$RNA,method = "pearson")


ggplot(xmodal,aes(x = ATAC,y = RNA))+
  geom_point()+
  geom_text_repel(aes(label = label))+
  geom_hline(yintercept = c(0.2,-0.2),linetype = 2)+
  geom_vline(xintercept = c(0.2,-0.2),linetype = 2)+
  geom_smooth(method = "lm")+labs(subtitle = lm_eqn(xmodal),title = paste0("spearman rho = ",round(cortest$estimate,2), ", p-value = ", scales::scientific(cortest$p.value, digits = 3)))+
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7))
# ggsave("figures/current_figure_drafts/ATAC_RNA_gene_score_correlation_plot.pdf", width = 8, height = 8)

ggplot(xmodal,aes(x = ATAC,y = RNA))+
  geom_point()+
  geom_hline(yintercept = c(0.2,-0.2),linetype = 2)+
  geom_vline(xintercept = c(0.2,-0.2),linetype = 2)+
  geom_smooth(method = "lm")+labs(subtitle = lm_eqn(xmodal),title = paste0("spearman rho = ",round(cortest$estimate,2), ", p-value = ", scales::scientific(cortest$p.value, digits = 3)))+
  coord_cartesian(xlim = c(-0.7,0.7), ylim = c(-0.7,0.7))
ggsave("figures/current_figure_drafts/ATAC_RNA_gene_score_correlation_plot_no_labels.pdf", width = 8, height = 8)

#### Plot by rank
lm_eqn <- function(df){
  m <- lm(RNA_rank ~ ATAC_rank, df);
  eq <- substitute(italic(RNA_rank) == a + b %.% italic(ATAC_rank)*","~~italic(r)^2~"="~r2, 
                   list(a = format(coef(m)[1], digits = 2), 
                        b = format(coef(m)[2], digits = 2), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.expression(eq);                 
}
cortest <- cor.test(xmodal$ATAC_rank,xmodal$RNA_rank,method = "spearman")
ggplot(xmodal,aes(x = ATAC_rank,y = RNA_rank))+
  geom_point()+
  geom_text_repel(aes(label = label))+
  geom_smooth(method = "lm")+labs(subtitle = lm_eqn(xmodal),title = paste0("spearman rho = ",round(cortest$estimate,2), ", p-value = ", scales::scientific(cortest$p.value, digits = 3)))
# ggsave("figures/current_figure_drafts/ATAC_RNA_gene_score_correlation_plot.pdf", width = 8, height = 8)
