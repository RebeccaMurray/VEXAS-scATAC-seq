library(ggplot2)
library(patchwork)
library(Signac)
library(Seurat)
library(future)
library(GenomeInfoDb)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(dsb)
library(JASPAR2020)
library(TFBSTools)
library(parallel)
library(tximport)
library(Biostrings)
set.seed(1234)
theme_set(theme_classic())

## Define input/output paths
plots_path = "/gpfs/commons/home/rmurray/rscripts/VEXAS_GoTChA_signac/output/integrated_plots/"
cellranger_outs_path = "/gpfs/commons/home/rmurray/VEXAS_GoTChA/"

## Define which samples to process
sample.list <- list(
  `BM5_LLL_CD34_plus`="BM5_LLL_CD34_plus",
  `BM6_LLL_CD34_plus`="BM6_LLL_CD34_plus",
  `BM7_LLL_CD34_plus`="BM7_LLL_CD34_plus",
  `VEX_BM8_POS_ATAC`="VEX_BM8_POS_ATAC",
  `VEX_BM10_POS_ATAC`="VEX_BM10_POS_ATAC",
  `VEX_BM11_LIN_NEG_ATAC` = "VEX_BM11_LIN_NEG_ATAC",
  `SG_VX16_ATAC` = "SG_VX16_ATAC",
  `SG_VX17_ATAC` = "SG_VX17_ATAC",
  `SG_VX18_ATAC` = "SG_VX18_ATAC"
)


######################## Begin analysis - Create Seurat objects for each sample with merged peaks set ##################
# 
# # Read in the peaks for each sample
# read_peaks <- function(sample_name) {
#   df <- read.table(
#     file = paste0(cellranger_outs_path, sample_name, "/outs/peaks.bed"),
#     col.names = c("chr", "start", "end")
#   )
#   df <- makeGRangesFromDataFrame(df)
#   return(df)
# }
# 
# sample.peaks <- mclapply(sample.list, read_peaks, mc.cores =  detectCores())
# 
# # Reduce the genomic ranges for each sample
# combined.peaks <- GenomicRanges::reduce(x = do.call(c, unlist(sample.peaks, use.names = FALSE)))
# 
# # Very minimal QC filtering for peaks
# filter_bad_peaks <- function(combined.peaks.obj) {
#   peakwidths <- width(combined.peaks.obj)
#   combined.peaks.obj <- combined.peaks.obj[peakwidths < 15000 & peakwidths > 20]
#   return(combined.peaks.obj)
# }
# 
# combined.peaks <- filter_bad_peaks(combined.peaks)
# 
# 
# # Build an individual Seurat object for each sample
# create_fragment_object <- function(sample_name, combined_peaks_obj) {
#   md <- read.table(
#     file = paste0(cellranger_outs_path, sample_name, "/outs/singlecell.csv"),
#     stringsAsFactors = FALSE,
#     sep = ",",
#     header = TRUE,
#     row.names = 1
#   )[-1, ] # remove the first row
# 
#   md <- md[md$passed_filters > 500, ]
# 
#   frags <- CreateFragmentObject(
#     path = paste0(cellranger_outs_path, sample_name, "/outs/fragments.tsv.gz"),
#     cells = rownames(md)
#   )
# 
#   counts <- FeatureMatrix(
#     fragments = frags,
#     features = combined_peaks_obj,
#     cells = rownames(md)
#   )
# 
#   atac.obj <- CreateSeuratObject(CreateChromatinAssay(counts, fragments = frags), project = sample_name, assay = "ATAC", meta.data=md)
#   atac.obj$origin <- sample_name
# 
#   return(atac.obj)
# }
# 
# atac.objs <- mclapply(sample.list, function(sample_name) {
#   create_fragment_object(sample_name = sample_name, combined_peaks_obj = combined.peaks)
# }, mc.cores = detectCores())
# 
# # ## Load all samples
# # atac.objs <- lapply(sample.list, function(x) {
# #   readRDS(paste0("data/seurat_objects/individual_samples/", x, ".rds"))
# # })
# 
# 
# ###################################### Perform individual sample QC ###################################
# 
# # extract gene annotations from EnsDb
# annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# 
# # Change to UCSC style, hg38
# seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
# genome(annotations) <- "hg38"
# 
# 
# print('Adding QC metrics...')
# atac.objs <- mclapply(atac.objs, function(atac.obj) {
#     Annotation(atac.obj) <- annotations
# 
#     # compute nucleosome signal score per cell
#     atac.obj <- NucleosomeSignal(object = atac.obj)
# 
#     # compute TSS enrichment score per cell
#     atac.obj <- TSSEnrichment(object = atac.obj, fast = T)
# 
#     # add blacklist ratio and fraction of reads in peaks
#     atac.obj$pct_reads_in_peaks <- atac.obj$peak_region_fragments / atac.obj$passed_filters * 100
#     atac.obj$blacklist_ratio <- atac.obj$blacklist_region_fragments / atac.obj$peak_region_fragments
#     return(atac.obj)
#   }
# , mc.cores = detectCores())
# 
# ## Save all samples
# mclapply(names(atac.objs), function(x) {
#   se <- atac.objs[[x]]
#   saveRDS(se, paste0("data/seurat_objects/individual_samples/", x, ".rds"))
# }, mc.cores = detectCores())
# 
# 
# print('Plotting QC metrics...')
# mclapply(atac.objs, function(atac.obj) {
#     p <- VlnPlot(atac.obj, features = c('pct_reads_in_peaks', 'peak_region_fragments',
#                                         'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
#                  pt.size = 0.1,
#                  ncol = 5,
#     )
#     ggsave(paste0(plots_path, atac.obj@project.name, "_violin_plots.pdf"), device = "pdf", width = 15, height = 6)
#     p <- VlnPlot(atac.obj, features = c('nCount_ATAC', 'pct_reads_in_peaks',
#                                         'TSS.enrichment', 'blacklist_ratio', 'nucleosome_signal'),
#                  pt.size = 0.1,
#                  ncol = 5,
#     )
#     ggsave(paste0(plots_path, atac.obj@project.name, "_violin_plots_with_nCount.pdf"), device = "pdf", width = 15, height = 6)
#     return(p)
#   }
# , mc.cores = detectCores())
# 

################################## QC filtering ###############################

## Load all samples
atac.objs <- lapply(sample.list, function(x) {
  readRDS(paste0("data/seurat_objects/individual_samples/", x, ".rds"))
})

print_cell_count_info <- function(atac.obj) {
  cat(atac.obj@project.name)
  cat('\n')
  cat(nrow(atac.obj@meta.data))
  cat('\n\n')
}
cat('Before filtering:\n')
tmp <- lapply(atac.objs, print_cell_count_info)

## Rename with more formal names
meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_2024-03-29.csv")
names.mapping.list <- meta.data.mapping$updated_name
names(names.mapping.list) <- meta.data.mapping$orig.ident

## Print QC scatter plot
nCount.ATAC.limits.list <- list(
  SG_VX18_ATAC = c(500, 10**3.75),
  SG_VX17_ATAC = c(500, 10**3.75),
  SG_VX16_ATAC = c(500, 10**3.75),
  VEX_BM11_LIN_NEG_ATAC = c(500, 10**3.75),
  VEX_BM10_POS_ATAC = c(1000, 10**4.75),
  VEX_BM8_POS_ATAC = c(500, 10**4.25),
  BM7_LLL_CD34_plus = c(500, 10**4.5), 
  BM6_LLL_CD34_plus = c(500, 10**4.25),
  BM5_LLL_CD34_plus = c(1500, 30000)
)

plist <- lapply(atac.objs, function(x) {
  limits <- nCount.ATAC.limits.list[[x@project.name]]
  if (x@project.name == "SG_VX17_ATAC") {
    x@meta.data %>% 
      filter(TSS.enrichment < 70) %>%  ## Removes 1 outlier cell that skews the plot
      ggplot(aes(x = log10(nCount_ATAC), y = TSS.enrichment)) +
      geom_bin2d(bins = 70) +
      scale_fill_continuous(type = "viridis") +
      ggtitle(names.mapping.list[[x@project.name]]) +
      geom_vline(xintercept = log10(limits[[1]]), linetype = "longdash", color = "red") +
      geom_vline(xintercept = log10(limits[[2]]), linetype = "longdash", color = "red") +
      geom_hline(yintercept = 2.5, linetype = "longdash", color = "red")
  } else {
    x@meta.data %>% 
      ggplot(aes(x = log10(nCount_ATAC), y = TSS.enrichment)) +
      geom_bin2d(bins = 70) +
      scale_fill_continuous(type = "viridis") +
      ggtitle(names.mapping.list[[x@project.name]]) +
      geom_vline(xintercept = log10(limits[[1]]), linetype = "longdash", color = "red") +
      geom_vline(xintercept = log10(limits[[2]]), linetype = "longdash", color = "red") +
      geom_hline(yintercept = 2.5, linetype = "longdash", color = "red")
  }
  
})
names(plist) <- names(atac.objs)
names(plist) <- plyr::mapvalues(names(plist), from = names(names.mapping.list), names.mapping.list)
p.combined <- wrap_plots(plist[names(plist) %>% sort()], ncol = 3)
ggsave(filename = "figures/current_figure_drafts/pre_filter_QC_heatmaps_with_limits_20240328.pdf", plot = p.combined, width = 10, height = 8)


## Print QC violin plots
plist <- lapply(atac.objs, function(x) {
  md <- x@meta.data
  p1 <- md %>% 
    ggplot(aes(x = "Sample", y = pct_reads_in_peaks)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Percent reads in peaks")
  p2 <- md %>% 
    ggplot(aes(x = "Sample", y = nucleosome_signal)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Nucleosome signal")
  p3 <- md %>% 
    ggplot(aes(x = "Sample", y = blacklist_ratio)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Blacklist ratio") 
  p.combined <- (p1 | p2 | p3) + plot_annotation(title = x@project.name)
  return(patchwork::wrap_elements(p.combined))
})
ggsave(filename = paste0("figures/current_figure_drafts/pre_filter_QC_violins.pdf"), plot = p.combined, width = 15, height = 15)




## CD34+
if (exists("SG_VX18_ATAC", atac.objs)) { ## checked
  atac.objs$`SG_VX18_ATAC` <- subset(
    atac.objs$`SG_VX18_ATAC`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**3.75 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("SG_VX17_ATAC", atac.objs)) { ## checked
  atac.objs$`SG_VX17_ATAC` <- subset(
    atac.objs$`SG_VX17_ATAC`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**3.75 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("SG_VX16_ATAC", atac.objs)) { ## checked
  atac.objs$`SG_VX16_ATAC` <- subset(
    atac.objs$`SG_VX16_ATAC`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**3.75 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}


if (exists("VEX_BM11_LIN_NEG_ATAC", atac.objs)) { ## checked
  atac.objs$`VEX_BM11_LIN_NEG_ATAC` <- subset(
    atac.objs$`VEX_BM11_LIN_NEG_ATAC`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**3.75 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 15 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("VEX_BM10_POS_ATAC", atac.objs)) { ## Checked
  atac.objs$`VEX_BM10_POS_ATAC` <- subset(
    atac.objs$`VEX_BM10_POS_ATAC`,
    subset = nCount_ATAC > 1000 &
      nCount_ATAC < 10**4.75 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("VEX_BM8_POS_ATAC", atac.objs)) { ## Checked
  atac.objs$`VEX_BM8_POS_ATAC` <- subset(
    atac.objs$`VEX_BM8_POS_ATAC`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**4.25 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("BM7_LLL_CD34_plus", atac.objs)) { ## Checked
  atac.objs$`BM7_LLL_CD34_plus` <- subset(
    atac.objs$`BM7_LLL_CD34_plus`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**4.5 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("BM6_LLL_CD34_plus", atac.objs)) { ## checked
  atac.objs$`BM6_LLL_CD34_plus` <- subset(
    atac.objs$`BM6_LLL_CD34_plus`,
    subset = nCount_ATAC > 500 &
      nCount_ATAC < 10**4.25 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

if (exists("BM5_LLL_CD34_plus", atac.objs)) { ## checked
  atac.objs$`BM5_LLL_CD34_plus` <- subset(
    atac.objs$`BM5_LLL_CD34_plus`,
    subset = nCount_ATAC > 1500 &
      nCount_ATAC < 30000 &
      TSS.enrichment > 2.5 &
      pct_reads_in_peaks > 30 &
      nucleosome_signal < 5 &
      blacklist_ratio < 0.05
  )
}

cat('After manual filtering:\n')
lapply(atac.objs, print_cell_count_info)

## Print QC scatter plot
plist <- lapply(atac.objs, function(x) {
  x@meta.data %>% 
    ggplot(aes(x = log10(nCount_ATAC), y = TSS.enrichment)) +
    geom_bin2d(bins = 70) +
    scale_fill_continuous(type = "viridis") +
    ggtitle(x@project.name)
  
})
p.combined <- wrap_plots(plist, ncol = 3)
ggsave(filename = "figures/current_figure_drafts/post_filter_QC_heatmaps.pdf", plot = p.combined)

## Print QC violin plots
plist <- lapply(atac.objs, function(x) {
  md <- x@meta.data
  p1 <- md %>% 
    ggplot(aes(x = "Sample", y = pct_reads_in_peaks)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Percent reads in peaks")
  p2 <- md %>% 
    ggplot(aes(x = "Sample", y = nucleosome_signal)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Nucleosome signal")
  p3 <- md %>% 
    ggplot(aes(x = "Sample", y = blacklist_ratio)) +
    geom_violin(fill = "salmon") + xlab(element_blank()) + ggtitle("Blacklist ratio") 
  p.combined <- (p1 | p2 | p3) + plot_annotation(title = x@project.name)
  return(patchwork::wrap_elements(p.combined))
})
ggsave(filename = paste0("figures/current_figure_drafts/post_filter_QC_violins.pdf"), plot = p.combined, width = 15, height = 15)



################################### Remove amulet-identified doublets ##################################

doublet.barcodes <- list.files("/gpfs/commons/groups/landau_lab/VEXAS/doublet_detection/amulet/output/", recursive = T, pattern = "MultipletBarcodes_01.txt")
## Make doublet QC plots
doublet.df.list <- lapply(atac.objs, function(atac.obj) {
  doublet.bcs <- read_csv(paste0("/gpfs/commons/groups/landau_lab/VEXAS/doublet_detection/amulet/output/", atac.obj@project.name, "/MultipletBarcodes_01.txt"), col_names = "barcode") %>% pull(barcode)
  md <- atac.obj@meta.data %>% 
    rownames_to_column("barcode") %>% 
    mutate(doublet_status = if_else(barcode %in% doublet.bcs, "Doublet", "Singlet")) %>% 
    mutate(doublet_status = factor(doublet_status, levels = c("Singlet", "Doublet"))) %>% 
    mutate(barcode = paste0(orig.ident, "_", barcode))
  return(md)
})
doublet.df <- do.call(rbind, doublet.df.list)
doublet.df %>% count(orig.ident, doublet_status)
doublet.df <- doublet.df %>% mutate(formatted_name = plyr::mapvalues(orig.ident, from = names(names.mapping.list), to = names.mapping.list)) %>% 
  mutate(formatted_name = factor(formatted_name, levels = names.mapping.list %>% sort()))

## nCount_ATAC boxplots
doublet.df %>% 
  ggplot(aes(x = doublet_status, y = nCount_ATAC)) +
  geom_boxplot() +
  xlab(element_blank()) +
  facet_wrap(~ formatted_name, scales = "free")
ggsave("figures/current_figure_drafts/ATAC_doublet_plots_nCount_ATAC_20240328.pdf", dpi = 300, width = 6, height = 6)

## nFeature_ATAC boxplots
doublet.df %>% 
  ggplot(aes(x = doublet_status, y = nFeature_ATAC)) +
  geom_boxplot() +
  xlab(element_blank()) +
  facet_wrap(~ formatted_name, scales = "free")
ggsave("figures/current_figure_drafts/ATAC_doublet_plots_nFeature_ATAC_20240328.pdf", dpi = 300, width = 6, height = 6)

## Percent detected doublets per sample
doublet.df %>% 
  group_by(formatted_name) %>% 
  summarize(doublet_fraction = sum(doublet_status == "Doublet") / n()) %>% 
  ggplot(aes(x = reorder(formatted_name, -doublet_fraction), y = doublet_fraction)) +
  geom_col() +
  xlab(element_blank()) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Detected doublets (%)")
ggsave("figures/current_figure_drafts/ATAC_doublet_plots_percents_ATAC_20240328.pdf", dpi = 300, width = 6, height = 4)

## Percent detected doublets per sample (color scheme v2)
doublet.df %>% 
  group_by(formatted_name) %>% 
  count(doublet_status) %>% 
  mutate(frac = n / sum(n)) %>% 
  # summarize(doublet_fraction = sum(doublet_status == "Doublet") / n(), singlet_fraction = sum(doublet_status == "Sing") / n()) %>% 
  ggplot(aes(x = formatted_name, y = frac, fill = doublet_status)) +
  geom_col() +
  xlab(element_blank()) +
  scale_fill_manual(values = c(Singlet = "grey80", Doublet = "black")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab("Percent of cells")
ggsave("figures/current_figure_drafts/ATAC_doublet_plots_percents_ATAC_color_scheme_v2.pdf", dpi = 300, width = 6, height = 4)


## Actually remove the doublets
atac.objs <- lapply(atac.objs, function(atac.obj) {
  doublet.bcs <- read_csv(paste0("/gpfs/commons/groups/landau_lab/VEXAS/doublet_detection/amulet/output/", atac.obj@project.name, "/MultipletBarcodes_01.txt"), col_names = "barcode") %>% pull(barcode)
  atac.obj <- subset(atac.obj, cells = setdiff(atac.obj@meta.data %>% rownames, doublet.bcs))
  return(atac.obj)
})

cat('After removing doublets:\n')
lapply(atac.objs, print_cell_count_info)

################################### Normalize samples individually ##################################

print("Starting normalization...")

## Rename the cell barcodes so they're unique
atac.objs <- lapply(atac.objs, function(atac.obj) {
  RenameCells(atac.obj, add.cell.id = atac.obj@project.name)
})

## Normalize, compute LSI
atac.objs <- mclapply(atac.objs, function(atac.obj) {
  atac.obj <- FindTopFeatures(atac.obj, min.cutoff = 10)
  atac.obj <- RunTFIDF(atac.obj)
  atac.obj <- RunSVD(atac.obj)
  return(atac.obj)
}, mc.cores = detectCores())


### Merge the samples (without integrating)
print('Performing merge...')
atac.combined <- merge(
  x = atac.objs[[1]],
  y = atac.objs[seq(2, length(atac.objs))]
)

## Normalize, compute LSI while merged
print('Normalizing and computing LSI for merged dataset...')
atac.combined <- FindTopFeatures(atac.combined, min.cutoff = 10)
atac.combined <- RunTFIDF(atac.combined)
atac.combined <- RunSVD(atac.combined)


##################################### Perform the integration ####################################

print('Finding integration anchors...')
integration.anchors <- FindIntegrationAnchors(
  object.list = atac.objs,
  anchor.features = rownames(atac.objs[[1]]),
  reduction = "rlsi",
  dims = 2:50
  # k.anchor = 2,
  # k.filter = 5,
  # k.score = 30
)

print('Integrating embeddings...')
integrated <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = atac.combined[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:50
  # k.weight = 30
)

print('Running UMAP...')
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50, return.model = TRUE, seed.use = 1)
integrated <- FindNeighbors(object = integrated, reduction = 'integrated_lsi', dims = 2:50)
integrated <- FindClusters(object = integrated, verbose = FALSE, algorithm = 3, resolution = 1.5, random.seed = 1)


################################ Adding cell types via bridge integration ##########################

## Add the bridge integration results
df <- read_csv("data/bridge_integration_metadata/all_samples_through_BM18_07_24_23.csv") %>% 
  select(predicted.l1, predicted.l1.score, predicted.l2, predicted.l2.score, bridge_integration_barcodes) %>% 
  column_to_rownames("bridge_integration_barcodes")
integrated <- AddMetaData(integrated, df)


################################ Adding genotyping info ############################################

## genotyping file to atac sample name mapping
geno.file.to.sample.name = list(
  `BM5_LLL_CD34_plus` = "BM5_LLL_CD34_plus",
  `BM6_LLL_CD34_plus` = "BM6_LLL_CD34_plus",
  `BM7_LLL_CD34_plus` = "BM7_LLL_CD34_plus",
  `VEX_BM8_POS_redone` = "VEX_BM8_POS_ATAC",
  `VEX_BM10_POS_redone` = "VEX_BM10_POS_ATAC",
  `VEX_BM11_LIN_NEG_redone` = "VEX_BM11_LIN_NEG_ATAC",
  `VEX_BM16` = "SG_VX16_ATAC",
  `SG_VX17_ATAC` = "SG_VX17_ATAC",
  `SG_VX18_ATAC` = "SG_VX18_ATAC"
)

## Load genotyping results
geno.files <- list.files("/gpfs/commons/home/rmurray/rscripts/running_knn_genotyping_VEXAS/data/gotcha_pipeline_outputs", full.names = F) %>% grep("VEX_BM11_LIN_NEG$|VEX_BM10_POS$|BM8_POS$|re_only$|second_variant", ., invert = T, value = T)
geno.dfs <- lapply(geno.files, function(x) { 
  sample.name = geno.file.to.sample.name[[x]]
  read_csv(paste0("/gpfs/commons/home/rmurray/rscripts/running_knn_genotyping_VEXAS/data/gotcha_pipeline_outputs/", x, "/", x, "_final_output_manual_adjustment.csv")) %>% 
    mutate(WhiteListMatch = paste0(sample.name, "_", WhiteListMatch, "-1")) %>% 
    column_to_rownames("WhiteListMatch")
})
geno.df <- do.call(rbind, geno.dfs)

## Make sure barcodes match
geno.bcs <- rownames(geno.df)
atac.bcs <- integrated@meta.data %>% rownames()
intersect.bcs <- intersect(geno.bcs, atac.bcs)
print(setdiff(geno.bcs, atac.bcs))
geno.df <- geno.df[atac.bcs, ] ## Adds NAs for missing barcodes

## Add genotyping to metadata
integrated <- AddMetaData(integrated, geno.df)

integrated@meta.data %>% count(orig.ident, genotype_pred_manual)

################################## Add ADT data + run DSB #########################################


pull.counts.and.run.DSB <- function(sample.name, dir.name) {
  # Load adt data for both real and empty droplets
  print(sample.name)
  base.asap.path.all.barcodes <- paste0("/gpfs/commons/groups/landau_lab/VEXAS/ADT_quantification/ASAP/methods/kallisto/output/", dir.name, "_empty_droplets/counts_unfiltered/")
  adt_data <- Matrix::readMM(paste0(base.asap.path.all.barcodes, "cells_x_features.mtx"))
  adt_bcs <- read.csv(paste0(base.asap.path.all.barcodes, "cells_x_features.barcodes.txt"), header = FALSE) %>% pull(V1) %>% paste0(sample.name, "_", ., "-1")
  adt_tag_names <- read.csv(paste0(base.asap.path.all.barcodes, "cells_x_features.genes.txt"),  header = FALSE) %>% pull(V1)
  rownames(adt_data) <- adt_bcs
  colnames(adt_data) <- adt_tag_names
  adt_data <- as.matrix(adt_data)
  
  # Load real cell barcodes as defined by cellranger
  real.bcs <- read_table(paste0("/gpfs/commons/home/rmurray/VEXAS_GoTChA/", sample.name, "/outs/filtered_peak_bc_matrix/barcodes.tsv"), col_names = "BC") %>% 
    mutate(BC = paste0(sample.name, "_", BC)) %>% 
    pull(BC)
  
  head(real.bcs)
  head(rownames(adt_data))
  length(real.bcs)
  length(rownames(adt_data))
  
  setdiff(real.bcs, rownames(adt_data))
  setdiff(rownames(adt_data), real.bcs)
  
  # Filter for only barcodes present in the ADT data, the cellranger output, and the QC-filtered seurat object
  real.bcs <- real.bcs[(real.bcs %in% adt_bcs)]
  real.bcs <- intersect(integrated@meta.data %>% filter(orig.ident == sample.name) %>% rownames(), real.bcs)
  
  # Separate real and empty cells
  adt.cells = adt_data[real.bcs,]
  adt.empty = adt_data[!(rownames(adt_data) %in% real.bcs),]
  Prot = t(adt.cells)
  empty = t(adt.empty)
  
  isotypes <- adt_tag_names[grepl("CTRL|Isotype", adt_tag_names)]
  print(isotypes)
  
  # Run DSB
  normProt = DSBNormalizeProtein(
    cell_protein_matrix = Prot,
    empty_drop_matrix = empty,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = isotypes
  )
  
  # If there are any atac barcodes not present in adt, add them to the DSB output
  integrated.bcs <- integrated@meta.data %>% dplyr::filter(orig.ident == sample.name) %>% rownames()
  bc_missing_from_adt <- integrated.bcs[!(integrated.bcs %in% colnames(normProt))]
  print("Number of barcodes missing from normProt data for:")
  print(sample.name)
  print(length(bc_missing_from_adt))
  missing_matrix <- matrix(NA_real_, nrow=nrow(normProt), ncol=length(bc_missing_from_adt), dimnames=list(rownames(normProt), bc_missing_from_adt))
  normProt <- cbind(normProt, missing_matrix)
  
  return(normProt)
}


## Add CITE counts for each sample
sample.list <- list(
  BM5_LLL_CD34_plus = "BM5_LLL_CD34_plus_cite_reformatted",
  BM6_LLL_CD34_plus = "BM6_LLL_CD34_plus_CITE_reformatted",
  BM7_LLL_CD34_plus = "BM7_LLL_CD34_plus_CITE_reformatted",
  VEX_BM8_POS_ATAC = "VEX_BM8_POS_CITE_reformatted",
  VEX_BM10_POS_ATAC = "VEX_BM10_POS_CITE_reformatted",
  VEX_BM11_LIN_NEG_ATAC = "VEX_BM11_LIN_NEG_CITE_reformatted"
  
)

normalized.counts.list <- list()
for (sample in names(sample.list)) {
  adt.dir.name = sample.list[[sample]]
  
  print("Running DSB...")
  
  normalized.adt.counts <- pull.counts.and.run.DSB(sample, adt.dir.name)
  
  print("Done!")
  normalized.counts.list[[sample]] <- normalized.adt.counts
  print("Saved!")
  
}

normalized.counts.compiled <- do.call(cbind, normalized.counts.list)


# Make a list of the barcodes present in atac, but missing in adt
integrated.bcs <- rownames(integrated@meta.data)
bc_missing_from_adt <- integrated.bcs[!(integrated.bcs %in% colnames(normalized.counts.compiled))]
missing_matrix <- matrix(NA_real_, nrow=nrow(normalized.counts.compiled), ncol=length(bc_missing_from_adt), dimnames=list(rownames(normalized.counts.compiled), bc_missing_from_adt))
normalized.counts.compiled.expanded <- cbind(normalized.counts.compiled, missing_matrix)


## Add as a new assay
normalized.counts.compiled.expanded <- normalized.counts.compiled.expanded[, integrated@meta.data %>% rownames]
integrated[["ADT_DSB"]] <- CreateAssayObject(data = normalized.counts.compiled.expanded)


################################### Nicely format some fields for figures ##########################

integrated$Sample <- plyr::mapvalues(
  integrated$orig.ident, 
  c("BM5_LLL_CD34_plus", "BM6_LLL_CD34_plus", "BM7_LLL_CD34_plus", "VEX_BM8_POS_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM11_LIN_NEG_ATAC", "SG_VX16_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC"),
  c("PT05", "PT06", "PT07", "PT08", "PT10", "PT11", "PT16", "PT17", "PT18")
)
integrated$Sample <- factor(integrated$Sample, levels = c("PT05", "PT06", "PT07", "PT08", "PT10", "PT11", "PT16", "PT17", "PT18"))
integrated$Donor <- integrated$Sample

# celltypes.plot <- integrated@meta.data %>%  count(predicted.l2) %>% filter(n > 200) %>% pull(predicted.l2)
# integrated@meta.data[["CellType"]] <- integrated$predicted.l2
# integrated$CellType[!(integrated$CellType %in% celltypes.plot)] <- "Other"
# integrated$CellType <- integrated$cluster_celltype

integrated$Genotype <- integrated$genotype_pred_manual
integrated$Genotype[is.na(integrated$Genotype)] <- "NA"
integrated$Genotype <- factor(integrated$Genotype, levels = c("WT", "MUT", "NA"))

################################# Adding TSB (phospho-seq) data #####################################

phospho.seq.samples <- c(
  "SG_VX16_ATAC",
  "SG_VX17_ATAC",
  "SG_VX18_ATAC"
)

tsb.dfs <- lapply(phospho.seq.samples, function(sample.name) {
  adt.files <- file.path(paste0("/gpfs/commons/groups/landau_lab/VEXAS/ADT_quantification/Phospho_Seq/antibodies_phospho_seq/", sample.name, "/TSB_alevin/alevin_output/alevin/quants_mat.gz"))
  adt.names <- read_tsv(paste0("/gpfs/commons/groups/landau_lab/VEXAS/ADT_quantification/Phospho_Seq/antibodies_phospho_seq/", sample.name, "/TSB_alevin/features_TSB_ID_first.tsv"), col_names = c("feature", "barcode"))
  
  ## Load Alevin outputs
  adt.txi <- tximport(files = adt.files, type = "alevin")$counts
  rownames(adt.txi) <- plyr::mapvalues(rownames(adt.txi), from = adt.names$barcode, to = adt.names$feature)
  reverse.complemented.barcodes <- lapply(colnames(adt.txi), function(x) {
    dna.string <- DNAString(x)
    rev.comp <- reverseComplement(dna.string)
    return(toString(rev.comp))
  })
  colnames(adt.txi) <- paste0(sample.name, "_", reverse.complemented.barcodes, "-1")
  return(adt.txi)
})

## Add to object
tsb.df <- do.call(cbind, tsb.dfs)
barcode.list <- intersect(colnames(tsb.df), integrated@meta.data %>% rownames())
tsb.df.formatted <- tsb.df[, barcode.list]
integrated[["ADT_TSB"]] <- CreateAssayObject(counts = tsb.df.formatted)

## Normalize with CLR
integrated <- NormalizeData(integrated, normalization.method = "CLR", margin = 2, assay = "ADT_TSB")

################################# Call peaks separately for each cluster ############################

## Recall peaks using macs2 - run this using gotcha_signac_env_3_17_23
DefaultAssay(integrated) <- "ATAC"
macs2.peaks <- CallPeaks(integrated, 
                         group_by = "seurat_clusters",
                         outdir = "~/rscripts/VEXAS_ATAC_signac/data", 
                         fragment.tempdir = "~/rscripts/VEXAS_ATAC_signac/data", 
                         macs2.path = "/nfs/sw/macs2/macs2-2.2.7.1/python/bin/macs2")
macs2.counts <- FeatureMatrix(fragments = Fragments(integrated), features = macs2.peaks, cells = colnames(integrated))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
integrated[["macs2_peaks"]] <- CreateChromatinAssay(counts = macs2.counts, annotation = annotations)


################################### Run ChromVar for all clusters ####################################

## Get a list of motif position frequency matrices from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)

# add motif information
DefaultAssay(integrated) <- "ATAC"
integrated <- AddMotifs(object = integrated, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "ATAC", pfm = pfm)

# Compute a per-cell motif activity score
DefaultAssay(integrated) <- "ATAC"
integrated <- RunChromVAR(
  object = integrated,
  assay = "ATAC",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

## Name to ID matrix
id.to.name <- data.frame()
for (item.name in names(pfm)) {
  print(item.name)
  id.to.name <- rbind(id.to.name, data.frame(ID = pfm[[item.name]]@ID, name = pfm[[item.name]]@name))
}

## Replace chromvar assay ID rownames with gene names
rownames(integrated@assays$chromvar@data) <- plyr::mapvalues(rownames(integrated@assays$chromvar@data), from = id.to.name$ID, to = id.to.name$name)


#################################### Save! #######################################


saveRDS(integrated, file = "data/seurat_objects/vexas_ATAC_final_custom_limits.RDS")
