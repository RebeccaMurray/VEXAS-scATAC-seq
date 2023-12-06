library(Signac)
library(Seurat)
library(tidyverse)
library(JASPAR2020)
library(TFBSTools)
# library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(patchwork)
# library(ggrepel)
set.seed(1234)

theme_set(theme_bw())

# source("~/rscripts/general_helpers/volcano_plots.R")
# source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")

plan("multicore", workers = 10)
# options(future.globals.maxSize = 180000 * 1024^2) # for 250 Gb RAM

## Specify cluster as a command-line arg
args = commandArgs(trailingOnly=TRUE)
cluster <- args[[1]]

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")

atac.obj.subset <- subset(atac.obj, subset = cluster_celltype == cluster)
print("ATAC object:")
print(atac.obj.subset)

## Get a list of motif position frequency matrices from JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', species = "Homo sapiens", all_versions = FALSE)
)

## Recall peaks using macs2 - run this using gotcha_signac_env_3_17_23
DefaultAssay(atac.obj.subset) <- "ATAC"
macs2.peaks <- CallPeaks(atac.obj.subset, outdir = "~/rscripts/VEXAS_ATAC_signac/data", fragment.tempdir = "~/rscripts/VEXAS_ATAC_signac/data", macs2.path = "/nfs/sw/macs2/macs2-2.2.7.1/python/bin/macs2")
macs2.counts <- FeatureMatrix(fragments = Fragments(atac.obj.subset), features = macs2.peaks, cells = colnames(atac.obj.subset))
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
atac.obj.subset[["macs2_peaks_HSC"]] <- CreateChromatinAssay(counts = macs2.counts, annotation = annotations)


################# Adding motif information (only need to do once) ##################

# Add motif information
DefaultAssay(atac.obj.subset) <- "macs2_peaks_HSC"
atac.obj.subset <- AddMotifs(object = atac.obj.subset, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "macs2_peaks_HSC", pfm = pfm)

# Compute a per-cell motif activity score
DefaultAssay(atac.obj.subset) <- "macs2_peaks_HSC"
atac.obj.subset <- RunChromVAR(
  object = atac.obj.subset,
  assay = "macs2_peaks_HSC",
  new.assay.name = "chromvar_HSC",
  genome = BSgenome.Hsapiens.UCSC.hg38
)

################## Create an ID to name matrix for motifs #############

## Name to ID matrix
id.to.name <- data.frame()
for (item.name in names(pfm)) {
  print(item.name)
  id.to.name <- rbind(id.to.name, data.frame(ID = pfm[[item.name]]@ID, name = pfm[[item.name]]@name))
}

## Replace chromvar assay ID rownames with gene names
rownames(atac.obj.subset@assays$chromvar_HSC@data) <- plyr::mapvalues(rownames(atac.obj.subset@assays$chromvar_HSC@data), from = id.to.name$ID, to = id.to.name$name)

#################### Save updated ATAC obj ####################

saveRDS(atac.obj.subset, file = paste0("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_", cluster, "_chromvar_20231205.rds"))
# atac.obj.subset <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/atac_obj_cluster_celltype_HSC_chromvar_08_04_23.rds")
# 
# #################### Load expressed TFs ####################
# 
# expressed.genes <- read_csv(paste0("~/rscripts/VEXAS_RNA_seurat/data/expressed_gene_lists/genes_expressed_5_cells_", cluster, ".csv")) %>% pull(feature)
# tf.list <- rownames(atac.obj.subset@assays$chromvar@data)
# names(tf.list) <- tf.list
# 
# decide_if_keep_tf <- function(tf) {
#   tf <- tf
#   tf <- gsub("\\(.*\\)", "", tf)
#   tf.list <- strsplit(tf, "::") %>% unlist()
#   for (tf.current in tf.list) {
#     if(tf.current %in% expressed.genes) {
#       return(TRUE)
#     }
#   }
#   print(paste0("Excluding TF: ", tf))
#   return(FALSE)
# }
# 
# tf.list.annotated <- lapply(tf.list, function(tf) {decide_if_keep_tf(tf)})
# tf.list.filtered <- tf.list.annotated[tf.list.annotated == T] %>% names(.)
# 
# 
# ############### Using LMM - compare all MUT and WT cells ####################
# 
# # For use with Gotcha helper function included in package
# # atac.obj.subset.donor <- subset(atac.obj.subset, orig.ident %in% c("SG_VX16_ATAC", "VEX_BM10_POS_ATAC", "VEX_BM8_POS_ATAC") & predicted.l2 == "HSC" & predicted.l2.score > 0.3)
# atac.obj.subset.donor <- subset(atac.obj.subset, orig.ident %in% c("VEX_BM8_POS_ATAC", "VEX_BM10_POS_ATAC", "SG_VX16_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC"))
# atac.obj.subset.donor$clustering_col <- "Constant"
# print(table(atac.obj.subset.donor$Sample))
# 
# # ## Randomly downsample to 50 MUT/WT cells per sample
# # downsampled.bcs <- atac.obj.subset.donor@meta.data %>%
# #   rownames_to_column("barcode") %>%
# #   filter(Genotype != "NA") %>%
# #   group_by(orig.ident, Genotype) %>%
# #   slice_sample(n=35) %>%
# #   pull(barcode)
# # atac.obj.subset.donor@meta.data[downsampled.bcs, ] %>% group_by(orig.ident, genotype_pred_manual) %>% count()
# #
# # atac.obj.subset.donor <- subset(atac.obj.subset.donor, cells = downsampled.bcs)
# 
# da_motifs.all.cells.genotyped.lmm <- DiffLMMPseudocount(metadata = atac.obj.subset.donor@meta.data,
#                                                  provided.matrix = atac.obj.subset.donor@assays$chromvar@data[tf.list.filtered, ],
#                                                  sample.column = "orig.ident",
#                                                  cluster.column = "clustering_col",
#                                                  selected.clusters = "Constant",
#                                                  treatment.column = "genotype_pred_manual",
#                                                  treatment.levels = c("WT","MUT"),
#                                                  ncores = 18)
# 
# print("Saving...")
# da_motifs.all.cells.genotyped.lmm %>%
#   write_csv(file = "data/chromVAR/HSC_chromvar_analysis_with_lmm_expressed_5_cells.csv")
# 
# print("Done!")
# 
# # ## For later - load result
# da_motifs.all.cells.genotyped.lmm <- read_csv("data/chromVAR/JASPAR2020/HSC_chromvar_analysis_with_lmm_expressed_5_cells.csv")
# 
# da_motifs.all.cells.genotyped.lmm %>%
#   plot_volcano(., metadata.df = atac.obj.subset@meta.data, effect_column = "delta", stat_column = "pval", effect_line = 0.2, stat_line = 0.05)