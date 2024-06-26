# Sys.setenv(RETICULATE_PYTHON="/gpfs/commons/home/rmurray/.virtualenvs/r-reticulate-gotcha/bin/python")

library(Seurat)
library(Signac)
library(tidyverse)
library(patchwork)
library(ggrepel)
# library(Gotcha)
theme_set(theme_bw())

# source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMMPseudocount.R")
source("/gpfs/commons/home/rmurray/rscripts/general_helpers/DiffLMM2.R")

atac.obj <- readRDS("~/rscripts/VEXAS_ATAC_signac/data/seurat_objects/vexas_ATAC_final_custom_limits_celltypes.RDS")

## Load the VEXAS RNA object
rna.obj <- readRDS("~/rscripts/VEXAS_RNA_seurat/data/seurat_objects/vexas_final_RNA_data_celltypes_umap_updated_20231218.rds")


############################# Running the DE analysis #################################

# Helper function to filter for minimum expression - ATAC assay version
FilterGenes_ArchR <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    # object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@assays$RNA_ArchR@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@assays$RNA@RNA_ArchR[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@assays$RNA_ArchR@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      # object@assays$RNA_ArchR@data <- object@assays$RNA_ArchR@data[genes.use, ] # original way of subsetting
      
      ## modified way of subsetting, removing from entire obj
      counts <- GetAssayData(object, assay = "RNA_ArchR")
      counts <- counts[(which(rownames(counts) %in% genes.use)),]
      object <- subset(object, features = rownames(counts))

      return(object)
    } else {
      return(object)
    }
  }

# Helper function to filter for minimum expression - RNA assay version
FilterGenes <-
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    # object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@assays$RNA@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@assays$RNA@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@assays$RNA@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      # object@assays$RNA@data <- object@assays$RNA@data[genes.use, ] # original way of subsetting
      
      ## modified way of subsetting, removing from entire obj
      counts <- GetAssayData(object, assay = "RNA")
      counts <- counts[(which(rownames(counts) %in% genes.use)),]
      DefaultAssay(object) <- "RNA"
      object <- subset(object, features = rownames(counts))
      return(object)
    } else {
      return(object)
    }
  }


mclapply(c("HSC", "EMP", "LMPP", "MkP", "GMP", "CD14 Mono"), function(celltype) {
  
  print("Running for cluster...")
  print(celltype)
  
  se <- subset(atac.obj, CellType == celltype & genotype_pred_manual %in% c("MUT", "WT") & orig.ident %in% c("VEX_BM8_POS_ATAC", "VEX_BM10_POS_ATAC", "SG_VX16_ATAC", "SG_VX17_ATAC", "SG_VX18_ATAC"))
  
  ## Add an extra QC filtering step
  se <- subset(se, nCount_ATAC > 1500 & pct_reads_in_peaks > 50)
  se@meta.data %>% count(Genotype) %>% print()
  se@meta.data %>% count(Sample, Genotype) %>% print()
  
  ## Formatting the sample and donor columns
  meta.data.mapping <- read_csv("data/formatting_patient_ids/orig_ident_to_study_id_mapping_2024-03-28.csv")
  se$Sample <- plyr::mapvalues(se$orig.ident, from = meta.data.mapping$orig.ident, to = meta.data.mapping$updated_name)
  se$Sample <- factor(se$Sample, levels = meta.data.mapping$updated_name %>% sort())
  se@meta.data %>% count(orig.ident, Sample)
  se$Donor <- gsub("_BMMC|_CD34.*", "", se$Sample) %>% gsub("_[A|B]", "", .)
  se$Donor <- factor(se$Donor, levels = unique(se$Donor) %>% sort())
  se@meta.data %>% count(orig.ident, Sample, Genotype) %>% 
    dplyr::filter(Genotype %in% c("MUT", "WT")) %>% 
    pivot_wider(names_from = Genotype, values_from = n) %>% 
    write_csv("data/metadata/cell_counts_used_gene_score_analysis.csv")

  DefaultAssay(se) <- "RNA_ArchR"
  print(table(se$Sample))
  
  ## Subset RNA object for cluster, keep genes expressed in 10%
  rna.obj.se <- subset(rna.obj, cluster_celltype == celltype)
  rna.obj.se <- FilterGenes(object = rna.obj.se, min.value = 0, min.cells = floor(nrow(rna.obj.se@meta.data)*0.20))
  genes.to.test <- rna.obj@assays$RNA@data %>% rownames()
  se <- subset(se, features = intersect(genes.to.test, rownames(se@assays$RNA_ArchR@data)))
  
  # Subset for only genes with nonzero scores in 20 cells
  se <- FilterGenes_ArchR(object = se, min.value = 0, min.cells = floor(nrow(se@meta.data)*0.20))
  
 
  m <- se@assays$RNA_ArchR@data
  diffLMM.results = DiffLMM2(metadata = se@meta.data,
                                       provided.matrix = m,
                                       sample.column = "orig.ident",
                                       cluster.column = "CellType",
                                       selected.clusters = celltype,
                                       treatment.column = "genotype_pred_manual",
                                       treatment.levels = c("WT","MUT"),
                                       ncores = 18)  
  
  print("Saving output...")
  write_delim(diffLMM.results, file = paste0("data/differential_gene_scores_ArchR/VEXAS_MUT_vs_WT/", celltype, "_genotyping_updated_DiffLMM2_20240304.csv"), delim = ",")
  print("Done!")
}, mc.cores = detectCores())



