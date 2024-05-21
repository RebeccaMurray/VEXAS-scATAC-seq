library(annotatr)
library(tidyverse)
library(GenomicRanges)
library(GenomeInfoDb)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

source("~/rscripts/general_helpers/custom_fgsea_helpers.R")

# # Another option - SEACR peaks
# ## peak set files
# input.bed.files <- list.files("/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/", recursive = T, pattern = ".bed", full.names = T) %>%
#   grep("SEACR", ., value = T) %>%
#   grep("Treated", ., value = T) %>%
#   grep("relaxed", ., value = T)
# 
# granges.objs <- lapply(input.bed.files, function(x) {
#   peaks <- read_table(x, col_names = c("chr", "start", "end", "", "total_signal", "extra")) %>%
#     dplyr::select(chr, start, end, total_signal)
#   gr <- with(peaks, GRanges(chr, IRanges(start, end), mcols = data.frame(total_signal = total_signal)))
# 
# })

input.bed.files <- list.files("/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run//Treated_ATF4/peakCalling/MACS2_controls/", recursive = T, pattern = ".bed", full.names = T) %>%
  grep("MACS2", ., value = T) %>% 
  grep("Treated", ., value = T)


granges.objs <- lapply(input.bed.files, function(x) {
  peaks <- read_table(x, col_names = c("chr", "start", "end", "path", "total_signal")) %>%
    dplyr::select(chr, start, end, total_signal)
  gr <- with(peaks, GRanges(chr, IRanges(start, end), mcols = data.frame(total_signal = total_signal)))
  
})

# Build the annotations (a single GRanges object)
builtin_annotations() %>% grep("hg38", ., value = T)
annotations = build_annotations(genome = 'hg38', annotations = 'hg38_genes_promoters')

# Intersect the regions we read in with the annotations
granges.anno = annotate_regions(
  regions = granges.objs[[1]],
  annotations = annotations,
  ignore.strand = TRUE,
  quiet = FALSE)

# A GRanges object is returned
print(granges.anno)

## Coerce to datafame
granges.anno <- data.frame(granges.anno)

granges.anno %>% 
  write_csv("data/atf4_cut_and_run/annotated_peak_sets/Treated_ATF4_MACS2_controls_annotations_peaks.csv")


####################################### Run hypergeometric test - GO:BP ###########################

gene.list <- unique(granges.anno$annot.symbol)
ego <- enrichGO(gene = gene.list,
                keyType = "SYMBOL",
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                minGSSize = 30,
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

ego@result %>% View()

################################## Run hypergeometric test - custom modules ###########################

gene.list <- unique(granges.anno$annot.symbol)

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
atf4.df$TermID = "ATF4 targets (Han et al.)"
atf4.df <- atf4.df %>% select(TermID, geneID)
chop.df <- data.frame(geneID = read_csv("~/rscripts/VEXAS_RNA_seurat/data/gene_module_lists/Han_et_al_references/CHOP_only.txt") %>% filter(Feature != "N/A") %>% dplyr::rename(geneID = Feature))
chop.df$TermID = "CHOP only (Han et al.)"
chop.df <- chop.df %>% select(TermID, geneID)
additional.modules <- rbind(additional.modules, atf4.df)
additional.modules <- rbind(additional.modules, chop.df)

hs <- org.Hs.eg.db
mapping <- AnnotationDbi::select(hs, 
                                 keys = unique(additional.modules$geneID),
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")
additional.modules <- merge(additional.modules, mapping, by.x = "geneID", by.y = "SYMBOL")
additional.modules <- additional.modules %>% select(TermID, ENTREZID) %>% dplyr::rename(geneID = ENTREZID) %>% filter(!is.na(geneID))
m_t2g <- rbind(m_t2g, additional.modules)

## Convert our gene list to gene IDs
gene.list.df <- AnnotationDbi::select(hs, 
                                 keys = gene.list,
                                 columns = c("ENTREZID", "SYMBOL"),
                                 keytype = "SYMBOL")


res.custom <- enricher(gene.list.df$ENTREZID, gson = NULL, TERM2GENE = m_t2g, pvalueCutoff = 1)
res.custom@result %>% View
res.custom@result %>% head
res.custom@result %>% 
  write_csv(file = "data/atf4_cut_and_run/annotated_peak_sets/Treated_ATF4_MACS2_controls_annotated_peaks_hypergeometric_results.csv")

p.atf4.gsea <- res.custom@result %>% 
  mutate(ID = gsub("_", " ", ID) %>% str_to_title(.)) %>% 
  slice_min(order_by = pvalue, n = 10) %>%
  arrange(pvalue) %>% 
  ggplot(aes(x = reorder(ID, -pvalue), y = -log10(p.adjust))) +
  geom_col(fill = "#F99C1C", color = "black") +
  coord_flip() +
  xlab(element_blank()) +
  theme_classic() +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash") +
  ggtitle("TAK-243 Treated\npromoter peaks")
p.atf4.gsea
ggsave("figures/current_figure_drafts/atf4_cut_and_run_gsea_bar_plot.pdf", plot = p.atf4.gsea, device = "pdf", dpi = 300, width = 6, height = 3)

complete.gene.list <- AnnotationDbi::select(hs, 
                                      keys = gene.list,
                                      columns = c("ENTREZID", "SYMBOL"),
                                      keytype = "SYMBOL")
print_genes_for_module <- function(module_id_list) {
  print(module_id_list)
  module_id_list <- strsplit(module_id_list, "/")[[1]] %>% as.numeric()
  complete.gene.list %>% 
    filter(ENTREZID %in% module_id_list) %>% 
    arrange(SYMBOL) %>% 
    pull(SYMBOL)
}

print_genes_for_module(res.custom@result %>% filter(ID == "GATA1_01") %>% pull(geneID))  
print_genes_for_module(res.custom@result %>% filter(ID == "NP priming module") %>% pull(geneID))  
print_genes_for_module(res.custom@result %>% filter(ID == "ATF4_Q2") %>% pull(geneID))  
print_genes_for_module(res.custom@result %>% filter(ID == "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS") %>% pull(geneID))  
print_genes_for_module(res.custom@result %>% filter(ID == "KEGG_LYSOSOME") %>% pull(geneID))  


####################################### Run fgsea - custom list ###########################

treated.df <- granges.anno %>% select(mcols.total_signal, annot.symbol) %>% 
  filter(!is.na(annot.symbol)) %>% 
  arrange(-mcols.total_signal) %>% 
  group_by(annot.symbol) %>% 
  slice_max(n = 1, order_by = mcols.total_signal, with_ties = F) %>% 
  arrange(-mcols.total_signal)

ranks.list <- treated.df$mcols.total_signal
names(ranks.list) <- treated.df$annot.symbol

run_fgsea_for_list(ranks.list = ranks.list)


