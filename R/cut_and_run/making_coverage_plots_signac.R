## SPI1 promoter
atac.obj$const_col = "const"
pdf("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/figures/current_figure_drafts/SPI1_Treated_ATF4_20240328.pdf", width = 6, height = 8)
CoveragePlot(atac.obj, 
             group.by = "const_col",
             bigwig = list(ATF4 = "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_ATF4/alignment/bam/Treated_ATF4.bw",
                           IgG = "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_IgG/alignment/bam/Treated_IgG.bw"), 
             bigwig.type = "coverage", 
             bigwig.scale = "common",
             region = "chr11-47374069-47379756", 
             assay = "ATAC", peaks = F)

dev.off()

## TNFRSF14 promoter
atac.obj$const_col = "const"
pdf("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/figures/current_figure_drafts/TNFRSF14_Treated_ATF4_20240328.pdf", width = 6, height = 8)
CoveragePlot(atac.obj, 
             group.by = "const_col",
             bigwig = list(ATF4 = "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_ATF4/alignment/bam/Treated_ATF4.bw",
                           IgG = "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_IgG/alignment/bam/Treated_IgG.bw"), 
             bigwig.type = "coverage", 
             bigwig.scale = "common",
             region = "chr1-2554064-2558586", 
             assay = "ATAC", peaks = F)

dev.off()


# TNFRSF14 promoter
pdf("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/figures/current_figure_drafts/TNFRSF14_Treated_ATF4_20240328.pdf", width = 4, height = 3)
bigWigs = read_coldata(bws = c(
  "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_ATF4/alignment/bam/Treated_ATF4.bw"
), build = "hg38")
t = track_extract(colData = bigWigs, loci = "chr1:2554064-2558586")
track_plot(summary_list = t, y_min = 0, y_max = 0.8)
dev.off()
pdf("/gpfs/commons/home/rmurray/rscripts/VEXAS_ATAC_signac/figures/current_figure_drafts/TNFRSF14_Treated_IgG_20240328.pdf", width = 4, height = 3)
bigWigs = read_coldata(bws = c(
  "/gpfs/commons/groups/landau_lab/VEXAS/cut_and_run/Treated_IgG/alignment/bam/Treated_IgG.bw"
), build = "hg38")
t = track_extract(colData = bigWigs, loci = "chr1:2554064-2558586")
track_plot(summary_list = t, y_min = 0, y_max = 0.8)
dev.off()