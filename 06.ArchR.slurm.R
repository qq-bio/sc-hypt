
################################################################################

dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

# addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 1) 


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")
################################################################################
### load dataset
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"

proj = readRDS(outfile)


################################################################################
### motif enrichment for peaks linked to DEG
# DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.multiomics.DEG.out", header=T)
DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", header=T)
DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.2 & DEG$treatment=="SS-LS" & DEG$cell_type=="EC", ]

#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_2_tracks.R
markerGenes = DEG_use$gene_name
p <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "cell_grp",
  geneSymbol = markerGenes,
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  sizes = c(7, 0.2, 1, 1)
)

grid::grid.newpage()
grid::grid.draw(p$Nox4)

saveRDS(p, "rat.ss.LK.strain_wise.EC_DEG.plotTrack.rds")
