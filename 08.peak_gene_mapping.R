library(Signac)
library(Seurat)
library(rtracklayer)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/")
################################################################################
### lift over genomic coordinates


# sci_peaks_mm9 <- StringToGRanges(regions = rownames(sci.counts), sep = c("_", "_"))
peak_mm10 = proj@peakSet
peak_mm10 = StringToGRanges(regions = unique(p2geneDF$peakName), sep = c(":", "-"))
mm10_hg38 <- rtracklayer::import.chain("/xdisk/mliang1/qqiu/reference/liftover/mm10ToHg38.over.chain")
peak_hg38 <- rtracklayer::liftOver(x = peak_mm10, chain = mm10_hg38)
names(peak_hg38) <- unique(p2geneDF$peakName)

# discard any peaks that were mapped to >1 region in hg38
correspondence <- elementNROWS(peak_hg38)
peak_hg38 <- peak_hg38[correspondence == 1]
peak_hg38 <- unlist(peak_hg38)
p2geneDF <- p2geneDF[p2geneDF$peakName %in% names(peak_hg38), ]

# rename peaks with mm10 coordinates
rownames(sci.counts) <- GRangesToString(grange = sci_peaks_mm10)







################################################################################
### peak-gene linkage data load
proj = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds")

p2g <- getPeak2GeneLinks(
  ArchRProj = proj,
  corCutOff = 0.45,
  resolution = 10000,
  returnLoops = FALSE
)

p2geneDF = metadata(proj@peakSet)$Peak2GeneLinks
p2geneDF$gene_symbol = mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName = (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
p2geneDF = as.data.frame(p2geneDF)

# p <- plotPeak2GeneHeatmap(proj, groupBy="subclass_level2")

DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.multiomics.DEG.out", header=T)

label_genes = intersect(DEG$gene_name, unique(p2geneDF[p2geneDF$FDR<0.05, ]$gene_symbol))


promoterGR <- promoters(getGenes(proj))
mPromoterGR <- promoterGR[promoterGR$symbol %in% label_genes]
mP2G_GR <- p2gGR[p2gGR$symbol %in% label_genes]

plotLoops <- getPeak2GeneLinks(proj, corCutOff=corrCutoff, resolution = 100)[[1]]
sol <- findOverlaps(resize(plotLoops, width=1, fix="start"), mPromoterGR)
eol <- findOverlaps(resize(plotLoops, width=1, fix="end"), mPromoterGR)
plotLoops <- c(plotLoops[from(sol)], plotLoops[from(eol)])
plotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
plotLoops <- plotLoops[width(plotLoops) > 100]

# sub_plot_loop_list <- list()
# for(pn in names(plot_loop_list)){
#   subPlotLoops <- plot_loop_list[[pn]]
#   sol <- findOverlaps(resize(subPlotLoops, width=1, fix="start"), mPromoterGR)
#   eol <- findOverlaps(resize(subPlotLoops, width=1, fix="end"), mPromoterGR)
#   subPlotLoops <- c(subPlotLoops[from(sol)], subPlotLoops[from(eol)])
#   subPlotLoops$symbol <- c(mPromoterGR[to(sol)], mPromoterGR[to(eol)])$symbol
#   sub_plot_loop_list[[pn]] <- subPlotLoops[width(subPlotLoops) > 100]
# }

# Bracket plot regions around loops
plotRegions <- lapply(label_genes, function(x){
  gr <- range(plotLoops[plotLoops$symbol == x])
  lims <- grLims(gr)
  gr <- GRanges(
    seqnames = seqnames(gr)[1],
    ranges = IRanges(start=lims[1], end=lims[2])
  )
  gr
}) %>% as(., "GRangesList") %>% unlist()
plotRegions <- resize(plotRegions, 
                      width=width(plotRegions) + 0.05*width(plotRegions), 
                      fix="center")



#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_2_tracks.R
markerGenes = DEG_filter[order(DEG_filter$p_val),]$Human.gene.name
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
grid::grid.draw(p$Pik3r3)

saveRDS(p, "mouse.LK.DEG_GWAS.plotTrack.rds")




################################################################################
### peak-gene linkage data load
























