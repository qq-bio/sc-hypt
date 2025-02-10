################################################################################
### https://github.com/GreenleafLab/ArchR/discussions/880
### https://github.com/GreenleafLab/ArchR/blob/c44b50c56323d5995fb3d2fb28a9072255fec644/R/AnnotationPeaks.R#L1053-L1146
computeEnrichment=function (matches = NULL, compare = NULL, background = NULL)
{
  matches <- .getAssay(matches, grep("matches", names(assays(matches)),
                                     value = TRUE, ignore.case = TRUE))
  matchCompare <- matches[compare, , drop = FALSE]
  matchBackground <- matches[background, , drop = FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)
  pOut <- data.frame(feature = colnames(matches), CompareFrequency = matchCompareTotal,
                     nCompare = nrow(matchCompare), CompareProportion = matchCompareTotal/nrow(matchCompare),
                     BackgroundFrequency = matchBackgroundTotal, nBackground = nrow(matchBackground),
                     BackgroundProporition = matchBackgroundTotal/nrow(matchBackground))
  pOut$Enrichment <- pOut$CompareProportion/pOut$BackgroundProporition
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x) {
    p <- -phyper(pOut$CompareFrequency[x] - 1, pOut$BackgroundFrequency[x],
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x],
                 pOut$nCompare[x], lower.tail = FALSE, log.p = TRUE)
    return(p/log(10))
  }) %>% unlist %>% round(4)
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)),
                          0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
  pOut
}

.getAssay=function (se = NULL, assayName = NULL)
{
  .assayNames <- function(se) {
    names(SummarizedExperiment::assays(se))
  }
  if (is.null(assayName)) {
    o <- SummarizedExperiment::assay(se)
  }
  else if (assayName %in% .assayNames(se)) {
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }
  else {
    stop(sprintf("assayName '%s' is not in assayNames of se : %s",
                 assayName, paste(.assayNames(se), collapse = ", ")))
  }
  return(o)
}



################################################################################

dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TFBSTools)
load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')

addArchRThreads(threads = 1) 

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")


################################################################################
### subset input data
# # outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
# proj = readRDS(outfile)
# subset_outfile = gsub("rds", "EC_sub.rds", outfile)
# 
# strain = gsub(".LK.archr.rds", "", basename(outfile))
# outputDirectory = paste0("ArchRSubset_", strain)
# 
# cell_list = rownames(proj@cellColData[proj@cellColData$subclass_level1=="EC" & 
#                                         !(proj@cellColData$subclass_level2 %in% c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")), ])
# proj_tmp = subsetArchRProject(proj, cells = cell_list, outputDirectory = outputDirectory, force = TRUE, dropCells=FALSE)
# proj = proj_tmp
# getAvailableMatrices(proj)
# 
# 
# ### calling peaks
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "subclass_level2", force = TRUE)
# proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "subclass_level2", 
#                                pathToMacs2 = "/home/u1/qqiu/.conda/envs/macs2/bin/macs2",
#                                genomeSize = 1.86e9)
# proj <- addPeakMatrix(proj)
# saveRDS(proj, subset_outfile)


outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
subset_outfile = gsub("rds", "EC_sub.rds", outfile)
proj = readRDS(subset_outfile)


# ### identifying markers
# markersPeaks <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "PeakMatrix", 
#   groupBy = "subclass_level2",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# marker_outfile = gsub("rds", "marker_peaks.rds", subset_outfile)
# saveRDS(markersPeaks, marker_outfile)
# 
# # markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05")
# # DEP_outfile = gsub("rds", "EC_sub.marker_peaks.out", outfile)
# # write.table(markerList, DEP_outfile, col.names = T, row.names = T, sep='\t', quote = F)
# 
# 
# ### motif enrichment in differential peaks
# if("Motif" %ni% names(proj@peakAnnotation)){
#   proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
# }
# 
# enrichMotifs <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# marker_outfile = gsub("rds", "marker_peaks.motif.rds", subset_outfile)
# saveRDS(enrichMotifs, marker_outfile)
# 

# df <- data.frame(TF = rownames(enrichMotifs), mlog10Padj = assay(enrichMotifs)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)


### chromvar deviation
proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj, 
  peakAnnotation = "Motif",
  force = TRUE
)
saveRDS(proj, subset_outfile)
