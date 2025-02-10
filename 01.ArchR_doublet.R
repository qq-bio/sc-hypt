dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(ArchR)
library(SeuratDisk)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)
library(parallel)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 8) 

## process the raw data & merge & primary QC
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")

# load genome and annotation
# load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
# addArchRGenome("mm10")

## process each sample
# para_list = data.frame(sample_ID = character(),
#                        n_cell_ATAC = numeric(),
#                        n_doublet_ATAC = numeric())

sample_list = c( # paste0("MLK", c(1:6)),  
                 # paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92",
                 # paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4))
                 "RLKW4")

for(i in sample_list){
  sample_id = i
  
  if(grepl("MLK", sample_id)){
    geneAnnotation = geneAnnoMm10
    genomeAnnotation = genomeAnnoMm10
  }else if(grepl("RLK", sample_id)){
    load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
  }
  
  fragment_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/", sample_id, "/outs/atac_fragments.tsv.gz")
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", sample_id, ".ArchR_doublet.rds")
  
  # using ArchR to detect doublet in ATAC data
  input_file = getInputFiles(paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/", sample_id, "/outs"))
  ArrowFiles = createArrowFiles(
    inputFiles = input_file,
    # outputNames = rep(sample_id, length(input_file)),
    sampleNames = names(input_file),
    minTSS = 2, #Dont set this too high because you can always increase later
    minFrags = 100,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    force = TRUE
  )
  
  archr_object <- ArchRProject(
    ArrowFiles = ArrowFiles,
    copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation)
  
  archr_object <- addDoubletScores(archr_object,
                                   k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
                                   knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
                                   LSIMethod = 1, #tf-logidf
                                   UMAPParams = list(n_neighbors = 30, min_dist = 0.3, metric = "cosine", verbose =FALSE),
                                   dimsToUse = 1:30)
  
  # if R^2 less than 0.9 -> little heterogeneity of cells -> worse doublet calling -> skip doublet prediction
  ## R^2 = 0.99582 for MLK1
  
  # archr_object2 = filterDoublets(archr_object)
  
  saveRDS(archr_object, outfile)
  
}






################################################################################

# archr_object <- addIterativeLSI(
#   ArchRProj = archr_object,
#   useMatrix = "TileMatrix", 
#   name = "IterativeLSI", 
#   iterations = 2, 
#   clusterParams = list( #See Seurat::FindClusters
#     resolution = c(0.2), 
#     sampleCells = 10000, 
#     n.start = 10
#   ), 
#   varFeatures = 25000, 
#   dimsToUse = 1:30
# )
# 
# archr_object <- addUMAP(
#   ArchRProj = archr_object, 
#   reducedDims = "IterativeLSI", 
#   name = "UMAP", 
#   nNeighbors = 30, 
#   minDist = 0.5, 
#   metric = "cosine"
# )
# 
# plotEmbedding(ArchRProj = archr_object, colorBy = "colData", name = "DoubletEnrichment", embedding = "UMAP")
