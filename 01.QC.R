dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 8) 


################################################################################
### generate arrow

## process the raw data & merge & primary QC
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")

sample_list = c( paste0("MLK", c(1:6)),  
                 paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92",
                 paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)))

for(i in sample_list){
  sample_id = i
  
  if(grepl("MLK", sample_id)){
    geneAnnotation = geneAnnoMm10
    genomeAnnotation = genomeAnnoMm10
  }else if(grepl("RLK", sample_id)){
    load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
  }
  
  fragment_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/", sample_id, "/outs/atac_fragments.tsv.gz")
  ArrowFiles = createArrowFiles(
    inputFiles = fragment_file,
    outputNames = sample_id,
    sampleNames = sample_id,
    minTSS = 2, #Dont set this too high because you can always increase later
    minFrags = 100,
    addTileMat = TRUE,
    addGeneScoreMat = TRUE,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    force = TRUE
  )
  
}



################################################################################
### generate proj, subset, transfer annotation

angii_sample_list = c( paste0("MLK", c(1:6)) )
ss_sample_list = c( paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92" )
shr_sample_list = c( paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)) )
model_list = list(angii_sample_list, ss_sample_list, shr_sample_list)


for(sample_list in model_list){
  
  if(grepl("MLK", sample_list[1])){
    geneAnnotation = geneAnnoMm10
    genomeAnnotation = genomeAnnoMm10
  }else if(grepl("RLK", sample_list[1])){
    load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
  }
  
  if(sample_list[1]=="MLK1"){
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
    outputDirectory = "angii"
  }else if(sample_list[1]=="RLK1"){
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
    outputDirectory = "ss"
  }else if(sample_list[1]=="RLKS1"){
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
    outputDirectory = "sp"
  }
  
  ArrowFiles = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", sample_list, ".arrow")
  
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles, 
    outputDirectory = outputDirectory,
    copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )
  
  seurat_object = readRDS(seurat_file)
  proj = subsetArchRProject(proj, cells = gsub("_", "#", colnames(seurat_object)),
                            outputDirectory = outputDirectory, force=T)
  proj$cellNames_mod = gsub("#", "_", proj$cellNames)
  
  seurat_object = RenameCells(seurat_object, new.names = gsub("_", "#", colnames(seurat_object)))
  DefaultAssay(seurat_object) = "RNA"
  seRNA = as.SingleCellExperiment(seurat_object, assay = "RNA")
  
  #### add scRNA info: https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
  proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)
  
  #### https://github.com/GreenleafLab/ArchR/discussions/1563
  seurat.umap <- seurat_object@reductions$wnn.umap.harmony@cell.embeddings
  df <- DataFrame(row.names=proj$cellNames,
                  "seurat#wnn.umap.harmony1" = seurat.umap[proj$cellNames_mod, "wnnUMAPHarmony_1"],
                  "seurat#wnn.umap.harmony2" = seurat.umap[proj$cellNames_mod, "wnnUMAPHarmony_2"],
                  check.names = FALSE)
  proj@embeddings$wnn.umap.harmony <- SimpleList(df = df, params = list())
  
  seurat.meta = seurat_object@meta.data[proj$cellNames_mod, ]
  for(j in colnames(seurat.meta)){
    proj@cellColData[, j] = seurat.meta[, j]
  }
  
  saveRDS(proj, outfile)
  
}


# plotEmbedding(ArchRProj = proj, name = "subclass_level2", embedding = "wnn.umap.harmony")



################################################################################
### call peaks with macs2
pathToMacs2 <- "/home/u1/qqiu/.conda/envs/macs2/bin/macs2"

proj <- addGroupCoverages(proj, 
                          groupBy = "subclass_level2",
                          minCells = 15)

proj <- addReproduciblePeakSet(
  ArchRProj = proj, 
  groupBy = "subclass_level2", 
  pathToMacs2 = pathToMacs2
)




### gene acticity score

markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)



### DE peaks




### peak to gene linkage





### motif activity
