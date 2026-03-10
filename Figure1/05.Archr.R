
library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn7)

addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 1) 


################################################################################
### ArchR preprocessing
################################################################################
## process the raw data & merge & primary QC
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")

angii_sample_list = c( paste0("MLK", c(1:6)) )
ss_sample_list = c( paste0("RLK", c(1:3, 5:7, 10:13)), "RLK82", "RLK92" )
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
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
    outputDirectory = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/angii"
    output_tmp = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/angii_tmp"
  }else if(sample_list[1]=="RLK1"){
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LK.multiomics.anno.L2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
    outputDirectory = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/ss"
    output_tmp = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/ss_tmp"
  }else if(sample_list[1]=="RLKS1"){
    seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.rds"
    outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
    outputDirectory = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/sp"
    output_tmp = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/sp_tmp"
  }
  
  ArrowFiles = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", sample_list, ".arrow")
  
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    # sampleNames = sample_list,
    outputDirectory = output_tmp,
    copyArrows = TRUE, 
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation
  )

  seurat_object = readRDS(seurat_file)
  proj = subsetArchRProject(proj, cells = gsub("_", "#", colnames(seurat_object)),
                            outputDirectory = outputDirectory, force=T)
  proj$cellNames_mod = gsub("#", "_", proj$cellNames)

  ### https://github.com/GreenleafLab/ArchR/discussions/1563
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

  proj <- addTileMatrix(proj)
  proj <- addGeneScoreMatrix(proj)
  
  proj <- readRDS(outfile)
  proj <- loadArchRProject(outputDirectory,showLogo=FALSE)
  proj$Cell_type_L1 <- make.names(proj$Cell_type_L1)
  
  pathToMacs2 <- "/home/u1/qqiu/.conda/envs/macs2/bin/macs2"
  proj <- addGroupCoverages(proj,
                            groupBy = "Cell_type_L1",
                            minCells = 50,
                            force = T)
  saveRDS(proj, outfile)
  
  proj <- addReproduciblePeakSet(
    ArchRProj = proj,
    groupBy = "Cell_type_L1",
    pathToMacs2 = pathToMacs2,
    reproducibility = "1",
    threads = 1,
    minCells = 50,
    genomeSize = 1.86e9 # for rat: 70% * genome size  = 0.7 * 2.65e9 = 1.86e9
  )
  proj <- addPeakMatrix(proj)
  proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", species = "mus musculus")

  counts = seurat_object@assays$RNA@counts
  #### keep the expression of genes overlapped with ArchR gene annotation (https://github.com/GreenleafLab/ArchR/issues/1084)
  counts = counts[rownames(counts) %in% geneAnnotation$genes$symbol, ]
  colnames(counts)=gsub("[SN]*_", "#", colnames(counts))
  metadata = seurat_object@meta.data
  rownames(metadata)=gsub("[SN]*_", "#", rownames(metadata))
  index = match(rownames(counts), geneAnnotation$genes$symbol)
  ordered_gr = geneAnnotation$genes[index]; names(ordered_gr) <- NULL
  seRNA = SummarizedExperiment(counts, metadata=metadata,
                               rowRanges = ordered_gr)
  seRNA$Group <- paste0(seurat_object@meta.data[gsub("#", "_", colnames(seRNA)), "Sample_ID"])

  proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)

  
  #### add reducedDim? (https://github.com/GreenleafLab/ArchR/discussions/1377)
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
  )

  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(
      resolution = 0.2,
      sampleCells = 10000,
      n.start = 10
    ),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
  )

  proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
  proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
  proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)

  saveRDS(proj, outfile)

  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "LSI_Combined",
    seRNA = seRNA,
    addToArrow = T,
    groupRNA = "Group",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un",
    force = TRUE
  )

  saveRDS(proj, outfile)
  
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "LSI_Combined" # IterativeLSI
  )

  saveArchRProject(ArchRProj = proj, outputDirectory = outputDirectory, load = FALSE)
  
  
}



################################################################################
### Reproducible peaks
################################################################################

peak_folder_m = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/angii/PeakCalls/Cell_type_L1/"
peak_folder_ss = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/ss/PeakCalls/Cell_type_L1/"
peak_folder_shr = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/sp/PeakCalls/Cell_type_L1/"

peak_folder_list = c(peak_folder_m, peak_folder_ss, peak_folder_shr)
project_name = c("AngII", "Salt-sensitive", "Spontaneous")

peak_summary <- list()

k <- 1

for (i in seq_along(peak_folder_list)) {
  
  pi <- peak_folder_list[i]
  proj <- project_name[i]
  
  rds_list <- list.files(
    path = pi,
    pattern = "reproduciblePeaks\\.gr\\.rds$",
    full.names = TRUE
  )
  
  for (ri in rds_list) {
    
    ct <- gsub("-reproduciblePeaks\\.gr\\.rds$", "", basename(ri))
    
    gr_tmp <- readRDS(ri)
    gr_sub <- gr_tmp[gr_tmp$Reproducibility > 1]
    n_peak  <- length(gr_sub)
    
    peak_summary[[k]] <- data.frame(
      project   = proj,
      cell_type = ct,
      n_peak    = n_peak,
      stringsAsFactors = FALSE
    )
    
    k <- k + 1
  }
}

peak_summary_df <- do.call(rbind, peak_summary)

peak_summary_df$cell_type = factor(peak_summary_df$cell_type, levels = c(cell_order, "Immune.cell"))
ggplot(peak_summary_df, aes(x = cell_type, y = n_peak, fill = project)) +
  geom_bar(stat = "identity", position=position_dodge()) +
  theme_classic() +
  labs(x = 'Cell type', y = "Number of reproducible peaks", fill = "Model") +
  theme(text = element_text(family="Arial"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(colour = "black")
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/atac.peak_num.barplot.png", width=512/96, height=300/96, dpi=300)





################################################################################
### Representative loci
################################################################################
proj_m = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds")
proj_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds")
proj_shr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds")

proj_m$Cell_type_L1 <- factor(
  proj_m$Cell_type_L1,
  levels = intersect(cell_order, unique(proj_m$Cell_type_L1))
)

gene_list = c("Lrp2", "Umod", "Slc12a1", "Slc12a3", "Aqp2", "Cdh5")
proj_list = c("proj_m", "proj_ss", "proj_shr")

embedding_list = list()
for(gi in gene_list){
  
  for(pi in proj_list){
    name_tmp = paste0(pi, "_", gi)
    
    proj_tmp = get(pi)
    p_tmp <- plotEmbedding(
      ArchRProj = proj_tmp,
      colorBy = "GeneScoreMatrix",
      name = gi,   
      embedding = "wnn.umap.harmony",  
      plotAs = "points",   
      imputeWeights = getImputeWeights(proj) 
    )
    
    embedding_list[name_tmp] = p_tmp
  }
  
}

embedding_list$proj_m_Lrp2 %>%
  ggplot(aes(x, y, color = color)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "log2(NormCounts+1)",
       title = "ATAC-derived gene activity: Lrp2\n(AngII)") +
  scale_color_gradientn(
    colours = ArchR::paletteContinuous("solarExtra")
  ) +
  theme_ArchR() +
  theme(text = element_text(family = "Arial"))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/atac.lrp2.angii.umap.png", width = 398/96, height = 442/96, dpi = 300)


embedding_list$proj_ss_Slc12a1 %>%
  ggplot(aes(x, y, color = color)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "log2(NormCounts+1)",
       title = "ATAC-derived gene activity: Slc12a1\n(Salt-sensitive)") +
  scale_color_gradientn(
    colours = ArchR::paletteContinuous("solarExtra")
  ) +
  theme_ArchR() +
  theme(text = element_text(family = "Arial"))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/atac.slc12a1.ss.umap.png", width = 398/96, height = 442/96, dpi = 300)


embedding_list$proj_shr_Slc12a3 %>%
  ggplot(aes(x, y, color = color)) +
  geom_point(size = 0.5) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "log2(NormCounts+1)",
       title = "ATAC-derived gene activity: Slc12a3\n(Spontaneous)") +
  scale_color_gradientn(
    colours = ArchR::paletteContinuous("solarExtra")
  ) +
  theme_ArchR() +
  theme(text = element_text(family = "Arial"))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/atac.Slc12a3.shr.umap.png", width = 398/96, height = 442/96, dpi = 300)



p <- plotBrowserTrack(
  ArchRProj = proj_m,
  groupBy = "Cell_type_L1",
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  sizes = c(7, 0.2, 1, 1)
)

grid::grid.newpage()
grid::grid.draw(p$Cdh5)
plotPDF(plotList = p, 
        name = "atac.marker.angii.browser_check.pdf", 
        ArchRProj = proj_m, 
        addDOC = FALSE, width = 5, height = 3.5)





























