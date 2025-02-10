### shell code of preprocessing
# for i in `ls -d */`; do
# file=$i/outs/atac_fragments.tsv.gz
# temp_file=$(mktemp)
# 
# gunzip -c "$file" | \
# sed -E 's/^(1[0-9]|[1-9]|20|X|Y|MT)\t/chr&/' | \
# gzip > "$temp_file"
# 
# mv "$temp_file" "$file"
# 
# echo "Processed and replaced $file"
# done
# 
# for i in `ls -d */`; do
# mv $i/outs/atac_fragments.tsv.gz.tbi $i/outs/atac_fragments.org.tsv.gz.tbi
# done
# 
# for dir in /xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/*/; do
# if [ "$(basename $dir)" != "RLK10" ]; then
# cd "$dir/outs/"
# gunzip -c atac_fragments.tsv.gz | bgzip > atac_fragments.tsv.bgz && tabix -p bed -b 2 -e 3 atac_fragments.tsv.bgz
# echo "Processed and replaced $dir"
# fi
# done




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

# addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 1) 

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")

################################################################################
### generate arrow

## process the raw data & merge & primary QC
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")
# 
# sample_list = c( paste0("MLK", c(1:6)),  
#   paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92",
#   paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)))
# 
# for(i in sample_list){
#   sample_id = i
#   
#   if(grepl("MLK", sample_id)){
#     geneAnnotation = geneAnnoMm10
#     genomeAnnotation = genomeAnnoMm10
#   }else if(grepl("RLK", sample_id)){
#     load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
#   }
#   
#   fragment_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/", sample_id, "/outs/atac_fragments.tsv.bgz")
#   ArrowFiles = createArrowFiles(
#     inputFiles = fragment_file,
#     outputNames = sample_id,
#     sampleNames = sample_id,
#     minTSS = 2, #Dont set this too high because you can always increase later
#     minFrags = 100,
#     addTileMat = TRUE,
#     addGeneScoreMat = TRUE,
#     geneAnnotation = geneAnnotation,
#     genomeAnnotation = genomeAnnotation,
#     force = TRUE
#   )
#   
# }



################################################################################
### generate proj, subset, transfer annotation
### multiome processing tutorial: https://greenleaflab.github.io/ArchR_2020/Ex-Analyze-Multiome.html
# 
# angii_sample_list = c( paste0("MLK", c(1:6)) )
# ss_sample_list = c( paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92" )
# shr_sample_list = c( paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)) )
# # model_list = list(angii_sample_list, ss_sample_list, shr_sample_list)
# model_list = list(ss_sample_list, shr_sample_list)
# 
# 
# for(sample_list in model_list){
#   
#   if(grepl("MLK", sample_list[1])){
#     geneAnnotation = geneAnnoMm10
#     genomeAnnotation = genomeAnnoMm10
#   }else if(grepl("RLK", sample_list[1])){
    # library(BSgenome.Rnorvegicus.UCSC.rn7)
    # load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
#   }
#   
#   if(sample_list[1]=="MLK1"){
#     seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds"
#     outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
#     outputDirectory = "angii"
#   }else if(sample_list[1]=="RLK1"){
#     seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds"
#     outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
#     outputDirectory = "ss"
#   }else if(sample_list[1]=="RLKS1"){
#     seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds"
#     outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
#     outputDirectory = "sp"
#   }
#   
  # ArrowFiles = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", sample_list, ".arrow")
  # proj <- ArchRProject(
  #   ArrowFiles = ArrowFiles,
  #   outputDirectory = outputDirectory,
  #   copyArrows = TRUE, #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
  #   geneAnnotation = geneAnnotation,
  #   genomeAnnotation = genomeAnnotation
  # )
  # 
  # seurat_object = readRDS(seurat_file)
  # proj = subsetArchRProject(proj, cells = gsub("_", "#", colnames(seurat_object)),
  #                           outputDirectory = outputDirectory, force=T)
  # proj$cellNames_mod = gsub("#", "_", proj$cellNames)
  # 
  # #### https://github.com/GreenleafLab/ArchR/discussions/1563
  # seurat.umap <- seurat_object@reductions$wnn.umap.harmony@cell.embeddings
  # df <- DataFrame(row.names=proj$cellNames,
  #                 "seurat#wnn.umap.harmony1" = seurat.umap[proj$cellNames_mod, "wnnUMAPHarmony_1"],
  #                 "seurat#wnn.umap.harmony2" = seurat.umap[proj$cellNames_mod, "wnnUMAPHarmony_2"],
  #                 check.names = FALSE)
  # proj@embeddings$wnn.umap.harmony <- SimpleList(df = df, params = list())
  # 
  # seurat.meta = seurat_object@meta.data[proj$cellNames_mod, ]
  # for(j in colnames(seurat.meta)){
  #   proj@cellColData[, j] = seurat.meta[, j]
  # }
  # 
  # pathToMacs2 <- "/home/u1/qqiu/.conda/envs/macs2/bin/macs2"
  # proj <- addGroupCoverages(proj,
  #                           groupBy = "subclass_level2",
  #                           minCells = 15, threads = 1)
  # proj <- addReproduciblePeakSet(
  #   ArchRProj = proj,
  #   groupBy = "subclass_level2",
  #   pathToMacs2 = pathToMacs2,
  #   genomeSize = 1.86e9 # for rat: 70% * genome size  = 0.7 * 2.65e9 = 1.86e9
  # )
  # proj <- addPeakMatrix(proj)
  # saveRDS(proj, outfile)
  # 
  # proj <- addMotifAnnotations(proj, motifSet="cisbp", name="Motif", force=TRUE, species = "mus musculus")
  # 
  # counts = seurat_object@assays$RNA@counts
  # #### keep the expression of genes overlapped with ArchR gene annotation (https://github.com/GreenleafLab/ArchR/issues/1084)
  # counts = counts[rownames(counts) %in% geneAnnotation$genes$symbol, ]
  # colnames(counts)=gsub("[SN]*_", "#", colnames(counts))
  # metadata = seurat_object@meta.data
  # rownames(metadata)=gsub("[SN]*_", "#", rownames(metadata))
  # index = match(rownames(counts), geneAnnotation$genes$symbol)
  # ordered_gr = geneAnnotation$genes[index]; names(ordered_gr) <- NULL
  # seRNA = SummarizedExperiment(counts, metadata=metadata,
  #                              rowRanges = ordered_gr)
  # seRNA$Group <- paste0(seurat_object@meta.data[gsub("#", "_", colnames(seRNA)), "orig.ident"])
  # 
  # proj <- addGeneExpressionMatrix(input = proj, seRNA = seRNA, force = TRUE)
  # 
  # #### add reducedDim? (https://github.com/GreenleafLab/ArchR/discussions/1377)
  # proj <- addIterativeLSI(
  #   ArchRProj = proj,
  #   clusterParams = list(
  #     resolution = 0.2,
  #     sampleCells = 10000,
  #     n.start = 10
  #   ),
  #   saveIterations = FALSE,
  #   useMatrix = "TileMatrix",
  #   depthCol = "nFrags",
  #   name = "LSI_ATAC"
  # )
  # 
  # proj <- addIterativeLSI(
  #   ArchRProj = proj,
  #   clusterParams = list(
  #     resolution = 0.2,
  #     sampleCells = 10000,
  #     n.start = 10
  #   ),
  #   saveIterations = FALSE,
  #   useMatrix = "GeneExpressionMatrix",
  #   depthCol = "Gex_nUMI",
  #   varFeatures = 2500,
  #   firstSelection = "variable",
  #   binarize = FALSE,
  #   name = "LSI_RNA"
  # )
  # 
  # proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")
  # proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
  # proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)
  # 
  # saveRDS(proj, outfile)
  # 
  # proj <- addGeneIntegrationMatrix(
  #   ArchRProj = proj, 
  #   useMatrix = "GeneScoreMatrix",
  #   matrixName = "GeneIntegrationMatrix",
  #   reducedDims = "LSI_Combined",
  #   seRNA = seRNA,
  #   addToArrow = T,
  #   force = TRUE,
  #   groupRNA = "Group",
  #   nameCell = "predictedCell_Un",
  #   nameGroup = "predictedGroup_Un",
  #   nameScore = "predictedScore_Un"
  # )
  
#   proj = readRDS(outfile)
#   
  # proj <- addPeak2GeneLinks(
  #   ArchRProj = proj,
  #   reducedDims = "LSI_Combined" # IterativeLSI
  # )
#   
#   saveArchRProject(ArchRProj = proj, outputDirectory = outputDirectory, load = FALSE)
#   
#   saveRDS(proj, outfile)
#   
# }




################################################################################
### load dataset
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"

proj = readRDS(outfile)

proj <- addPeak2GeneLinks(
  ArchRProj = proj,
  k = 100,  # Change to 50 for mouse-VSMC
  reducedDims = "LSI_Combined"  # IterativeLSI
)



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
p2geneDF$cell_type = "all"

P2G_outfile = gsub("rds", "all.p2g.out", outfile)
write.table(p2geneDF, P2G_outfile, col.names = T, row.names = F, sep='\t', quote = F)




################################################################################
### differential peaks (treatment vs control)
# cell_list = rownames(proj@cellColData[proj@cellColData$subclass_level1=="EC", ])
# proj_tmp = subsetArchRProject(proj, cells = cell_list, outputDirectory = outputDirectory, force = TRUE)
# proj = proj_tmp
proj$cell_grp = paste0(proj$subclass_level2, "_", proj$species, "_", proj$treatment)
cell_grp_count = table(proj$cell_grp)
cell_list = unique(proj$subclass_level2)
cell_list = setdiff(unique(proj@cellColData[proj@cellColData$subclass_level1=="EC", ]$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC"))
species_list = unique(proj$species)
treatment_list = intersect(c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"),
                           unique(proj$treatment))

marker_merge = c()
for(ci in cell_list){

  for(si in species_list){

    for(ti in treatment_list[-1]){

      bgdGroups = paste0(ci, "_", si, "_", treatment_list[1])
      useGroups = paste0(ci, "_", si, "_", ti)

      bgd_count <- ifelse(bgdGroups %in% names(cell_grp_count), cell_grp_count[bgdGroups], 0)
      use_count <- ifelse(useGroups %in% names(cell_grp_count), cell_grp_count[useGroups], 0)

      if(bgd_count > 10 & use_count > 10){

        markerTest <- getMarkerFeatures(
          ArchRProj = proj,
          useMatrix = "PeakMatrix",
          normBy = "ReadsInPromoter",
          groupBy = "cell_grp",
          testMethod = "wilcoxon",
          bias = c("TSSEnrichment", "log10(nFrags)"),
          useGroups = useGroups,
          bgdGroups = bgdGroups,
          threads = 1
        )

        markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.05")
        markerList = markerList[[1]]

        if(nrow(markerList)>0){

          markerList$cell = ci
          markerList$species = si
          markerList$control = treatment_list[1]
          markerList$treatment = ti

          marker_merge = rbind(marker_merge, markerList)

        }
      }

    }

  }

}

DEP_outfile = gsub("rds", "EC_sub.DEP.out", outfile)
write.table(marker_merge, DEP_outfile, col.names = T, row.names = T, sep='\t', quote = F)

# pma <- markerPlot(seMarker = markerTest, name = "PT_C57BL/6_AngII 28d", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")




# ### differential peaks (SS vs SD)
# proj$cell_grp = paste0(proj$subclass_level1, "_", proj$species, "_", proj$treatment)
# cell_grp_count = table(proj$cell_grp)
# cell_list = unique(proj$subclass_level1)
# 
# marker_merge = c()
# for(ci in cell_list){
# 
#       bgdGroups = paste0(ci, "_SD_LS")
#       useGroups = paste0(ci, "_SS_LS")
#       
#       bgd_count <- ifelse(bgdGroups %in% names(cell_grp_count), cell_grp_count[bgdGroups], 0)
#       use_count <- ifelse(useGroups %in% names(cell_grp_count), cell_grp_count[useGroups], 0)
#       
#       if(bgd_count > 10 & use_count > 10){
#         
#         markerTest <- getMarkerFeatures(
#           ArchRProj = proj,
#           useMatrix = "PeakMatrix",
#           normBy = "ReadsInPromoter",
#           groupBy = "cell_grp",
#           testMethod = "wilcoxon",
#           bias = c("TSSEnrichment", "log10(nFrags)"),
#           useGroups = useGroups,
#           bgdGroups = bgdGroups,
#           threads = 1
#         )
#         
#         markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.05")
#         markerList = markerList[[1]]
#         
#         if(nrow(markerList)>0){
#           
#           markerList$cell = ci
#           markerList$species = "SS"
#           markerList$control = "SD_LS"
#           markerList$treatment = "SS_LS"
#           
#           marker_merge = rbind(marker_merge, markerList)
#           
#         }
#         
#       }
#       
# }
# 
# DEP_outfile = gsub("rds", "DEP_strainwise.out", outfile)
# write.table(marker_merge, DEP_outfile, col.names = T, row.names = T, sep='\t', quote = F)
# 
# # pma <- markerPlot(seMarker = markerTest, name = "EC_SS_LS", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 0.5", plotAs = "MA")




################################################################################
### peak to gene linkage

strain = gsub(".LK.archr.rds", "", basename(outfile))
outputDirectory = paste0("ArchRSubset_", strain)
### peak to gene linkage at cell type specific level
# p2geneDF_merge = c()
cell_type_list = unique(proj@cellColData[proj@cellColData$subclass_level1=="EC",]$subclass_level2)
for( i in cell_type_list){

  cell_list = rownames(proj@cellColData[proj@cellColData$subclass_level2==i, ])
  proj_tmp = subsetArchRProject(proj, cells = cell_list, outputDirectory = outputDirectory, force = TRUE, dropCells=FALSE)
  # proj_tmp = addPeakMatrix(proj_tmp)
  result1 <- tryCatch({
    proj_tmp <- addPeak2GeneLinks(
      ArchRProj = proj_tmp,
      k = 50,  # Change to 50 for mouse-VSMC
      reducedDims = "LSI_Combined"  # IterativeLSI
    )
    TRUE  # Return TRUE if successful
  }, error = function(e) {
    message("addPeak2GeneLinks failed: ", e$message)
    FALSE  # Return FALSE if there was an error
  })
  
  if (result1) {
    result2 <- tryCatch({
      p2g <- getPeak2GeneLinks(
        ArchRProj = proj_tmp,
        corCutOff = 0.45,
        resolution = 10000,
        returnLoops = FALSE
      )
      TRUE  # Return TRUE if successful
    }, error = function(e) {
      message("getPeak2GeneLinks failed: ", e$message)
      FALSE  # Return FALSE if there was an error
    })
    
    # If the second step is successful, apply other code here
    if (result2) {
      
      p2geneDF = metadata(proj_tmp@peakSet)$Peak2GeneLinks
      p2geneDF$gene_symbol = mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
      p2geneDF$peakName = (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
      p2geneDF = as.data.frame(p2geneDF)
      p2geneDF$cell_type = i
      
      P2G_outfile = gsub("rds", paste0(i, ".EC_sub.p2g.out"), outfile)
      write.table(p2geneDF, P2G_outfile, col.names = T, row.names = F, sep='\t', quote = F)
      
      # p2geneDF_merge = rbind(p2geneDF_merge, p2geneDF)
      
    } 
  }

}

# P2G_outfile = gsub("rds", "EC.p2g.out", outfile)
# write.table(p2geneDF_merge, P2G_outfile, col.names = T, row.names = F, sep='\t', quote = F)





################################################################################
### motif enrichment for cell type specific p2g
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"

proj = readRDS(outfile)
matches = getMatches(proj, "Motif")

strain = gsub(".LK.archr.rds", "", basename(outfile))
p2g_list = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", paste0(strain, ".*EC_sub.p2g.out"))
for(i in p2g_list){
  
  dat = read.table(i, header=T, sep='\t')
  
  selected_peaks = dat[! is.na(dat$FDR) & dat$FDR<0.05, ]$idxATAC %>% unique()
  cor_tf = computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
  cor_tf$cell_type = i
  
  motif_outfile = gsub("out", "motif.out", i)
  write.table(cor_tf, motif_outfile, col.names = T, row.names = F, sep='\t', quote = F)
  
}


################################################################################
### motif enrichment for DEP
strain = gsub(".LK.archr.rds", "", basename(outfile))

DEP_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.DEP.out ")
DEP = read.table(DEP_file, header=T, sep='\t')

motif_merge = c()
matches = getMatches(proj, "Motif")
cell_type_list = unique(proj@cellColData$subclass_level1)
for( i in cell_type_list){
  
  selected_peaks = DEP[! is.na(DEP$FDR) & DEP$FDR<0.05 & DEP$cell==i, ]$idx %>% unique()
  cor_tf = computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
  cor_tf$cell_type = i
  
  motif_merge = rbind(motif_merge, cor_tf)
  
}

motif_outfile = gsub("rds", "DEP_motif.out", outfile)
write.table(motif_merge, motif_outfile, col.names = T, row.names = F, sep='\t', quote = F)


# df <- data.frame(TF = rownames(cor_tf), mlog10Padj = cor_tf$mlog10p)
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
#   geom_point(size = 1) +
#   ggrepel::geom_label_repel(
#     data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
#     size = 4,
#     nudge_x = 2,
#     color = "black"
#   ) + theme_ArchR() +
#   ylab("-log10(P-adj) Motif Enrichment") +
#   xlab("Rank Sorted TFs Enriched") +
#   scale_color_gradientn(colors = paletteContinuous(set = "comet"))





################################################################################
### motif enrichment for peaks linked to DEG
DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/EC.DEG.all.out", header=T)
# DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", header=T)
# DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.2 & DEG$treatment=="SS-LS" & DEG$cell_type=="EC", ]

strain = gsub(".LK.archr.rds", "", basename(outfile))
strain_list = c("mouse"="C57BL/6", "rat.ss"="SS", "rat.sp"="SHR")
strain_use = strain_list[strain]

p2g_list = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", paste0(strain, ".*EC_sub.p2g.out"))
p2geneDF = c()
for(i in p2g_list){
  dat = read.table(i, header=T, sep='\t')
  p2geneDF = rbind(p2geneDF, dat)
}

motif_merge = c()
matches = getMatches(proj, "Motif")
cell_type_list = unique(proj@cellColData[proj@cellColData$subclass_level1=="EC",]$subclass_level2)
for( i in cell_type_list){
  
  DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25 & DEG$tissue=="LK" & DEG$strain==strain_use & DEG$cell_type==i, ]
  selected_peaks = p2geneDF[! is.na(p2geneDF$FDR) & p2geneDF$FDR<0.05 & p2geneDF$cell_type==i & p2geneDF$gene_symbol %in% DEG_use$gene_name, ]$idxATAC %>% unique()
  
  if(length(selected_peaks)>0){
    cor_tf = computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
    cor_tf$cell_type = i
    
    motif_merge = rbind(motif_merge, cor_tf)
    
  }
  
}

motif_outfile = gsub("rds", "EC_sub.DEG_p2g_motif.out", outfile)
write.table(motif_merge, motif_outfile, col.names = T, row.names = F, sep='\t', quote = F)





################################################################################
### motif enrichment for peaks linked to marker genes
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.EC_sub.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj = readRDS(outfile)

strain = gsub(".LK.archr.EC_sub.rds", "", basename(outfile))

p2g_list = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", paste0(strain, ".*p2g.out"))
p2geneDF = c()
for(i in p2g_list){
  dat = read.table(i, header=T, sep='\t')
  p2geneDF = rbind(p2geneDF, dat)
}

DEG_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", strain, ".LK.EC.allmarker.0.25.long.txt")
DEG = read.table(DEG_file, header = T, sep = "\t")
DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25, ]

motif_merge = c()
matches = getMatches(proj, "Motif")
cell_list = unique(DEG_use$cluster)
for( i in cell_list){
  
  DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25 & DEG$cluster==i, ]
  selected_peaks = p2geneDF[! is.na(p2geneDF$FDR) & p2geneDF$FDR<0.05 & p2geneDF$gene_symbol %in% DEG_use$gene, ]$idxATAC %>% unique()
  
  if(length(selected_peaks)>0){
    cor_tf = computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
    cor_tf$cell_type = i
    
    motif_merge = rbind(motif_merge, cor_tf)
    
  }
  
}

motif_outfile = gsub("rds", "marker_gene_p2g_motif.out", outfile)
write.table(motif_merge, motif_outfile, col.names = T, row.names = F, sep='\t', quote = F)






################################################################################
### motif enrichment for peaks linked to marker genes (using p2g from ECs)
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.EC_sub.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj = readRDS(outfile)

strain = gsub(".LK.archr.EC_sub.rds", "", basename(outfile))

p2g_list = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", paste0(strain, ".*EC.*p2g.out"))
p2geneDF = c()
for(i in p2g_list){
  dat = read.table(i, header=T, sep='\t')
  p2geneDF = rbind(p2geneDF, dat)
}

DEG_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", strain, ".LK.EC.allmarker.0.25.long.txt")
DEG = read.table(DEG_file, header = T, sep = "\t")
DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25, ]

motif_merge = c()
matches = getMatches(proj, "Motif")
cell_list = unique(DEG_use$cluster)
for( i in cell_list){
  
  DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25 & DEG$cluster==i, ]
  selected_peaks = p2geneDF[! is.na(p2geneDF$FDR) & p2geneDF$FDR<0.05 & p2geneDF$gene_symbol %in% DEG_use$gene, ]$idxATAC %>% unique()
  
  if(length(selected_peaks)>0){
    cor_tf = computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
    cor_tf$cell_type = i
    
    motif_merge = rbind(motif_merge, cor_tf)
    
  }
  
}

motif_outfile = gsub("rds", "marker_gene_EC_p2g_motif.out", outfile)
write.table(motif_merge, motif_outfile, col.names = T, row.names = F, sep='\t', quote = F)

head(motif_merge[motif_merge$cell_type=="EC_Mecom",])






################################################################################
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
# proj = readRDS(outfile)
# 
# strain = gsub(".LK.archr.rds", "", basename(outfile))
# outputDirectory = paste0("ArchRSubset_", strain)
# subset_outfile = gsub("rds", "EC_sub.rds", outfile)
# 
# cell_list = rownames(proj@cellColData[proj@cellColData$subclass_level1=="EC" & 
#                                       !(proj@cellColData$subclass_level2 %in% c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")), ])
# proj_tmp = subsetArchRProject(proj, cells = cell_list, outputDirectory = outputDirectory, force = TRUE, dropCells=FALSE)
# proj = proj_tmp
# getAvailableMatrices(proj)
# proj <- addGroupCoverages(ArchRProj = proj, groupBy = "subclass_level2", force = TRUE)
# proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "subclass_level2", pathToMacs2 = "/home/u1/qqiu/.conda/envs/macs2/bin/macs2")
# proj <- addPeakMatrix(proj)
# subset_outfile = gsub("rds", "EC_sub.rds", outfile)
# saveRDS(proj, subset_outfile)

# markersPeaks <- getMarkerFeatures(
#   ArchRProj = proj, 
#   useMatrix = "PeakMatrix", 
#   groupBy = "subclass_level2",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.05")
# 
# 
# motifsUp <- peakAnnoEnrichment(
#   seMarker = markersPeaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
# 
# df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
# df <- df[order(df$mlog10Padj, decreasing = TRUE),]
# df$rank <- seq_len(nrow(df))
# heatmapEM <- plotEnrichHeatmap(motifsUp, n = 7, transpose = TRUE)

# if("Motif" %ni% names(proj@peakAnnotation)){
#   proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
# }
# 
# proj <- addBgdPeaks(proj)
# 
# proj <- addDeviationsMatrix(
#   ArchRProj = proj, 
#   peakAnnotation = "Motif",
#   force = TRUE
# )
# saveRDS(proj, subset_outfile)




################################################################################
### compare different motif enrichment results
library(dplyr)

ss_motif_ap = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_p2g_motif.out", header=T, sep='\t')
ss_motif_ec = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')
shr_motif_ap = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_p2g_motif.out", header=T, sep='\t')
shr_motif_ec = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')

for(i in c("ss_motif_ap", "ss_motif_ec", "shr_motif_ap", "shr_motif_ec")){
  
  p = get(i) %>%
    filter(mlog10Padj > 2) %>%
    group_by(cell_type) %>%
    slice_max(order_by = mlog10Padj, n = 10) %>%  # Select top 10 features by mlog10Padj within each cell_type
    ungroup() %>%  # Ungroup after selection
    mutate(feature = factor(feature, levels = unique(feature))) %>%  # Convert feature to a factor
    ggplot(aes(x = cell_type, y = feature, fill = mlog10Padj)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(motif_merge$mlog10Padj, na.rm = TRUE)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # Add a black border
      panel.background = element_rect(fill = NA, color = NA),  # Ensure no background color
      panel.grid = element_blank()  # Remove grid lines for a cleaner look
    ) +
    labs(title = i,
         subtitle = "Top 10 motifs enriched in each cell type",
         x = "", y = "",
         fill = "-log10(p-adj)")
  print(p)
  
}

length(intersect(ss_motif_ap[ss_motif_ap$mlog10Padj>2 & ss_motif_ap$cell_type=="EC_Mecom",]$feature, 
                 ss_motif_ec[ss_motif_ec$mlog10Padj>2 & ss_motif_ec$cell_type=="EC_Mecom",]$feature))
# 67
length(intersect(shr_motif_ap[shr_motif_ap$mlog10Padj>2 & shr_motif_ap$cell_type=="EC_Mecom",]$feature, 
                 shr_motif_ec[shr_motif_ec$mlog10Padj>2 & shr_motif_ec$cell_type=="EC_Mecom",]$feature))
# 111
length(intersect(ss_motif_ec[ss_motif_ec$mlog10Padj>2 & ss_motif_ec$cell_type=="EC_Mecom",]$feature, 
                 shr_motif_ec[shr_motif_ec$mlog10Padj>2 & shr_motif_ec$cell_type=="EC_Mecom",]$feature))
# 43
length(intersect(ss_motif_ap[ss_motif_ap$mlog10Padj>2 & ss_motif_ap$cell_type=="EC_Mecom",]$feature, 
                 shr_motif_ap[shr_motif_ap$mlog10Padj>2 & shr_motif_ap$cell_type=="EC_Mecom",]$feature))
# 73


################################################################################
### subtype-specific motifs
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj = readRDS(outfile)

strain = gsub(".LK.archr.EC_sub.rds", "", basename(outfile))

motifMatrix <- getMatrixFromProject(proj, useMatrix = "MotifMatrix", binarize = FALSE)
motifZScores <- assay(motifMatrix, "z")
motifZScores <- assay(motifMatrix)
cellTypes <- proj$subclass_level2

motif_results <- lapply(unique(cellTypes), function(subtype) {
  targetCells <- which(cellTypes == subtype)
  otherCells <- which(cellTypes != subtype)
  
  p_values <- apply(motifZScores, 1, function(z) {
    wilcox.test(z[targetCells], z[otherCells])$p.value
  })
  
  p_adj <- p.adjust(p_values, method = "BH")
  
  mean_diff <- rowMeans(motifZScores[, targetCells]) - rowMeans(motifZScores[, otherCells])
  
  data.frame(
    motif = rownames(motifZScores),
    p_value = p_values,
    p_adj = p_adj,
    mean_diff = mean_diff,
    cell_type = subtype,
    higher_in_subtype = mean_diff > 0  # TRUE if higher in the target subtype
  )
})
motif_results_df <- bind_rows(motif_results)
motif_results_outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.EC.motif_dev.specificity.out")
write.table(motif_results_df, motif_results_outfile, col.names = T, row.names = F, sep='\t', quote = F)

table(motif_results_df[motif_results_df$p_value<0.05 & motif_results_df$higher_in_subtype==TRUE,]$cell_type)


motif_diff_ss = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.EC.motif_specificity.out", header=T, sep="\t")
intersect(motif_diff_ss[motif_diff_ss$p_value<0.05 & motif_diff_ss$higher_in_subtype==TRUE & motif_diff_ss$cell_type=="EC_Mecom",]$motif, 
          ss_motif_ec[ss_motif_ec$mlog10Padj>1.3 & ss_motif_ec$cell_type=="EC_Mecom",]$feature)
# [1] "Sox15_737"
intersect(motif_diff_ss[motif_diff_ss$p_value<0.05 & motif_diff_ss$higher_in_subtype==TRUE & motif_diff_ss$cell_type=="EC_Mecom",]$motif, 
          ss_motif_ap[ss_motif_ap$mlog10Padj>1.3 & ss_motif_ap$cell_type=="EC_Mecom",]$feature)
# [1] "Sox15_737"

motif_diff_ss = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.EC.motif_dev.specificity.out", header=T, sep="\t")
intersect(motif_diff_ss[motif_diff_ss$p_value<0.05 & motif_diff_ss$higher_in_subtype==TRUE & motif_diff_ss$cell_type=="EC_Mecom",]$motif, 
          ss_motif_ec[ss_motif_ec$mlog10Padj>1.3 & ss_motif_ec$cell_type=="EC_Mecom",]$feature)
intersect(motif_diff_ss[motif_diff_ss$p_value<0.05 & motif_diff_ss$higher_in_subtype==TRUE & motif_diff_ss$cell_type=="EC_Mecom",]$motif, 
          ss_motif_ap[ss_motif_ap$mlog10Padj>1.3 & ss_motif_ap$cell_type=="EC_Mecom",]$feature)

motif_diff_shr = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.EC.motif_dev.specificity.out", header=T, sep="\t")
intersect(motif_diff_shr[motif_diff_shr$p_value<0.05 & motif_diff_shr$higher_in_subtype==TRUE & motif_diff_shr$cell_type=="EC_Mecom",]$motif, 
          shr_motif_ec[shr_motif_ec$mlog10Padj>1.3 & shr_motif_ec$cell_type=="EC_Mecom",]$feature)
intersect(motif_diff_shr[motif_diff_shr$p_value<0.05 & motif_diff_shr$higher_in_subtype==TRUE & motif_diff_shr$cell_type=="EC_Mecom",]$motif, 
          shr_motif_ap[shr_motif_ap$mlog10Padj>1.3 & shr_motif_ap$cell_type=="EC_Mecom",]$feature)
# [1] "Foxi1_317"  "Foxi2_318"  "Foxi3_323"  "Gm5294_340" "Foxl1_362" 

motif_diff_shr = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.EC.motif_dev.specificity.out", header=T, sep="\t")
intersect(motif_diff_shr[motif_diff_shr$p_value<0.05 & motif_diff_shr$higher_in_subtype==TRUE & motif_diff_shr$cell_type=="EC_Mecom",]$motif, 
          shr_motif_ec[shr_motif_ec$mlog10Padj>1.3 & shr_motif_ec$cell_type=="EC_Mecom",]$feature)
intersect(motif_diff_shr[motif_diff_shr$p_value<0.05 & motif_diff_shr$higher_in_subtype==TRUE & motif_diff_shr$cell_type=="EC_Mecom",]$motif, 
          shr_motif_ap[shr_motif_ap$mlog10Padj>1.3 & shr_motif_ap$cell_type=="EC_Mecom",]$feature)
# [1] "Zfp300_180"




################################################################################
### subtype-specific motif through max deviation
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj = readRDS(outfile)

strain = gsub(".LK.archr.EC_sub.rds", "", basename(outfile))

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "subclass_level2")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

delta_list <- list()
for (i in seq_len(ncol(seZ))) {
  subtype_name <- colnames(seZ)[i]
  delta_values <- assay(seZ)[, i] - rowMeans(assay(seZ)[, -i])
  delta_df <- data.frame(
    feature = rowData(seZ)$name,
    delta = delta_values,
    cell_type = subtype_name
  )
  delta_list[[subtype_name]] <- delta_df
}
delta_df_combined <- do.call(rbind, delta_list)
delta_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.EC_sub.motif_delta.out")
write.table(delta_df_combined, delta_file, col.names = T, row.names = F, sep='\t', quote = F)

# head(delta_df_combined)

strain = "rat.ss"
strain = "rat.sp"

delta_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.EC_sub.motif_delta.out")
delta_df = read.table(delta_file, header=T, sep='\t')
enrich_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.EC_sub.marker_gene_EC_p2g_motif.out")
enrich_res = read.table(enrich_file, header=T, sep='\t')
# enrich_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.EC_sub.marker_gene_p2g_motif.out")
# enrich_res = read.table(enrich_file, header=T, sep='\t')

merged_data <- merge(delta_df, enrich_res, by = c("feature", "cell_type")) %>%
  dplyr::filter(cell_type == "EC_Mecom") %>%
  dplyr::mutate(TFRegulator = ifelse(mlog10Padj > 1.3 & delta > 0.35, "Yes", "No")) %>%
  arrange(TFRegulator)

ggplot(merged_data, aes(mlog10Padj, delta, color = TFRegulator)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = merged_data[merged_data$TFRegulator=="Yes",],
    aes(x = mlog10Padj, y = delta, label = feature ),
    size = 4,
    nudge_x = 0,
    color = "black",
    max.overlaps = 20
  )  +
  geom_vline(xintercept = 1.3, lty = "dashed") +
  geom_hline(yintercept = 0.35, lty = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("No"="darkgrey", "Yes"="firebrick3")) +
  xlab("-log10(p-adj) of peak enrichment analysis") +
  ylab("Z-score difference between subtypes") +
  labs(color = "Subtype specific\nTF regulator") +
  facet_wrap(~ cell_type)




ss_motif_ap = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_p2g_motif.out", header=T, sep='\t')
ss_motif_ec = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')
shr_motif_ap = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_p2g_motif.out", header=T, sep='\t')
shr_motif_ec = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')

delta_ss = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.motif_delta.out", header=T, sep='\t')
delta_shr = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.motif_delta.out", header=T, sep='\t')

overlap_motif = intersect(ss_motif_ap[ss_motif_ap$mlog10Padj>1.3 & ss_motif_ap$cell_type=="EC_Mecom",]$feature, 
                          shr_motif_ap[shr_motif_ap$mlog10Padj>1.3 & shr_motif_ap$cell_type=="EC_Mecom",]$feature)
overlap_motif = intersect(ss_motif_ec[ss_motif_ec$mlog10Padj>1.3 & ss_motif_ec$cell_type=="EC_Mecom",]$feature, 
                          shr_motif_ec[shr_motif_ec$mlog10Padj>1.3 & shr_motif_ec$cell_type=="EC_Mecom",]$feature)
overlap_motif = intersect(ss_motif_ec[ss_motif_ec$mlog10Padj>1.3 & ss_motif_ec$cell_type=="Ang EC",]$feature, 
                          shr_motif_ec[shr_motif_ec$mlog10Padj>1.3 & shr_motif_ec$cell_type=="Ang EC",]$feature)

merged_data <- merge(delta_ss, delta_shr, by = c("feature", "cell_type")) %>%
  # dplyr::filter(cell_type == "EC_Mecom") %>%
  dplyr::filter(cell_type == "Ang EC") %>%
  dplyr::filter(feature %in% overlap_motif) %>%
  dplyr::mutate(TFRegulator = ifelse(delta.x > quantile(delta.x, 0.75) & delta.y > quantile(delta.y, 0.75), "Yes", "No")) %>%
  arrange(TFRegulator)

merged_data[merged_data$TFRegulator=="Yes", ]$feature

ggplot(merged_data, aes(delta.x, delta.y, color = TFRegulator)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = merged_data[merged_data$TFRegulator=="Yes",],
    aes(x = delta.x, y = delta.y, label = feature ),
    size = 4,
    nudge_x = 0,
    color = "black",
    max.overlaps = Inf
  )  +
  geom_vline(xintercept = quantile(merged_data$delta.x, 0.75), lty = "dashed") +
  geom_hline(yintercept = quantile(merged_data$delta.y, 0.75), lty = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("No"="darkgrey", "Yes"="firebrick3")) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(min(merged_data$delta.y),
               max(merged_data$delta.y) * 2)) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(min(merged_data$delta.x),
               max(merged_data$delta.x) * 2)) +
  xlab("Z-score difference between subtypes\n(in salt-sensitive model)") +
  ylab("Z-score difference between subtypes\n(in spontaneous model)") +
  labs(color = "Subtype specific\nTF regulator") +
  facet_wrap(~ cell_type)






################################################################################
### positive TF
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
# outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj = readRDS(outfile)

strain = gsub(".LK.archr.EC_sub.rds", "", basename(outfile))

seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "subclass_level2")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = proj,
  useMatrix1 = "GeneIntegrationMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "LSI_Combined"
)

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])


p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = data.frame(corGSM_MM)[corGSM_MM$TFRegulator=="YES",],
    aes(x = cor, y = maxDelta, label = GeneIntegrationMatrix_name ),
    size = 4,
    nudge_x = 0,
    color = "black"
  ) + theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") +
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0),
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p





################################################################################
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds"
proj_ss = readRDS(outfile)
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
proj_shr = readRDS(outfile)

seurat_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
seurat_object_ss = readRDS(seurat_file)
seurat_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")
seurat_object_shr = readRDS(seurat_file)

cell_list = gsub("#", "_", rownames(proj_ss@embeddings$wnn.umap.harmony$df))
proj_ss@embeddings$wnn.umap.harmony$df[,1] = seurat_object_ss@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 1]
proj_ss@embeddings$wnn.umap.harmony$df[,2] = seurat_object_ss@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 2]

cell_list = gsub("#", "_", rownames(proj_shr@embeddings$wnn.umap.harmony$df))
proj_shr@embeddings$wnn.umap.harmony$df[,1] = seurat_object_shr@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 1]
proj_shr@embeddings$wnn.umap.harmony$df[,2] = seurat_object_shr@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 2]

# ss_motif_list = c("Myc_37", "Figla_54")
# shr_motif_list = c("Alx1_505", "Dbx1_480")
# motifs = ss_motif_list

proj = proj_ss
proj = proj_shr

motifs = c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497" )
motifs = c("Alx3_419", "En1_563", "Evx1_410", "Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc5_441", "Hoxc6_406", "Lhx3_459", "Lhx4_454", "Lmx1b_515", "Otp_436", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497")
motifs = c("Gm4881_296", "Etv2_270", "Gata4_386")

markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = sort(markerMotifs),
  embedding = "wnn.umap.harmony",
  imputeWeights = getImputeWeights(proj),
  size = 1.5,
  plotAs = "points"
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 12) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 4),p2))




strain = "rat.ss"
strain = "rat.sp"

if(strain == "rat.ss"){
  proj = proj_ss
}else{
  proj = proj_shr
}

motifs_list = c("Hoxb4_511", "Vax1_413", "Hoxa1_479")

for(motifs in motifs_list){
  
  markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
  markerMotifs
  
  markerMotifs <- grep("z:", markerMotifs, value = TRUE)
  
  p <- plotEmbedding(
    ArchRProj = proj,
    colorBy = "MotifMatrix",
    name = sort(markerMotifs),
    embedding = "wnn.umap.harmony",
    imputeWeights = getImputeWeights(proj),
    size = 2,
    plotAs = "points"
  ) + labs(x="UMAP 1", y="UMAP 2", title = markerMotifs, color = "Motif activity")
  print(p)
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/TFregulatorPlots/", strain, ".", motifs, ".umap.png")
  ggsave(outfile, width=523/96, height=377/96, dpi=300)
  
}

### identify peaks linked to specific motif: https://github.com/GreenleafLab/ArchR/discussions/707




################################################################################
#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_2_tracks.R

proj = proj_ss

p2g_list = Sys.glob("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss*EC*p2g.out")
p2geneDF = c()
for(i in p2g_list){
  dat = read.table(i, header=T, sep='\t')
  p2geneDF = rbind(p2geneDF, dat)
}
metadata(proj@peakSet)$Peak2GeneLinks = p2geneDF

markerGenes = c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
                "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd", "Vegfa", "Mecom")
p <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "subclass_level2",
  # groupBy = "cell_grp",
  # useGroups = c("EC_SS_LS","EC_SD_LS","EC_SD_HS 3d","EC_SS_HS 3d","EC_SS_HS 21d"),
  geneSymbol = markerGenes,
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  sizes = c(7, 0.2, 1, 1)
)

grid::grid.newpage()
grid::grid.draw(p$Mecom)

# p2geneDF[p2geneDF$gene_symbol=="Mecom", ]

saveRDS(p, "rat.ss.LK.strain_wise.EC_DEG.plotTrack.rds")



proj = proj_shr
markerGenes = c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
                "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd", "Vegfa", "Mecom")
p <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "subclass_level2",
  # groupBy = "cell_grp",
  # useGroups = c("EC_SS_LS","EC_SD_LS","EC_SD_HS 3d","EC_SS_HS 3d","EC_SS_HS 21d"),
  geneSymbol = markerGenes,
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  sizes = c(7, 0.2, 1, 1)
)

grid::grid.newpage()
grid::grid.draw(p$Mecom)

# p2geneDF[p2geneDF$gene_symbol=="Mecom", ]

saveRDS(p, "rat.shr.LK.strain_wise.EC_DEG.plotTrack.rds")










################################################################################
#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_3_TFregulators.R
##########################################################################################
# Identify regulatory targets of TFs 
##########################################################################################
# library(RcppAlgos)
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/archr_helpers.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/matrix_helpers.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/misc_helpers.R")


strain = "rat.ss"
strain = "rat.sp"

if(strain=="rat.ss"){
  atac_proj = proj_ss
  p2g_list = Sys.glob("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.*EC.*p2g.out")
  p2geneDF = c()
  for(i in p2g_list){
    dat = read.table(i, header=T, sep='\t')
    p2geneDF = rbind(p2geneDF, dat)
  }
}else{
  atac_proj = proj_shr
  p2g_list = Sys.glob("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.*EC.*p2g.out")
  p2geneDF = c()
  for(i in p2g_list){
    dat = read.table(i, header=T, sep='\t')
    p2geneDF = rbind(p2geneDF, dat)
  }
  
}

# ChromVAR deviations matrix: (rows motif names x cols cell names)
motifMatrix <- getMatrixFromProject(atac_proj, useMatrix="MotifMatrix")
deviationsMatrix <- assays(motifMatrix)$deviations

# GeneIntegration Matrix: (rows gene names x cols cell names)
GIMatrix <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
GImat <- assays(GIMatrix)$GeneIntegrationMatrix
rownames(GImat) <- rowData(GIMatrix)$name
GImat <- as(GImat[Matrix::rowSums(GImat) > 0,], "sparseMatrix") # Remove unexpressed genes

regulators <- c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497" )
deviationsMatrix <- deviationsMatrix[regulators,]

# Identify pseudobulks for performing matrix correlations
knn_groups <- getLowOverlapAggregates(atac_proj, target.agg=500, k=100, 
                                      overlapCutoff=0.8, dimReduc="LSI_Combined")

kgrps <- unique(knn_groups$group)

# GeneIntegrationMatrix
GIMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(GImat[,use_cells])
}) %>% do.call(cbind,.)
colnames(GIMatPsB) <- kgrps

# In rare instances, we can get pseudo-bulked genes that have zero averages
GIMatPsB <- GIMatPsB[Matrix::rowSums(GIMatPsB) > 0,]

# DeviationsMatrix
DevMatPsB <- lapply(kgrps, function(x){
  use_cells <- knn_groups[knn_groups$group==x,]$cell_name
  Matrix::rowMeans(deviationsMatrix[,use_cells])
}) %>% do.call(cbind,.)
colnames(DevMatPsB) <- kgrps

# Perform chromVAR deviations to Integrated RNA correlation analysis:
start <- Sys.time()
geneCorMat <- cor2Matrices(DevMatPsB, GIMatPsB)
colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
end <- Sys.time()
message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))

allGenes <- rownames(GIMatPsB) %>% sort() # Already filtered to only expressed genes


# regulators <- c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497" )
# markerGenes = c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
#                 "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd", "Vegfa", "Mecom")
# geneCorMat[geneCorMat$FDR<0.05 & geneCorMat$motifName %in% regulators & geneCorMat$symbol %in% markerGenes,]
# table(geneCorMat[abs(geneCorMat$Correlation)>0.25 & geneCorMat$motifName %in% regulators & geneCorMat$symbol %in% markerGenes,]$motifName)


# Get locations of motifs of interest:
motifPositions <- getPositions(atac_proj, name="Motif")
motifGR <- stack(motifPositions, index.var="motifName")

# Get peak to gene GR
corrCutoff <- 0.45
metadata(atac_proj@peakSet)$Peak2GeneLinks = p2geneDF
p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
p2gDF <- p2gDF[which(p2gDF$VarQATAC > 0.25 & p2gDF$VarQRNA > 0.25),]
p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
mcols(p2gGR) <- p2gDF
p2gGR

calculateLinkageScore <- function(motifLocs, p2gGR){
  # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
  ###################################
  # For a given gene, the LS = sum(corr peak R2 * motifScore)
  ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
  olGenes <- p2gGR[to(ol)]
  olGenes$motifScore <- motifLocs[from(ol)]$score
  olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
  LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
  LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
  LSdf$rank <- 1:nrow(LSdf)
  return(LSdf)
}

calculateMotifEnrichment <- function(motifLocs, p2gGR){
  # Calculate Motif enrichment per gene
  ###################################
  # For a given gene, calculate the hypergeometric enrichment of motifs in 
  # linked peaks (generally will be underpowered)
  motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
  m <- length(motifP2G) # Number of possible successes in background
  n <- length(p2gGR) - m # Number of non-successes in background
  
  motifLinks <- motifP2G$symbol %>% getFreqs()
  allLinks <- p2gGR$symbol %>% getFreqs()
  df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
  df$motifLinks[is.na(df$motifLinks)] <- 0
  df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
  df <- df[order(df$mLog10pval, decreasing=TRUE),]
  df$symbol <- rownames(df)
  return(df)
}

regPlotDir <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/TFregulatorPlots"
dir.create(regPlotDir, showWarnings = FALSE, recursive = TRUE)

# plot all TF regulators
markerGenes  <- c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
                                  "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd", "Vegfa", "Mecom") %>% unique()

res_list <- list()
##########################################################################################
# regulators <- c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497" )
regulators <- c("Hoxb4_511")

for(motif in regulators){
  motif_short <- strsplit(motif,"_")[[1]][1]
  # First get motif positions
  motifLocs <- motifGR[motifGR$motifName == motif]
  # Calculate Linkage Score for motif
  LS <- calculateLinkageScore(motifLocs, p2gGR)
  # Get just genes correlated to motif
  motifGeneCorDF <- unique(geneCorMat[geneCorMat$motifName == motif,])
  plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
  # Calculate motif enrichment per gene
  ME <- calculateMotifEnrichment(motifLocs, p2gGR)
  plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
  plot_df = plot_df[rowSums(is.na(plot_df))==0,]
  plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
  plot_df$toLabel <- "NO"
  topN <- 5
  plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
  plot_df$rank_LS <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
  plot_df$rank_Corr <- 1:nrow(plot_df)
  plot_df$toLabel[1:topN] <- "YES"
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
  plot_df$rank_Pval <- 1:nrow(plot_df)
  plot_df$toLabel[1:10] <- "YES"
  # plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
  # plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
  # plot_df$toLabel[1:topN] <- "YES"
  # Label any marker genes in window of interest
  LS_window <- quantile(plot_df$LS, 0.75)
  corr_window <- 0.25
  pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
  neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
  if(nrow(plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]) > 0){
    plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
  }
  plot_df[plot_df$LS < LS_window | abs(plot_df$Correlation) < corr_window,]$toLabel <- "NO"
  res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
  # Save dataframe of results
  save_df <- plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes),c(1:5)]
  save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
  saveRDS(save_df, paste0(regPlotDir, sprintf("/%s.regulatory_targets_%s.rds", strain, motif_short)))
  plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
  # Label motif as well
  plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
  plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
  # Threshold pvalue for plotting
  maxPval <- 5
  plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
  #Plot results
  p <- (
    ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
    #+ geom_point(size = 2)
    + ggrastr::geom_point_rast(size=2)
    + ggrepel::geom_text_repel(
      data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
      #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
      size=4,
      point.padding=0, # additional pading around each point
      box.padding=0.5,
      min.segment.length=0, # draw all line segments
      max.overlaps=Inf, # draw all labels
      #nudge_x = 2,
      color="black"
    ) 
    + geom_vline(xintercept=0, lty="dashed") 
    + geom_vline(xintercept=corr_window, lty="dashed", color="red")
    + geom_vline(xintercept=-corr_window, lty="dashed", color="red")
    + geom_hline(yintercept=LS_window, lty="dashed", color="red")
    + theme_classic()
    + theme(panel.grid.major=element_blank(), 
            panel.grid.minor= element_blank(), 
            plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
            aspect.ratio=1.0,
            #legend.position = "none", # Remove legend
            axis.text.x = element_text(angle=90, hjust=1))
    + ylab("Linkage Score") 
    + xlab("Motif Correlation to Gene") 
    + scale_color_gradientn(colors=wesanderson::wes_palette("Zissou1", type = "continuous"), limits=c(0, maxPval), name = "-log10(p-val)")
    + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
    + scale_x_continuous(limits = c(-0.85, 0.955)) # Force plot limits
    + ggtitle(sprintf("%s putative targets", motif_short))
  )
  
  pdf(paste0(regPlotDir, sprintf("/%s.%s_LS.pdf", strain, motif_short)), width=8, height=6)
  print(p)
  dev.off()
  
  outfile = paste0(regPlotDir, sprintf("/%s.%s_LS.png", strain, motif_short))
  print(p)
  ggsave(outfile, width=678/96, height=420/96, dpi=300)
  
  print(c(motif, length(unique(plot_df[plot_df$toLabel=="YES",]$symbol)), 
          length(unique(plot_df[plot_df$toLabel=="YES" & plot_df$symbol %in% markerGenes,]$symbol))))
}


