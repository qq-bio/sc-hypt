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

dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")
dyn.load("/usr/lib64/atlas/libsatlas.so.3")

library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)

# addArchRChrPrefix(chrPrefix = FALSE)
addArchRThreads(threads = 1) 


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
#     library(BSgenome.Rnorvegicus.UCSC.rn7)
#     load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
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
#   proj <- addPeak2GeneLinks(
#     ArchRProj = proj,
#     reducedDims = "LSI_Combined" # IterativeLSI
#   )
#   
#   saveArchRProject(ArchRProj = proj, outputDirectory = outputDirectory, load = FALSE)
#   
#   saveRDS(proj, outfile)
#   
# }




################################################################################
### load dataset
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"

proj = readRDS(outfile)

### dimensionality reduction
# plotEmbedding(ArchRProj = proj, name = "subclass_level2", embedding = "wnn.umap.harmony")


################################################################################
### differential peaks (treatment vs control)
proj$cell_grp = paste0(proj$subclass_level1, "_", proj$species, "_", proj$treatment)
cell_grp_count = table(proj$cell_grp)
cell_list = unique(proj$subclass_level1)
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

DEP_outfile = gsub("rds", "DEP.out", outfile)
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


### peak to gene linkage at cell type specific level
cell_type_list = unique(proj@cellColData$subclass_level1)
for( i in cell_type_list){

  cell_list = rownames(proj@cellColData[proj@cellColData$subclass_level1==i, ])
  proj_tmp = subsetArchRProject(proj, cells = cell_list, outputDirectory = "ArchRSubset_test")
  # proj_tmp = addPeakMatrix(proj_tmp)
  proj_tmp <- addPeak2GeneLinks(
      ArchRProj = proj_tmp,
      reducedDims = "LSI_Combined" # IterativeLSI
    )

  p2g <- getPeak2GeneLinks(
    ArchRProj = proj_tmp,
    corCutOff = 0.45,
    resolution = 10000,
    returnLoops = FALSE
  )

  p2geneDF = metadata(proj_tmp@peakSet)$Peak2GeneLinks
  p2geneDF$gene_symbol = mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
  p2geneDF$peakName = (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
  p2geneDF = as.data.frame(p2geneDF)
  p2geneDF$cell_type = i

}

P2G_outfile = gsub("rds", "p2g.out", outfile)
write.table(p2geneDF, P2G_outfile, col.names = T, row.names = F, sep='\t', quote = F)





################################################################################
### motif enrichment for DEP

motif_merge = c()
cell_type_list = unique(proj@cellColData$subclass_level1)
for( i in cell_type_list){
  
  selected_peaks = p2geneDF[! is.na(p2geneDF$FDR) & p2geneDF$FDR<0.05 & p2geneDF$cell_type==i, ]$idxATAC %>% unique()
  
  enrichRegions <- peakAnnoEnrichment(
    seMarker = getPeakSet(proj)[selected_peaks,],
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05"
  )
  
}









################################################################################
### motif enrichment for peaks linked to DEG
DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header=T)
# DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", header=T)
# DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.2 & DEG$treatment=="SS-LS" & DEG$cell_type=="EC", ]

strain = gsub(".LK.archr.rds", "", basename(outfile))
strain_list = c("mouse"="C57BL/6", "rat.ss"="SS", "rat.sp"="SHR")
strain_use = strain_list[strain]

marker_merge = c()
cell_type_list = unique(proj@cellColData$subclass_level1)
for( i in cell_type_list){
  
  DEG_use = DEG[DEG$p_val_adj<0.05 & abs(DEG$avg_log2FC)>0.25 & DEG$tissue=="LK" & DEG$strain==strain_use & DEG$cell_type==i, ]
  selected_peaks = p2geneDF[! is.na(p2geneDF$FDR) & p2geneDF$FDR<0.05 & p2geneDF$cell_type==i & p2geneDF$gene_symbol %in% DEG_use$gene_name, ]$idxATAC %>% unique()
  
  customEnrichment( getPeakSet(proj)[selected_peaks,], getMatches(proj), bgdPeaks = getBgdPeaks(proj, force=TRUE) )
  
  enrichRegions <- peakAnnoEnrichment(
    seMarker = getPeakSet(proj)[selected_peaks,],
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05"
  )
  
}



#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_2_tracks.R
markerGenes = DEG_use$gene_name
markerGenes = c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
                "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd")
p <- plotBrowserTrack(
  ArchRProj = proj,
  groupBy = "cell_grp",
  useGroups = c("EC_SS_LS","EC_SD_LS","EC_SD_HS 3d","EC_SS_HS 3d","EC_SS_HS 21d"),
  geneSymbol = markerGenes,
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(proj),
  sizes = c(7, 0.2, 1, 1)
)

grid::grid.newpage()
grid::grid.draw(p$Lrp2)

p2geneDF[p2geneDF$gene_symbol=="Lrp2", ]

saveRDS(p, "rat.ss.LK.strain_wise.EC_DEG.plotTrack.rds")













label_genes = intersect(DEG_use$gene_name, unique(p2geneDF[p2geneDF$FDR<0.05, ]$gene_symbol))

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





### motif activity
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
# 
# proj@cellColData$subclass_level2 = as.character(proj@cellColData$subclass_level2)
# 
# proj$subclass_level2.1 = sub(" ", ".", proj$subclass_level2)
# markersPeaks <- getMarkerFeatures(
#   ArchRProj = proj,
#   useMatrix = "PeakMatrix",
#   groupBy = "subclass_level2.1",
#   bias = c("TSSEnrichment", "log10(nFrags)"),
#   testMethod = "wilcoxon"
# )
# 
# p2g <- getPeak2GeneLinks(
#   ArchRProj = proj,
#   corCutOff = 0.45,
#   resolution = 10000,
#   returnLoops = FALSE
# )
# 
# p2geneDF = metadata(proj@peakSet)$Peak2GeneLinks
# p2geneDF$gene_symbol = mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
# p2geneDF$peakName = (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2geneDF$idxATAC]
# p2geneDF = as.data.frame(p2geneDF)

# p <- plotPeak2GeneHeatmap(proj, groupBy="subclass_level2")


DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.multiomics.DEG.out", header=T)
# peaks = getPeakSet(proj)
# cor_peaks = peaks[na.omit(unique(p2geneDF[(p2geneDF$FDR<0.05) & (p2geneDF$gene_symbol %in% DEG$gene_name), ]$idxATAC)), ]
# cor_peaks = SummarizedExperiment(cor_peaks)

# peaks = getMatrixFromProject(proj, useMatrix = "PeakMatrix")

proj$subclass_level2.1 = sub(" ", ".", proj$subclass_level2)
markersPeaks <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "subclass_level2.1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

cor_peaks = subset(markersPeaks, rownames(markersPeaks) %in% unique(p2geneDF[(p2geneDF$FDR<0.05) & (p2geneDF$gene_symbol %in% DEG$gene_name), ]$idxATAC))
# metadata(cor_peaks)$Params$useMatrix = "PeakMatrix"

# motifsUp <- peakAnnoEnrichment(
#   seMarker = cor_peaks,
#   ArchRProj = proj,
#   peakAnnotation = "Motif",
#   cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
# )
#

matches = getMatches(proj, "Motif")
idx = na.omit(unique(p2geneDF[(p2geneDF$FDR<0.05) & (p2geneDF$gene_symbol %in% DEG$gene_name), ]$idxATAC))

cor_tf = computeEnrichment(matches, idx, seq_len(nrow(matches)))

df <- data.frame(TF = rownames(cor_tf), mlog10Padj = cor_tf$mlog10p)
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
    size = 4,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() +
  ylab("-log10(P-adj) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))



customEnrichment=function (matches = NULL, compare = NULL, background = NULL)
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

### https://github.com/GreenleafLab/ArchR/discussions/880
### https://github.com/GreenleafLab/ArchR/blob/c44b50c56323d5995fb3d2fb28a9072255fec644/R/AnnotationPeaks.R#L1053-L1146
customEnrichment( ranges = peaks1, matches = getMatches(proj) )






### motif enrichment in differential peaks
# proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif")
proj <- addBgdPeaks(proj)

proj <- addDeviationsMatrix(
  ArchRProj = proj,
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(proj, name = "MotifMatrix", plot = TRUE)
plotVarDev

motifs <- rownames(cor_tf)[1:10]
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs


p <- plotEmbedding(
  ArchRProj = proj,
  colorBy = "MotifMatrix",
  name = sort(markerMotifs),
  embedding = "wnn.umap.harmony",
  imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))


dep_tf = as.character(lapply(strsplit(rownames(cor_tf[cor_tf$mlog10Padj>2, ]), "_"), function(x) x[[1]][1]))
DEG[DEG$gene_name %in% dep_tf, ]

################################################################################
# seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "subclass_level2.1")
seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "cell_grp")

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
### gene acticity score

markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8",
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

heatmapGS <- markerHeatmap()


