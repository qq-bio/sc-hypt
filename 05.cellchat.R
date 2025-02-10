dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(CellChat)
library(Seurat)
# library(ggsci)
# library(ggupset)
# library(VennDiagram)
# library(ggvenn)
# library(ggtern)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))



setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/")

################################################################################
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

CellChatDB = CellChatDB.mouse
future::plan("multicore", workers = 4) 


e <- new.env()
for(i in input_file){
  
  tissue = gsub("\\.(RNA|multiomics)\\.anno\\.v2\\.rds", "", basename(i))
  
  seurat_object = readRDS(i)
  Idents(seurat_object) <- "subclass_level2"
  seurat_object@active.ident = factor(seurat_object@meta.data[, "subclass_level2"])
  Idents(seurat_object) = "subclass_level2"
  
  strain_list = unique(seurat_object$strain)
  for(si in strain_list){
    
    treatment_list = unique(seurat_object@meta.data[seurat_object$strain==si, ]$treatment)
    
    for(ti in treatment_list){
      
      seurat_object_tmp = subset(seurat_object, strain==si & treatment==ti)
      cellchat = createCellChat(seurat_object_tmp, group.by = "subclass_level2")
      cellchat@DB = CellChatDB
      cellchat = subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      
      output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
      
      with(e, {
        assign(output_name, cellchat)
      })
      
    }
    
  }
  
}



saveRDS(e, "cellchat.rds")


