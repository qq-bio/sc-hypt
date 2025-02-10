#!/usr/bin/env Rscript
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

setwd("/scratch/u/qqiu/project/multiomics/stereo_seq/result/")

file_list = list.files("./","h5ad")



# library(SeuratDisk)
# 
# Convert("952BR_E4.bin50.h5ad", dest = "h5seurat", overwrite = TRUE)
# 
# seuratObject <- ReadH5AD("849BL_A4.bin50.h5ad")
# 
# h5ls("849BL_A4.bin50.h5ad")
# bin = h5read("849BL_A4.bin50.h5ad", "bpMatrix_1", bit64conversion='bit64')
# dim(bin)
# [1]     1 26462 26462



seurat_object = readRDS("952BR_E4.bin50.rds")
seurat_object = readRDS("952BR_B4.bin50.rds")
SpatialFeaturePlot(seurat_object, features = "Slc12a1")
SpatialFeaturePlot(seurat_object, features = "Nphs1")


for(i in file_list){
  
  Load10X_Spatial()
  plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
  plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
  wrap_plots(plot1, plot2)
  
  SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
  
}
