library(Seurat)
library(harmony)
library(stringr)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")


merge_seurat_files <- function(files) {
  obj_list <- lapply(files, function(f) {
    so <- readRDS(f)
    sample_id <- str_remove(basename(f), "_doubletfinder\\.rds$")
    so$sample <- sample_id
    so <- RenameCells(so, new.names = paste0(sample_id, "_", colnames(so)))
    return(so)
  })
  
  merged <- Reduce(function(x, y) merge(x, y), obj_list)
  merged <- JoinLayers(merged)
  return(merged)
}


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder/")
so_list <- c("MLV1SN_doubletfinder.rds", "MLV2SN_doubletfinder.rds", "MLV3SN_doubletfinder.rds",
             "MLV4SN_doubletfinder.rds", "MLV5SN_doubletfinder.rds", "MLV6SN_doubletfinder.rds")

so_merged <- merge_seurat_files(so_list)
so_merged <- QC_harmony(so_merged)
saveRDS(so_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.cluster.v2.rds")


so_list <- c("RLV1SN_doubletfinder.rds", "RLV2SN_doubletfinder.rds", "RLV3SN_doubletfinder.rds", "RLV4SN_doubletfinder.rds", 
             "RLVW1SN_doubletfinder.rds", "RLVW2SN_doubletfinder.rds", "RLVW3SN_doubletfinder.rds", "RLVW4SN_doubletfinder.rds")
so_merged <- merge_seurat_files(so_list)
so_merged <- QC_harmony(so_merged)
saveRDS(so_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.cluster.v2.rds")


