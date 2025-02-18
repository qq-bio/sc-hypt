.libPaths("/xdisk/mliang1/qqiu/R/library_Rv4.2_puma9/", include.site = F)

library(scrat)
# library(Seurat)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scrat/")
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")

so_tmp <- FindVariableFeatures(seurat_object, selection.method = "vst")



env <- scrat.new(list(dataset.name="Cross-organ EC"))
# recommend to use expression values in logarithmic scale
env$indata <- seurat_object@assays$RNA@data
env$group.labels <- seurat_object$seurat_clusters

scrat.run(env)
saveRDS(env, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/scrat/ec.scvi.scvi.rds")
