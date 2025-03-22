library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(sceasy)


seurat = sceasy::convertFormat("/xdisk/mliang1/qqiu/data/HeartCellAtlas/visium-OCT_LV_lognormalised.updated.h5ad", from = "anndata", to = "seurat")

FeaturePlot(seurat, "IL1R1", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "EC6_ven", reduction = "spatial", split.by = "sangerID")

VlnPlot(seurat, "IL1R1", group.by = "spatial", split.by = "sangerID")


ec_list = colnames(seurat@meta.data)[grepl("EC", colnames(seurat@meta.data))]
for(i in ec_list){
  tryCatch({
    result <- cor.test(seurat@meta.data[, i], seurat@assays$RNA@data["IL1R1", ])
    print(c(i, result$estimate, result$p.value))
  }, error = function(e) {
  })
}



