dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(SeuratDisk)
library(DropletUtils)

input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds"
seurat_object = readRDS(input_file)
seurat_object$ATAC <- NULL
i <- sapply(seurat_object@meta.data, is.factor)
seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)
SaveH5Seurat(seurat_object, filename = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.h5Seurat", overwrite = T)
Convert("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.h5Seurat", dest = "h5ad", overwrite = T)

write10xCounts(x = seurat_object@assays$RNA@counts, path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.count.10x.h5",
               gene.symbol=rownames(seurat_object), barcodes=colnames(seurat_object), version='3', overwrite = T)


input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds"
seurat_object = readRDS(input_file)
seurat_object$ATAC <- NULL
# i <- sapply(seurat_object@meta.data, is.factor)
# seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)
# SaveH5Seurat(seurat_object, filename = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.h5Seurat", overwrite = T)
# Convert("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.h5Seurat", dest = "h5ad", overwrite = T)

write10xCounts(x = seurat_object@assays$RNA@counts, path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.count.10x.h5",
               gene.symbol=rownames(seurat_object), barcodes=colnames(seurat_object), version='3', overwrite = T)


input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds"
seurat_object = readRDS(input_file)
# seurat_object$ATAC <- NULL
# i <- sapply(seurat_object@meta.data, is.factor)
# seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)
# SaveH5Seurat(seurat_object, filename = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.h5Seurat", overwrite = T)
# Convert("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.h5Seurat", dest = "h5ad", overwrite = T)

write10xCounts(x = seurat_object@assays$RNA@counts, path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.count.10x.h5",
               gene.symbol=rownames(seurat_object), barcodes=colnames(seurat_object), version='3', overwrite = T)



input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds"
seurat_object = readRDS(input_file)
i <- sapply(seurat_object@meta.data, is.factor)
seurat_object@meta.data[i] <- lapply(seurat_object@meta.data[i], as.character)
SaveH5Seurat(seurat_object, filename = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.multiomics.anno.h5Seurat")
Convert("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.multiomics.anno.h5Seurat", dest = "h5ad")

write10xCounts(x = seurat_object@assays$RNA@counts, path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.multiomics.anno.count.10x.h5",
               gene.symbol=rownames(seurat_object), barcodes=colnames(seurat_object), version='3')

