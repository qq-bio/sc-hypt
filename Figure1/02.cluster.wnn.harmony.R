
library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(clustree)


################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")

angii_infile <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/mouse.LK.multiomics.merged.rds"
ss_infile <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.LK.multiomics.merged.rds"
sp_infile <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.LK.multiomics.merged.rds"

infile <- sp_infile
seurat_object <- readRDS(infile)

# Load data, perform quality control, and filter cells based on features
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- subset(seurat_object, nFeature_RNA > 200 & 
                          nCount_RNA < 20000 & 
                          doublet_pc.2 == "Singlet" & 
                          nCount_ATAC < 100000 & 
                          nCount_ATAC > 400 & 
                          TSSEnrichment > 2)

# Cell cycle scoring
s.genes <- intersect(c(str_to_title(cc.genes$s.genes), "Cenpu"), rownames(seurat_object))
g2m.genes <- intersect(c(str_to_title(cc.genes$g2m.genes), "Pimreg", "Jpt1"), rownames(seurat_object))
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score

# Process RNA data
DefaultAssay(seurat_object) <- "RNA"
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object)
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "pca", 
                            reduction.save = "harmony.rna", assay.use = "RNA", project.dim = FALSE)

# Process ATAC data
DefaultAssay(seurat_object) <- "ATAC"
seurat_object <- RunTFIDF(seurat_object)
seurat_object <- FindTopFeatures(seurat_object, min.cutoff = 'q0')
seurat_object <- RunSVD(seurat_object)
seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "lsi", 
                            reduction.save = "harmony.atac", assay.use = "ATAC", project.dim = FALSE)

# Integration and clustering
seurat_object <- FindMultiModalNeighbors(seurat_object, reduction.list = list("harmony.rna", "harmony.atac"), 
                                         dims.list = list(1:30, 2:30))
seurat_object <- RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", 
                         reduction.key = "wnnUMAPHarmony_")
seurat_object <- FindClusters(seurat_object, graph.name = "wsnn", algorithm = 3, 
                              resolution = seq(0.1, 1, 0.1))

seurat_meta <- seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

reso_para <- paste0('wsnn_res.', 0.4)
seurat_object@active.ident <- seurat_object@meta.data[, reso_para]
seurat_object$seurat_clusters <- seurat_object@meta.data[, reso_para]

saveRDS(seurat_object, gsub(".merged.rds", ".cluster.rds", infile))





