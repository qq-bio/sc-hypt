
library(Seurat)
library(harmony)
library(clustree)




################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")

input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/mouse.HYP.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/mouse.LV.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.HYP.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.LV.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.MSA.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.HYP.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.LV.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.MSA.RNA.merged.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.MCA.RNA.merged.rds"
)


for(i in input_file){
  
  seurat_object <- readRDS(i)
  
  seurat_object <- subset(seurat_object, nFeature_RNA>200 &
                           nCount_RNA<20000 &
                           percent.mt<20 &
                           doublet_pc.2=="Singlet")
  
  # Cell cycle scoring
  s.genes <- intersect(c(str_to_title(cc.genes$s.genes), "Cenpu"), rownames(seurat_object))
  g2m.genes <- intersect(c(str_to_title(cc.genes$g2m.genes), "Pimreg", "Jpt1"), rownames(seurat_object))
  seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
  
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))
  seurat_object <- RunHarmony(seurat_object, group.by.vars = "orig.ident")
  
  seurat_object <- RunUMAP(seurat_object, reduction = "harmony", dims = 1:30)
  seurat_object <- FindNeighbors(seurat_object, reduction = "harmony", dims = 1:30)
  seurat_object <- FindClusters(seurat_object, resolution = seq(0.1, 1, 0.1))
  
  seurat_meta <- seurat_object@meta.data
  clustree(seurat_meta, prefix = "RNA_snn_res.")
  
  reso_para <- paste0('wsnn_res.', 0.4)
  seurat_object@active.ident <- seurat_object@meta.data[, reso_para]
  seurat_object$seurat_clusters <- seurat_object@meta.data[, reso_para]
  
  saveRDS(seurat_object, gsub(".merged.rds", ".cluster.rds", i))
  
}


