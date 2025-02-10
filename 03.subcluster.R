source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DE_enrichment.R")

library(Seurat)
# library(SeuratDisk)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster")





input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

project_list = c("AngII", "Salt-sensitive", "Spontaneous"); names(project_list)=c("mouse", "rat.ss", "rat.sp")

for( i in input_file){
  
  seurat_object = readRDS(i)
  
  dataset = gsub("\\.(RNA|multiomics)+.*rds", "", basename(i), perl = T)
  project = gsub("\\.[A-Z]+", "", dataset)
  tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
  # title = paste0(project_list[project], "-", tissue, "(N=",
  #                prettyNum(ncol(seurat_object), big.mark = ',', scientific=F), ")")
  
  if(tissue=="LK"){
    reduction = "wnn.umap.harmony"
  }else{
    reduction = "umap"
  }
  if(tissue=="HYP"){
    seurat_object$subclass_level2 = seurat_object$subclass_level1
  }
  
  seurat_object_tmp = subset(seurat_object, subclass_level2 %in% c("Astrocytes"))
  seurat_object_tmp = QC_harmony(seurat_object_tmp)
  
  DE_analysis(seurat_object_tmp, "RNA_snn_res.0.1", "mouse.HYP.Astrocytes")
  
  tmp = table(seurat_object_tmp$RNA_snn_res.0.1, seurat_object_tmp$treatment)
  t(tmp)/colSums(tmp)
  
  
  Idents(seurat_object_tmp) = "RNA_snn_res.0.1"
  UMAP_plot1 = DimPlot(seurat_object_tmp, label = F, reduction = reduction) + theme(legend.position = "none")
  UMAP_plot2 = DimPlot(seurat_object_tmp, label = T, reduction = reduction) + theme(legend.position = "none")
  UMAP_plot2
  
  
  DA_plot(seurat_object_tmp@meta.data, outfile="rat.HYP.major_cluster.proportion_change.pdf", project = "project")
  
  enrichment_analysis(seurat_object_mouse, cell="Astrocytes", outfile="mouse.HYP.astrocytes")
  
  
  
}










infile="../QC/mouse.HYP.doubletfinder.merged.rds"
seurat_object = readRDS(infile)

seurat_object = QC_harmony(seurat_object)

enrichment_analysis(seurat_object_mouse, cell="Astrocytes", outfile="mouse/mouse.HYP.astrocytes")
