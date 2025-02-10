library(Matrix)
library(Seurat)
library(sctransform)
library(harmony)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(R.utils)
library(ggsci)
library(presto)
library(clustree)
library(ROGUE)
library(UpSetR)
library(SingleR)
library(celldex)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/00.initial_setting.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/panglao_anno.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/cluster_sample_sum.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")

input_file = c(
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds",

"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds",

"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds",
"/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds"
)


anno_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/major_cluster.anno.txt", header=T, sep="\t", fill=T)

sample_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_sample_info.txt", header=T, sep='\t')
rownames(sample_info) = sample_info$seqID

project_list = c("AngII", "Salt-sensitive", "Spontaneous"); names(project_list)=c("mouse", "rat.ss", "rat.sp")

for(i in input_file){
  
  dataset = gsub("\\.(RNA|multiomics)+.cluster.rds", "", basename(i), perl = T)
  project = gsub("\\.[A-Z]+", "", dataset)
  outfile = gsub("cluster.rds", "anno.rds", i)
  
  seurat_object = readRDS(i)
  
  seurat_object$cell_id = colnames(seurat_object)
  seurat_object$seqID2 = sample_info[as.character(seurat_object$orig.ident),]$seqID2
  seurat_object$sample = sample_info[as.character(seurat_object$orig.ident),]$sample
  seurat_object$strain = factor(sample_info[as.character(seurat_object$orig.ident),]$strain,
                                 levels=c("C57BL/6", "WKY", "SHR", "SD", "SS"))
  seurat_object$tissue = factor(sample_info[as.character(seurat_object$orig.ident),]$tissue.abbr,
                                levels=c("HYP", "MCA", "LV", "LK", "MSA"))
  seurat_object$treatment = factor(sample_info[as.character(seurat_object$orig.ident),]$treatment,
                                   levels=c("Saline 3d", "AngII 3d", "AngII 28d",
                                        "10w", "26w", "LS", "HS 3d", "HS 21d"))
  seurat_object$project = factor(project_list[project], levels=c("AngII", "Salt-sensitive", "Spontaneous"))
  
  
  if(grepl("HYP", dataset)){
    HYP_levels = c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
                   "Astrocyte", "Microglia", "Activated microglia", 
                   "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
                   "Tanycyte", "Ependymal cell", "Pars tuberalis cell", 
                   "EC", "E/P transition cell", "Pericyte", "Fibroblast")
    merged_levels = HYP_levels
  }else if(grepl("LV", dataset)){
    LV_levels = c("CM", "EC", "Pericyte", "Fibroblast", "IMM")
    merged_levels = LV_levels
  }else if(grepl("LK", dataset)){
    LK_levels = c("POD", "PT", "TL", "TAL", "DCT", "DCT/CT", "CT", "CD", "IC", "EC", "VSMC", "Fibroblast", "IMM")
    merged_levels = LK_levels
  }else if(grepl("MCA", dataset)){
    MCA_levels = c("Microglia", #"Neuron", "Myelinating OL", "OPC", "Astrocyte"
                   "EC", "VSMC", "Pericyte", "Fibroblast", "IMM")
    merged_levels = MCA_levels
  }else if(grepl("MSA", dataset)){
    MSA_levels = c("EC", "VSMC", "Fibroblast", "IMM", "Adipocyte")
    merged_levels = MSA_levels
  }
  
  class_levels = c("neurons", "glial cells", "muscle cells", "epithelial cells", "endothelial cells", 
                   "stromal cells", "immune cells", "adipocytes", "endocrine cells")
  
  anno_info_use = anno_info[anno_info$dataset == dataset, ]
  reso = unique(anno_info_use$resolution)
  cluster_use = anno_info_use[!(grepl("remove", anno_info_use$note)),]$cluster
  subclass_level1 = anno_info_use[!(grepl("remove", anno_info_use$note)),]$subclass_level1
  class = anno_info_use[!(grepl("remove", anno_info_use$note)),]$class
  
  Idents(seurat_object) = reso
  seurat_object$seurat_clusters = seurat_object@meta.data[, reso]
  seurat_object = subset(seurat_object, seurat_clusters %in% cluster_use)
  seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)
  
  seurat_object@meta.data$"subclass_level1" = factor(subclass_level1[match(seurat_object$seurat_clusters, cluster_use)],
                                                          levels=merged_levels)
  seurat_object = subset(seurat_object, subclass_level1 %in% merged_levels)
  
  seurat_object@meta.data$"class" = factor(class[match(seurat_object$seurat_clusters, cluster_use)],
                                                     levels=class_levels)
  
  saveRDS(seurat_object, outfile)
  print(c(dataset, all(unique(subclass_level1) %in% merged_levels)))
}




# 
# ################################################################################
# ### immune cell annotation using singleR
# ################################################################################
# # QC & harmony - mouse
# LK_infile="mouse.LK.multiomics.anno.rds"
# MCA_infile="mouse.MCA.RNA.anno.rds"
# LV_infile="mouse.LV.RNA.anno.rds"
# 
# LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
# MCA_object = readRDS(MCA_infile)
# LV_object = readRDS(LV_infile)
# 
# DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.5")
# FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
# DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.5")
# FeaturePlot(MCA_object, features = "Ptprc")
# DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
# FeaturePlot(LV_object, features = "Ptprc")
# 
# LK_object@active.ident <- factor(LK_object@active.ident)
# names(LK_object@active.ident) <- colnames(LK_object)
# LK_immune = subset(LK_object, wsnn_res.0.5 %in% c(12, 17, 18))
# LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(2))
# MCA_object@active.ident <- factor(MCA_object@active.ident)
# names(MCA_object@active.ident) <- colnames(MCA_object)
# MCA_immune = subset(MCA_object, RNA_snn_res.0.5 %in% c(18))
# 
# ### prepare ref data
# ref_object = merge(MCA_immune, y = LV_immune,
#                    project = "mouse.immune"
# )
# ref_object = NormalizeData(ref_object)
# ref_object = FindVariableFeatures(ref_object)
# ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
# ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
# ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
# ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
# ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)
# 
# ### prepare query data
# query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
#                                   min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
# query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
# query_object = NormalizeData(query_object)
# 
# transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
#                                        k.filter = NA, reference.reduction = 'pca', dims = 1:30)
# 
# query_object = MapQuery(anchorset = transfer_anchors,
#                         query = query_object, reference = ref_object,
#                         reference.reduction = "pca",
#                         reduction.model = "umap")
# 
# merged_object = merge(ref_object, query_object)
# merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
# merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
# merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
# merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
# merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))
# 
# saveRDS(merged_object, "mouse.immune_cell.cluster.rds")
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# # QC & harmony - rat.ss
# LK_infile="rat.ss.LK.multiomics.anno.rds"
# MCA_infile="rat.ss.MCA.RNA.anno.rds"
# MSA_infile="rat.ss.MSA.RNA.anno.rds"
# LV_infile="rat.ss.LV.RNA.anno.rds"
# 
# LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
# MCA_object = readRDS(MCA_infile)
# MSA_object = readRDS(MSA_infile)
# LV_object = readRDS(LV_infile)
# 
# DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.4")
# FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
# DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.4")
# FeaturePlot(MCA_object, features = "Ptprc")
# DimPlot(MSA_object, label = T, group.by = "RNA_snn_res.0.4")
# FeaturePlot(MSA_object, features = "Ptprc")
# DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
# FeaturePlot(LV_object, features = "Ptprc")
# 
# LK_object@active.ident <- factor(LK_object@active.ident)
# names(LK_object@active.ident) <- colnames(LK_object)
# LK_immune = subset(LK_object, wsnn_res.0.4 %in% c(4, 13, 15))
# LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(4, 5))
# MCA_immune = subset(MCA_object, RNA_snn_res.0.4 %in% c(2))
# MSA_immune = subset(MSA_object, RNA_snn_res.0.4 %in% c(7))
# 
# ### prepare ref data
# ref_object = merge(MCA_immune, y = list(LV_immune, MSA_immune),
#                    project = "rat.ss.immune"
# )
# ref_object = NormalizeData(ref_object)
# ref_object = FindVariableFeatures(ref_object)
# ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
# ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
# ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
# ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
# ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)
# 
# 
# ### prepare query data
# query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
#                                   min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
# query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
# query_object = NormalizeData(query_object)
# 
# transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
#                                        k.filter = NA, reference.reduction = 'pca', dims = 1:30)
# 
# query_object = MapQuery(anchorset = transfer_anchors,
#                         query = query_object, reference = ref_object,
#                         reference.reduction = "pca",
#                         reduction.model = "umap")
# 
# merged_object = merge(ref_object, query_object)
# merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
# merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
# merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
# merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
# merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))
# 
# saveRDS(merged_object, "rat.ss.immune_cell.cluster.rds")
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# LK_infile="rat.sp.LK.multiomics.anno.rds"
# MCA_infile="rat.sp.MCA.RNA.anno.rds"
# MSA_infile="rat.sp.MSA.RNA.anno.rds"
# LV_infile="rat.sp.LV.RNA.anno.rds"
# 
# LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
# MCA_object = readRDS(MCA_infile)
# MSA_object = readRDS(MSA_infile)
# LV_object = readRDS(LV_infile)
# 
# DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.4")
# FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
# DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.5")
# FeaturePlot(MCA_object, features = "Ptprc")
# DimPlot(MSA_object, label = T, group.by = "RNA_snn_res.0.4")
# FeaturePlot(MSA_object, features = "Ptprc")
# DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
# FeaturePlot(LV_object, features = "Ptprc")
# 
# LK_object@active.ident <- factor(LK_object@active.ident)
# names(LK_object@active.ident) <- colnames(LK_object)
# LK_immune = subset(LK_object, wsnn_res.0.4 %in% c(12))
# LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(4))
# MCA_immune = subset(MCA_object, RNA_snn_res.0.5 %in% c(5, 18))
# MSA_immune = subset(MSA_object, RNA_snn_res.0.4 %in% c(7))
# 
# ### prepare ref data
# ref_object = merge(MCA_immune, y = list(LV_immune, MSA_immune),
#                    project = "rat.sp.immune"
# )
# ref_object = NormalizeData(ref_object)
# ref_object = FindVariableFeatures(ref_object)
# ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
# ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
# ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
# ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
# ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)
# 
# 
# ### prepare query data
# query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
#                                   min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
# query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
# query_object = NormalizeData(query_object)
# 
# transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
#                                        k.filter = NA, reference.reduction = 'pca', dims = 1:30)
# 
# query_object = MapQuery(anchorset = transfer_anchors,
#                         query = query_object, reference = ref_object,
#                         reference.reduction = "pca",
#                         reduction.model = "umap")
# 
# merged_object = merge(ref_object, query_object)
# merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
# merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
# merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
# merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
# merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))
# 
# saveRDS(merged_object, "rat.sp.immune_cell.cluster.rds")
# 
# 
# 
# 
# 
# ################################################################################
# ### immune cell annotation using singleR
# 
# input_file = c("mouse.immune_cell.cluster.rds",
#                "rat.ss.immune_cell.cluster.rds",
#                "rat.sp.immune_cell.cluster.rds")
# 
# mimd.sc <- celldex::ImmGenData()
# # merged_levels = c("Microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
# #                                       "NK cells", "NKT", "T cells", "B cells")
# 
# i = input_file[1]
# seurat_object = readRDS(i)
# outfile = gsub("cluster.rds", "anno.rds", i)
# DimPlot(seurat_object, label = T, group.by = c("RNA_snn_res.1"))
# seurat_object$seurat_clusters = seurat_object$RNA_snn_res.1
# main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
# main.group$labels = c(main.group$labels)
# seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
# DimPlot(seurat_object, label = T, group.by = "subclass_level1")
# FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
# saveRDS(seurat_object, outfile)
# 
# 
# i = input_file[2]
# seurat_object = readRDS(i)
# outfile = gsub("cluster.rds", "anno.rds", i)
# DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
# seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
# main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
# seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
# DimPlot(seurat_object, label = T, group.by = "subclass_level1")
# FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
# seurat_object@meta.data[seurat_object$subclass_level1=="Endothelial cells", ]$subclass_level1 = "Macrophages"
# saveRDS(seurat_object, outfile)
# 
# 
# i = input_file[3]
# seurat_object = readRDS(i)
# outfile = gsub("cluster.rds", "anno.rds", i)
# DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
# seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
# main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
# seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
# DimPlot(seurat_object, label = T, group.by = "subclass_level1")
# FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
# seurat_object@meta.data[seurat_object$subclass_level1=="Endothelial cells", ]$subclass_level1 = "Macrophages"
# saveRDS(seurat_object, outfile)



################################################################################
### map EC results back

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.cluster.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.cluster.rds"
)

anno_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/EC.anno.txt", header=T, sep='\t')
for(i in input_file){
  
  dataset = gsub(".cluster.rds", "", basename(i))
  EC_object = readRDS(i)
  outfile = gsub("cluster.rds", "anno.rds", i)
  
  anno_info_use = anno_info[anno_info$dataset == dataset, ]
  reso = unique(anno_info_use$resolution)
  cluster_use = anno_info_use$cluster
  subclass_level2 = anno_info_use$subclass_level2
  
  Idents(EC_object) = reso
  EC_object$seurat_clusters = EC_object@meta.data[, reso]
  
  EC_object@meta.data$"subclass_level2" = subclass_level2[match(EC_object$seurat_clusters, cluster_use)]
  
  saveRDS(EC_object, outfile)
  
  # seurat_file = gsub("subcluster", "cluster", i)
  # if(grepl("LK", i)){
  #   seurat_file = gsub("EC.cluster.rds", "multiomics.anno.rds", seurat_file)
  # }else{
  #   seurat_file = gsub("EC.cluster.rds", "RNA.anno.rds", seurat_file)
  # }
  # outfile = gsub("anno", "anno.v2", seurat_file)
  # 
  # seurat_object = readRDS(seurat_file)
  # seurat_object$subclass_level2 = as.character(seurat_object$subclass_level1)
  # seurat_object@meta.data[colnames(EC_object), ]$subclass_level2 = EC_object$subclass_level2
  # 
  # saveRDS(seurat_object, outfile)
}






# ################################################################################
# ### map immune cell results back
# ang_immune_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.immune_cell.anno.rds")
# ss_immune_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.immune_cell.anno.rds")
# sp_immune_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.immune_cell.anno.rds")
# subcluster_replace_list = c(as.character(ang_immune_object$subclass_level1), 
#                             as.character(ss_immune_object$subclass_level1), 
#                             as.character(sp_immune_object$subclass_level1))
# names(subcluster_replace_list) = c(colnames(ang_immune_object), colnames(ss_immune_object), colnames(sp_immune_object))
# 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
# )
# 
# 
# project_list = c("AngII", "Salt-sensitive", "Spontaneous"); names(project_list)=c("mouse", "rat.ss", "rat.sp")
# 
# for(i in input_file){
#   
#   # dataset = gsub("\\.[RNA|multiomics]+.anno.rds", "", basename(i), perl = T)
#   # project = gsub("\\.[A-Z]+", "", dataset)
#   
#   # outfile = gsub("anno", "anno.v2", i)
#   
#   seurat_object = readRDS(i)
#   
#   subclass_level1 = as.character(seurat_object$subclass_level1)
#   names(subclass_level1) = colnames(seurat_object)
#   replace_list = intersect(names(subclass_level1), names(subcluster_replace_list))
#   subclass_level1[replace_list] = subcluster_replace_list[replace_list]
#   seurat_object$subclass_level1 = factor(subclass_level1, levels = cell_order)
#   seurat_object = subset(seurat_object, subclass_level1 %in% cell_order)
#   
#   print(c(i, all(unique(subclass_level1) %in% cell_order)))
#   
#   saveRDS(seurat_object, i)
# }
# 
# 



################################################################################

blank_theme = theme(axis.line=element_blank(),axis.text.x=element_blank(),
                    axis.text.y=element_blank(),axis.ticks=element_blank(),
                    axis.title.x=element_blank(),
                    axis.title.y=element_blank(),legend.position="none",
                    panel.background=element_blank(),panel.border=element_rect(colour = "black"),panel.grid.major=element_blank(),
                    panel.grid.minor=element_blank(),plot.background=element_blank())

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

project_list = c("AngII", "Salt-sensitive", "Spontaneous"); names(project_list)=c("mouse", "rat.ss", "rat.sp")

pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.umap.v2.pdf", width=350/96, height=350/96)

for( i in input_file){
  
  filename = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/",
                    gsub("anno.rds", "umap.wo_label.png", basename(i)))
  
  seurat_object = readRDS(i)
  
  dataset = gsub("\\.(RNA|multiomics)+.*rds", "", basename(i), perl = T)
  project = gsub("\\.[A-Z]+", "", dataset)
  tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
  title = paste0(project_list[project], "-", tissue, "(N=",
                 prettyNum(ncol(seurat_object), big.mark = ',', scientific=F), ")")
  
  if(tissue=="LK"){
    reduction = "wnn.umap.harmony"
  }else{
    reduction = "umap"
  }
  if(tissue=="HYP"){
    seurat_object$subclass_level2 = seurat_object$subclass_level1
  }
  
  p = DimPlot(seurat_object, reduction = reduction, group.by = "subclass_level1", label = F, repel = F) + 
    blank_theme +
    scale_color_manual(values = cell_col) + 
    ggtitle(title)
  
  print(p)
  print(i)
  # ggsave(filename, width=450/96, height=450/96, dpi=300)
  
}

dev.off()


## individual legend
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_type.legend.png", width=2247, height=1896, res=300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(cell_col), pch=19, pt.cex=2, cex=1, bty='n',
       col = cell_col, ncol=2)
mtext("Cell type", at=0.1, cex=1.2)
dev.off()






meta_table = c()
selected_col = c("orig.ident", "seqID2", "sample", "strain", "tissue", "treatment",
                 "seurat_clusters", "cell_id", "project", "subclass_level1", "class")
for(i in input_file){

  # dataset = gsub("\\.(RNA|multiomics)+.*rds", "", basename(i), perl = T)
  # tmp = get(dataset)
  seurat_object = readRDS(i)
  meta_table = rbind(meta_table, seurat_object@meta.data[selected_col])
}

write.table(meta_table, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/multi-HYP.meta.out")




################################################################################

meta_table = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/multi-HYP.meta.out")
meta_table = meta_table[!(meta_table$tissue=="MCA" & meta_table$project %in% c("AngII", "Salt-sensitive")), ]

meta_table$tissue = factor(meta_table$tissue, levels = tissue_order)
meta_table$subclass_level1 = factor(meta_table$subclass_level1, levels = cell_order)
meta_table$strain = factor(meta_table$strain, levels = strain_order)
meta_table$treatment = factor(meta_table$treatment, levels = treatment_order)

cluster = "subclass_level1"
strain = "strain"
project = "project"
treatment = "treatment"
tissue = "tissue"

# meta_table[meta_table$subclass_level2=="Pars tuberalis cells", ]$class = "endocrine cells"
# meta_table[meta_table$subclass_level2=="ABC", ]$class = "endothelial cells"

dat = meta_table[, c(cluster, strain, project, treatment, tissue)]
pct = dat %>% group_by(strain, treatment, tissue) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  group_by_all() %>% 
  summarise(Count = n()) %>% ungroup() %>% 
  mutate(percentage = Count/n, id = row_number()) %>% as.data.frame()

colnames(pct)[1] = c("Cell_type")
pct$Cell_type = as.factor(pct$Cell_type)
# pct$treatment = as.factor(pct$treatment)

# pdf("stackedbarplot.pdf", width=4, height = 4)
pctBar = ggplot(data = pct, aes(y = percentage, x = treatment, fill = Cell_type, alluvium = Cell_type)) + 
  geom_bar(stat="identity", width = .7) +
  labs(title="", fill='',
       x='', y='Cell composition') +
  theme(legend.position = "None",
        plot.title = element_text(size=15, colour = 'black'),
        axis.text.y = element_text(size=10, colour = 'black'),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
        axis.title = element_text(size=15, colour = 'black'),
        strip.text = element_text(size=10)) +
  scale_x_discrete(limits = rev(levels(pct$Sample))) +
  geom_flow(alpha= 0.3, lty = 1, color = "white",
            curve_type = "linear", width = .7) +
  coord_flip() +
  scale_fill_manual(values = cell_col)

pctBar = pctBar + facet_nested(project + strain ~ tissue, scales = "free_y", space = "free_y")
print(pctBar)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.stacked_barplot.png", width=755/96, height=461/96, dpi=300)


cluster = "class"
strain = "strain"
project = "project"
treatment = "treatment"
tissue = "tissue"

# meta_table[meta_table$subclass_level2=="Pars tuberalis cells", ]$class = "endocrine cells"
# meta_table[meta_table$subclass_level2=="ABC", ]$class = "endothelial cells"

dat = meta_table[, c(cluster, strain, project, treatment, tissue)]
pct = dat %>% group_by(strain, treatment, tissue) %>% 
  mutate(n=n()) %>% 
  ungroup() %>% 
  group_by_all() %>% 
  summarise(Count = n()) %>% ungroup() %>% 
  mutate(percentage = Count/n, id = row_number()) %>% as.data.frame()

colnames(pct)[1] = c("Cell_type")
pct$Cell_type = as.factor(pct$Cell_type)
# pct$treatment = as.factor(pct$treatment)

# pdf("stackedbarplot.pdf", width=4, height = 4)
pctBar = ggplot(data = pct, aes(y = percentage, x = treatment, fill = Cell_type, alluvium = Cell_type)) + 
  geom_bar(stat="identity", width = .7) +
  labs(title="", fill='',
       x='', y='Cell composition') +
  theme(legend.position = "right",
        plot.title = element_text(size=15, colour = 'black'),
        axis.text.y = element_text(size=10, colour = 'black'),
        axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
        axis.title = element_text(size=15, colour = 'black'),
        strip.text = element_text(size=10)) +
  # scale_x_discrete(limits = rev(levels(pct$Sample))) + 
  geom_flow(alpha= 0.3, lty = 1, color = "white",
            curve_type = "linear", width = .7) +
  coord_flip() +
  scale_fill_manual(values = class_col)

pctBar = pctBar + facet_nested(project + strain ~ tissue, scales = "free_y", space = "free_y")
print(pctBar)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_class.stacked_barplot.png", width=947/96, height=461/96, dpi=300)


# 
# input_file = c("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.immune_cell.anno.rds",
#                "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.immune_cell.anno.rds",
#                "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.immune_cell.anno.rds")
# 
# for( i in input_file[3]){
#   
#   seurat_object = readRDS(i)
#   
#   dataset = gsub("\\.(RNA|multiomics)+.*rds", "", basename(i), perl = T)
#   project = gsub("\\.[A-Z]+", "", dataset)
#   tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
#   title = paste0(project_list[project], "-", tissue, "(N=",
#                  prettyNum(ncol(seurat_object), big.mark = ',', scientific=F), ")")
#   
#   if(tissue=="LK"){
#     reduction = "wnn.umap.harmony"
#   }else{
#     reduction = "umap"
#   }
#   if(tissue=="HYP"){
#     seurat_object$subclass_level2 = seurat_object$subclass_level1
#   }
#   
#   p = DimPlot(seurat_object, reduction = reduction, group.by = "subclass_level2", 
#               label = T, repel = T, split.by = "species") + 
#     blank_theme +
#     scale_color_manual(values = cell_col) + 
#     ggtitle("")
#   print(p) # width=358&height=354
# }
# 


input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

immune_list =   c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
                  "NK cells", "NKT", "T cells", "B cells")

### HYP
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds")
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", "Avp",
                "Slc1a2", "Cx3cr1", "P2ry12", "Tgfbr1", "Cspg4", "Pdgfra","Mbp", "St18", 
                "Col23a1", "Tmem212", "Flt1", "Pecam1", "Ebf1",
                "Atp13a5", "Ptgds")
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds")
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Cx3cr1", "P2ry12", "Tgfbr1", "Arhgap15", "Cd74",
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Lama2")
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds")
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Cx3cr1", "P2ry12", "Tgfbr1", "Arhgap15",
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Atp13a5", "Lama2", "Ptgds")
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


### LV
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Fhl2", "Ttn", "Flt1", "Pecam1", "Pdgfrb", "Rgs5", "Dcn", "Ptprc"
                # c("Ccr2"), # Monocytes
                # c("Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("Cd3e", "Klrb1c", "Ccl5"), # NKT
                # c("Cd19", "Ms4a1", "Cd79a") # B cells
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Fhl2", "Ttn", "Vwf", "Pecam1", "Pdgfrb", "Rgs5", "Dcn", "Ptprc"
                # c("Ly6c2", "Ccr2", "Sell"), # Monocytes
                # c("Adgre1", "Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("S100a8"), # Neutrophils
                # c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                # c("Cd3e", "Klrb1c", "Ccl5"), # NKT
                # c("Cd3d", "Cd4", "Cd8a"), # T cells
                # c("Cd19", "Ms4a1", "Cd79a") # B cells
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Fhl2", "Ttn", "Vwf", "Pecam1", "Pdgfrb", "Rgs5", "Dcn", "Ptprc"
                # c("Ly6c2", "Ccr2", "Sell"), # Monocytes
                # c("Adgre1", "Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                # c("Cd3e", "Klrb1c", "Ccl5", "Cxcr6"), # NKT
                # c("Cd19", "Ms4a1", "Cd79a") # B cells
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


### LK
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Nphs1", "Nphs2",
                "Lrp2", "Slc5a12", # "Slc22a6", "Slc7a13", "Slc13a3", "Slc16a9",
                "Cryab", 
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Flt1", 
                "Myh11", "Notch3", "Fbln5", "Pdgfra", # "Ren1", "Robo1", 
                # "Ncam1", 
                "Ptprc"
                # c("Ccr2"), # Monocytes
                # c("Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("Cd3e", "Klrb1c", "Ccl5"), # NKT
                # c("Cd19", "Ms4a1", "Cd79a") # B cells
                
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Nphs1", "Nphs2",
                "Lrp2", "Slc5a12", # "Slc22a6", "Slc7a13", "Slc13a3", "Slc16a9",
                "Cryab", 
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Pecam1", 
                "Myh11", "Notch3", "Fbln5", "Pdgfra", # "Ren1", "Robo1", 
                # "Ncam1", 
                "Ptprc"
                # c("Ly6c2", "Ccr2", "Sell"), # Monocytes
                # c("Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("S100a8", "Ly6g", "Mpo"), # Neutrophils
                # c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                # c("Cd3e", "Klrb1c", "Ccl5", "Cxcr6"), # NKT
                # c("Cd3d", "Cd4", "Cd8a"), # T cells
                # c("Cd19", "Ms4a1", "Cd79a") # B cells
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Nphs1", "Nphs2",
                "Lrp2", "Slc5a12", # "Slc22a6", "Slc7a13", "Slc13a3", "Slc16a9",
                "Cryab", 
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Pecam1", 
                # "Myh11", "Notch3", 
                "Fbln5", "Pdgfra", # "Ren1", "Robo1", 
                # "Ncam1", 
                "Ptprc"
                # c("Ly6c2", "Ccr2", "Sell"), # Monocytes
                # c("Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("Ncr1", "Klrb1c", "Gzmb", "Ccl5"), # NK cells
                # c("Cd3e", "Cxcr6") # NKT
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)



### MCA
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Vwf", "Mecom",
                "Myh11", 
                "Pdgfrb", "Rgs5", 
                "Ranbp3l", "Lama2",
                c("Cx3cr1", "P2ry12"), # Microglia
                "Ptprc"
                # c("Mrc1", "Cd163"), # Macrophages
                # c("Flt3", "Cd74"), # DC
                # c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                # c("Klrb1c", "Ccl5", "Il2rb"), 
                # c("Cd161", "Cxcr6") # NKT
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


### MSA
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Vwf", "Mecom",
                "Myh11", 
                "Dcn",
                "Adipoq",
                "Ptprc"
                # c("Mrc1", "Cd163")
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds")
seurat_object$subclass_level1 = factor(seurat_object$subclass_level1, levels = c(levels(seurat_object$subclass_level1), "Immune cell"))
seurat_object@meta.data[seurat_object$subclass_level1 %in% immune_list, ]$subclass_level1 = "Immune cell"
marker_list = c("Vwf", "Mecom",
                "Myh11", 
                "Dcn",
                "Adipoq",
                "Ptprc"
                # c("Mrc1", "Cd163"),
                # c("Flt3", "Cd74"), # DC
                # c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                # c("Cd3e", "Cd161", "Il2rb", "Cxcr6") # NKT
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - ", unique(seurat_object$tissue))
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)



################################################################################
mouse_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(mouse_hyp$project), " - ", unique(mouse_hyp$tissue), ")")
FeaturePlot(mouse_hyp, features = "Cd74") + labs(title = title)

ss_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(ss_hyp$project), " - ", unique(ss_hyp$tissue), ")")
FeaturePlot(ss_hyp, features = "Cd74")

sp_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(sp_hyp$project), " - ", unique(sp_hyp$tissue), ")")
FeaturePlot(ss_hyp, features = "Cd74")


mouse_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(mouse_hyp$project), " - ", unique(mouse_hyp$tissue), ")")
FeaturePlot(mouse_hyp, features = "Cd74") + labs(title = title)

ss_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(ss_hyp$project), " - ", unique(ss_hyp$tissue), ")")
FeaturePlot(ss_hyp, features = "Cd74") + labs(title = title)

sp_hyp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds")
title = paste0("Cd74 (", unique(sp_hyp$project), " - ", unique(sp_hyp$tissue), ")")
FeaturePlot(sp_hyp, features = "Cd74") + labs(title = title)
