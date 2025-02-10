dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")
# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DE_enrichment.R")

library(Seurat)
library(SeuratDisk)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")


sample_info = read.table("/xdisk/mliang1/qqiu/project/multiomics/data/Multiomics_sample_info.txt", header=T, sep='\t')
rownames(sample_info) = sample_info$seqID

# QC & harmony - mouse
LK_infile="mouse.LK.multiomics.cluster.rds"
MCA_infile="mouse.MCA.RNA.cluster.rds"
LV_infile="mouse.LV.RNA.cluster.rds"
HYP_infile="mouse.HYP.RNA.cluster.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, wsnn_res.0.5 %in% c(4))
LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, RNA_snn_res.0.1 %in% c(0))
MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, RNA_snn_res.0.5 %in% c(16))
HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, RNA_snn_res.0.9 %in% c(16))

### prepare ref data
ref_object = merge(MCA_EC, y = list(LV_EC, HYP_EC),
                      project = "mouse.EC"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)

### prepare query data
query_object = CreateSeuratObject(counts = LK_EC@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_EC@meta.data)
query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
query_object = NormalizeData(query_object)
# query_object = FindVariableFeatures(query_object)
# query_object = ScaleData(query_object, vars.to.regress = c("percent.mt"))
# query_object = RunPCA(query_object, features = VariableFeatures(object = query_object))
# query_object = RunHarmony(query_object, group.by.vars = "orig.ident")
# query_object = RunUMAP(query_object, reduction = "harmony", dims=1:30)
# query_object = FindNeighbors(query_object, reduction = "harmony", dims=1:30)


transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
  k.filter = NA, reference.reduction = 'pca', dims = 1:30)

query_object = MapQuery(anchorset = transfer_anchors,
                             query = query_object, reference = ref_object,
                             reference.reduction = "pca",
                             reduction.model = "umap")

merged_object = merge(ref_object, query_object)
merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
merged_object = FindClusters(merged_object, resolution = seq(0.1, 1, 0.1))

saveRDS(merged_object, "mouse.EC.cluster.rds")


DimPlot(merged_object, label = T, group.by = c("tissue", "treatment", "seqID2"))
FeaturePlot(merged_object, features = marker_list)
DotPlot(merged_object, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")


marker_list = c("Pecam1", "Kdr", "Cldn5", "Emcn", "Cdh5",
                "Apj", "Nr2f2", "Vwf",
                "Gja4", "Mecom")



################################################################################
# QC & harmony - mouse

LK_infile="rat.ss.LK.multiomics.cluster.rds"
MCA_infile="rat.ss.MCA.RNA.cluster.rds"
MSA_infile="rat.ss.MSA.RNA.cluster.rds"
LV_infile="rat.ss.LV.RNA.cluster.rds"
HYP_infile="rat.ss.HYP.RNA.cluster.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, wsnn_res.0.4 %in% c(7, 24))
LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, RNA_snn_res.0.1 %in% c(0, 6, 8))
MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, RNA_snn_res.0.4 %in% c(12))
MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, RNA_snn_res.0.4 %in% c(4, 6))
HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, RNA_snn_res.1.2 %in% c(28, 33))

### prepare ref data
ref_object = merge(MCA_EC, y = list(LV_EC, HYP_EC, MSA_EC),
                   project = "rat.ss.EC"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)

### prepare query data
query_object = CreateSeuratObject(counts = LK_EC@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_EC@meta.data)
query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
query_object = NormalizeData(query_object)

transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
                                       k.filter = NA, reference.reduction = 'pca', dims = 1:30)

query_object = MapQuery(anchorset = transfer_anchors,
                        query = query_object, reference = ref_object,
                        reference.reduction = "pca",
                        reduction.model = "umap")

merged_object = merge(ref_object, query_object)
merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
merged_object = FindClusters(merged_object, resolution = seq(0.1, 1, 0.1))

saveRDS(merged_object, "rat.ss.EC.cluster.rds")


DimPlot(merged_object, label = T, group.by = c("tissue", "treatment", "seqID2"))
FeaturePlot(merged_object, features = c("Ptprc"))









################################################################################
LK_infile="rat.sp.LK.multiomics.cluster.rds"
MCA_infile="rat.sp.MCA.RNA.cluster.rds"
MSA_infile="rat.sp.MSA.RNA.cluster.rds"
LV_infile="rat.sp.LV.RNA.cluster.rds"
HYP_infile="rat.sp.HYP.RNA.cluster.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, wsnn_res.0.4 %in% c(9))
LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, RNA_snn_res.0.1 %in% c(0, 6))
MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, RNA_snn_res.0.5 %in% c(4, 9))
MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, RNA_snn_res.0.4 %in% c(3, 5, 8))
HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, RNA_snn_res.1 %in% c(24))

### prepare ref data
ref_object = merge(MCA_EC, y = list(LV_EC, HYP_EC, MSA_EC),
                   project = "rat.ss.EC"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)

### prepare query data
query_object = CreateSeuratObject(counts = LK_EC@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_EC@meta.data)
query_object[["percent.mt"]] = PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
query_object = NormalizeData(query_object)

transfer_anchors = FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
                                       k.filter = NA, reference.reduction = 'pca', dims = 1:30)

query_object = MapQuery(anchorset = transfer_anchors,
                        query = query_object, reference = ref_object,
                        reference.reduction = "pca",
                        reduction.model = "umap")

merged_object = merge(ref_object, query_object)
merged_object[["pca"]] = merge(ref_object[["pca"]], query_object[["ref.pca"]])
merged_object = RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = F)
merged_object = RunUMAP(object = merged_object, reduction = "harmony", dims = 1:30)
merged_object = FindNeighbors(merged_object, reduction = "harmony", dims=1:30)
merged_object = FindClusters(merged_object, resolution = seq(0.1, 1, 0.1))

saveRDS(merged_object, "rat.sp.EC.cluster.rds")


DimPlot(merged_object, label = T, group.by = c("tissue", "treatment", "seqID2"))
FeaturePlot(merged_object, features = c("Mrc1"))








################################################################################
ang_immune = readRDS("mouse.immune_cell.cluster.rds")
ss_immune = readRDS("rat.ss.immune_cell.cluster.rds")
sp_immune = readRDS("rat.sp.immune_cell.cluster.rds")

for(i in c(0.1)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("../DEG/mouse.immune_cell.reso_", i)
  DE_analysis(ang_immune, cluster = cluster, outfile = outfile)
}

for(i in c(0.3)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("../DEG/rat.ss.immune_cell.reso_", i)
  DE_analysis(ss_immune, cluster = cluster, outfile = outfile)
}

for(i in c(0.3)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("../DEG/rat.sp.immune_cell.reso_", i)
  DE_analysis(sp_immune, cluster = cluster, outfile = outfile)
}



ang_immune_mkr = read.table("../DEG/mouse.immune_cell.reso_0.1.allmarker.0.25.wide.txt", header=T)
ss_immune_mkr = read.table("../DEG/rat.ss.immune_cell.reso_0.3.allmarker.0.25.wide.txt", header=T)
sp_immune_mkr = read.table("../DEG/rat.sp.immune_cell.reso_0.3.allmarker.0.25.wide.txt", header=T)

seurat_object = ang_immune
marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd79a", "Meg3", "Pid1", "Cdk14", "Cd36", "Csmd1",
                "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
                "Flt1", "Mecom", "Myh11", "Ncam1", "Syt1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.1", reduction = "umap")


seurat_object = ss_immune
marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", "Cdk14", "Cd36", "Csmd1",
                "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
                "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.3", reduction = "umap")


seurat_object = sp_immune
marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", "Cdk14", "Cd36", "Csmd1",
                "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
                "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.3", reduction = "umap")













