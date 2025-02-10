dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(SeuratDisk)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)
library(SingleR)
library(celldex)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")


################################################################################
# QC & harmony - mouse
LK_infile="mouse.LK.multiomics.anno.rds"
MCA_infile="mouse.MCA.RNA.anno.rds"
LV_infile="mouse.LV.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
LV_object = readRDS(LV_infile)

DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.5")
FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.5")
FeaturePlot(MCA_object, features = "Ptprc")
DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
FeaturePlot(LV_object, features = "Ptprc")

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_immune = subset(LK_object, wsnn_res.0.5 %in% c(12, 17, 18))
LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(2))
MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_immune = subset(MCA_object, RNA_snn_res.0.5 %in% c(18))

### prepare ref data
ref_object = merge(MCA_immune, y = LV_immune,
                      project = "mouse.immune"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)

### prepare query data
query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
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
merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))

saveRDS(merged_object, "mouse.immune_cell.cluster.rds")







################################################################################
# QC & harmony - rat.ss
LK_infile="rat.ss.LK.multiomics.anno.rds"
MCA_infile="rat.ss.MCA.RNA.anno.rds"
MSA_infile="rat.ss.MSA.RNA.anno.rds"
LV_infile="rat.ss.LV.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)

DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.4")
FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.4")
FeaturePlot(MCA_object, features = "Ptprc")
DimPlot(MSA_object, label = T, group.by = "RNA_snn_res.0.4")
FeaturePlot(MSA_object, features = "Ptprc")
DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
FeaturePlot(LV_object, features = "Ptprc")

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_immune = subset(LK_object, wsnn_res.0.4 %in% c(4, 13, 15))
LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(4, 5))
MCA_immune = subset(MCA_object, RNA_snn_res.0.4 %in% c(2))
MSA_immune = subset(MSA_object, RNA_snn_res.0.4 %in% c(7))

### prepare ref data
ref_object = merge(MCA_immune, y = list(LV_immune, MSA_immune),
                   project = "rat.ss.immune"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)


### prepare query data
query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
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
merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))

saveRDS(merged_object, "rat.ss.immune_cell.cluster.rds")







################################################################################
LK_infile="rat.sp.LK.multiomics.anno.rds"
MCA_infile="rat.sp.MCA.RNA.anno.rds"
MSA_infile="rat.sp.MSA.RNA.anno.rds"
LV_infile="rat.sp.LV.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)

DimPlot(LK_object, label = T, reduction = "wnn.umap.harmony", group.by = "wsnn_res.0.4")
FeaturePlot(LK_object, reduction = "wnn.umap.harmony", features = "Ptprc")
DimPlot(MCA_object, label = T, group.by = "RNA_snn_res.0.5")
FeaturePlot(MCA_object, features = "Ptprc")
DimPlot(MSA_object, label = T, group.by = "RNA_snn_res.0.4")
FeaturePlot(MSA_object, features = "Ptprc")
DimPlot(LV_object, label = T, group.by = "RNA_snn_res.0.1")
FeaturePlot(LV_object, features = "Ptprc")

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_immune = subset(LK_object, wsnn_res.0.4 %in% c(12))
LV_immune = subset(LV_object, RNA_snn_res.0.1 %in% c(4))
MCA_immune = subset(MCA_object, RNA_snn_res.0.5 %in% c(5, 18))
MSA_immune = subset(MSA_object, RNA_snn_res.0.4 %in% c(7))

### prepare ref data
ref_object = merge(MCA_immune, y = list(LV_immune, MSA_immune),
                   project = "rat.sp.immune"
)
ref_object = NormalizeData(ref_object)
ref_object = FindVariableFeatures(ref_object)
ref_object = ScaleData(ref_object, vars.to.regress = c("percent.mt"))
ref_object = RunPCA(ref_object, features = VariableFeatures(object = ref_object), npcs = 30)
ref_object = RunHarmony(ref_object, group.by.vars = "orig.ident")
ref_object = RunUMAP(ref_object, reduction = "harmony", dims=1:30, return.model = T)
ref_object = FindNeighbors(ref_object, reduction = "harmony", dims=1:30)


### prepare query data
query_object = CreateSeuratObject(counts = LK_immune@assays$RNA@counts, 
                                  min.cells = 3, min.features = 200, meta.data = LK_immune@meta.data)
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
merged_object = FindClusters(merged_object, resolution = seq(0.5, 3, 0.5))

saveRDS(merged_object, "rat.sp.immune_cell.cluster.rds")





################################################################################
### immune cell annotation using singleR

input_file = c("mouse.immune_cell.cluster.rds",
               "rat.ss.immune_cell.cluster.rds",
               "rat.sp.immune_cell.cluster.rds")

mimd.sc <- celldex::ImmGenData()
# merged_levels = c("Microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
#                                       "NK cells", "NKT", "T cells", "B cells")
                  
i = input_file[1]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = c("RNA_snn_res.1"))
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.1
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
main.group$labels = c(main.group$labels)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
DimPlot(seurat_object, label = T, group.by = "subclass_level1")
FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
saveRDS(seurat_object, outfile)


i = input_file[2]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
DimPlot(seurat_object, label = T, group.by = "subclass_level1")
FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
seurat_object@meta.data[seurat_object$subclass_level1=="Endothelial cells", ]$subclass_level1 = "Macrophages"
saveRDS(seurat_object, outfile)


i = input_file[3]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
DimPlot(seurat_object, label = T, group.by = "subclass_level1")
FeaturePlot(seurat_object, c("Ptprc", "Mrc1"))
seurat_object@meta.data[seurat_object$subclass_level1=="Endothelial cells", ]$subclass_level1 = "Macrophages"
saveRDS(seurat_object, outfile)





### visualize singleR results
input_file = c("mouse.immune_cell.anno.rds",
               "rat.ss.immune_cell.anno.rds",
               "rat.sp.immune_cell.anno.rds")

mimd.sc <- celldex::ImmGenData()

i = input_file[1]
seurat_object = readRDS(i)
DimPlot(seurat_object, label = T, group.by = c("subclass_level1"))
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=as.factor(seurat_object$subclass_level1))

score_mtx <- main.group$scores
score_mtx_melted <- melt(score_mtx) %>%
  group_by(Var1) %>%
  mutate(annotation = if_else(value == max(value), "*", ""))

ggplot(score_mtx_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = annotation)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.3, 
                       limits = c(min(score_mtx_melted$value), max(score_mtx_melted$value))) +
  theme_minimal() +
  labs(title = "Immune cell annotation (AngII)", 
       x = "ImmGen reference cell type", 
       y = "Final annotation", 
       fill = "Score") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))



i = input_file[2]
seurat_object = readRDS(i)
DimPlot(seurat_object, label = T, group.by = c("subclass_level1"))
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=as.factor(seurat_object$subclass_level1))

score_mtx <- main.group$scores
score_mtx_melted <- melt(score_mtx) %>%
  group_by(Var1) %>%
  mutate(annotation = if_else(value == max(value), "*", ""))

ggplot(score_mtx_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = annotation)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.3, 
                       limits = c(min(score_mtx_melted$value), max(score_mtx_melted$value))) +
  theme_minimal() +
  labs(title = "Immune cell annotation (Salt-sensitive)", 
       x = "ImmGen reference cell type", 
       y = "Final annotation", 
       fill = "Score") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))



i = input_file[3]
seurat_object = readRDS(i)
DimPlot(seurat_object, label = T, group.by = c("subclass_level1"))
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=as.factor(seurat_object$subclass_level1))

score_mtx <- main.group$scores
score_mtx_melted <- melt(score_mtx) %>%
  group_by(Var1) %>%
  mutate(annotation = if_else(value == max(value), "*", ""))

ggplot(score_mtx_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  geom_text(aes(label = annotation)) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.3, 
                       limits = c(min(score_mtx_melted$value), max(score_mtx_melted$value))) +
  theme_minimal() +
  labs(title = "Immune cell annotation (Spontaneous)", 
       x = "ImmGen reference cell type", 
       y = "Final annotation", 
       fill = "Score") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/")
mouse_so = readRDS("mouse.immune_cell.anno.rds")
ss_so = readRDS("rat.ss.immune_cell.anno.rds")
sp_so = readRDS("rat.sp.immune_cell.anno.rds")

# mouse_so$sxt = paste0(mouse_so$strain, "-", mouse_so$treatment)
# mouse_so$subclass_level2 = paste0(mouse_so$subclass_level1, ":C", mouse_so$RNA_snn_res.1)
# markers = DE_analysis(mouse_so, "subclass_level2", outfile = "mouse.immune.")
# DimPlot(mouse_so, label = T, group.by = c("subclass_level1", "tissue", "sxt"))
# 
# ss_so$sxt = paste0(ss_so$strain, "-", ss_so$treatment)
# ss_so$subclass_level2 = paste0(ss_so$subclass_level1, ":C", ss_so$RNA_snn_res.1)
# markers = DE_analysis(ss_so, "subclass_level2", outfile = "ss.immune.")
# DimPlot(ss_so, label = T, group.by = c("subclass_level2", "tissue", "sxt"))
# 
# sp_so$sxt = paste0(sp_so$strain, "-", sp_so$treatment)
# sp_so$subclass_level2 = paste0(sp_so$subclass_level1, ":C", sp_so$RNA_snn_res.1)
# markers = DE_analysis(sp_so, "subclass_level2", outfile = "sp.immune.")
# DimPlot(sp_so, label = T, group.by = c("subclass_level2"))
# 
# 
# barplot = stacked_barplot(mouse_so@meta.data, 
#                           cluster = "subclass_level2", 
#                           tissue = "tissue",
#                           treatment = "treatment")
# wrap_plots(barplot, nrow = 1, widths = c(2,2,1))



seurat_object = mouse_so
cell_list = table(seurat_object$subclass_level1)
cell_list = names(cell_list)[cell_list>=10]
seurat_object = subset(seurat_object, subclass_level1 %in% cell_list)
seurat_object = subset(seurat_object, subclass_level1 != "Microglia")
marker_list = c("Ptprc", 
                c("Cd19", "Ms4a1", "Cd79a"), # B cells
                c("Itgax", "Flt3", "Cd74"), # DC
                c("Adgre1", "Mrc1", "Cd163"), # Macrophages
                c("Ly6c2", "Ccr2"), # Monocytes
                c("Cd3e", "Klrb1c", "Ccl5") # NKT
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - Immune cells")
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = ss_so
cell_list = table(seurat_object$subclass_level1)
cell_list = names(cell_list)[cell_list>=10]
seurat_object = subset(seurat_object, subclass_level1 %in% cell_list)
seurat_object = subset(seurat_object, subclass_level1 != "Microglia")
marker_list = c("Ptprc", 
                c("Cd19", "Ms4a1", "Cd79a"), # B cells
                c("Itgax", "Flt3", "Cd74"), # DC
                c("Adgre1", "Mrc1", "Cd163"), # Macrophages
                c("Cd14", "Fcgr3", "Itgam", "Cd68", "Ccr2", "Ly6c1", "Emr1"), # Monocytes
                c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                c("Cd3e", "Klrb1c", "Ccl5"), # NKT
                c("S100a8", "Ly6g", "Mpo"), # Neutrophils
                c("Cd3d", "Cd4", "Cd8a") # T cells
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - Immune cells")
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)


seurat_object = sp_so
cell_list = table(seurat_object$subclass_level1)
cell_list = names(cell_list)[cell_list>=10]
seurat_object = subset(seurat_object, subclass_level1 %in% cell_list)
seurat_object = subset(seurat_object, subclass_level1 != "Microglia")
marker_list = c("Ptprc", 
                c("Cd19", "Ms4a1", "Cd79a"), # B cells
                c("Itgax", "Flt3", "Cd74"), # DC
                c("Adgre1", "Mrc1", "Cd163"), # Macrophages
                c("Cd14", "Fcgr3", "Itgam", "Cd68", "Ccr2", "Ly6c1", "Emr1"), # Monocytes
                c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
                c("Cd3e", "Klrb1c", "Ccl5") # NKT
)
marker_list = unique(marker_list)
title = paste0(unique(seurat_object$project), " - Immune cells")
DotPlot(seurat_object, features = marker_list, group.by = "subclass_level1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = title)




################################################################################
# 
# 
# input_file = c("mouse.immune_cell.cluster.rds",
#                "rat.ss.immune_cell.cluster.rds",
#                "rat.sp.immune_cell.cluster.rds")
# 
# anno_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/immune_cell.anno.txt", header=T, sep="\t", fill=T)
# 
# sample_info = read.table("/xdisk/mliang1/qqiu/project/multiomics/data/Multiomics_sample_info.txt", header=T, sep='\t')
# rownames(sample_info) = sample_info$seqID
# 
# project_list = c("AngII", "Salt-sensitive", "Spontaneous"); names(project_list)=c("mouse", "rat.ss", "rat.sp")
# tissue_list = c("3rd mesenteric artery"="MSA", "Cardiac left ventricle"="LV", "Hypothalamus"="HYP", "Kidney"="LK", "Middle cerebral artery"="MCA")
# 
# for(i in input_file){
#   
#   dataset = gsub(".cluster.rds", "", basename(i), perl = T)
#   project = gsub(".immune_cell", "", dataset)
#   outfile = gsub("cluster.rds", "anno.rds", i)
#   
#   seurat_object = readRDS(i)
#   
#   seurat_object = FindClusters(seurat_object, resolution = c(2, 3))
#   
#   merged_levels = c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
#                     "NK cells", "NKT", "T cells", "B cells")
# 
#   
#   anno_info_use = anno_info[anno_info$dataset == dataset, ]
#   reso = unique(anno_info_use$resolution)
#   cluster_use = anno_info_use[!(grepl("remove", anno_info_use$note)),]$cluster
#   subclass_level1 = anno_info_use[!(grepl("remove", anno_info_use$note)),]$subclass_level1
#   class = anno_info_use[!(grepl("remove", anno_info_use$note)),]$class
#   
#   Idents(seurat_object) = reso
#   seurat_object$seurat_clusters = seurat_object@meta.data[, reso]
#   seurat_object = subset(seurat_object, seurat_clusters %in% cluster_use)
#   seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)
#   
#   seurat_object@meta.data$"subclass_level2" = factor(subclass_level1[match(seurat_object$seurat_clusters, cluster_use)],
#                                                      levels=merged_levels)
#   
#   HYP_object = readRDS(paste0(project, ".HYP.RNA.anno.rds"))
#   HYP_mg = subset(HYP_object, subclass_level1 %in% c("Microglia",  "Activated microglia"))
#   HYP_mg$subclass_level2 = HYP_mg$subclass_level1
#   
#   seurat_object = merge(seurat_object, HYP_mg)
#   
#   seurat_object = NormalizeData(seurat_object)
#   seurat_object = FindVariableFeatures(seurat_object)
#   seurat_object = ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
#   seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 30)
#   seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident")
#   seurat_object = RunUMAP(seurat_object, reduction = "harmony", dims=1:30, return.model = T)
#   seurat_object = FindNeighbors(seurat_object, reduction = "harmony", dims=1:30)
#   
# 
#   seurat_object$cell_id = colnames(seurat_object)
#   seurat_object$seqID2 = sample_info[as.character(seurat_object$orig.ident),]$seqID2
#   seurat_object$sample = sample_info[as.character(seurat_object$orig.ident),]$sample
#   seurat_object$species = factor(sample_info[as.character(seurat_object$orig.ident),]$strain,
#                                  levels=c("C57BL/6", "WKY", "SHR", "SD", "SS"))
#   seurat_object$tissue = factor(tissue_list[sample_info[as.character(seurat_object$orig.ident),]$tissue],
#                                 levels=c("HYP", "MCA", "LV", "LK", "MSA"))
#   seurat_object$treatment = factor(sample_info[as.character(seurat_object$orig.ident),]$treatment,
#                                    levels=c("Saline 3d", "AngII 3d", "AngII 28d",
#                                             "10w", "26w", "LS", "HS 3d", "HS 21d"))
#   seurat_object$project = factor(project_list[project], levels=c("AngII", "Salt-sensitive", "Spontaneous"))
#   
#   
#   saveRDS(seurat_object, outfile)
#   print(c(dataset, all(unique(subclass_level1) %in% merged_levels)))
#   
#   # p = DimPlot(seurat_object, reduction = "umap", group.by = "subclass_level2",
#   #             pt.size = 0.2, raster = F, label = T, ncol=5) 
#   # print(p)
#   # 
# }
# 
# 
# # 
# 
# 
# ################################################################################
# 
# 
# input_file = c("mouse.immune_cell.anno.rds",
#                "rat.ss.immune_cell.anno.rds",
#                "rat.sp.immune_cell.anno.rds")
# 
# angii = readRDS(input_file[1])
# ss = readRDS(input_file[2])
# sp = readRDS(input_file[3])
# 
# 
# 
# pseudo_angii = AggregateExpression(angii, assays = "RNA", return.seurat = T, group.by = c("subclass_level2", "species", "tissue", "treatment"))
# pseudo_ss = AggregateExpression(ss, assays = "RNA", return.seurat = T, group.by = c("subclass_level2", "species", "tissue", "treatment"))
# pseudo_sp = AggregateExpression(sp, assays = "RNA", return.seurat = T, group.by = c("subclass_level2", "species", "tissue", "treatment"))
# 
# pseudo_angii = NormalizeData(pseudo_angii, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
# pseudo_ss = NormalizeData(pseudo_ss, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
# pseudo_sp = NormalizeData(pseudo_sp, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
# 
# cell_count = c(table(apply(angii@meta.data[c("subclass_level2", "species", "tissue", "treatment")], 1, function(x) paste0(x, collapse = "_"))),
#                table(apply(ss@meta.data[c("subclass_level2", "species", "tissue", "treatment")], 1, function(x) paste0(x, collapse = "_"))),
#                table(apply(sp@meta.data[c("subclass_level2", "species", "tissue", "treatment")], 1, function(x) paste0(x, collapse = "_"))))
# names(cell_count) = gsub(" ", ".", names(cell_count))
# 
# exclude_list = c("Microglia_LK", "Microglia_LV", "Microglia_MSA", "Monocytes_MCA", "DC_MCA")
# 
# expr_list = list(pseudo_angii@assays$RNA@data, pseudo_ss@assays$RNA@data, pseudo_sp@assays$RNA@data)
# expr_mat = Reduce(merge, lapply(expr_list, function(x) data.frame(x, rn = row.names(x))))
# rownames(expr_mat) = expr_mat[,1]; expr_mat = expr_mat[,-1]
# expr_mat = expr_mat[,!(sapply(colnames(expr_mat), function(item) any(sapply(exclude_list, grepl, item))))]
# expr_mat = expr_mat[rowSums(expr_mat==0)<ncol(expr_mat), ]
# cell_type = unlist(lapply(strsplit(colnames(expr_mat), "_"), function(x) x[1]))
# strain = unlist(lapply(strsplit(colnames(expr_mat), "_"), function(x) x[2]))
# organ_list = unlist(lapply(strsplit(colnames(expr_mat), "_"), function(x) x[3]))
# treatment = unlist(lapply(strsplit(colnames(expr_mat), "_"), function(x) x[4]))
# 
# 
# # deg_merged = read.table("../DEG/DEG.merged.out")
# 
# 
# breaks = c(10, 100, 1000, 2000)
# 
# cell_list = c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
#   "NK cells", "NKT", "T cells", "B cells")
# for(cell in cell_list){
#   
#   selected_col = grepl(gsub(" ", ".", cell), colnames(expr_mat))
#   
#   gene_used = unique(deg_merged[deg_merged$cell_type %in% cell, ]$gene_name)
#   expr_mat_use = expr_mat[intersect(gene_used, rownames(expr_mat)), selected_col]
#   expr_mat_scale = t(scale(t(expr_mat_use)))
#   
#   top_Annotation = HeatmapAnnotation(Cell_type = cell_type[selected_col],
#                                      Strain = strain[selected_col],
#                                      Tissue = organ_list[selected_col],
#                                      Treatment = treatment[selected_col],
#                                      Cell_count = anno_barplot(as.numeric(cell_count[colnames(expr_mat)[selected_col]]))
#                                      )
#   p = Heatmap(expr_mat_scale, name = "Contribution",
#               top_annotation = top_Annotation,
#               show_row_names = F, show_column_names = T, 
#               cluster_rows = T ,cluster_columns = T, 
#               # row_names_side = 'left',col = col_map,
#               row_title_rot = 0, 
#               show_heatmap_legend = T,
#               column_split = cell_type[selected_col],
#               heatmap_legend_param = list(color_bar='continuous'))
#   
#   outfile = paste0(cell, ".pseudobulk.heatmap.pdf")
#   pdf(outfile, width=12, height=5)
#   print(p)
#   dev.off()
# }
# 
# 
# 
# 
# 
# hc = hclust(dist(expr_mat))
# Heatmap(matrix(nc = 0, nr = 10), cluster_rows = hc, 
#         right_annotation = rowAnnotation(
#           "Cell type" = cell_type,
#           Tissue = organ_list,
#           Treatment = treatment,
#           "Cell count" = anno_barplot(cell_count[colnames(expr_mat)])),
#         row_split = cell_type)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################
# sample_info = read.table("/xdisk/mliang1/qqiu/project/multiomics/data/Multiomics_sample_info.txt", header=T, sep='\t')
# rownames(sample_info) = sample_info$seqID
# 
# 
# ang_immune = readRDS("mouse.immune_cell.cluster.rds")
# ss_immune = readRDS("rat.ss.immune_cell.cluster.rds")
# sp_immune = readRDS("rat.sp.immune_cell.cluster.rds")
# 
# tmp=as.data.frame(table(sp_immune$RNA_snn_res.2))
# tmp[order(as.numeric(as.character(tmp$Var1))),]
# 
# # for(i in c(0.4)){
# #   cluster = paste0("RNA_snn_res.", i)
# #   outfile = paste0("../DEG/mouse.immune_cell.reso_", i)
# #   DE_analysis(ang_immune, cluster = cluster, outfile = outfile)
# # }
# # 
# # for(i in c(0.3)){
# #   cluster = paste0("RNA_snn_res.", i)
# #   outfile = paste0("../DEG/rat.ss.immune_cell.reso_", i)
# #   DE_analysis(ss_immune, cluster = cluster, outfile = outfile)
# # }
# # 
# # for(i in c(0.3)){
# #   cluster = paste0("RNA_snn_res.", i)
# #   outfile = paste0("../DEG/rat.sp.immune_cell.reso_", i)
# #   DE_analysis(sp_immune, cluster = cluster, outfile = outfile)
# # }
# 
# ang_immune = FindClusters(ang_immune, resolution = 2)
# 
# ang_immune$cell_id = colnames(ang_immune)
# ang_immune$seqID2 = sample_info[as.character(ang_immune$orig.ident),]$seqID2
# ang_immune$sample = sample_info[as.character(ang_immune$orig.ident),]$sample
# ang_immune$species = sample_info[as.character(ang_immune$orig.ident),]$strain
# ang_immune$tissue = sample_info[as.character(ang_immune$orig.ident),]$tissue
# ang_immune$treatment = sample_info[as.character(ang_immune$orig.ident),]$treatment
# 
# ang_immune_anno = read.table("auto_annotation/mouse.immune_cell.csv", sep=',', header=T, row.names = 1)
# ang_immune$auto.main = ang_immune_anno[colnames(ang_immune),]$main.group
# ang_immune$auto.fine = ang_immune_anno[colnames(ang_immune),]$fine.group
# 
# # cell_remain = colnames(ang_immune)[!(ang_immune$RNA_snn_res.2 %in% c(3, 4, 6, 13, 18, 20))]
# # ang_immune@active.ident <- factor(ang_immune@active.ident)
# # names(ang_immune@active.ident) <- colnames(ang_immune)
# # ang_immune = subset(ang_immune, cell_id %in% cell_remain)
# 
# tmp = table(ang_immune$auto.main, ang_immune$RNA_snn_res.2)
# t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# 
# # tmp = table(ang_immune$auto.main, paste(ang_immune$tissue))
# # t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# # 
# # tmp = table(ang_immune$auto.main, paste(ang_immune$tissue, ang_immune$treatment))
# # t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# # 
# # tmp = table(ang_immune$auto.main, paste(ang_immune$RNA_snn_res.2, ang_immune$tissue))
# # tmp / colSums(tmp) * 100
# 
# DimPlot(ang_immune, reduction = "umap", split.by = c("auto.main"), group.by = "tissue",
#          pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# DimPlot(ang_immune, reduction = "umap", split.by = c("auto.main"), group.by = "treatment",
#         pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# 
# DotPlot(seurat_object, group.by = "RNA_snn_res.2", 
#         split.by = "tissue", cols = c("blue", "green", "red"), 
#         features = marker_list) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   # scale_color_gradientn(colors = colGEX) +
#   labs(x="", y="")
# 
# 
# 
# 
# 
# 
# 
# ss_immune = FindClusters(ss_immune, resolution = 3)
# ss_immune$cell_id = colnames(ss_immune)
# ss_immune$seqID2 = sample_info[as.character(ss_immune$orig.ident),]$seqID2
# ss_immune$sample = sample_info[as.character(ss_immune$orig.ident),]$sample
# ss_immune$species = sample_info[as.character(ss_immune$orig.ident),]$strain
# ss_immune$tissue = sample_info[as.character(ss_immune$orig.ident),]$tissue
# ss_immune$treatment = sample_info[as.character(ss_immune$orig.ident),]$treatment
# 
# ss_immune_anno = read.table("auto_annotation/rat.ss.immune_cell.csv", sep=',', header=T, row.names = 1)
# ss_immune$auto.main = ss_immune_anno[colnames(ss_immune),]$main.group
# ss_immune$auto.fine = ss_immune_anno[colnames(ss_immune),]$fine.group
# 
# seurat_object = ss_immune
# marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", 
#                 "Cdk14", "Cd36", "Csmd1",
#                 "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4", 
#                 "Cd3e", "Cd8a", "Klrb1c",
#                 "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
# umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.3", reduction = "umap")
# 
# tmp = table(ss_immune$auto.main, ss_immune$RNA_snn_res.3)
# t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# 
# 
# DimPlot(ss_immune, reduction = "umap", split.by = c("auto.main"), group.by = "tissue",
#         pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# cell_remain = colnames(ss_immune)[!(ss_immune$RNA_snn_res.2 %in% c(7, 10, 15, 18, 21, 31, 32, 34, 36))]
# ss_immune = subset(ss_immune, cell_id %in% cell_remain)
# 
# 
# DimPlot(ss_immune, reduction = "umap", split.by = c("auto.main"), group.by = "treatment",
#         pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# 
# DotPlot(seurat_object, group.by = "RNA_snn_res.1", 
#         split.by = "tissue", cols = c("blue", "green", "red"), 
#         features = marker_list) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   # scale_color_gradientn(colors = colGEX) +
#   labs(x="", y="")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# sp_immune = FindClusters(sp_immune, resolution = c(2, 3))
# sp_immune$cell_id = colnames(sp_immune)
# sp_immune$seqID2 = sample_info[as.character(sp_immune$orig.ident),]$seqID2
# sp_immune$sample = sample_info[as.character(sp_immune$orig.ident),]$sample
# sp_immune$species = sample_info[as.character(sp_immune$orig.ident),]$strain
# sp_immune$tissue = sample_info[as.character(sp_immune$orig.ident),]$tissue
# sp_immune$treatment = sample_info[as.character(sp_immune$orig.ident),]$treatment
# 
# sp_immune_anno = read.table("auto_annotation/rat.sp.immune_cell.csv", sep=',', header=T, row.names = 1)
# sp_immune$auto.main = sp_immune_anno[colnames(sp_immune),]$main.group
# sp_immune$auto.fine = sp_immune_anno[colnames(sp_immune),]$fine.group
# 
# seurat_object = sp_immune
# marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", "Cdk14", "Cd36", "Csmd1",
#                 "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
#                 "Cd3e", "Cd8a", "Klrb1c",
#                 "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
# marker_list = c("Cald1", "Fbln5", "Pdgfra", "Dcn", 
#                 "Myh11", "Notch3", 
#                 "Ptgds",
#                 "Ebf1", "Atp13a5", "Rgs5", "Lama2",
#                 "Ranbp3l",
#                 "Pecam1", "Mecom", "vwf")
# umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.2", reduction = "umap")
# 
# cell_remain = colnames(sp_immune)[!(sp_immune$auto.main %in% c("Stromal cells", "Fibroblasts") & 
#                                       sp_immune$RNA_snn_res.2 %in% c(6, 12, 13, 14, 22))]
# sp_immune = subset(sp_immune, cell_id %in% cell_remain)
# 
# tmp = table(sp_immune$auto.main, sp_immune$RNA_snn_res.2)
# t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# 
# tmp = table(sp_immune$auto.main, sp_immune$RNA_snn_res.1)
# t(apply(tmp, 1, function(x) x / colSums(tmp) * 100))
# 
# DimPlot(sp_immune, reduction = "umap", split.by = c("auto.main"), group.by = "tissue",
#         pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# DimPlot(sp_immune, reduction = "umap", split.by = c("auto.main"), group.by = "treatment",
#         pt.size = 0.2, raster = F, label = F, ncol=5) 
# 
# 
# DotPlot(seurat_object, group.by = "RNA_snn_res.1", 
#         split.by = "tissue", cols = c("blue", "green", "red"), 
#         features = marker_list) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   # scale_color_gradientn(colors = colGEX) +
#   labs(x="", y="")
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# sp_immune_anno = read.table("auto_annotation/rat.sp.immune_cell.csv", sep=',', header=T, row.names = 1)
# 
# 
# 
# sp_immune$auto.main = sp_immune_anno[colnames(sp_immune),]$main.group
# sp_immune$auto.fine = sp_immune_anno[colnames(sp_immune),]$fine.group
# 
# 
# 
# 
# 
# ang_immune_mkr = read.table("../DEG/mouse.immune_cell.reso_0.1.allmarker.0.25.wide.txt", header=T)
# ss_immune_mkr = read.table("../DEG/rat.ss.immune_cell.reso_0.3.allmarker.0.25.wide.txt", header=T)
# sp_immune_mkr = read.table("../DEG/rat.sp.immune_cell.reso_0.3.allmarker.0.25.wide.txt", header=T)
# 
# seurat_object = ang_immune
# marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd79a", "Meg3", "Pid1", "Cdk14", "Cd36", "Csmd1",
#                 "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
#                 "Flt1", "Mecom", "Myh11", "Ncam1", "Syt1")
# umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.2", reduction = "umap")
# 
# 
# seurat_object = ss_immune
# marker_list = c("Ptprc", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", "Cdk14", "Cd36", "Csmd1",
#                 "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
#                 "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
# umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.3", reduction = "umap")
# 
# 
# seurat_object = sp_immune
# marker_list = c("Ptprc", "Il3ra", "Thbd", "Cd1c", 
#                 "Fcgr3a", "Cd14", "Mrc1", "F13a1", "Rbpj", "Bank1", "Cd74", "Mki67", "Top2a", "Cdk14", "Cd36", "Csmd1",
#                 "Gm2682","Skap1","Cd226","Themis", "Cd247", "Camk4",
#                 "Pecam1", "Mecom", "Myh11", "Ncam1", "Syt1")
# umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.3", reduction = "umap")
# 
# 
# 
# 
# 
# DotPlot(seurat_object, group.by = "RNA_snn_res.2", 
#         split.by = "tissue", cols = c("blue", "green", "red"), 
#         features = marker_list) +
#   theme(axis.text.x = element_text(angle = 45, hjust=1)) +
#   # scale_color_gradientn(colors = colGEX) +
#   labs(x="", y="")
# 
# 
# 
# DimPlot(seurat_object, reduction = "umap", split.by = c("auto.main"), 
#         pt.size = 0.2, raster = F, label = T) 
# 


