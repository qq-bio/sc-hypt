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

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster")


# sample_info = read.table("/xdisk/mliang1/qqiu/project/multiomics/data/Multiomics_sample_info.txt", header=T, sep='\t')
# rownames(sample_info) = sample_info$seqID

################################################################################

LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds"
MCA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
saveRDS(LK_EC, "mouse.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
saveRDS(LV_EC, "mouse.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
saveRDS(MCA_EC, "mouse.MCA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC"))
saveRDS(HYP_EC, "mouse.HYP.EC.cluster.rds")



LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds"
MCA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds"
MSA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
saveRDS(LK_EC, "rat.ss.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
saveRDS(LV_EC, "rat.ss.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
saveRDS(MCA_EC, "rat.ss.MCA.EC.cluster.rds")

MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, subclass_level1 %in% c("EC"))
saveRDS(MSA_EC, "rat.ss.MSA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC"))
saveRDS(HYP_EC, "rat.ss.HYP.EC.cluster.rds")




# LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds"
MCA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
MSA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
saveRDS(LK_EC, "rat.sp.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
saveRDS(LV_EC, "rat.sp.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
saveRDS(MCA_EC, "rat.sp.MCA.EC.cluster.rds")

MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, subclass_level1 %in% c("EC"))
saveRDS(MSA_EC, "rat.sp.MSA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC"))
saveRDS(HYP_EC, "rat.sp.HYP.EC.cluster.rds")




snRNA_list = c("mouse.LV.EC.cluster.rds", "mouse.MCA.EC.cluster.rds",
               "mouse.HYP.EC.cluster.rds",
               "rat.ss.LV.EC.cluster.rds", "rat.ss.MCA.EC.cluster.rds",
               "rat.ss.MSA.EC.cluster.rds", "rat.ss.HYP.EC.cluster.rds",
               "rat.sp.LV.EC.cluster.rds", "rat.sp.MCA.EC.cluster.rds",
               "rat.sp.MSA.EC.cluster.rds", "rat.sp.HYP.EC.cluster.rds")

snMulti_list = c("mouse.LK.EC.cluster.rds",
                 "rat.ss.LK.EC.cluster.rds",
                 "rat.sp.LK.EC.cluster.rds")

for( i in snRNA_list){

  seurat_object = readRDS(i)
  seurat_object = NormalizeData(seurat_object)
  seurat_object = FindVariableFeatures(seurat_object)
  seurat_object = ScaleData(seurat_object)
  seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object), npcs = 20)
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident",project.dim = F)
  seurat_object = RunUMAP(seurat_object, reduction = "harmony", dims=1:20, return.model = T)
  seurat_object = FindNeighbors(seurat_object, reduction = "harmony", dims=1:20)
  seurat_object = FindClusters(seurat_object, resolution = c(seq(0.1, 1, 0.1), 1, 1.5, 2))

  saveRDS(seurat_object, i)
}



for( i in snMulti_list){

  seurat_object = readRDS(i)
  ### process the RNA data
  DefaultAssay(seurat_object) = "RNA"
  seurat_object = NormalizeData(seurat_object)
  seurat_object = FindVariableFeatures(seurat_object)
  seurat_object = ScaleData(seurat_object)
  seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "pca",
                             reduction.save = "harmony.rna",assay.use = "RNA", project.dim = F)

  ## process the ATAC data
  DefaultAssay(seurat_object) = "ATAC"
  seurat_object@assays$ATAC@counts
  seurat_object= RunTFIDF(seurat_object)
  seurat_object = FindTopFeatures(seurat_object, min.cutoff='q0')
  seurat_object = RunSVD(seurat_object)
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "lsi",
                             reduction.save = "harmony.atac", assay.use = "ATAC", project.dim = F)

  ### integration & cluster
  seurat_object = FindMultiModalNeighbors(seurat_object, reduction.list = list("harmony.rna", "harmony.atac"),
                                          dims.list = list(1:20, 2:20))
  seurat_object = RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", reduction.key = "wnnUMAPHarmony_")
  seurat_object = FindClusters(seurat_object, graph.name = "wsnn", algorithm = 3,
                               resolution = c(seq(0.1, 1, 0.1), 1, 1.5, 2))

  DefaultAssay(seurat_object) = "RNA"
  saveRDS(seurat_object, i)
}

# DimPlot(seurat_object, label = T, group.by = c("tissue", "treatment", "seqID2", "wsnn_res.1"), reduction = "wnn.umap.harmony")









for(i in c(snRNA_list, snMulti_list)){
  
  seurat_object = readRDS(i)
  
  if(grepl("LK", i)){
    p = DimPlot(seurat_object, label = T, group.by = c("tissue", "treatment", "seqID2", "wsnn_res.1"), reduction = "wnn.umap.harmony")
  }else{
    p = DimPlot(seurat_object, label = T, group.by = c("tissue", "treatment", "seqID2", "RNA_snn_res.1"), reduction = "umap")
  }
  print(p)
  
  for(j in c(0.1)){
    
    if(grepl("LK", i)){
      cluster = paste0("wsnn_res.", j)
    }else{
      cluster = paste0("RNA_snn_res.", j)
    }
    
    outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", gsub("cluster.rds", paste0("reso_", j), i))
    
    seurat_object@active.ident = seurat_object@meta.data[, cluster]
    Idents(seurat_object) = cluster
    
    DefaultAssay(seurat_object) = "RNA"
    
    # cluster_sample_tbl = table(seurat_object$seqID2, seurat_object@meta.data[, cluster])
    # cluster_sample_sum(cluster_sample_tbl, paste0(outfile, ".cluster_sample_sum.xlsx"))
    
    markers = FindAllMarkers(seurat_object, group.by=cluster, only.pos = TRUE, min.pct = 0.25)
    markers$pct.diff = markers$pct.1 - markers$pct.2
    markers = markers %>% group_by(cluster) %>%
      dplyr::arrange(desc(pct.diff), .by_group=TRUE)
    write.table(markers,file=paste0(outfile, ".allmarker.0.25.long.txt"), sep="\t")
    
    marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
      dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
      dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
      dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                                   v.names="gene", sep=" ", direction = "wide")
    colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
    write.table(marker_tbl,file=paste0(outfile, ".allmarker.0.25.wide.txt"), sep="\t", row.names = F)
    
    panglao_anno(marker_tbl, paste0(outfile, ".allmarker.0.25.wide.xlsx"), anno_file="/xdisk/mliang1/qqiu/reference/single_cell_marker_ref/Cell-EC-2020.txt")
    
  }
  
}









snRNA_list = c("mouse.LV.EC.cluster.rds", "mouse.MCA.EC.cluster.rds", 
               "mouse.HYP.EC.cluster.rds", 
               "rat.ss.LV.EC.cluster.rds", "rat.ss.MCA.EC.cluster.rds",
               "rat.ss.MSA.EC.cluster.rds", "rat.ss.HYP.EC.cluster.rds",
               "rat.sp.LV.EC.cluster.rds", "rat.sp.MCA.EC.cluster.rds",
               "rat.sp.MSA.EC.cluster.rds", "rat.sp.HYP.EC.cluster.rds")

snMulti_list = c("mouse.LK.EC.cluster.rds", 
                 "rat.ss.LK.EC.cluster.rds", 
                 "rat.sp.LK.EC.cluster.rds")

angii_LV = readRDS("mouse.LV.EC.cluster.rds")
ss_LV = readRDS("rat.ss.LV.EC.cluster.rds")
sp_LV = readRDS("rat.sp.LV.EC.cluster.rds")

EC = read.table("/xdisk/mliang1/qqiu/reference/single_cell_marker_ref/Cell-EC-2020.txt", header=T, sep='\t')
for(i in unique(EC[grepl("heart", EC$cell.type),]$cell.type)){
  marker_list = EC[EC$cell.type==i,]$official.gene.symbol
  p=DotPlot(sp_LV, features = marker_list, group.by = "RNA_snn_res.1") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title=i)
  print(p)
}

marker_list = c("Mecom", "Kdr", "Myl2", "Mb", "Tnnt2", 
                "Fbln5", "Hey1", "Vwf", "Reln", "Flt4", "Mmrn1", "Nrp2", "Trp53i11")
DotPlot(angii_LV, features = marker_list, group.by = "RNA_snn_res.0.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
marker_list = c("Mecom", "Kdr", "Cxcl12", "Btnl9", "Nrp2", "Vwf", "Hmcn1", "Fbln5", "Sox17",
                "Cd274", "Lfit3", "Rsad", "Gbp2", "Reln", "Flt4", "Mmrn1", "Ptprc")
DotPlot(ss_LV, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
marker_list = c("Mecom", "Fbln5", "Cxcl12", "Btnl9", "Kdr", 
                "Nrp2", "Hmcn1", "Vwf", "Tnnt2", "Mki67", "Reln", "Abcc9")
marker_list = c("Mecom", "Kdr", "Cxcl12", "Btnl9", "Nrp2", "Vwf", "Hmcn1", 
                "Fbln5", "Reln", "Flt4", "Mmrn1", "Mki67", "Abcc9", "Myh6", "Rora", "Lama2")
DotPlot(sp_LV, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
# DimPlot(sp_LV)


angii_HYP = readRDS("mouse.HYP.EC.cluster.rds")
ss_HYP = readRDS("rat.ss.HYP.EC.cluster.rds")
sp_HYP = readRDS("rat.sp.HYP.EC.cluster.rds")

marker_list = c("Mgp", "Cytl1", "Fbln5", "Bmx", "S100a6", "Azin1", "Pi16",
                "Gkn3", "Hey1", "Edn3", "Glul", "Slc26a10", "Cxcl12", "Spock2",
                "Car4", "Itm2a", "Tmbs10", "Icam1", "Vwf", "Ackr1", "Isg15",
                "Ifit1", "Ifit3", "Ifit3b", "Plvap", "Plpp3", "Kdr")
FeaturePlot(angii_HYP, features = marker_list)
DotPlot(angii_HYP, features = marker_list, group.by = "RNA_snn_res.0.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
marker_list = c("Wnt5b", "Arhgap10", "Glul", "Zeb2", 
                "Uchl1", "Clip1", "Spock2", "Vwf", "Slc7a1",
                "Apod", "Hs3st1", "Ptn", "Fth1",
                "Kdr", "Plpp1", "Plpp3")
FeaturePlot(ss_HYP, features = marker_list)
DotPlot(ss_HYP, features = marker_list, group.by = "RNA_snn_res.0.3") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")

marker_list = c("Wnt5b", "Arhgap10", "Glul", "Zeb2", 
                "Uchl1", "Clip1", "Spock2", "Vwf", "Kdr", "Plpp1", "Plpp3")
FeaturePlot(sp_HYP, features = marker_list)
DotPlot(sp_HYP, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")






angii_LK = readRDS("mouse.LK.EC.cluster.rds")
ss_LK = readRDS("rat.ss.LK.EC.cluster.rds")
sp_LK = readRDS("rat.sp.LK.EC.cluster.rds")

EC = read.table("/xdisk/mliang1/qqiu/reference/single_cell_marker_ref/Cell-EC-2020.txt", header=T, sep='\t')
for(i in unique(EC[grepl("kidney", EC$cell.type),]$cell.type)){
  marker_list = EC[EC$cell.type==i,]$official.gene.symbol
  p=DotPlot(sp_LK, features = marker_list, group.by = "wsnn_res.1") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title=i)
  print(p)
}



marker_list = c("Mecom", "Kdr", "Rapgef4", "Plat", "Aqp1",  "Mmrn1", "Ptprc")
DimPlot(sp_LK, label = T, group.by = "wsnn_res.1", reduction="wnn.umap.harmony") +
FeaturePlot(sp_LK, features = marker_list, reduction = "wnn.umap.harmony")
DotPlot(ss_LK, features = marker_list, group.by = "wsnn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(ss_LK, features = marker_list, group.by = "wsnn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(sp_LK, features = marker_list, group.by = "wsnn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")


markers = FindAllMarkers(sp_LK, group.by="wsnn_res.1", only.pos = TRUE, min.pct = 0.25)
marker_list = markers %>%
  mutate(pct.diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice_head(n=5) %>% ungroup() %>%
  select(gene) %>% unique()

umap_dotplot(sp_LK, marker_list, cluster="wsnn_res.1", reduction = "wnn.umap.harmony")





angii_LK = readRDS("mouse.LK.EC.anno.rds")
ss_LK = readRDS("rat.ss.LK.EC.anno.rds")
sp_LK = readRDS("rat.sp.LK.EC.abbo.rds")

marker_list = c("Kdr", "Rapgef4", "Nlgn1", "Tmtc1", "Sema3d", "Unc5c", "Sema3a", "Sgcd", "Mecom", "Aqp1", "Scin", "Sox17", 
                "Plat", "Ehd3", "Nrp1", "Smad6", "Nostrin", "Mmrn1", "Tbx1", "Pkhd1l1", "Prox1",
                "Rbpms", "Plscr2.1", "Nf1", "Mt-nd4l", "Mt-nd5", "Mt-nd1", "Mt-co2")
DotPlot(angii_LK, features = marker_list, group.by = "subclass_level2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = "AngII (LK)") 
DotPlot(ss_LK, features = marker_list, group.by = "subclass_level2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = "Salt-sensitive (LK)")
DotPlot(sp_LK, features = marker_list, group.by = "subclass_level2") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = "Spontaneous (LK)")







angii_HYP = readRDS("mouse.HYP.EC.cluster.rds")
ss_HYP = readRDS("rat.ss.HYP.EC.cluster.rds")
sp_HYP = readRDS("rat.sp.HYP.EC.cluster.rds")

EC = read.table("/xdisk/mliang1/qqiu/reference/single_cell_marker_ref/Cell-EC-2020.txt", header=T, sep='\t')
for(i in unique(EC[grepl("brain", EC$cell.type),]$cell.type)){
  marker_list = EC[EC$cell.type==i,]$official.gene.symbol
  p=DotPlot(sp_HYP, features = marker_list, group.by = "RNA_snn_res.1") + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title=i)
  print(p)
}

DimPlot(sp_HYP, label = T, group.by = "RNA_snn_res.0.3") +
  FeaturePlot(sp_HYP, features = c("Mecom", "Kdr"))

DimPlot(ss_HYP, label = T, group.by = "RNA_snn_res.1") +
  # FeaturePlot(ss_HYP, features = "Emcn") +
  DotPlot(ss_HYP, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")



marker_list = c("Syt1", "Mecom", "Kdr", "Vwf", "Flt1", "Pecam1", "Lgr5",
                "Lama2", "Zeb2")
DotPlot(angii_HYP, features = marker_list, group.by = "RNA_snn_res.0.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(ss_HYP, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(sp_HYP, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")









angii_MCA = readRDS("mouse.MCA.EC.cluster.rds")
ss_MCA = readRDS("rat.ss.MCA.EC.cluster.rds")
sp_MCA = readRDS("rat.sp.MCA.EC.cluster.rds")

marker_list = c("Syt1", "Mecom", "Flt1", "Pecam1", "Kdr", "Vwf", "Lgr5",
                "Bmx", "Smad6", "Flnb", "Pak7", "Zeb2", "Lama2", "Acta2", "Col8a1", "Dkk2", "Prdm16",
                "Plp1", "Il1r1", "Spock2", "Myh11", "Ranbp3l", "Ptn")
DimPlot(angii_MCA, label = T, group.by = "RNA_snn_res.0.1") +
  FeaturePlot(angii_MCA, features = c("Vwf", "Acta2"))

DimPlot(ss_MCA, label = T, group.by = c("RNA_snn_res.1", "treatment")) +
  FeaturePlot(ss_MCA, features = c("Vwf", "Acta2"))

DimPlot(sp_MCA, label = T, group.by = "RNA_snn_res.1") +
  FeaturePlot(sp_MCA, features = c("Vwf", "Acta2"))

DotPlot(angii_MCA, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(ss_MCA, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(sp_MCA, features = marker_list, group.by = "RNA_snn_res.0.5") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")




ss_MSA = readRDS("rat.ss.MSA.EC.cluster.rds")
sp_MSA = readRDS("rat.sp.MSA.EC.cluster.rds")

marker_list = c("Syt1", "Mecom", "Flt1", "Pecam1", "Kdr", "Vwf", "Lgr5",
                "Bmx", "Smad6", "Flnb", "Pak7", "Zeb2", "Lama2", "Acta2", "Col8a1", "Dkk2", "Prdm16",
                "Plp1", "Il1r1", "Spock2", "Myh11", "Ranbp3l", "Ptn")
marker_list = c("Syt1", "Mecom", "Flt1", "Pecam1", "Kdr", "Vwf", "Lgr5", "Rnf220", "Apba1", "Nrxn3", "Adgrg6", "Emcn", "Ccdc80",
                "Bmx", "Smad6", "Flnb", "Pak7", "Zeb2", "Lama2", "Acta2", "Col8a1", "Dkk2", "Fbln5", "Prdm16", "Chrm3",
                "Plp1", "Il1r1", "Spock2", "Myh11", "Ranbp3l", "Ptn", "Plin1", "Pparg", "Daam1", "Tmod2", "Gab1")
DimPlot(ss_MSA, label = T, group.by = "RNA_snn_res.1") +
  FeaturePlot(ss_MSA, features = c("Vwf", "Mecom"))

DimPlot(sp_MSA, label = T, group.by = "RNA_snn_res.0.5") +
  FeaturePlot(sp_MSA, features = c("Vwf", "Mecom"))

DotPlot(ss_MSA, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
DotPlot(sp_MSA, features = marker_list, group.by = "RNA_snn_res.0.5") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")







################################################################################

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

}







################################################################################
### merge DEG results
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*EC.DEG_all*")

deg_merged = c()

for(i in input_file){
  
  dat_tmp = try(read.table(i, header = T), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/EC.DEG.all.out", sep='\t', col.names = T, row.names = F)



### merge DEG results for strain-wise results
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*EC.strain_wise.DEG_all*")

deg_merged = c()

for(i in input_file){
  
  dat_tmp = try(read.table(i, header = T), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/EC.strain_wise.DEG.all.out", sep='\t', col.names = T, row.names = F)






################################################################################

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

so_list = c()
for( i in input_file[1:3] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

ec_mouse_bulk <- AggregateExpression(seurat_object, assays = "RNA", 
                               group.by = c("strain", "treatment", "tissue", "subclass_level2"), 
                               slot = 'counts', return.seurat = T)
ec_mouse_bulk <- NormalizeData(ec_mouse_bulk, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
ec_mouse_bulk <- FindVariableFeatures(ec_mouse_bulk)


so_list = c()
for( i in input_file[4:7] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

ec_ss_bulk <- AggregateExpression(seurat_object, assays = "RNA", 
                                     group.by = c("strain", "treatment", "tissue", "subclass_level2"), 
                                     slot = 'counts', return.seurat = T)

ec_ss_bulk <- NormalizeData(ec_ss_bulk, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
ec_ss_bulk <- FindVariableFeatures(ec_ss_bulk)


so_list = c()
for( i in input_file[8:12] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

ec_sp_bulk <- AggregateExpression(seurat_object, assays = "RNA", 
                                   group.by = c("strain", "treatment", "tissue", "subclass_level2"), 
                                   slot = 'counts', return.seurat = T)

ec_sp_bulk <- NormalizeData(ec_sp_bulk, assay='RNA',normalization.method = "RC",scale.factor = 1000000)
ec_sp_bulk <- FindVariableFeatures(ec_sp_bulk)

# mouse2rat = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2rat.out.txt", header = T, sep = "\t")

gene_use = intersect(intersect(VariableFeatures(ec_mouse_bulk), VariableFeatures(ec_ss_bulk)), VariableFeatures(ec_sp_bulk))

mouse_df = ec_mouse_bulk@assays$RNA@data[gene_use, ]
colnames(mouse_df) = paste0("C57BL/6_", colnames(mouse_df))
ec_bulk_merged = cbind(mouse_df,
                       ec_ss_bulk@assays$RNA@data[gene_use, ],
                       ec_sp_bulk@assays$RNA@data[gene_use, ])


cor_matrix = cor(as.matrix(ec_bulk_merged))
cor_matrix[is.na(cor_matrix)] <- 0

anno_df = data.frame(cell_group = sapply(colnames(cor_matrix), as.character))
anno_df[, c("strain", "treatment", "tissue", "cell_type1", "cell_type2")] = do.call(rbind, strsplit(anno_df$cell_group, "_"))
anno_df$cell_type = ifelse(anno_df$strain==anno_df$cell_type2, anno_df$cell_type1, paste0(anno_df$cell_type1, "_", anno_df$cell_type2))

top_annotation <- HeatmapAnnotation(
  Strain = anno_df$strain,
  Tissue = anno_df$tissue,
  Treatment = anno_df$treatment,
  col = list(
    Strain = species_col,
    Tissue = tissue_col
  )
)

heatmap <- Heatmap(cor_matrix, 
                   name = "Similarity",
                   # col = colorRamp2(c(0, 1), c("white", "red")),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_columns = "euclidean",
                   clustering_method_rows = "complete",
                   clustering_method_columns = "complete",
                   top_annotation = top_annotation,
                   column_title_rot = 45,
                   row_title_rot = 360,
                   heatmap_legend_param = list(title = "Pearson\ncorrelation", legend_direction = "vertical"))

# Draw heatmap
draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")







################################################################################


input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

so_list = c()
for( i in input_file[1:3] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
vg_mouse <- FindVariableFeatures(seurat_object)


so_list = c()
for( i in input_file[4:7] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
vg_ss <- FindVariableFeatures(seurat_object)

so_list = c()
for( i in input_file[8:12] ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

seurat_object = merge(so_list[[1]], y = so_list[2:length(so_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
vg_sp <- FindVariableFeatures(seurat_object)

vg_use = intersect(intersect(VariableFeatures(vg_mouse), VariableFeatures(vg_ss)), VariableFeatures(vg_sp))


ec_list = c()
for( i in input_file[1:3]){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  assign(name_tmp, so_tmp)
  ec_list = c(ec_list, get(name_tmp))
  
}

seurat_object = merge(ec_list[[1]], y = ec_list[2:length(ec_list)])
seurat_object$sxtxt = paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)
seurat_object = subset(seurat_object, subclass_level2 %in% setdiff(unique(seurat_object$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
seurat_object = NormalizeData(seurat_object)
seurat_object = ScaleData(seurat_object)
seurat_object = RunPCA(seurat_object, features =vg_use, npcs = 20)
seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident",project.dim = F)
seurat_object = RunUMAP(seurat_object, reduction = "harmony", dims=1:20, return.model = T)


DimPlot(seurat_object, label = T, group.by = "sxtxt")
  FeaturePlot(angii_MCA, features = c("Vwf", "Acta2"))

marker_list = c("Syt1", "Mecom", "Flt1", "Pecam1", "Kdr", "Vwf", "Lgr5",
                "Bmx", "Smad6", "Flnb", "Pak7", "Zeb2", "Lama2", "Acta2", "Col8a1", "Dkk2", "Prdm16",
                "Plp1", "Il1r1", "Spock2", "Myh11", "Ranbp3l", "Ptn")











input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

ec_damage_genes <- c(  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", # circulatory system process
                       "Arhgap24", "Igf1r", "Pkhd1", # negative regulation of intracellular signal transduction (GO)
                       "Ank3","Cacnb4","Efna5","Erbb4","Ptprd","Sdk1" # cell junction organization
)

so_list = c()
for( i in input_file ){
  
  name_tmp = gsub(".EC.anno.rds", "", basename(i))
  so_tmp = readRDS(i)
  
  so_tmp = subset(so_tmp, subclass_level2 %in% setdiff(unique(so_tmp$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
  so_tmp$sxt = paste0(so_tmp$strain, " - ", so_tmp$treatment)
  
  so_tmp <- AddModuleScore(
    object = so_tmp,
    features = list(ec_damage_genes),
    name = "EC_Deficiency_Score"
  )
  
  reduction = ifelse("wnn.umap.harmony" %in% names(so_tmp@reductions), "wnn.umap.harmony", "umap")
  p1 = DimPlot(so_tmp, label = T, group.by = "subclass_level2", reduction = reduction) + labs(title = unique(name_tmp))
  p2 = FeaturePlot(so_tmp, features = "EC_Deficiency_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = reduction)
  p3 = RidgePlot(so_tmp, features = "EC_Deficiency_Score1", group.by = "subclass_level2")
  
  # p4 = so_tmp %>% subset(., subclass_level2 %in% c("EC_Mecom")) %>% 
  #   VlnPlot(features = "EC_Damage_Score1", group.by = "sxt", pt.size = 0.5)
  
  p5 = so_tmp@meta.data %>%
    group_by(sxt, subclass_level2) %>%               # Group by group and cell type
    dplyr::summarise(cell_count = n()) %>%              # Count the number of cells in each group-cell type combination
    dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each group
    ggplot(aes(x = sxt, y = proportion, fill = subclass_level2)) +  # Plot the results
    geom_bar(stat = "identity", position = "fill") +  # Create a stacked bar plot with proportions
    scale_y_continuous(labels = scales::percent) +    # Format the y-axis as percentages
    labs(x = "Group", y = "Proportion of Cell Types", fill = "Cell Type",
         title = "Proportion of Cell Types in Different Groups") +  # Add labels and title
    theme_minimal()
  
  combined_plot <- p1 + p2 + p3 + p5 +
    patchwork::plot_layout(ncol = 1)
  
  print(combined_plot)
  
  assign(name_tmp, so_tmp)
  so_list = c(so_list, get(name_tmp))
  
}

for( i in so_list){
  
  so_tmp = i
  name = paste(unique(so_tmp$strain), unique(so_tmp$tissue))
  reduction = ifelse("wnn.umap.harmony" %in% names(so_tmp@reductions), "wnn.umap.harmony", "umap")
  p1 = FeaturePlot(so_tmp, features = "Mecom", cols = c("lightblue", "red"), pt.size = 0.5, reduction = reduction) + labs(title = paste0(name, "\nMecom"))
  p2 = FeaturePlot(so_tmp, features = "EC_Deficiency_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = reduction)
  
  print(p1+p2)
  
}





