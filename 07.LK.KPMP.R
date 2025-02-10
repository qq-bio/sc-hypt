library(Seurat)
library(sctransform)
library(ggplot2)
library(cowplot)


################################################################################
### integration of reference data
# sc_ref = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/data/KPMP.ref.sc.rds")
# sn_ref = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/data/KPMP.ref.sn.rds")
# 
# for(i in c("sc_ref", "sn_ref")){
#   
#   seurat_object = get(i)
#   seurat_object = NormalizeData(seurat_object, verbose = FALSE)
#   seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst", 
#                                        nfeatures = 2000, verbose = FALSE)
#   
#   assign(i, seurat_object)
#   
# }
# 
# 
# reference_list = list(sc_ref, sn_ref)
# reference_anchors = FindIntegrationAnchors(object.list = reference_list, dims = 1:30)
# 
# reference_integrated = IntegrateData(anchorset = reference_anchors, dims = 1:30)
# DefaultAssay(reference_integrated) = "integrated"
# reference_integrated = ScaleData(reference_integrated, verbose = FALSE)
# reference_integrated = RunPCA(reference_integrated, npcs = 30, verbose = FALSE)
# reference_integrated = RunUMAP(reference_integrated, reduction = "pca", dims = 1:30)
# # p1 = DimPlot(reference_integrated, reduction = "umap", group.by = "tech")
# # p2 = DimPlot(reference_integrated, reduction = "umap", group.by = "sc_meta$subclass.l1", label = TRUE, 
# #              repel = TRUE) + NoLegend()
# # plot_grid(p1, p2)
# 
# saveRDS(reference_integrated, "/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.ref.integrated.RDS")

# rm(c(sc_ref, sn_ref))


# ################################################################################
# ### cell classification using integrated reference
# #### reference: https://satijalab.org/seurat/archive/v3.0/integration.html
# 
# # reference_integrated = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.ref.integrated.rds")
# reference_integrated = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/data/KPMP.ref.sc.rds")
# reference_integrated = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/data/KPMP.ref.sn.rds")
# 
# 
# for(i in c("sn", "mo")){ # c("sc", "sn", "mo")
# 
#   input_file = paste0("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.", i, ".merged.cluster.rds")
#   seurat_object = readRDS(input_file)
# 
#   seurat_object = NormalizeData(seurat_object, verbose = FALSE)
#   seurat_object = FindVariableFeatures(seurat_object, selection.method = "vst",
#                                        nfeatures = 2000, verbose = FALSE)
# 
#   reference_anchors = FindTransferAnchors(reference = reference_integrated, query = seurat_object,
#                                           normalization.method = "SCT",
#                                           dims = 1:30, recompute.residuals = FALSE)
#   predictions = TransferData(anchorset = reference_anchors, refdata = reference_integrated$subclass.l3,
#                              dims = 1:30)
#   seurat_object = AddMetaData(seurat_object, metadata = predictions)
# 
#   output = paste0("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.", i, ".query.l3.rds")
#   saveRDS(seurat_object, output)
# 
# }









# ################################################################################
# ### process transferred dataset
# 
# setwd("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster")
# clinic_info = read.table("/xdisk/mliang1/qqiu/data/KPMP/KPMP.Clinical.Tissue_type.csv", sep=',', header=T)
# rownames(clinic_info) = clinic_info$Participant.ID
# 
# for(i in c("sn")){
# 
#   cluster_file = paste0("KPMP.", i, ".merged.cluster.rds")
#   label_file = paste0("KPMP.", i, ".query.l3.label.rds")
# 
#   cluster_object = readRDS(cluster_file)
#   label_object = readRDS(label_file)
# 
#   all(colnames(cluster_object)==rownames(label_object))
#   cluster_object$predicted.id = label_object$predicted.id
# 
#   cluster_object$orig.ident[cluster_object$orig.ident=="Feb-32"] = "32-2"
#   cluster_object$tissue.type = clinic_info[cluster_object$orig.ident, ]$Tissue.Type
# 
# 
# 
#   reso = "SCT_snn_res.1"
# 
#   p1 = DimPlot(cluster_object, reduction = "umap", group.by = reso,
#                pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") +
#     labs(x="UMAP 1", y="UMAP 2")
#   p2 = DimPlot(cluster_object, reduction = "umap", group.by = "predicted.id",
#                pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") +
#     labs(x="UMAP 1", y="UMAP 2")
# 
#   p1 + p2
# 
# 
#   tbl = table(cluster_object@meta.data[, reso], cluster_object$predicted.id)
#   tbl_sum = tbl/rowSums(tbl)*100
#   sum_df = data.frame( cluster = rownames(tbl_sum),
#                        prop = apply(tbl_sum, 1, max),
#                        cell_type = colnames(tbl)[apply(tbl_sum, 1, which.max)] )
# 
#   cluster_object@meta.data$"class.l3" = sum_df[cluster_object@meta.data[, reso], ]$cell_type
#   DimPlot(cluster_object, reduction = "umap", group.by = "class.l3",
#           pt.size = 0.2, raster = F, label = T, repel = T) + theme(legend.position = "none") +
#     labs(x="UMAP 1", y="UMAP 2")
# 
#   DotPlot(cluster_object,
#           features = unique(marker_list), group.by = reso) +
#     theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
# 
#   mod = read.table(paste0("KPMP.", i, ".l3.final.txt"), sep='\t', header=T)
#   rownames(mod) = mod$cluster
#   cluster_object@meta.data$"class.l3" = mod[cluster_object@meta.data[, reso], ]$cell_type_mod
#   cluster_object = subset(cluster_object, class.l3 %in% na.omit(mod$cell_type_mod))
# 
#   marker_list = as.character(unlist(sapply(mod[order(mod$cell_type_mod),]$markers, function(x) strsplit(x, ", "))))
# 
#   DotPlot(cluster_object,
#           features = na.omit(unique(marker_list)), group.by = "class.l3") +
#     theme(axis.text.x = element_text(angle = 45, hjust=1))
# 
#   DimPlot(cluster_object, reduction = "umap", group.by = "class.l3",
#           pt.size = 0.2, raster = F, label = T, repel = T) + theme(legend.position = "none") +
#     labs(x="UMAP 1", y="UMAP 2")
#   
# 
#   # sample_info = cluster_object@meta.data[,c("orig.ident", "tissue.type", "class.l3")]
#   # 
#   # outfile = paste0("KPMP.", i, ".anno.l3.out")
#   # write.table(sample_info, outfile, col.names = T, row.names = T, sep = '\t', quote = F)
#   # # system("sed -i '1s/^/cell.id\t/' KPMP.sn.anno.l3.out")
# }
# 
# 
# 
# 
# 
# 
# 
# marker_list = c("PTPRQ", "WT1", "NTNG1", "NPHS1", "NPHS2", "CLIC5", "PODXL",
#                 "CLDN1", "VCAM1", "CFH", "RBFOX1", "ALDH1A2",
#                 "LRP2", "CUBN", "SLC13A1",
#                 "SLC5A12", "SLC13A3", "SLC22A6", "PRODH2", "SLC5A2", "SLC22A8",
#                 "SLC5A12", "SLC13A3", "SLC22A6", "SLC34A1", "SLC22A7",
#                 "SLC22A7", "MOGAT1", "SLC5A11", "SLC22A24", "SLC7A13", "SLC5A8", "ABCC3", "SATB2",
#                 "CRYAB", "TACSTD2", "SLC44A5", "KLRG2", "COL26A1", "BOC",
#                 "VCAM1", "SLC39A8", "AQP1", "LRRC4C", "LRP2", "UNC5D", "SATB2",
#                 "SATB2", "JAG1", "ADGRL3", "ID1",
#                 "CLDN1", "AKR1B1", "CLDN4", "BCL6", "SH3GL3", "SLC14A2", "SMOC2",
#                 "CLDN1", "AKR1B1", "CLDN4", "BCL6", "SH3GL3", "BCAS1", "CLCNKA", "CLDN10", "PROX1",
#                 "CASR", "SLC12A1", "UMOD",
#                 "NELL1", "ESRRB", "EGF", "CLDN14", "PROX1", "MFSD4A", "KCTD16", "RAP1GAP", "ANK2", "CYFIP2",
#                 "NELL1", "ESRRB", "EGF", "PPM1E", "GP2", "ENOX1", "TMEM207", "TMEM52B", "CLDN16", "WNK1",
#                 "NOS1", "ROBO2", "CALCR", "PPFIA2", "PAPPA2",
#                 "SLC12A3", "CNNM2", "FGF13", "KLHL3", "LHX1", "TRPM6",
#                 "TRPM7", "ADAMTS17", "ITPKB", "ZNF385D", "HS6ST2",
#                 "TRPV5", "SLC8A1", "SCN2A", "HSD11B2", "CALB1",
#                 "SLC8A1", "SCN2A", "HSD11B2", "CALB1",
#                 "KITLG", "PCDH7",
#                 "RALYL", "TOX", "SGPP1", "SCNN1G", "SCNN1B", "KCNIP1",
#                 "GATA3", "AQP2", "AQP3",
#                 "SCNN1G", "SCNN1B", "FXYD4", "SOX5", "PDE10A", "SLC25A29", "ST6GAL1", "PAPPA",
#                 "SCNN1G", "SCNN1B", "FXYD4", "SOX5", "SYK", "FAM81A", "PROM1", "KCNK13",
#                 "FXYD4", "SOX5", "PHACTR1", "PCDH7", "SLC14A2", "HS3ST5",
#                 "TACSTD2", "TP63", "GPX2", "FXYD3", "KRT5",
#                 "ATP6V0D2", "ATP6V1C2", "TMEM213", "CLNK",
#                 "SLC4A1", "SLC26A7", "HS6ST3", "NXPH2", "LEF1", "ADGRF5",
#                 "SLC4A1", "SLC26A7", "SLC8A1", "SCN2A", "CALB1",
#                 "SLC4A1", "SLC26A7", "KIT", "AQP6", "STAP1", "FAM184B", "CALCA",
#                 "SLC4A9", "SLC35F3", "SLC26A4", "INSRR", "TLDC2",
#                 "CD34", "PECAM1", "PTPRB", "MEIS2", "FLT1", "EMCN",
#                 "EMCN", "HECW2", "PLAT", "ITGA8", "EHD3", "KDR", "SOST",
#                 "BTNL9", "ADAMTS6", "PALMD", "AQP1", "TM4SF1", "VEGFC", "CCDC3", "CDH5", "SERPINE2", "FBLN5", "CXCL12", "SOX17",
#                 "BTNL9", "ADAMTS6", "PALMD", "AQP1", "TM4SF1", "MCTP1", "SLC14A1", "ENPP2", "LYPD6B",
#                 "CEACAM1", "DNASE1L3", "PLVAP", "PITPNC1", "GRB10", "SLCO2A1", "RAPGEF4",
#                 "CEACAM1", "DNASE1L3", "PLVAP", "GPM6A", "EDIL3", "TLL1", "ZNF385D", "NR2F2",
#                 "MMRN1", "CD36", "TBX1", "PKHD1L1", "PROX1",
#                 "NOTCH3", "PDGFRB", "ITGA8",
#                 "PIP5K1B", "ROBO1", "PIEZO2", "DAAM2", "PHTF2", "GATA3", "POSTN",
#                 "PIP5K1B", "ROBO1", "REN", "PDE10A", "ABCC8", "COL13A1", "GRID2",
#                 "NTRK3", "MYH11", "RGS6", "ADRA1A", "LDB3", "MCAM",
#                 "NTRK3", "CCDC102B", "RGS5", "ABCC9", "ADCY3", "ADGRB3",
#                 "COL1A1", "COL1A2", "C7", "NEGR1", "FBLN5", "DCN", "CDH11",
#                 "LAMA2", "GGT5", "LUM", "AEBP1", "C1S", "SFRP1", "MEG3", "CXCL12",
#                 "SYT1", "TNC", "PLCXD3", "GABRG3", "GREB1L", "KCNK2",
#                 "SYT1", "TNC", "ADAMTSL1", "CLMP", "NCAM1", "FREM1",
#                 "SYT1", "TNC", "IGFBP5", "MGP", "BGN", "IGFBP2",
#                 "SYNPO2", "PCDH7", "KCNMA1", "LMOD1", "TTLL7", "DTNA", "COL14A1",
#                 "LAMA2", "COL16A1", "SULF1", "COL6A3", "NTM", "GLI2", "COL5A1", "FN1",
#                 "COL16A1", "SULF1", "SLC24A3", "CDH13", "SMOC2", "ANO3", "RXFP1",
#                 "PTPRC",
#                 "BANK1", "BLK", "MS4A1", "BACH2",
#                 "IGKC", "TENT5C", "MZB1", "FCRL5", "CD38", "JCHAIN",
#                 "CD96", "CD247", "THEMIS", "BCL11B", "CAMK4", "IL7R",
#                 "CD96", "CD247", "RUNX3", "GNLY", "NKG7", "CCL5", "KLRF1", "CCL4", "GZMA",
#                 "MS4A2", "CPA3", "KIT",
#                 "F13A1", "MRC1", "CD163", "STAB1", "SLC1A3", "CD14", "FOLR2",
#                 "MSR1", "ITGAX", "HLA-DQA1", "HLA_DRB1", "CSF2RA", "CD14", "TRPM2",
#                 "ITGAX", "HLA-DQA1", "HLA-DRA", "CSF2RA", "CIITA", "WDFY4", "FLT3", "ZNF366", "CADM1", "ZBTB46", "CLEC9A",
#                 "IRF8", "CUX2", "P2RY14", "IL3RA", "CLEC4C",
#                 "CTSS", "IRAK3", "TCF7L2", "TNFRSF1B", "FCN1", "HLA-DRA", "FCGR3A",
#                 "S100A9", "S100A8", "IFITM2", "FCGR3B", "CD1C",
#                 "CDH19", "NRXN1", "GINS3"
# )
# 
# 
# 
# 
# setwd("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster")
# 
# for(i in c("sc", "mo")){ # c("sc", "sn", "mo")
#   
#   cluster_file = paste0("KPMP.", i, ".merged.cluster.rds")
#   label_file = paste0("KPMP.", i, ".anno.l3.out")
#   
#   cluster_object = readRDS(cluster_file)
#   label_object = read.table(label_file, header=T, sep='\t')
#   rownames(label_object) = label_object$cell.id
#   
#   all(colnames(cluster_object)==rownames(label_object))
#   cluster_object$class.l3 = label_object[colnames(cluster_object),]$class.l3
#   cluster_object$tissue.type = label_object[colnames(cluster_object),]$tissue.type
#   
#   reso = "SCT_snn_res.1"
#   
#   p1 = DimPlot(cluster_object, reduction = "umap", group.by = reso,
#                pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") +
#     labs(x="UMAP 1", y="UMAP 2")
#   p2 = DimPlot(cluster_object, reduction = "umap", group.by = "class.l3",
#                pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") +
#     labs(x="UMAP 1", y="UMAP 2")
#   
#   p1 + p2
#   
#   output = paste0("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.", i, ".query.l3.rds")
#   saveRDS(cluster_object, output)
#   
# }



################################################################################
### EC subcluster
dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Matrix)
library(Seurat)
library(sctransform)
library(harmony)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(R.utils)
library(ggsci)
library(presto)
library(ggh4x)

source("/xdisk/mliang1/qqiu/project/multiomics-CKD/src/function/panglao_anno.R")

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

setwd("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster")


for(i in c("sc", "sn", "mo")){
  
  
  infile = paste0("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.", i, ".query.l3.rds")
  seurat_object_orig = readRDS(infile)
  
  DefaultAssay(seurat_object_orig) <- "SCT"
  
  ec_list = unique(seurat_object_orig$class.l3)[grepl("^EC", unique(seurat_object_orig$class.l3))]
  seurat_object = subset(seurat_object_orig, class.l3 %in% ec_list)
  rm(seurat_object_orig)

  # seurat_object = NormalizeData(seurat_object)
  seurat_object = FindVariableFeatures(seurat_object)
  seurat_object = ScaleData(seurat_object)
  seurat_object = RunPCA(seurat_object)
  # seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", project.dim = F)
  
  ElbowPlot(seurat_object, ndims = 50, reduction = "pca") + labs(x="PCA")
  
  seurat_object = RunUMAP(seurat_object, reduction = "pca", dims=1:30)
  seurat_object = FindNeighbors(seurat_object, reduction = "pca", dims=1:30)
  seurat_object = FindClusters(seurat_object, resolution = c(0.1, 0.2, 0.5, 0.8, 1, 2, 4, 5, 10))
  
  
  ## resolution selection
  assay="SCT"
  for(resolution in c(0.5)){
    reso_para = paste0(assay, '_snn_res.', resolution)
    seurat_object@active.ident = seurat_object@meta.data[, reso_para]
    seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]

    p = DimPlot(seurat_object, reduction = "umap", group.by = c(reso_para),
                pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") +
      labs(x="UMAP 1", y="UMAP 2", title = paste0("reso=", resolution), color="Cluster")
    print(p)
  }
  
  # reso_para = paste0('RNA_snn_res.', 0.5)
  # seurat_object@active.ident = seurat_object@meta.data[, reso_para]
  # seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]
  
  outfile = paste0("KPMP.", i, ".EC.cluster.rds")
  saveRDS(seurat_object, outfile)
  
}



for(i in c("sc", "sn", "mo")){
  
  infile = paste0("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.", i, ".EC.cluster.rds")
  seurat_object = readRDS(infile)
  
  DefaultAssay(seurat_object) <- "SCT"
  
  DimPlot(seurat_object, group.by = "class.l3") + FeaturePlot(seurat_object, "MECOM")
  reso = 'SCT_snn_res.1'
  DimPlot(seurat_object, group.by = reso, label = T) + VlnPlot(seurat_object, "MECOM", group.by = reso)
  VlnPlot(seurat_object, "MECOM", group.by = reso, split.by = "tissue.type")
  
  seurat_object@meta.data %>%
    group_by(orig.ident, !!sym(reso), tissue.type) %>%  # Group by sample, cell type, and tissue type
    dplyr::summarise(cell_count = n(), .groups = 'drop') %>%  # Count the number of cells in each group-cell type combination
    ungroup() %>%
    group_by(orig.ident) %>%
    dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each sample
    filter(tissue.type %in% c("Healthy Reference", "AKI", "CKD")) %>%  # Filter relevant tissue types
    dplyr::mutate(tissue.type = factor(tissue.type, levels = c("Healthy Reference", "AKI", "CKD"))) %>%
    group_by(!!sym(reso), tissue.type) %>%  # Group by cell type and tissue type
    dplyr::summarise(
      mean_proportion = mean(proportion),
      sd_proportion = sd(proportion),
      .groups = 'drop'  # Explicitly drop the last grouping
    ) %>% 
    ggplot(aes(x = tissue.type, y = mean_proportion, fill = get(reso))) +  # Plot the results
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Create a bar plot with means
    geom_errorbar(aes(ymin = mean_proportion - sd_proportion, ymax = mean_proportion + sd_proportion),
                  position = position_dodge(width = 0.8), width = 0.2) +  # Add error bars
    scale_y_continuous(labels = scales::percent) +  # Format the y-axis as percentages
    # scale_fill_manual(values = lk_ec_col) +
    labs(x = "Treatment", y = "Mean proportion of cell types") +  # Add labels and title
    theme(    axis.text.y = element_text(colour = 'black'),
              axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
              legend.position = "None"
    ) +
    facet_nested(get(reso) ~ ., scales = "free") 
  
}




seurat_object_sc = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.sc.EC.cluster.rds")
seurat_object_sn = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.sn.EC.cluster.rds")
seurat_object_mo = readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.mo.EC.cluster.rds")

table(unique(seurat_object_sc@meta.data[,c("orig.ident", "tissue.type")])$tissue.type)
table(unique(seurat_object_sn@meta.data[,c("orig.ident", "tissue.type")])$tissue.type)
table(unique(seurat_object_mo@meta.data[,c("orig.ident", "tissue.type")])$tissue.type)

seurat_object = seurat_object_sc
Idents(seurat_object) = "SCT_snn_res.0.5"
seurat_object$cluster = seurat_object$SCT_snn_res.0.5
markers = FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE)
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC.allmarker.0.25.long.txt", sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.25,] %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                               v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC.allmarker.0.25.wide.txt", sep="\t", row.names = F)


FeaturePlot(seurat_object, "MECOM") + theme(plot.title = element_text(face = "plain")) +
  labs(x="UMAP 1", y="UMAP 2", title="MECOM expression\nin human kidney EC")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3f.kpmp.ec.mecom.png", width=338/96, height=285/96, dpi=300)



pathway_res = read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC_Mecom.metascape/metascape_result.xlsx", 2)
pathway_res = pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category_mod = lapply(strsplit(pathway_res$Category, " "), function(x) x[1]) %>% as.character()
pathway_res$pathway = paste0(pathway_res$Description, " (", pathway_res$category_mod, ")")
ggplot(pathway_res, aes(x = -1*`Log(q-value)`, y = fct_reorder(pathway, -1*`Log(q-value)`), fill = -1*`Log(q-value)`)) +  # Plot the results
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  # Add labels and title
  theme(    axis.text = element_text(colour = 'black'),
            # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
            legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_gradient(low = "white", high="red", limits = c(0, 16))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3g.kpmp.mecom_ec.pathway_barplot.png", width=706/96, height=285/96, dpi=300)











# ec_damage_genes <- c(  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", # circulatory system process
#                        "Arhgap24", "Igf1r", "Pkhd1", # negative regulation of intracellular signal transduction (GO)
#                        "Ank3","Cacnb4","Efna5","Erbb4","Ptprd","Sdk1" # cell junction organization
# )
# ec_damage_genes = toupper(ec_damage_genes)
# 
# seurat_object <- AddModuleScore(
#   object = seurat_object,
#   features = list(ec_damage_genes),
#   name = "EC_Deficiency_Score"
# )
# p1 = DimPlot(seurat_object, label = T)
# p2 = FeaturePlot(seurat_object, features = "EC_Deficiency_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(title = "EC deficiency score", x="UMAP 1", y="UMAP 2", color="Score")
# p1 + p2
# 
# VlnPlot(seurat_object, ec_damage_genes)


seurat_object@meta.data %>%
  group_by(orig.ident, cluster, tissue.type) %>%  # Group by sample, cell type, and tissue type
  dplyr::summarise(cell_count = n(), .groups = 'drop') %>%  # Count the number of cells in each group-cell type combination
  ungroup() %>%
  group_by(orig.ident) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each sample
  filter(tissue.type %in% c("Healthy Reference", "AKI", "CKD")) %>%  # Filter relevant tissue types
  dplyr::mutate(tissue.type = factor(tissue.type, levels = c("Healthy Reference", "AKI", "CKD"))) %>%
  group_by(cluster, tissue.type) %>%  # Group by cell type and tissue type
  dplyr::summarise(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion),
    .groups = 'drop'  # Explicitly drop the last grouping
  ) %>% 
  filter(cluster %in% c(10)) %>% 
  ggplot(aes(x = tissue.type, y = mean_proportion, fill = cluster)) +  # Plot the results
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Create a bar plot with means
  geom_errorbar(aes(ymin = mean_proportion - sd_proportion, ymax = mean_proportion + sd_proportion),
                position = position_dodge(width = 0.8), width = 0.2) +  # Add error bars
  scale_y_continuous(labels = scales::percent) +  # Format the y-axis as percentages
  # scale_fill_manual(values = lk_ec_col) +
  labs(x = "", y = "Mean proportion of cell types", title = "MECOM+ EC") +  # Add labels and title
  theme(    axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
            legend.position = "None"
  )



seurat_object@meta.data %>%
  group_by(orig.ident, cluster, tissue.type) %>%  # Group by sample, cell type, and tissue type
  dplyr::summarise(cell_count = n(), .groups = 'drop') %>%  # Count the number of cells in each group-cell type combination
  ungroup() %>%
  group_by(orig.ident) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each sample
  filter(tissue.type %in% c("Healthy Reference", "AKI", "CKD")) %>%  # Filter relevant tissue types
  dplyr::mutate(tissue.type = factor(tissue.type, levels = c("Healthy Reference", "AKI", "CKD"))) %>%
  filter(cluster %in% c(10)) %>% 
  ggplot(aes(x = tissue.type, y = proportion, fill = cluster)) +  # Plot the results
  geom_boxplot() +  # Create a bar plot with means
  scale_y_continuous(labels = scales::percent) +  # Format the y-axis as percentages
  scale_fill_manual(values = "#B1C968") +
  labs(x = "", y = "Proportion", title = "MECOM+ EC\n(KPMP)") +  # Add labels and title
  theme(    axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
            legend.position = "None"
  ) +
  stat_compare_means(method = "wilcox.test", label = "p.format", 
                     comparisons = list(c("Healthy Reference", "AKI"), 
                                        c("Healthy Reference", "CKD"), 
                                        c("AKI", "CKD")))  # Pairwise comparison

