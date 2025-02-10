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

# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/panglao_anno.R")
# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DE_enrichment.R")
# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")



## process for VLMC
mouse.HYP = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.rds")
rat.sp.HYP = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds")
mouse.MCA = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds")
rat.ss.MCA = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds")
rat.sp.MCA = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds")


mouse.MCA.deg = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.MCA.reso_0.5.allmarker.0.25.long.txt")
rat.ss.MCA.deg = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.MCA.reso_0.4.allmarker.0.25.long.txt")
rat.sp.MCA.deg = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.MCA.reso_0.5.allmarker.0.25.long.txt")

mouse.MCA.c9 = mouse.MCA.deg[mouse.MCA.deg$cluster==9 & mouse.MCA.deg$p_val_adj<0.05 & mouse.MCA.deg$avg_log2FC>1, ]$gene
mouse.MCA.c20 = mouse.MCA.deg[mouse.MCA.deg$cluster==20 & mouse.MCA.deg$p_val_adj<0.05 & mouse.MCA.deg$avg_log2FC>1, ]$gene
rat.ss.MCA.c4 = rat.ss.MCA.deg[rat.ss.MCA.deg$cluster==4 & rat.ss.MCA.deg$p_val_adj<0.05 & rat.ss.MCA.deg$avg_log2FC>1, ]$gene
rat.ss.MCA.c5 = rat.ss.MCA.deg[rat.ss.MCA.deg$cluster==5 & rat.ss.MCA.deg$p_val_adj<0.05 & rat.ss.MCA.deg$avg_log2FC>1, ]$gene
rat.ss.MCA.c14 = rat.ss.MCA.deg[rat.ss.MCA.deg$cluster==14 & rat.ss.MCA.deg$p_val_adj<0.05 & rat.ss.MCA.deg$avg_log2FC>1, ]$gene
rat.ss.MCA.c13 = rat.ss.MCA.deg[rat.ss.MCA.deg$cluster==13 & rat.ss.MCA.deg$p_val_adj<0.05 & rat.ss.MCA.deg$avg_log2FC>1, ]$gene
rat.sp.MCA.c0 = rat.sp.MCA.deg[rat.sp.MCA.deg$cluster==0 & rat.sp.MCA.deg$p_val_adj<0.05 & rat.sp.MCA.deg$avg_log2FC>1, ]$gene
rat.sp.MCA.c4 = rat.sp.MCA.deg[rat.sp.MCA.deg$cluster==4 & rat.sp.MCA.deg$p_val_adj<0.05 & rat.sp.MCA.deg$avg_log2FC>1, ]$gene
rat.sp.MCA.c6 = rat.sp.MCA.deg[rat.sp.MCA.deg$cluster==6 & rat.sp.MCA.deg$p_val_adj<0.05 & rat.sp.MCA.deg$avg_log2FC>1, ]$gene
rat.sp.MCA.c7 = rat.sp.MCA.deg[rat.sp.MCA.deg$cluster==7 & rat.sp.MCA.deg$p_val_adj<0.05 & rat.sp.MCA.deg$avg_log2FC>1, ]$gene

marker_list = c("Bnc2", "Fbxl7", "Slc4a10", "Eya2", "Slc47a1", "Slc38a2", "Adamtsl3", "Col25a1"
)


umap_dotplot(mouse.HYP, c("Dcn", "Vtn", "Lama2", "Ptgds", "Ranbp3l"), cluster="RNA_snn_res.0.9", reduction = "umap")
umap_dotplot(rat.sp.HYP, c("Dcn", "Vtn", "Lama2", "Ptgds", "Ranbp3l"), cluster="RNA_snn_res.1", reduction = "umap")
umap_dotplot(mouse.MCA, c("Dcn", "Vtn", "Lama2", "Ptgds", "Ranbp3l"), cluster="RNA_snn_res.0.5", reduction = "umap")
umap_dotplot(rat.ss.MCA, c("Dcn", "Vtn", "Lama2", "Ptgds", "Ranbp3l"), cluster="RNA_snn_res.0.4", reduction = "umap")
umap_dotplot(rat.sp.MCA, c("Dcn", "Vtn", "Lama2", "Ptgds", "Ranbp3l"), cluster="RNA_snn_res.0.5", reduction = "umap")



## process for MSA
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.MSA.merged.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.MSA.merged.rds"


infile = ss_infile
seurat_ss = readRDS(infile)

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.4)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("rat.ss.MSA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = paste0("../DEG/", outfile))
}

saveRDS(seurat_object, paste0("../cluster/", "rat.ss.MSA.RNA.cluster.rds"))

seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds")
marker_list = c("Myh11", "Adipoq", "Vwf", "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4")
stacked_barplot(seurat_object@meta.data, species="species", cluster="RNA_snn_res.0.4")
DA_plot(seurat_object@meta.data, outfile="rat.ss.MSA.proportion_change.pdf", cluster="RNA_snn_res.0.4")




infile = sp_infile
seurat_sp = readRDS(infile)

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.4)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("rat.sp.MSA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}

saveRDS(seurat_object, "rat.sp.MSA.RNA.cluster.rds")

seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds")
marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4")
stacked_barplot(seurat_object@meta.data, species="species", cluster="RNA_snn_res.0.4")
DA_plot(seurat_object@meta.data, outfile="rat.sp.MSA.proportion_change.pdf", cluster="RNA_snn_res.0.4")





tmp=as.data.frame(table(seurat_sp$RNA_snn_res.0.4))
tmp[order(as.numeric(as.character(tmp$Var1))),]






################################################################################
merged_levels = c("ECs", "VSMCs", "FIBs", "Adipocytes", 
                  "Immune cells", "Neurons")



seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds")
seurat_object = subset(seurat_object, orig.ident!="RMSA4SN")

Idents(seurat_object) = "RNA_snn_res.0.4"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.4
seurat_object = subset(seurat_object, seurat_clusters %in% c(0:8))
seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)

new.cluster.ids_umap <- c("VSMCs", "VSMCs", "Adipocytes", "ECs", "FIBs", 
                          "ECs", "Neurons", "Immune cells", "ECs")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$seurat_clusters, c(0:8))],
                                                        levels=merged_levels)
# seurat_object$project = factor(c("AngII", "Spontaneous", "Spontaneous", "Salt-sensitive", "Salt-sensitive")[
#   match(seurat_object$species, c("C57BL/6", "WKY", "SHR", "SD", "SS"))
# ], levels=c("AngII", "Salt-sensitive", "Spontaneous"))

marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster = "new.cluster.ids_umap")
stacked_barplot(seurat_object@meta.data, species="species")

saveRDS(seurat_object, "rat.sp.MSA.RNA.cluster.final.rds")
write.table(seurat_object@meta.data, "rat.sp.MSA.RNA.cluster.final.meta.out")



seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds")
Idents(seurat_object) = "RNA_snn_res.0.4"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.4
seurat_object = subset(seurat_object, seurat_clusters %in% c(0:8, 11))
seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)

new.cluster.ids_umap <- c("VSMCs", "VSMCs", "Adipocytes", "FIBs", "ECs", 
                          "Neurons", "ECs", "Immune cells", "VSMCs", "VSMCs")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$seurat_clusters, c(0:8, 11))],
                                                        levels=merged_levels)
# seurat_object$project = factor(c("AngII", "Spontaneous", "Spontaneous", "Salt-sensitive", "Salt-sensitive")[
#   match(seurat_object$species, c("C57BL/6", "WKY", "SHR", "SD", "SS"))
# ], levels=c("AngII", "Salt-sensitive", "Spontaneous"))

marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster = "new.cluster.ids_umap")
stacked_barplot(seurat_object@meta.data, species="species")

saveRDS(seurat_object, "rat.ss.MSA.RNA.cluster.final.rds")
write.table(seurat_object@meta.data, "rat.ss.MSA.RNA.cluster.final.meta.out")

















seurat_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds")
seurat_sp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds")


seurat_object = seurat_ss
marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4", reduction = "umap")


seurat_object = seurat_sp
marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4", reduction = "umap")





seurat_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.ss.HYP.RNA.cluster.rds")
seurat_ss = FindClusters(seurat_ss, resolution = c(1.1, 1.2))
saveRDS(seurat_ss, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds")


seurat_ang = readRDS("/xdisk/mliang1/qqiu/project/multiomics/cluster/mouse/mouse.HYP.RNA.cluster.rds")
seurat_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds")
seurat_sp = readRDS("/xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.sp.HYP.RNA.cluster.rds")


seurat_object = seurat_ang
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", "Avp",
                "Slc1a2", "Tgfbr1", "Cspg4", "Pdgfra","Mbp", "St18", 
                "Col23a1", "Tmem212", "Flt1", "Ebf1",
                "Atp13a5", "Rgs5", "Ptgds", "Slc47a1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.9", reduction = "umap")

seurat_object = seurat_ss
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Tgfbr1", "Arhgap15", "Lyn", "Cd74",
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Lama2")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.1.2", reduction = "umap")


seurat_object = seurat_sp
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Tgfbr1", "Arhgap15", "Lyn", 
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Atp13a5", "Lama2", "Ptgds")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.1", reduction = "umap")







ang_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds"
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds"


seurat_ang = readRDS(ang_infile)
seurat_ss = readRDS(ss_infile)
seurat_sp = readRDS(sp_infile)

seurat_object = seurat_ang
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5", "Ebf1", "Atp13a5", "Lama2",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4", "Ptgds",
                "Ptprc", "Mrc1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.5", reduction = "umap")


seurat_object = seurat_ss
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5", "Grm8",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4",
                "Ptprc", "Mrc1",
                "Lama2", "Ptgds")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4", reduction = "umap")



seurat_object = seurat_sp
marker_list = c("Syt1", "Slc1a2", "Mbp", "Plp1", "St18","Pdgfra", "Cspg4",
                "Flt1", "Vwf", "Mecom", "Pecam1", 
                "Myh11", 
                "Fbln5", "Pdgfrb", "Rgs5", "Ptgds", "Dcn", 
                # "Reln", "Prox1", 
                "Ptprc", "Tgfbr1", "Cx3cr1", "Mrc1", "Grm8")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.5", reduction = "umap")















################################################################################
## process for MCA
ang_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds"
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds"


infile = ang_infile
seurat_ang = readRDS(infile)
seurat_object = seurat_ang

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.5)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("rat.ss.MCA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = paste0("../DEG/", outfile))
}

saveRDS(seurat_object, paste0("../cluster/", "rat.ss.MCA.RNA.cluster.rds"))

marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Ptgds", "Lama2",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1",
                "Cspg4", "Pdgfra", "Pdgfrb", "Ebf1", "Atp13a5",
                "Ptprc", "Mrc1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.5", reduction = "umap")
FeaturePlot(seurat_object, reduction = "umap", features = "percent.mt")
stacked_barplot(seurat_object@meta.data, species="species", cluster="RNA_snn_res.0.4")
DA_plot(seurat_object@meta.data, outfile="rat.ss.MCA.proportion_change.pdf", cluster="RNA_snn_res.0.4")





infile = ss_infile
seurat_ss = readRDS(infile)
seurat_object = seurat_ss

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.4)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("rat.ss.MCA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = paste0("../DEG/", outfile))
}

saveRDS(seurat_object, paste0("../cluster/", "rat.ss.MCA.RNA.cluster.rds"))

seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds")
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", "Ptgds", "Lama2", "Ranbp3l",
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4",
                "Ptprc", "Mrc1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4", reduction = "umap")
stacked_barplot(seurat_object@meta.data, species="species", cluster="RNA_snn_res.0.4")
DA_plot(seurat_object@meta.data, outfile="rat.ss.MCA.proportion_change.pdf", cluster="RNA_snn_res.0.4")






ang_marker = read.table("../DEG/mouse.MCA.reso_0.5.allmarker.0.25.wide.txt", header = T)
ss_marker = read.table("../DEG/rat.ss.MCA.reso_0.4.allmarker.0.25.wide.txt", header = T)
sp_marker = read.table("../DEG/rat.sp.MCA.reso_0.5.allmarker.0.25.wide.txt", header = T)




infile = sp_infile
seurat_sp = readRDS(infile)
seurat_object = seurat_sp

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.5)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("rat.sp.MCA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}

saveRDS(seurat_object, "rat.sp.MCA.RNA.cluster.rds")




tmp=as.data.frame(table(seurat_object$RNA_snn_res.0.5))
tmp[order(as.numeric(as.character(tmp$Var1))),]







seurat_object = seurat_ang
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4",
                "Lama2", "Ptgds",
                "Ptprc", "Mrc1")
marker_list = head(ang_marker$Cluster.12, n=20)
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.5", reduction = "umap")
FeaturePlot(seurat_object, reduction = "umap", features = "percent.mt")
stacked_barplot(seurat_object@meta.data, species="species", cluster="wsnn_res.0.5")
DA_plot(seurat_object@meta.data, outfile="rat.sp.MSA.proportion_change.pdf", cluster="wsnn_res.0.6")





seurat_object = seurat_ss
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5", "Atp13a5", "Ranbp3l",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4",
                "Lama2", "Ptgds",
                "Ptprc", "Mrc1")
marker_list = head(ss_marker$Cluster.15, n=20)
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4", reduction = "umap")





seurat_object = seurat_sp
marker_list = c("Syt1", "Mbp", "Plp1", "St18", "Slc1a2", "Tgfbr1", 
                "Cx3cr1", "Flt1", "Vwf", "Myh11", "Fbln5", "Pdgfrb", "Rgs5", "Atp13a5", "Ranbp3l",
                "Mecom", "Pecam1", "Dcn", "Reln", "Prox1", "Pdgfra", "Cspg4",
                "Col1a1", "Col1a2", "Lum", "Fbln1",
                "Lama2", "Ptgds",
                "Ptprc", "Mrc1", "Grm8")
marker_list = head(sp_marker$Cluster.14, n=20)
marker_list = Reduce(intersect, list(sp_marker$Cluster.0, sp_marker$Cluster.3, 
                                     sp_marker$Cluster.6, sp_marker$Cluster.7))
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.5", reduction = "umap")












infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds"
seurat_object = readRDS(infile)

seurat_object = QC_harmony(seurat_object)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.3)){
  cluster = paste0("RNA_snn_res.", i)
  outfile = paste0("mouse.MCA.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}

saveRDS(seurat_object, "mouse.MCA.RNA.cluster.rds")









seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds")
marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.4")
stacked_barplot(seurat_object@meta.data, species="species", cluster="RNA_snn_res.0.4")
DA_plot(seurat_object@meta.data, outfile="rat.sp.MCA.proportion_change.pdf", cluster="RNA_snn_res.0.4")
















infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.cluster.rds"
seurat_ang = readRDS(infile)

infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.cluster.rds"
seurat_ss = readRDS(infile)

infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.cluster.rds"
seurat_sp = readRDS(infile)


tmp=as.data.frame(table(seurat_sp$RNA_snn_res.0.1))
tmp[order(as.numeric(as.character(tmp$Var1))),]

seurat_object = seurat_ang
DefaultAssay(seurat_object) = "RNA"
marker_list = c("Flt1", "Pecam1", "Fhl2", "Ttn", "Ptprc", "Dcn", "Pdgfrb", "Rgs5")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.1", reduction = "umap")


seurat_object = seurat_sp
DefaultAssay(seurat_object) = "RNA"
marker_list = c("Vwf", "Pecam1", "Dcn", "Fhl2", "Ttn", "Pdgfrb", "Rgs5", "Ptprc", "Nrxn1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.1", reduction = "umap")


seurat_object = seurat_ss
DefaultAssay(seurat_object) = "RNA"
marker_list = c("Vwf", "Pecam1", "Dcn", "Fhl2", "Ttn", "Pdgfrb", "Rgs5", "Ptprc", "Flt4", "Prox1", "Nrxn1")
umap_dotplot(seurat_object, marker_list, cluster="RNA_snn_res.0.1", reduction = "umap")



































































################################################################################
merged_levels = c("ECs", "VSMCs", "MSCs", "Adipocytes", 
                  "Immune cells", "Neurons")



seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds")
seurat_object = subset(seurat_object, orig.ident!="RMSA4SN")

Idents(seurat_object) = "RNA_snn_res.0.4"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.4
seurat_object = subset(seurat_object, seurat_clusters %in% c(0:8))
seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)

new.cluster.ids_umap <- c("VSMCs", "VSMCs", "Adipocytes", "ECs", "MSCs", 
                          "ECs", "Neurons", "Immune cells", "ECs")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$seurat_clusters, c(0:8))],
                                                        levels=merged_levels)
# seurat_object$project = factor(c("AngII", "Spontaneous", "Spontaneous", "Salt-sensitive", "Salt-sensitive")[
#   match(seurat_object$species, c("C57BL/6", "WKY", "SHR", "SD", "SS"))
# ], levels=c("AngII", "Salt-sensitive", "Spontaneous"))

marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster = "new.cluster.ids_umap")
stacked_barplot(seurat_object@meta.data, species="species")

saveRDS(seurat_object, "rat.sp.MCA.RNA.cluster.final.rds")
write.table(seurat_object@meta.data, "rat.sp.MSA.RNA.cluster.final.meta.out")








seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds")
Idents(seurat_object) = "RNA_snn_res.0.4"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.4
seurat_object = subset(seurat_object, seurat_clusters %in% c(0:8, 11))
seurat_object$seurat_clusters = droplevels(seurat_object$seurat_clusters)

new.cluster.ids_umap <- c("VSMCs", "VSMCs", "Adipocytes", "MSCs", "ECs", 
                          "Neurons", "ECs", "Immune cells", "VSMCs", "VSMCs")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$seurat_clusters, c(0:8, 11))],
                                                        levels=merged_levels)
# seurat_object$project = factor(c("AngII", "Spontaneous", "Spontaneous", "Salt-sensitive", "Salt-sensitive")[
#   match(seurat_object$species, c("C57BL/6", "WKY", "SHR", "SD", "SS"))
# ], levels=c("AngII", "Salt-sensitive", "Spontaneous"))

marker_list = c("Myh11", # SMC
                "Adipoq", # adipocyte
                "Vwf", "Mecom", "Pecam1", # endothelial cells
                "Dcn", # mesenchymal stromal cells
                "Reln", "Prox1", # neuron
                "Ptprc")
umap_dotplot(seurat_object, marker_list, cluster = "new.cluster.ids_umap")
stacked_barplot(seurat_object@meta.data, species="species")

saveRDS(seurat_object, "rat.ss.MSA.RNA.cluster.final.rds")
write.table(seurat_object@meta.data, "rat.ss.MSA.RNA.cluster.final.meta.out")













































################################################################################
# cluster resolution evaluation by clustree and ROGUE
angii_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds"
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds"

infile = angii_infile
seurat_object = readRDS(infile)

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.4)){
  cluster = paste0("wsnn_res.", i)
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}



seurat_expr = seurat_obj@assays$RNA@data


seurat_expr = seurat_expr
ent.res = SE_fun(seurat_expr)
rogue.value = CalculateRogue(ent.res, platform = "UMI")
rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)




































DefaultAssay(seurat_object) = "ATAC"
gene.activities = GeneActivity(seurat_object)
seurat_object[['gene.activity']] = CreateAssayObject(counts = gene.activities)
seurat_object = NormalizeData(
  object = seurat_object,
  assay = 'gene.activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat_object$nCount_)
)


seurat_object = readRDS(paste0(species, ".multiomics.rna.cluster.rds"))
markers = FindAllMarkers(seurat_object, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
write.table(markers,file=paste0(outfolder, outfile, ".multiomics.allmarker.0.25.long.txt"), sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                        v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file=paste0(outfolder, outfile, ".multiomics.allmarker.0.25.wide.txt"), sep="\t", row.names = F)

panglao_anno(marker_tbl, paste0(outfolder, outfile, ".multiomics.allmarker.0.25.wide.xlsx"))






















### generate gene activity score using cicero
library(cicero)
library(SeuratWrappers)
library(rtracklayer)

annotations = import("/scratch/u/qqiu/refdata/cellranger-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz")
gene_anno = unique(annotations[annotations$type=="gene" &
                                seqnames(annotations) %in% standardChromosomes(annotations) , c("gene_name")])

promoter = promoters(gene_anno, upstream=1000, downstream=1000)
promoterDF = data.frame(chromosome=seqnames(promoter),
                        start=start(promoter),
                        end=end(promoter),
                        gene=promoter$gene_name) %>% unique()

mm10 = genomeAnnoMm10$chromSizes
mm10 = as.data.frame(mm10)[,c(1,4)]


seurat_object = readRDS("MLK1_4.multiomics.rna.cluster.rds")
DefaultAssay(seurat_object) = "ATAC"

ds_cds = as.cell_data_set(x = seurat_object, assay="ATAC", default.reduction = "wnn.umap.harmony")
ds_cicero = make_cicero_cds(ds_cds, reduced_coordinates = reducedDims(ds_cds)$WNN.UMAP.HARMONY)

ds_conns = run_cicero(ds_cicero, mm10, sample_num = 100) 



input_cds = annotate_cds_by_site(input_cds, gene_annotation_sub)

# generate unnormalized gene activity matrix
unnorm_ga = build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga = unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes = pData(input_cds)$num_genes_expressed
names(num_genes) = row.names(pData(input_cds))

# normalize
cicero_gene_activities = normalize_gene_activities(unnorm_ga, num_genes)

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
unnorm_ga2 = unnorm_ga
cicero_gene_activities = normalize_gene_activities(list(unnorm_ga, unnorm_ga2), 
                                                    num_genes)