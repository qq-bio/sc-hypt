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
library(R.utils)
library(ggsci)
library(presto)
library(clustree)
library(ROGUE)
library(UpSetR)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/panglao_anno.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/cluster_sample_sum.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")

# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.cluster.rds

# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds


infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds"
seurat_ang = readRDS(infile)
seurat_object = seurat_ang
DefaultAssay(seurat_object) = "RNA"

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
for(i in c(0.5)){
  cluster = paste0("wsnn_res.", i)
  outfile = paste0("mouse.LK.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}



infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds"
seurat_ss = readRDS(infile)
seurat_object = seurat_ss
DefaultAssay(seurat_object) = "RNA"

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
for(i in c(0.4)){
  cluster = paste0("wsnn_res.", i)
  outfile = paste0("rat.ss.LK.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}



infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds"
seurat_sp = readRDS(infile)
seurat_object = seurat_sp
DefaultAssay(seurat_object) = "RNA"

seurat_meta = seurat_object@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
for(i in c(0.4)){
  cluster = paste0("wsnn_res.", i)
  outfile = paste0("rat.sp.LK.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}

# seurat_expr = seurat_object@assays$RNA@data
# ent.res = SE_fun(seurat_expr)
# rogue.value = CalculateRogue(ent.res, platform = "UMI")
# rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)






tmp=as.data.frame(table(seurat_object$wsnn_res.0.4))
tmp[order(as.numeric(as.character(tmp$Var1))),]














seurat_object = seurat_ang
DefaultAssay(seurat_object) = "RNA"
marker_list = c("Nphs1", "Nphs2",
                "Cryab", 
                "Lrp2", "Slc5a12", "Slc22a6", "Slc7a13", "Slc13a3", "Slc16a9",
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Flt1", 
                "Fbln5", "Pdgfra", "Myh11", "Notch3", "Ren1", "Robo1", 
                "Ncam1", 
                "Ptprc" 
)
umap_dotplot(seurat_object, marker_list, cluster="wsnn_res.0.5", reduction = "wnn.umap.harmony")
FeaturePlot(seurat_object, reduction = "wnn.umap.harmony", features = "percent.mt")
stacked_barplot(seurat_object@meta.data, species="species", cluster="wsnn_res.0.5")
DA_plot(seurat_object@meta.data, outfile="rat.sp.MSA.proportion_change.pdf", cluster="wsnn_res.0.6")





seurat_object = seurat_sp
DefaultAssay(seurat_object) = "RNA"
marker_list = c("Nphs1", "Nphs2",
                "Cryab", # "Rbfox1", # "Cldn1", "Vcam1", 
                "Lrp2", "Slc5a12", "Slc22a6", "Slc7a13", # "Slc7a12", "Havcr1", "Vcam1", "Mki67",
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Pecam1", 
                "Cald1", "Fbln5", "Pdgfra", "Myh11", "Notch3", "Robo1", "Renin", # "Notch3", "Postn",
                "Ncam1", 
                "Ptprc" # "Grip1", "Dcn"
)
umap_dotplot(seurat_object, marker_list, cluster="wsnn_res.0.4", reduction = "wnn.umap.harmony")







seurat_object = seurat_ss
DefaultAssay(seurat_object) = "RNA"
# seurat_object@active.ident <- factor(seurat_object@active.ident)
# names(seurat_object@active.ident) <- colnames(seurat_object)
# seurat_object_tmp = subset(seurat_object, DoubletScore==0 & doublet_pc.optm=="Singlet" & DoubletEnrichment<1)
marker_list = c("Nphs1", "Nphs2",
                "Cryab", # "Rbfox1", # "Cldn1", "Vcam1", 
                "Lrp2", "Slc5a12", "Slc22a6", "Slc7a13", # "Slc7a12", "Havcr1", "Vcam1", "Mki67",
                "Slc12a1", "Slc12a3", "Slc8a1", "Aqp2", 
                "Atp6v0d2", "Kit", "Slc26a4", "Pecam1", 
                "Cald1", "Fbln5", "Pdgfra", "Myh11", "Robo1", "Renin", # "Notch3", "Postn",
                "Ncam1", 
                "Ptprc" 
)
umap_dotplot(seurat_object, marker_list, cluster="wsnn_res.0.4", reduction = "wnn.umap.harmony")





gene = "Ncam1"
p1=FeaturePlot(seurat_ang, reduction = "wnn.umap.harmony", features = gene)
p2=FeaturePlot(seurat_ss, reduction = "wnn.umap.harmony", features = gene)
p3=FeaturePlot(seurat_sp, reduction = "wnn.umap.harmony", features = gene)
p1 + p2 + p3


fib_list = c("C7", "Fbln5", "Mfap4", "Lum", "Fmo2", "Col5a1")
vsmc_list = c("Myh11", "Robo1", "Notch3", "Itga8", "Daam2", "Gata3", "Piezo2", "Pip5k1b")
pec_list = c("Cldn1", "Vcam1", "Cfh", "Rbfox1", "Aldh1a2")
dtl_list = c("Cryab", "Tacstd2", "Slc44a5", "Klrg2", "Col26a1", "Boc")
pc_list = c("Gata3", "Aqp2", "Scnn1g", "Scnn1b" ,"Fxyd4", "Sox5")
select_list = c("Tafa1", "Nav3", "Epha7")
gene=pc_list
p1=DotPlot(seurat_ang, features = gene, group.by = "wsnn_res.0.5") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p2=DotPlot(seurat_ss, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p3=DotPlot(seurat_sp, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p1 + p2 + p3




ang_marker = read.table("../DEG/mouse.LK.reso_0.5.allmarker.0.25.wide.txt", header = T)
ss_marker = read.table("../DEG/rat.ss.LK.reso_0.4.allmarker.0.25.wide.txt", header = T)
sp_marker = read.table("../DEG/rat.sp.LK.reso_0.4.allmarker.0.25.wide.txt", header = T)


input = list(ang_c15=unique(ang_marker$Cluster.15), ang_c19=unique(ang_marker$Cluster.19),
             sp_c13=unique(sp_marker$Cluster.13), sp_c14=unique(sp_marker$Cluster.14))

upset(fromList(input), order.by = "freq")

l1 = intersect(ang_marker$Cluster.15, sp_marker$Cluster.14)
l2 = intersect(ang_marker$Cluster.19, sp_marker$Cluster.13)
ang19_sp13 = setdiff(l2, l1) # tafa1 and nav3
ang15_sp14 = setdiff(l1, l2) # epha7
gene=ang15_sp14
p1=DotPlot(seurat_ang, features = gene, group.by = "wsnn_res.0.5") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p2=DotPlot(seurat_ss, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p3=DotPlot(seurat_sp, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p1 + p3


input = list(ang_c15=unique(ang_marker$Cluster.15), ang_c19=unique(ang_marker$Cluster.19),
             ss_c16=unique(ss_marker$Cluster.16), ss_c20=unique(ss_marker$Cluster.20), 
             ss_c21=unique(ss_marker$Cluster.21), ss_c22=unique(ss_marker$Cluster.22), 
             ss_c23=unique(ss_marker$Cluster.23))

upset(fromList(input), order.by = "freq")

ang19_sp13 = setdiff(l2, l1) # tafa1 and nav3
ang15_sp14 = setdiff(l1, l2) # epha7
gene=ang15_sp14
p1=DotPlot(seurat_ang, features = gene, group.by = "wsnn_res.0.5") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p2=DotPlot(seurat_ss, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p3=DotPlot(seurat_sp, features = gene, group.by = "wsnn_res.0.4") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
p1 + p3



prop=table(seurat_ang$wsnn_res.0.5, seurat_ang$treatment)/colSums(table(seurat_ang$wsnn_res.0.5, seurat_ang$treatment))





outfolder = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/"




resolution = 0.5
outfile = paste0("rat.sp.LK.", "reso_", resolution)

reso_para = paste0('wsnn_res.', resolution)
seurat_object@active.ident = seurat_object@meta.data[, reso_para]
seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]
DefaultAssay(seurat_object) = "RNA"

p = DimPlot(seurat_object, reduction = "wnn.umap.harmony", group.by = c(reso_para), 
            pt.size = 0.2, raster = F, label = T) + theme(legend.position = "right") + 
  labs(x="UMAP 1", y="UMAP 2", title = paste0("reso=", resolution), color="Cluster") 
print(p)

### generate plots & tables for first-step cluster result
# cluster_sample_tbl = table(seurat_object$orig.ident, seurat_object$seurat_clusters)
# cluster_sample_sum(cluster_sample_tbl, paste0(outfolder, outfile, ".cluster_sample_sum.xlsx"))

markers = FindAllMarkers(seurat_object, group.by=reso_para, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
write.table(markers,file=paste0(outfolder, outfile, ".allmarker.0.25.long.txt"), sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                        v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file=paste0(outfolder, outfile, ".allmarker.0.25.wide.txt"), sep="\t", row.names = F)

panglao_anno(marker_tbl, paste0(outfolder, outfile, ".allmarker.0.25.wide.xlsx"))






marker_list = c("Nphs2", "Nphs1", # Podocytes
                "Aldh1a2","Rbfox1", # Parietal epithelial cell
                "Lrp2", # Proximal tubules
                "Slc22a8", "Slc22a7",
                "Havcr1", "Vcam1", # Injured proximal tubules
                "Slc44a5",
                "Slc12a1", "Umod", # Thick ascending limbs
                "Nos1", # MD
                "Sls12a3", "Trpm6", # Distal convoluted tubules
                "Slc8a1", "Calb1", # Connecting tubules
                "Scnn1g", "Aqp2", # Collecting ducts
                "Atp6v0d2", "Atp6v1c2", 
                "Slc26a7", # Type A intercalated cells
                "Slc26a4", # Type B intercalated cells
                
                "Flt1", # Endothelial cells
                "Dnase1l3", "Plvap", # EC-PTC
                "Hecw2", "Itga8", # EC-GC
                "Adamts6", "Palmd", # EC-AEA
                "Cd36", "Mmrn1", # EC-LYM
                
                "Pdgfrb", "Notch3", "Postn", # pericytes/Mesangial cells
                "C7", "Negr1", 
                "Col1a1", 
                
                "Ptprc", "Ms4a1", "Bank1", # B cellS
                "Cd74", 
                "Nrxn1")


DotPlot(seurat_object, group.by = "seurat_clusters",
        features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  # scale_color_gradientn(colors = colGEX) +
  labs(x="", y="") +
  guides(color=guide_colorbar("Average\nExpression"),
         size=guide_legend("Percent\nExpressed"))


