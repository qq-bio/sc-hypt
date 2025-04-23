library(Seurat)
library(Signac)
library(dplyr)



###
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.merged.rds")
cluster_order = c("M0610", "M24", "C18", "M5813",
                  "C1", "C7", "C9", 
                  "C12", "C21", "C15", "C19",
                  "C20", "C14", "C23", 
                  "C22", "C16", 
                  "C3", "C11", "C17")

seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels=cluster_order)
Idents(seurat_object) = "seurat_clusters"






###
ec_fibro <- c("Tgfb1", "Il1b", "Il6", "Edn1", "Notch1", "Snai1", "Snai2", "Twist1", "Vegfa", "Ccn2", "Mmp2", "Mmp9")

DotPlot(seurat_object, features = unique(ec_fibro)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(title = "Fibrosis")






###
seurat_lv = subset(seurat_object, seurat_clusters %in% c("M0610", "M24", "C7", "M5813"))
cluster_order = c("M24", "M0610", "C7", "M5813")
seurat_lv$seurat_clusters <- factor(seurat_lv$seurat_clusters, levels=cluster_order)
Idents(seurat_lv) = "seurat_clusters"


markers = FindAllMarkers(seurat_lv, group.by="seurat_clusters", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(rank = row_number()) %>% ungroup()
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.lv.DEG.out", sep="\t")






deg_ec <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")


lv_sub_genes <- c("Nrp2", "Col4a1", "Col4a2", "Nrp1", "Kdr", "Epas1", "Eng", "Smad6", "Efnb2", "Notch1", "Il1r1", "Calcrl", "Vwf", "Vegfc") # "Piezo2"
angio_in_m24 <- c("Cd34","Cd36","Sp1","Stab1","Tspan18") #Cd34 up in ss, Tsapn18 down in shr
angio_in_c7 <- c("Cdh13","Col4a1","Col4a2","Col15a1","Notch1","Sat1","Nrp2","Rnf213","Emc10") # most are down in shr
angio_in_m0610 <- c("Angpt2","Atp2b4","Cd36","Cdh5","Col4a2","Col4a3","Fyn","Pkm","Plcg1","Prkca","Slc12a2","Sp1","Sp100","Tgfbr2","Tie1","Nrp1","Emilin1","Plxnd1","Ppp1r16b","Adamts9","Rhoj","Tspan18","Emc10") # most are down in shr
angio_in_m5813 <- c("Col4a2","Ednra","Rnf213","Emc10") #Rnf213 down in shr, up in SD

gene_list <- c(angio_in_m24, angio_in_c7, angio_in_m0610, angio_in_m5813)
gene_list <- angio_in_m24
deg_ec[deg_ec$gene_name %in% gene_list & deg_ec$p_val_adj<0.05, ]
DotPlot(seurat_lv, features = unique(gene_list)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")




DotPlot(seurat_lv, features = unique(lv_sub_genes), split.by = "strain", cols = species_col) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")


