library(Seurat)
library(dplyr)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/")


i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/",
                 gsub("rds", "DEG_all.out", basename(i)))

seurat_object <- readRDS(i)
cluster = "seurat_clusters"
Idents(seurat_object) = cluster
meta_table = seurat_object@meta.data

project_list = unique(meta_table$project)

deg_merged = c()
for(pi in project_list){
  
  strain_list = unique(meta_table[meta_table$project==pi, ]$strain)
  for(si in strain_list){
    
    cell_list = unique(meta_table[meta_table$project==pi &
                                    meta_table$strain==si, cluster])
    
    treatment = meta_table[meta_table$project==pi &
                             meta_table$strain==si, ]$treatment
    treatment = intersect(treatment_order, unique(treatment))
    control = treatment[1]
    
    for(cell in cell_list){
      
      cell.1 = rownames(meta_table[meta_table$project==pi &
                                     meta_table[,cluster]==cell &
                                     meta_table$treatment==control &
                                     meta_table$strain==si, ])
      
      if(length(cell.1)>3){
        
        for(ti in treatment[-1]){
          
          cell.2 = rownames(meta_table[meta_table$project==pi &
                                         meta_table[,cluster]==cell &
                                         meta_table$treatment==ti &
                                         meta_table$strain==si, ])
          
          if(length(cell.2)>3){
            
            deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)
            
            # deg$gene_orig = rownames(deg)
            # deg$gene = convert_to_h_gene(rownames(deg))
            
            deg$gene_name = rownames(deg)
            deg$cell_type = cell
            deg$pct.diff = deg$pct.1 - deg$pct.2
            deg$project = pi
            deg$strain = si
            deg$control = control
            deg$treatment = ti
            deg$control_size = length(cell.1)
            deg$treatment_size = length(cell.2)
            deg$test_gene = nrow(deg)
            
            # deg = deg[deg$p_val_adj<0.05, ]
            
            deg_merged = rbind(deg_merged, deg)
            
            
          }
          
        }
        
      }
      
    }
    
  }
  
}

write.table(deg_merged, outfile)




################################################################################
i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
seurat_object <- readRDS(i)
cluster = "seurat_clusters"
Idents(seurat_object) = cluster

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")
# deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type, deg_merged[deg_merged$p_val_adj<0.05,]$strain)
table(deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5,]$cell_type, deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5,]$strain)

tbl = table(deg_merged[deg_merged$p_val_adj<0.05,]$gene_name)
tbl[order(tbl, decreasing = T)]
# ENSRNOG00000065867 ENSRNOG00000062930     AABR07000398.1 ENSRNOG00000070284 ENSRNOG00000068039             Zbtb16             Camk1d              Cmss1             Cacnb2              Lars2               Pbx1 
# 31                 24                 22                 22                 18                 16                 10                 10                  9                  8                  8 
# Ptprg             Rnf213            Auts2l1                Dst           Hsp90ab1     AABR07044900.1 ENSRNOG00000056487              Herc6              Plcl1                Ttn                B2m 
# 8                  8                  7                  7                  7                  6                  6                  6                  6                  6                  5 
# Cdh13 ENSRNOG00000067673              Grik2              Hmcn1              Mecom              Meis2             Mt-co3             Mt-nd4         RGD1565355              Srrm2             Zfp958 
# 5                  5                  5                  5                  5                  5                  5                  5                  5                  5                  5 
# AABR07007032.1     AABR07040864.1               Ano4             Arglu1             Col4a2               Dlc1 ENSRNOG00000055562               Etl4              Fhod3               Gbp1              Hdac9 
# 4                  4                  4                  4                  4                  4                  4                  4                  4                  4                  4 
# Id1             Luc7l3            Mt-atp6             Mt-co1             Mt-co2             Mt-cyb                Mx1               Nfia              Plcb4            Prpf38b              Rasa4 
# 4                  4                  4                  4                  4                  4                  4                  4                  4                  4                  4 
# RT1-T24-1          RT1-T24-4              Stab1            Tspan18             Zbtb20              Ahnak           AY036118               Bst2            Ccdc85a              Dach1              Ddx46 
# 4                  4                  4                  4                  4

tbl = table(deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5,]$gene_name)
tbl[order(tbl, decreasing = T)]
# ENSRNOG00000065867 ENSRNOG00000062930 ENSRNOG00000070284     AABR07000398.1 ENSRNOG00000068039             Zbtb16             Cacnb2             Camk1d              Cmss1              Lars2             Rnf213 
# 28                 22                 22                 19                 18                 16                  9                  9                  8                  8                  8 
# Hsp90ab1               Pbx1     AABR07044900.1            Auts2l1 ENSRNOG00000056487              Plcl1                Ttn                B2m                Dst              Grik2              Herc6 
# 7                  7                  6                  6                  6                  6                  6                  5                  5                  5                  5 
# Hmcn1             Mt-nd4             Zfp958     AABR07040864.1             Arglu1             Col4a2               Dlc1 ENSRNOG00000055562 ENSRNOG00000067673               Etl4              Fhod3 
# 5                  5                  5                  4                  4                  4                  4                  4                  4                  4                  4 
# Gbp1                Id1             Luc7l3              Meis2             Mt-co2             Mt-co3             Mt-cyb                Mx1            Prpf38b              Rasa4         RGD1565355 
# 4                  4                  4                  4                  4                  4                  4                  4                  4                  4                  4 
# RT1-T24-1          RT1-T24-4              Srrm2              Stab1            Tspan18              Ahnak               Ano4           AY036118               Bst2              Cdh13              Ddx46 
# 4                  4                  4                  4                  4 
for(i in names(tbl[tbl>=4])){
  print(deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$gene_name==i,])
}

# Il1r1+ EC
tbl = table(deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type %in% c("C23", "C20", "M5813"),]$gene_name)
tbl[order(tbl, decreasing = T)]
# ENSRNOG00000068039 ENSRNOG00000062930 ENSRNOG00000065867 ENSRNOG00000070284     AABR07000398.1     AABR07044900.1            Auts2l1              Ednra              Herc6             Inpp4b                Me3 
# 5                  4                  4                  4                  3                  2                  2                  2                  2                  2                  2 
# Rnf213              Stab1             Zbtb16
# 2                  2                  2
for(i in names(tbl[tbl>=2])){
  print(deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type %in% c("C23", "C20", "M5813") & deg_merged$gene_name==i,])
}


tbl = table(deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type %in% c("C14"),]$gene_name)
tbl[order(tbl, decreasing = T)]
# Etl4             Cadps2             Camk1d               Cdk8             Cep112              Cmss1               Gnas               Gphn           Hsp90ab1               Ica1             Jarid2 
# 4                  2                  2                  2                  2                  2                  2                  2                  2                  2                  2 
# Lars2                Mbp              Mcf2l               Mobp               Nav1              Nova2              Prkce              Ptgds               Rtn3             Slc6a6             Zbtb16 
# 2                  2                  2                  2                  2                  2                  2                  2                  2                  2                  2 
for(i in names(tbl[tbl>=2])){
  print(deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type %in% c("C14") & deg_merged$gene_name==i,])
}

deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="Cdk8",]



deg_merged[deg_merged$gene_name %in% c("Klf2", "Klf4", "Nos3", "Cdh5") & deg_merged$p_val_adj<0.05, ]
# p_val avg_log2FC pct.1 pct.2    p_val_adj gene_name cell_type pct.diff     project strain control treatment control_size treatment_size test_gene
# Cdh590  5.071221e-08  0.2449179 0.305 0.196 2.149691e-03      Cdh5         6    0.109 Spontaneous    SHR     10w       26w          997           2053      6868
# Klf2110 2.843230e-18 -2.3522984 0.058 0.485 1.205245e-13      Klf2        19   -0.427 Spontaneous    WKY     10w       26w          173            200      9524
# Klf4108 6.499390e-07 -1.6300259 0.075 0.285 2.755092e-02      Klf4        19   -0.210 Spontaneous    WKY     10w       26w          173            200      9524
deg_merged[deg_merged$gene_name %in% c("Prg5", "Myo10", "Tip60", "Ef1a2") & deg_merged$p_val_adj<0.05, ]
# p_val avg_log2FC pct.1 pct.2   p_val_adj gene_name cell_type pct.diff        project strain control treatment control_size treatment_size test_gene
# Myo1074 1.080692e-07   0.432542 0.615 0.515 0.004581053     Myo10         0      0.1 Salt-sensitive     SD      LS     HS 3d         1224            796      7166

deg_merged[deg_merged$cell_type==7 & deg_merged$p_val_adj<0.05,]
table(deg_merged[deg_merged$cell_type==7 & deg_merged$p_val_adj<0.05,]$gene_name)
# AABR07000398.1     AABR07044900.1             Cacnb2             Camk1d               Chd3              Cmss1               Dlc1 ENSRNOG00000062930 
# 4                  1                  1                  1                  1                  1                  1                  1 
# ENSRNOG00000065867 ENSRNOG00000068039 ENSRNOG00000069480 ENSRNOG00000070284             Mt-co1             Mt-nd4              Myrip              Plcb4 
# 2                  2                  1                  1                  1                  1                  1                  1 
# Ptk2              Ptprg          RT1-T24-4              Tmtc2                Ttn             Zbtb16 
# 1                  1                  1                  1                  1                  2 

deg_merged[deg_merged$cell_type==7 & deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="AABR07000398.1",]
FeaturePlot(seurat_object, "AABR07000398.1", split.by = "sxtxt")
deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="AABR07000398.1",]
deg_merged[deg_merged$cell_type==7 & deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="Zbtb16",]
deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="Zbtb16",]
deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$gene_name=="Reln",]


sort(table(deg_merged[deg_merged$cell_type==4 & deg_merged$p_val_adj<0.05,]$gene_name))
# Acta1             Atp2a2              Cdh13 
# 2                  2                  2 
# Cox7b              Cryab              Csrp3                Lpl                 Mb               Myl2              Tnnc1 
# 2                  2                  2                  2                  2                  2                  2 
# Tnni3              Tshz2              Uqcrq 
# 2                  2                  2 
for(i in c("Acta1", "Atp2a2", "Cdh13", "Cox7b", "Cryab", "Csrp3", "Lpl", 
           "Mb", "Myl2", "Tnnc1", "Tnni3", "Tshz2", "Uqcrq")){
  print(deg_merged[deg_merged$cell_type==4 & deg_merged$p_val_adj<0.05 & deg_merged$gene_name==i,])
}

# "Acta1", "Cox7b", "Cryab", "Csrp3", "Tnnc1", "Uqcrq"


################################################################################
i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds"
seurat_object <- readRDS(i)
cluster = "seurat_clusters"
Idents(seurat_object) = cluster

seurat_object$condition = ifelse(grepl("(C57BL/6-Ang|SS-HS|SHR-26w)", seurat_object$sxtxt), "Hypertension", "Normotension")
seurat_object$strain_cat = ifelse(seurat_object$strain %in% c(""), "Hypertensive", "Normotensive")

FeaturePlot(seurat_object, "Lpl", split.by = "condition")
FeaturePlot(seurat_object, "Lpl", split.by = "strain")
DotPlot(seurat_object, features = "Lpl", split.by = "condition")


################################################################################
library(metap)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.DEG_all.out")
deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)


deg_merged$sample_size <- deg_merged$control_size + deg_merged$treatment_size

deg_wide <- deg_merged %>%
  select(gene_name, cell_type, strain, treatment, p_val, avg_log2FC, sample_size) %>%
  pivot_wider(names_from = c(strain, treatment),
              values_from = c(p_val, avg_log2FC, sample_size),
              names_sep = "_") %>%
  as.data.frame()

# deg_wide <- deg_wide[rowSums(is.na(deg_wide[, c("p_val_C57BL/6_AngII 28d", "p_val_C57BL/6_AngII 3d", 
#                                                 "p_val_SS_HS 3d", "p_val_SS_HS 21d", "p_val_SHR_26w")])) < 3, ]


deg_wide$hypt_sqrt_weights <- apply(deg_wide[, c("sample_size_C57BL/6_AngII 28d", "sample_size_C57BL/6_AngII 3d", 
                                                 "sample_size_SS_HS 3d", "sample_size_SS_HS 21d", "sample_size_SHR_26w")], 
                                    1, function(sizes) {
                                      valid_sizes <- na.omit(sizes)  # Exclude NA values
                                      if (length(valid_sizes) > 0) {
                                        sum(sqrt(valid_sizes))  # Compute sum of sqrt(sample size)
                                      } else {
                                        NA  # Return NA if all sample sizes are missing
                                      }
                                    })

deg_wide$hypt_combined_p <- apply(deg_wide[, c("p_val_C57BL/6_AngII 28d", "p_val_C57BL/6_AngII 3d", 
                                               "p_val_SS_HS 3d", "p_val_SS_HS 21d", "p_val_SHR_26w")], 1, function(p) {
                                                 valid_p <- na.omit(as.numeric(p))  # Remove NA values
                                                 if (length(valid_p) > 0) {
                                                   sumz(valid_p, weights = hypt_sqrt_weights)$p  # Compute combined p-value
                                                 } else {
                                                   NA  # Return NA if all p-values are missing
                                                 }
                                               })

deg_wide$hypt_weighted_FC <- apply(deg_wide[, c("fc_C57BL/6_AngII 28d", "fc_C57BL/6_AngII 3d", 
                                                "fc_SS_HS 3d", "fc_SS_HS 21d", "fc_SHR_26w")], 1, function(fc) {
                                                  valid_hypt_fc <- na.omit(as.numeric(fc))  # Remove NA values
                                                  valid_hypt_weights <- deg_wide$hypt_sqrt_weights[!is.na(fc)]
                                                  if (length(valid_fc) > 0) {
                                                    sum(valid_fc * valid_weights) / sum(valid_weights)
                                                    
                                                })

deg_wide$hypt_meta_score <- -log10(deg_wide$hypt_combined_p) * deg_wide$hypt_weighted_FC

deg_wide$hypt_meta_pval_adj <- p.adjust(deg_wide$hypt_combined_p, method = "fdr")



deg_wide$nt_sqrt_weights <- sqrt(deg_wide$`sample_size_SD_HS 3d`) + sqrt(deg_wide$`sample_size_WKY_26w`)

deg_wide$nt_combined_p <- apply(deg_wide[, c("p_val_SD_HS 3d", "p_val_WKY_26w")], 1, function(p) {
  sumz(as.numeric(p), weights = deg_wide$nt_sqrt_weights)$p
})

deg_wide$nt_weighted_FC <- apply(deg_wide[, c("fc_SD_HS 3d", "fc_WKY_26w")], 1, function(fc) {
  sum(as.numeric(fc) * deg_wide$nt_sqrt_weights) / sum(deg_wide$nt_sqrt_weights)
})

deg_wide$nt_meta_score <- -log10(deg_wide$nt_combined_p) * deg_wide$nt_weighted_FC

deg_wide$nt_meta_pval_adj <- p.adjust(deg_wide$nt_combined_p, method = "fdr")













