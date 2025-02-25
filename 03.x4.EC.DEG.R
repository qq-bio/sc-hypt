library(Seurat)
library(dplyr)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/")


i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds"
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
i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds"
seurat_object <- readRDS(i)
cluster = "seurat_clusters"
Idents(seurat_object) = cluster

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.DEG_all.out")
# deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type, deg_merged[deg_merged$p_val_adj<0.05,]$strain)

tbl = table(deg_merged[deg_merged$p_val_adj<0.05,]$gene_name)
tbl[order(tbl, decreasing = T)]
# ENSRNOG00000065867     AABR07000398.1 ENSRNOG00000062930 ENSRNOG00000070284 ENSRNOG00000068039             Zbtb16             Camk1d              Cmss1 
# 31                 29                 28                 24                 22                 18                 12                 11 
# Lars2              Ptprg                Ttn                Dst             Rnf213              Srrm2     AABR07044900.1             Cacnb2 
# 10                 10                  9                  8                  7                  7                  6                  6 
# ENSRNOG00000067673             Mt-co3               Pbx1              Plcb4             Arglu1            Auts2l1              Dach1             Luc7l3 
# 6                  6                  6                  6                  5                  5                  5                  5 
# Mt-co1             Mt-nd4          RT1-T24-4               Chd7              Ddx46               Ddx5               Dlc1              Eif3a 
# 5                  5                  5                  4                  4                  4                  4                  4 
# ENSRNOG00000055562               Etl4              Fabp4            Filip1l            Gm19951              Herc6              Hmcn1           Hsp90ab1 
# 4                  4                  4                  4                  4                  4                  4                  4 
# Mb              Meis2              Myrip              Plcl1              Pnisr            Prpf38b              Rasa4         RGD1565355 
# 4                  4                  4                  4                  4                  4                  4                  4 
# RT1-T24-1            Snrnp70               Uaca             Zbtb20     AABR07007032.1              Acin1              Ahnak              Akap9 
# 4                  4                  4                  4


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













