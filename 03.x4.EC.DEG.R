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
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.DEG_all.out")
deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type)
table(deg_merged[deg_merged$p_val_adj<0.05,]$cell_type, deg_merged[deg_merged$p_val_adj<0.05,]$strain)








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













