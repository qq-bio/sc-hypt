
library(GenomicRanges)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(patchwork)
library(stringr)
library(leiden)
library(lsa)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")



################################################################################
### load reference files
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) <- ensembl$Gene.stable.ID

r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) <- c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) <- c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")





################################################################################
### Cell Type Enrichment Analysis od SNP-related DEG
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")

gwas_merge <- read.table("gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t')
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged$pct_avg <- (deg_merged$pct.1*deg_merged$control_size + deg_merged$pct.2*deg_merged$treatment_size)/rowSums(deg_merged[, c("control_size", "treatment_size")])

fisher_test_df <- data.frame()
num_permutations <- 1000

# Perform analysis for each trait
for (trait in unique(gwas_merge$trait)) {
  
  # Filter SNPs for the current trait
  SNPS <- unique(gwas_merge[gwas_merge$trait == trait, ]$SNPS)
  
  # Subset merged data based on selected SNPs
  eqtl_merge_use <- eqtl_merge[eqtl_merge$SNP %in% SNPS, ]
  e2g_merge_use <- e2g_merge[e2g_merge$SNP %in% SNPS, ]
  proximal_merge_use <- proximal_merge[proximal_merge$SNP %in% SNPS, ]
  
  # Combine SNP-gene mappings
  snp_gene_prox <- data.frame(SNP = proximal_merge_use$SNP, gene = proximal_merge_use$Gene_ID)
  snp_gene_eqtl <- data.frame(SNP = eqtl_merge_use$SNP, gene = eqtl_merge_use$gene_id_mod)
  snp_gene_e2g <- data.frame(SNP = e2g_merge_use$SNP, gene = e2g_merge_use$gene)
  snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))
  # snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")
  
  # Merge gene names with SNP-gene pairs
  snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
  snp_gene_df$pair <- ifelse(snp_gene_df$Gene.name=="", paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-"), paste(snp_gene_df$SNP, snp_gene_df$Gene.name, sep = "-"))
  snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
  snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
  snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat <- snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_id_rat
  snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse <- snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_id_mouse
  
  for (strain in c("C57BL/6", "SHR", "SS")) {
    treatment_list <- unique(deg_merged[deg_merged$strain == strain, ]$treatment)
    
    for (treatment in treatment_list) {
      tissue_list <- unique(deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment, ]$tissue)
      
      for (tissue in tissue_list) {
        cell_list <- unique(deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment & deg_merged$tissue == tissue, ]$cell_type)
        
        for (cell_type in cell_list) {
          deg_merged_use <- deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment & deg_merged$tissue == tissue & deg_merged$cell_type == cell_type, ]
          # deg_list <- unique(deg_merged_use[deg_merged_use$p_val_adj < 0.05 & abs(deg_merged_use$avg_log2FC) > 0.25, ]$gene_name)
          deg_list <- unique(deg_merged_use[deg_merged$pct_avg>0.1, ]$gene_name)
          
          # # Define SNP-related and all genes based on strain
          # if (strain == "C57BL/6") {
          #   snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse)
          #   all_gene_list <- unique(c(m2h[m2h$gene_id != "", ]$gene_name_mouse, unique(deg_merged_use$gene_name)))
          # } else {
          #   snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat)
          #   all_gene_list <- unique(c(r2h[r2h$gene_id != "", ]$gene_name_rat, unique(deg_merged_use$gene_name)))
          # }
          # 
          if (strain == "C57BL/6") {
            snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse)
            all_gene_list <- unique(c(m2h[m2h$gene_id != "", ]$gene_name_mouse, unique(deg_merged_use$gene_short_name)))
            snp_covered_list <- paste(unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse) & snp_gene_df$gene_name_mouse %in% deg_list, ]$pair), collapse = ",")
          } else {
            snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat)
            all_gene_list <- unique(c(r2h[r2h$gene_id != "", ]$gene_name_rat, unique(deg_merged_use$gene_short_name)))
            snp_covered_list <- paste(unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat) & snp_gene_df$gene_name_rat %in% deg_list, ]$pair), collapse = ",")
          }
          
          gene_covered_list = paste0(intersect(snp_gene_list, deg_list), collapse = ", ")
          
          # Prepare contingency table for Fisher's exact test
          a <- length(intersect(snp_gene_list, deg_list))
          b <- length(setdiff(snp_gene_list, deg_list))
          c <- length(setdiff(deg_list, snp_gene_list))
          d <- length(setdiff(all_gene_list, union(snp_gene_list, deg_list)))
          
          contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, 
                                      dimnames = list(c("SNP-related", "Not SNP-related"), c("DEG", "Not DEG")))
          
          fisher_test_result <- fisher.test(contingency_table, alternative = "greater")
          observed_log_p_value <- -log10(fisher_test_result$p.value)
          
          # Permutation test
          permutation_log_p_values <- numeric(num_permutations)
          for (i in 1:num_permutations) {
            permuted_snp_list <- sample(all_gene_list, length(snp_gene_list))
            permuted_overlap_genes <- intersect(permuted_snp_list, deg_list)
            a_perm <- length(permuted_overlap_genes)
            b_perm <- length(setdiff(permuted_snp_list, deg_list))
            c_perm <- length(setdiff(deg_list, permuted_snp_list))
            d_perm <- length(setdiff(all_gene_list, union(permuted_snp_list, deg_list)))
            
            perm_contingency_table <- matrix(c(a_perm, b_perm, c_perm, d_perm), nrow = 2, byrow = TRUE,
                                             dimnames = list(c("SNP-related", "Not SNP-related"), c("DEG", "Not DEG")))
            
            perm_fisher_test_result <- fisher.test(perm_contingency_table, alternative = "greater")
            permutation_log_p_values[i] <- -log10(perm_fisher_test_result$p.value)
          }
          
          # Calculate normalized enrichment score (NES)
          mean_permutation_log_p_value <- mean(permutation_log_p_values)
          NES <- observed_log_p_value / mean_permutation_log_p_value
          
          # Append results to the dataframe
          fisher_test_df <- rbind(fisher_test_df, 
                                  data.frame(trait = trait, strain = strain, treatment = treatment, tissue = tissue, cell_type = cell_type, 
                                             SNP_genes = length(snp_gene_list), DEG = length(deg_list), expr_genes = length(unique(deg_merged_use$gene_name)), 
                                             SNP_DEG = a, SNP_not_DEG = b, DEG_not_SNP = c, neither = d, 
                                             p_value = fisher_test_result$p.value, gene_covered = gene_covered_list, SNP_covered = snp_covered_list, NES = NES))
        }
      }
    }
  }
}

# Adjust p-values for multiple testing and write results to file
fisher_test_df$p.adj <- p.adjust(fisher_test_df$p_value, method = "BH")
write.table(fisher_test_df, "trait.fisher.snp_expr_gene.enrichment_results.w_SNP.out", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)




# 
# ################################################################################
# ### Cell Type Enrichment Analysis od SNP-related DEG
# ################################################################################
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")
# 
# gwas_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
# proximal_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_snp_gene_merged.txt", header = TRUE, sep = '\t')
# eqtl_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/eqtl_merged.txt", header = TRUE, sep = '\t')
# e2g_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/e2g_merged.txt", header = TRUE, sep = '\t')
# 
# deg_LK <- read.table("/xdisk/mliang1/qqiu/project/sex-difference/DEG/RLK.DEG.out", sep = '\t', header = TRUE)
# deg_LV <- read.table("/xdisk/mliang1/qqiu/project/sex-difference/DEG/RLV.DEG.out", sep = '\t', header = TRUE)
# deg_LK$tissue <- "LK"
# deg_LV$tissue <- "LV"
# colnames(deg_LK) <- gsub("HS.*d", "HS", colnames(deg_LK))
# colnames(deg_LV) <- gsub("HS.*d", "HS", colnames(deg_LV))
# 
# deg_merged <- rbind(deg_LK, deg_LV)
# deg_merged$strain <- "SS"
# deg_merged$treatment <- "HS"
# 
# fisher_test_df <- data.frame()
# num_permutations <- 1000
# 
# # Perform analysis for each trait
# for (trait in unique(gwas_merge$trait)) {
#   
#   # Filter SNPs for the current trait
#   SNPS <- unique(gwas_merge[gwas_merge$trait == trait, ]$SNPS)
#   
#   # Subset merged data based on selected SNPs
#   eqtl_merge_use <- eqtl_merge[eqtl_merge$SNP %in% SNPS, ]
#   e2g_merge_use <- e2g_merge[e2g_merge$SNP %in% SNPS, ]
#   proximal_merge_use <- proximal_merge[proximal_merge$SNP %in% SNPS, ]
#   
#   # Combine SNP-gene mappings
#   snp_gene_prox <- data.frame(SNP = proximal_merge_use$SNP, gene = proximal_merge_use$Gene_ID)
#   snp_gene_eqtl <- data.frame(SNP = eqtl_merge_use$SNP, gene = eqtl_merge_use$gene_id_mod)
#   snp_gene_e2g <- data.frame(SNP = e2g_merge_use$SNP, gene = e2g_merge_use$gene)
#   snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))
#   
#   # Merge gene names with SNP-gene pairs
#   snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
#   snp_gene_df$pair <- ifelse(snp_gene_df$Gene.name=="", paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-"), paste(snp_gene_df$SNP, snp_gene_df$Gene.name, sep = "-"))
#   snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
#   snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
#   snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat <- snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_id_rat
#   snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse <- snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_id_mouse
#   
#   for (strain in c("SS")) {
#     treatment_list <- unique(deg_merged[deg_merged$strain == strain, ]$treatment)
#     
#     for (treatment in treatment_list) {
#       tissue_list <- unique(deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment, ]$tissue)
#       
#       for (tissue in tissue_list) {
#         cell_list <- unique(deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment & deg_merged$tissue == tissue, ]$cell_type)
#         
#         for (cell_type in cell_list) {
#           deg_merged_use <- deg_merged[deg_merged$strain == strain & deg_merged$treatment == treatment & deg_merged$tissue == tissue & deg_merged$cell_type == cell_type, ]
#           
#           for (sex in c("f", "m")) {
#             
#             pval_col <- paste0("p_val_adj.", sex)
#             logfc_col <- paste0("avg_log2FC.", sex)
#             
#             sex <- ifelse(sex=="f", "Female", "Male")
#             
#             deg_list <- unique(deg_merged_use[deg_merged_use[[pval_col]] < 0.05 & abs(deg_merged_use[[logfc_col]]) > 0.25, ]$gene_short_name)
#             
#             # Define SNP-related and all genes based on strain
#             if (strain == "C57BL/6") {
#               snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse)
#               all_gene_list <- unique(c(m2h[m2h$gene_id != "", ]$gene_name_mouse, unique(deg_merged_use$gene_short_name)))
#               snp_covered_list <- paste(unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse) & snp_gene_df$gene_name_mouse %in% deg_list, ]$pair), collapse = ",")
#             } else {
#               snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat)
#               all_gene_list <- unique(c(r2h[r2h$gene_id != "", ]$gene_name_rat, unique(deg_merged_use$gene_short_name)))
#               snp_covered_list <- paste(unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat) & snp_gene_df$gene_name_rat %in% deg_list, ]$pair), collapse = ",")
#             }
#             
#             gene_covered_list = paste0(intersect(snp_gene_list, deg_list), collapse = ", ")
#             
#             # Prepare contingency table for Fisher's exact test
#             a <- length(intersect(snp_gene_list, deg_list))
#             b <- length(setdiff(snp_gene_list, deg_list))
#             c <- length(setdiff(deg_list, snp_gene_list))
#             d <- length(setdiff(all_gene_list, union(snp_gene_list, deg_list)))
#             
#             contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE, 
#                                         dimnames = list(c("SNP-related", "Not SNP-related"), c("DEG", "Not DEG")))
#             
#             fisher_test_result <- fisher.test(contingency_table, alternative = "greater")
#             observed_log_p_value <- -log10(fisher_test_result$p.value)
#             
#             # Permutation test
#             permutation_log_p_values <- numeric(num_permutations)
#             for (i in 1:num_permutations) {
#               permuted_snp_list <- sample(all_gene_list, length(snp_gene_list))
#               permuted_overlap_genes <- intersect(permuted_snp_list, deg_list)
#               a_perm <- length(permuted_overlap_genes)
#               b_perm <- length(setdiff(permuted_snp_list, deg_list))
#               c_perm <- length(setdiff(deg_list, permuted_snp_list))
#               d_perm <- length(setdiff(all_gene_list, union(permuted_snp_list, deg_list)))
#               
#               perm_contingency_table <- matrix(c(a_perm, b_perm, c_perm, d_perm), nrow = 2, byrow = TRUE,
#                                                dimnames = list(c("SNP-related", "Not SNP-related"), c("DEG", "Not DEG")))
#               
#               perm_fisher_test_result <- fisher.test(perm_contingency_table, alternative = "greater")
#               permutation_log_p_values[i] <- -log10(perm_fisher_test_result$p.value)
#             }
#             
#             # Calculate normalized enrichment score (NES)
#             mean_permutation_log_p_value <- mean(permutation_log_p_values)
#             NES <- observed_log_p_value / mean_permutation_log_p_value
#             
#             # Append results to the dataframe
#             fisher_test_df <- rbind(fisher_test_df, 
#                                     data.frame(trait = trait, strain = strain, sex = sex, treatment = treatment, tissue = tissue, cell_type = cell_type, 
#                                                SNP_genes = length(snp_gene_list), DEG = length(deg_list), expr_genes = length(unique(deg_merged_use$gene_name)), 
#                                                SNP_DEG = a, SNP_not_DEG = b, DEG_not_SNP = c, neither = d, 
#                                                p_value = fisher_test_result$p.value, gene_covered = gene_covered_list, SNP_covered = snp_covered_list, NES = NES)) 
#             
#           }
#         }
#       }
#     }
#   }
# }
# 
# # Adjust p-values for multiple testing and write results to file
# fisher_test_df$p.adj <- p.adjust(fisher_test_df$p_value, method = "BH")
# write.table(fisher_test_df, "trait.fisher.sex.snp_expr.gene.enrichment_results.out", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
