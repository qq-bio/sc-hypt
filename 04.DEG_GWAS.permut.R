library(tidyr)
library(dplyr)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")


# ################################################################################
# ### load reference files
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

### h2m
r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) = c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) = c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")
# 
# 
# 
# gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
# # proximal_merge = read.table("proximal_merged.txt", header=T, sep='\t', quote = "")
# proximal_merge = read.table("gwas_snp_gene_merged.txt", header=T, sep='\t')
# eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
# e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')
# 
# deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
# deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]
# 
# 
# 
# # Initialize dataframe to store results
# permutation_test_df <- data.frame()
# 
# # Number of permutations
# num_permutations <- 100000
# 
# # Perform the analysis
# for (trait in unique(gwas_merge$trait)) {
#   
#   SNPS <- unique(gwas_merge[gwas_merge$trait == trait, ]$SNPS)
#   eqtl_merge_use <- eqtl_merge[eqtl_merge$SNP %in% SNPS, ]
#   e2g_merge_use <- e2g_merge[e2g_merge$SNP %in% SNPS, ]
#   proximal_merge_use <- proximal_merge[proximal_merge$SNP %in% SNPS, ]
#   
#   snp_gene_prox <- data.frame(SNP = proximal_merge_use$SNP, gene = proximal_merge_use$Gene_ID)
#   snp_gene_eqtl <- data.frame(SNP = eqtl_merge_use$SNP, gene = eqtl_merge_use$gene_id_mod)
#   snp_gene_e2g <- data.frame(SNP = e2g_merge_use$SNP, gene = e2g_merge_use$gene)
#   snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))
#   
#   snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")
#   
#   snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
#   snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
#   snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
#   snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat <- snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_id_rat
#   snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse <- snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_id_mouse
#   
#   for (si in c("C57BL/6", "SHR", "SS")) {
#     treatment_list <- unique(deg_merged[deg_merged$strain == si, ]$treatment)
#     for (ti in treatment_list) {
#       tissue_list <- unique(deg_merged[deg_merged$strain == si & deg_merged$treatment == ti, ]$tissue)
#       for (tissue in tissue_list) {
#         cell_list <- unique(deg_merged[deg_merged$strain == si & deg_merged$treatment == ti & deg_merged$tissue == tissue, ]$cell_type)
#         for (ci in cell_list) {
#           deg_merged_use <- deg_merged[deg_merged$strain == si & deg_merged$treatment == ti & deg_merged$tissue == tissue & deg_merged$cell_type == ci, ]
#           expr_gene_list <- deg_merged_use$gene_name
#           deg_list <- deg_merged_use[deg_merged_use$p_val_adj < 0.05 & abs(deg_merged_use$avg_log2FC) > 0.5, ]$gene_name
#           
#           if (length(deg_list) > 0) {
#             if (si == "C57BL/6") {
#               snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse)
#               all_gene_list <- unique(m2h[m2h$gene_id != "", ]$gene_name_mouse)
#             } else {
#               snp_gene_list <- unique(snp_gene_df[!is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat)
#               all_gene_list <- unique(r2h[r2h$gene_id != "", ]$gene_name_rat)
#             }
#             
#             overlap_genes <- intersect(snp_gene_list, deg_list)
#             observed_overlap <- length(overlap_genes)
#             
#             # Permutation test
#             permutation_overlaps <- replicate(num_permutations, {
#               random_genes <- sample(all_gene_list, length(snp_gene_list))
#               length(intersect(random_genes, deg_list))
#             })
#             
#             p_value <- mean(permutation_overlaps >= observed_overlap)
#             
#             permutation_test_df <- rbind(permutation_test_df,
#                                          c(trait, si, ti, tissue, ci, length(snp_gene_list), length(deg_list), length(expr_gene_list), 
#                                            observed_overlap, p_value, paste0(overlap_genes, collapse = ", ")))
#           }
#         }
#       }
#     }
#   }
# }
# 
# # Convert to dataframe and set column names
# permutation_test_df <- as.data.frame(permutation_test_df)
# colnames(permutation_test_df) <- c("trait", "strain", "treatment", "tissue", "cell_type", "#SNP genes", "#DEG", "expr genes", 
#                                    "#SNP-DEG", "p.value", "gene_list")
# 
# # Adjust p-values
# permutation_test_df$p.adj <- p.adjust(permutation_test_df$p.value, method = "BH")
# write.table(permutation_test_df, "trait.permut.snp_gene.0712.out", col.names = T, row.names = F, sep = "\t", quote=F)









################################################################################
gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
# proximal_merge = read.table("proximal_merged.100k.txt", header=T, sep='\t', quote = "")
proximal_merge = read.table("gwas_snp_gene_merged.txt", header=T, sep='\t')
eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

fisher_test_df = c()
num_permutations = 1000
for(trait in unique(gwas_merge$trait)){
  
  SNPS = unique(gsub(" ", "", gwas_merge[gwas_merge$trait==trait, ]$SNPS))
  eqtl_merge_use = eqtl_merge[eqtl_merge$SNP %in% SNPS, ]
  e2g_merge_use = e2g_merge[e2g_merge$SNP %in% SNPS, ]
  proximal_merge_use = proximal_merge[proximal_merge$SNP %in% SNPS, ]
  
  snp_gene_prox = data.frame(SNP=proximal_merge_use$SNP, gene=proximal_merge_use$Gene_ID)
  snp_gene_eqtl = data.frame(SNP=eqtl_merge_use$SNP, gene=eqtl_merge_use$gene_id_mod)
  snp_gene_e2g = data.frame(SNP=e2g_merge_use$SNP, gene=e2g_merge_use$gene)
  snp_gene_df = unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))
  
  snp_gene_df$pair = paste(snp_gene_df$SNP, snp_gene_df$gene, sep="-")
  
  snp_gene_df = base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x="gene", by.y="Gene.stable.ID", all.x=T)
  snp_gene_df = base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
  snp_gene_df = base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
  snp_gene_df[snp_gene_df$gene_name_rat=="" & !(is.na(snp_gene_df$gene_name_rat)), ]$gene_name_rat = snp_gene_df[snp_gene_df$gene_name_rat=="" & !(is.na(snp_gene_df$gene_name_rat)), ]$gene_id_rat
  snp_gene_df[snp_gene_df$gene_name_mouse=="" & !(is.na(snp_gene_df$gene_name_mouse)), ]$gene_name_mouse = snp_gene_df[snp_gene_df$gene_name_mouse=="" & !(is.na(snp_gene_df$gene_name_mouse)), ]$gene_id_mouse
  
  for(si in c("C57BL/6", "SHR", "SS")){
    treatment_list = unique(deg_merged[deg_merged$strain==si, ]$treatment)
    for(ti in treatment_list){
      tissue_list = unique(deg_merged[deg_merged$strain==si & deg_merged$treatment==ti, ]$tissue)
      for(tissue in tissue_list){
        cell_list = unique(deg_merged[deg_merged$strain==si & deg_merged$treatment==ti & deg_merged$tissue==tissue, ]$cell_type)
        for(ci in cell_list){
          deg_merged_use = deg_merged[deg_merged$strain==si & deg_merged$treatment==ti & deg_merged$tissue==tissue & deg_merged$cell_type==ci, ]
          expr_gene_list = unique(deg_merged_use$gene_name)
          # deg_list = unique(deg_merged_use[deg_merged_use$p_val_adj<0.05 & abs(deg_merged_use$avg_log2FC)>0.5, ]$gene_name)
          deg_list = unique(deg_merged_use[deg_merged_use$p_val_adj<0.05 & abs(deg_merged_use$avg_log2FC)>0.25, ]$gene_name)
          
          # if(length(deg_list)>10){
          if(si=="C57BL/6"){
            snp_gene_list = unique(snp_gene_df[! is.na(snp_gene_df$gene_name_mouse),]$gene_name_mouse)
            all_gene_list = unique(c(m2h[m2h$gene_id!="",]$gene_name_mouse, expr_gene_list))
          }else{
            snp_gene_list = unique(snp_gene_df[! is.na(snp_gene_df$gene_name_rat),]$gene_name_rat)
            all_gene_list = unique(c(r2h[r2h$gene_id!="",]$gene_name_rat, expr_gene_list))
          }
          
          overlap_genes = intersect(snp_gene_list, deg_list)
          a <- length(overlap_genes)  # SNP-related genes that are DEGs
          b <- length(setdiff(snp_gene_list, deg_list))    # SNP-related genes that are not DEGs
          c <- length(setdiff(deg_list, snp_gene_list))    # DEGs that are not SNP-related
          d <- length(setdiff(all_gene_list, union(snp_gene_list, deg_list)))  # Genes that are neither
          
          contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE,
                                      dimnames = list(c("SNP-related", "Not SNP-related"),
                                                      c("DEG", "Not DEG")))
          
          fisher_test_result <- fisher.test(contingency_table, alternative = "greater")
          observed_log_p_value <- -log10(fisher_test_result$p.value)
          
          # Permutation test
          permutation_log_p_values <- numeric(num_permutations)
          for (i in 1:num_permutations) {
            permuted_snp_list <- sample(all_gene_list, length(snp_gene_list))
            permuted_overlap_genes = intersect(permuted_snp_list, deg_list)
            a_perm <- length(permuted_overlap_genes)
            b_perm <- length(setdiff(permuted_snp_list, deg_list))
            c_perm <- length(setdiff(deg_list, permuted_snp_list))
            d_perm <- length(setdiff(all_gene_list, union(permuted_snp_list, deg_list)))
            
            perm_contingency_table <- matrix(c(a_perm, b_perm, c_perm, d_perm), nrow = 2, byrow = TRUE,
                                             dimnames = list(c("SNP-related", "Not SNP-related"),
                                                             c("DEG", "Not DEG")))
            
            perm_fisher_test_result <- fisher.test(perm_contingency_table, alternative = "greater")
            permutation_log_p_values[i] <- -log10(perm_fisher_test_result$p.value)
          }
          
          mean_permutation_log_p_value <- mean(permutation_log_p_values, )
          NES <- observed_log_p_value / mean_permutation_log_p_value
          
          fisher_test_df = rbind(fisher_test_df,
                                 c(trait, si, ti, tissue, ci, length(snp_gene_list), length(deg_list), length(expr_gene_list), 
                                   a, b, c, d, fisher_test_result$p.value, paste0(overlap_genes, collapse = ", "), NES))
          
          
          # }
        }
      }
    }
  }
  
}

fisher_test_df = as.data.frame(fisher_test_df)
colnames(fisher_test_df) = c("trait", "strain", "treatment", "tissue", "cell_type", "#SNP genes", "#DEG", "expr genes", 
                             "#SNP-DEG", "#SNP-not-DEG", "#DEG-not-SNP", "#neither", "p.value", "gene_list", "NES")
fisher_test_df$p.adj = p.adjust(fisher_test_df$p.value, method = "BH")
write.table(fisher_test_df, "trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", col.names = T, row.names = F, sep = "\t", quote=F)













