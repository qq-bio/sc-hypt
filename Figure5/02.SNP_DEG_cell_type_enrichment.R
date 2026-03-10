
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
### DEG Analysis for Each Trait (Figure 4b)
################################################################################
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
snp_gene_df <- read.table("gwas_snp_gene.summary.out", header = T, sep = "\t")
snp_gene_df <- snp_gene_df %>% tidyr::separate_rows(., Gene_symbol_Mouse_Rat, sep=";")

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]

thresholds <- list(
  "p.adj < 0.05" = deg_merged %>% filter(p_val_adj < 0.05),
  "p.adj < 0.05 & |log2(FC)| > 0.25" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25),
  "p.adj < 0.05 & |log2(FC)| > 0.5" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
  "p.adj < 0.05 & |log2(FC)| > 1" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >1)
)

results <- data.frame()
for (trait in unique(gwas_merge$trait)) {
  SNPS <- unique(gwas_merge[gwas_merge$trait == trait, ]$SNPS)
  
  snp_gene_df_use = unique(snp_gene_df[snp_gene_df$SNP %in% SNPS,])
  snp_gene_list = unique(snp_gene_df_use$Gene_ID_Human)
  
  for (threshold_name in names(thresholds)) {
    deg_filtered <- thresholds[[threshold_name]]
    mapped_gene_list <- unique(c(snp_gene_df_use$Gene_symbol_Mouse_Rat))
    mapped_gene_list = setdiff(mapped_gene_list, "")
    deg_list <- unique(deg_filtered$gene_name)
    
    hit_genes <- intersect(mapped_gene_list, deg_list)
    num_hit_genes <- length(hit_genes)
    proportion_hit_genes <- num_hit_genes / length(snp_gene_list)
    num_DEG_snps = length(unique(snp_gene_df_use[snp_gene_df_use$Gene_symbol_Mouse_Rat %in% deg_list, ]$SNP))
    
    results <- rbind(results, data.frame(
      trait = trait,
      SNPS = length(unique(snp_gene_df_use$SNP)),
      SNP_genes = length(snp_gene_list),
      mapped_genes = length(mapped_gene_list),
      threshold = threshold_name,
      DEGs = num_hit_genes,
      DEG_SNPS = num_DEG_snps,
      proportion_DEG = proportion_hit_genes,
      proportion_SNP = num_DEG_snps/length(unique(snp_gene_df_use$SNP))
    ))
    
  }
}

results$label <- paste(results$DEGs, 
                       "(", round(results$proportion_DEG * 100, 2), "%)", 
                       ifelse(results$p.adj < 0.05, "*", ""), sep = "")

write.table(results, "trait.deg_prop.out", sep = ",", col.names = T, row.names = F, quote = F)


results_use = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/trait.deg_prop.out", sep = ",", header=T)
results_use$threshold = gsub("& ", "&\n", results_use$threshold)
threshold_col = RColorBrewer::brewer.pal(5, "Blues")[2:5]
names(threshold_col) = unique(results_use$threshold)

results_use = results_use[results_use$threshold=="p.adj < 0.05 &\n|log2(FC)| > 0.5",]
results_use$trait <- factor(results_use$trait, levels = results_use %>% group_by(trait) %>% dplyr::summarise(total_hits = sum(proportion_SNP)) %>% arrange(total_hits) %>% pull(trait))

p = ggplot(results_use) +
  geom_segment(aes(x = trait, xend = trait, y = proportion_DEG, yend = proportion_SNP), color = "gray", size = 1) +
  geom_point(aes(x = trait, y = proportion_DEG, label = DEGs, color = "SNP-related genes"), size = 3) +
  geom_point(aes(x = trait, y = proportion_SNP, label = DEG_SNPS, color = "SNPs"), size = 3) +
  geom_text(aes(x = trait, y = proportion_DEG-0.12, label = DEGs), color = "cadetblue4", size = 3) +
  geom_text(aes(x = trait, y = proportion_SNP+0.12, label = DEG_SNPS), color = "chocolate4", size = 3) +
  scale_color_manual(name = "", values = c("SNP-related genes" = "cadetblue4", "SNPs" = "chocolate4")) +
  scale_y_continuous(limits = c(0.2, 1.1), breaks = c(0.25, 0.5, 0.75, 1)) +
  labs(x = "Trait", y = "Proportion hit by DEGs") +
  theme_classic(base_family = "Arial") +
  theme(
    # legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(colour = 'black')
  ) +
  coord_flip()

print(p)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4b.deg_snp_number.png", width = 460 / 96, height = 178 / 96, dpi = 300)






################################################################################
### Cell Type Enrichment Analysis od SNP-related DEG
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")

gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
snp_gene_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_snp_gene.summary.out", header = T, sep = "\t")

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header=T)

fisher_test_df = c()
num_permutations = 1000
logFC_threshold = 0.25
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/trait.snp_gene.L1.fc_0.25.permut.out"

for(trait in unique(gwas_merge$trait)){
  
  SNPS = unique(gsub(" ", "", gwas_merge[gwas_merge$trait==trait, ]$SNPS))
  snp_gene_use = snp_gene_df[snp_gene_df$SNP %in% SNPS, ] %>%
    filter(!(is.na(Gene_symbol_Mouse_Rat) | Gene_symbol_Mouse_Rat=="")) %>%
    separate_longer_delim(
      col = Gene_symbol_Mouse_Rat,
      delim = ";"
    )
  
  for(si in c("C57BL/6", "SHR", "SS")){
    treatment_list = unique(deg_merged[deg_merged$strain==si, ]$treatment)
    for(ti in treatment_list){
      tissue_list = unique(deg_merged[deg_merged$strain==si & deg_merged$treatment==ti, ]$tissue)
      for(tissue in tissue_list){
        cell_list = unique(deg_merged[deg_merged$strain==si & deg_merged$treatment==ti & deg_merged$tissue==tissue, ]$cell_type)
        for(ci in cell_list){
          deg_merged_use = deg_merged[deg_merged$strain==si & deg_merged$treatment==ti & deg_merged$tissue==tissue & deg_merged$cell_type==ci, ]
          expr_gene_list = unique(deg_merged_use$gene_name)
          deg_list = unique(deg_merged_use[deg_merged_use$p_val_adj<0.05 & abs(deg_merged_use$avg_log2FC)>logFC_threshold, ]$gene_name)
          
          # if(length(deg_list)>10){
          if(si=="C57BL/6"){
            snp_gene_list = unique(snp_gene_use$Gene_symbol_Mouse_Rat)
            all_gene_list = unique(c(m2h[m2h$gene_id!="",]$gene_name_mouse, expr_gene_list))
          }else{
            snp_gene_list = unique(snp_gene_use$Gene_symbol_Mouse_Rat)
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
write.table(fisher_test_df, outfile, col.names = T, row.names = F, sep = "\t", quote=F)





################################################################################
### Visualize Fisher Test Results Using Dot Plot (Figure 4e)
################################################################################

# Load the Fisher test results
fisher_test_df <- read.table("trait.snp_gene.L1.fc_0.25.permut.out", header = TRUE, sep = "\t", comment.char = "")

# Set factor levels for tissue, strain, cell_type, and trait
fisher_test_df$tissue <- factor(fisher_test_df$tissue, levels = tissue_order)
fisher_test_df$strain <- factor(fisher_test_df$strain, levels = c("C57BL/6", "SS", "SHR"))
fisher_test_df$cell_type <- factor(fisher_test_df$cell_type, levels = cell_order)
fisher_test_df$trait <- factor(fisher_test_df$trait, 
                               levels = c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", 
                                          "stroke", "CAD", "eGFR", "BUN", "albuminuria"))

# Calculate DEG proportion and adjust p-values for zero entries
fisher_test_df$DEG_prop <- fisher_test_df$X.SNP.DEG / fisher_test_df$X.DEG
epsilon <- 1e-6
fisher_test_df$p.value <- ifelse(fisher_test_df$p.value == 0, epsilon, fisher_test_df$p.value)

# Filter significant results and those with NES values
test_df_use <- fisher_test_df[fisher_test_df$p.adj < 0.05 & fisher_test_df$X.SNP.DEG >= 5, ]
test_df_use <- test_df_use[!is.na(test_df_use$NES), ]

# Plot the dot plot
p <- ggplot() +
  geom_point(data = test_df_use[test_df_use$p.adj >= 0.01 & test_df_use$X.SNP.DEG >= 1, ],
             aes(x = strain, y = cell_type, size = X.SNP.DEG / 10, fill = NES),
             shape = 21, stroke = 0, color = "black") +
  geom_point(data = test_df_use[test_df_use$p.adj < 0.01 & test_df_use$X.SNP.DEG >= 1, ],
             aes(x = strain, y = cell_type, size = X.SNP.DEG / 10, fill = NES),
             shape = 21, stroke = 0.3, color = "black") +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks = c(5, 10, 15), labels = c(50, 100, 150)) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.spacing.y = unit(0.2, "lines"),
    panel.spacing.x = unit(0.1, "lines"),
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
    strip.text.x = element_text(color = 'black', angle = 90, hjust = 0, margin = margin(l = 0))
  ) +
  labs(fill = "NES", size = "Number of\nDE SNP-gene", y = "", x = "Model") +
  facet_nested(tissue ~ trait, scales = "free", space = "free")

# Display the plot and save as PNG
print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4e.cell_type.gwas.fisher.png", width = 693 / 96, height = 859 / 96, dpi = 300)








################################################################################
### identify cell community based on DE SNP-gene (Figure 4f-g)
################################################################################
fisher_test_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/tmp/trait.snp_gene.L1.fc_0.25.permut.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[fisher_test_df$strain %in% c("C57BL/6", "SHR", "SS"), ]
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]
fisher_test_df <- fisher_test_df %>% mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep = "_"))

snp_gene_df <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_snp_gene.summary.out", sep = "\t", header = T)
snp_gene_df_mod <- snp_gene_df %>%
  separate_longer_delim(
    col = Gene_symbol_Mouse_Rat,
    delim = ";"
  )

### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >= 0.25)
deg_use = deg_use = deg_filtered[deg_filtered$gene_name %in% snp_gene_df_mod$Gene_symbol_Mouse_Rat, ]

binary_matrix <- deg_use %>%
  mutate(condition = paste(strain, treatment, tissue, cell_type, sep = "_"),
         weight = sign(avg_log2FC) * -log10(p_val_adj)) %>%
  dplyr::select(gene_name, condition, weight) %>%
  dplyr::group_by(gene_name, condition) %>%
  dplyr::summarise(
    value = mean(weight, na.rm = TRUE),  # or sum(), max()
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = condition,
    values_from = value,
    values_fill = 0
  ) %>%
  as.data.frame()

row.names(binary_matrix) <- binary_matrix$gene_name
binary_matrix <- binary_matrix %>% dplyr::select(-gene_name) %>% as.matrix()

X <- t(as.matrix(binary_matrix))
sim <- proxy::simil(X, method = "cosine")
similarity_matrix <- as.matrix(sim)
similarity_matrix[is.na(similarity_matrix)] <- 0
similarity_matrix <- pmax(pmin(similarity_matrix, 1), 0)

use_knn <- TRUE
k <- 20

if (use_knn) {
  n <- nrow(similarity_matrix)
  knn_mask <- matrix(0, n, n)
  for (i in seq_len(n)) {
    ord <- order(similarity_matrix[i, ], decreasing = TRUE)
    ord <- ord[ord != i]
    keep <- head(ord, k)
    knn_mask[i, keep] <- 1
  }
  knn_mask <- (knn_mask + t(knn_mask)) > 0
  similarity_matrix <- similarity_matrix * knn_mask
}

diag(similarity_matrix) <- 0  # remove self-loops

similarity_graph <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "undirected", weighted = TRUE)

### consensus clustering
n_runs <- 1000
all_clusterings <- matrix(0, nrow = length(V(similarity_graph)), ncol = n_runs)

for (i in 1:n_runs) {
  all_clusterings[, i] <- leiden(similarity_graph, resolution_parameter = 1)
}

consensus_matrix <- matrix(0, nrow = length(V(similarity_graph)), ncol = length(V(similarity_graph)))
for (i in 1:n_runs) {
  clustering <- all_clusterings[, i]
  for (j in 1:length(clustering)) {
    for (k in j:length(clustering)) {
      if (clustering[j] == clustering[k]) {
        consensus_matrix[j, k] <- consensus_matrix[j, k] + 1
        consensus_matrix[k, j] <- consensus_matrix[k, j] + 1
      }
    }
  }
}
consensus_matrix <- consensus_matrix / n_runs
hclust_result <- hclust(as.dist(1 - consensus_matrix), method = "ward.D2")
consensus_clusters <- cutree(hclust_result, k = 6) 
names(consensus_clusters) = colnames(binary_matrix)

fisher_test_df <- fisher_test_df %>%
  mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep = "_")) %>%
  filter(cell_group %in% names(consensus_clusters)) %>%
  mutate(consensus_cluster = consensus_clusters[cell_group])

write.table(fisher_test_df, "trait.fisher.cluster.out", col.names = T, row.names = F, sep = "\t", quote=F)





# Summarize the NES scores for each community-trait combination
fisher_test_df = read.table("trait.fisher.cluster.out", sep="\t", header=T)

trait_order = c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", "CAD", "eGFR", "BUN", "albuminuria")

fisher_test_df = fisher_test_df %>%
  mutate(cluster_use = consensus_cluster)

nes_summary <- fisher_test_df %>%
  group_by(cluster_use, trait) %>%
  summarise(average_NES = mean(NES, na.rm = TRUE), 
            community_size = n(),
            NES_values = list(NES), .groups = 'drop') %>%
  ungroup()

nes_summary <- nes_summary %>%
  rowwise() %>%
  mutate(
    p_value = {
      community_nes <- unlist(NES_values)
      current_trait <- trait
      current_community <- cluster_use
      other_nes <- fisher_test_df %>%
        filter(trait == current_trait, cluster_use != current_community) %>%
        pull(NES)
      
      if (length(community_nes) > 1 && length(other_nes) > 1) {
        wilcox.test(community_nes, other_nes, alternative = "greater")$p.value
      } else if (length(community_nes) > 1 && length(other_nes) == 0) {
        NA  
      } else {
        1  
      }
    }
  ) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

nes_summary <- nes_summary %>%
  mutate(significant = ifelse(p_adj < 0.05, "*", ""),
         community = cluster_use,
         trait = factor(trait, levels = rev(trait_order)))

p1 <- ggplot(nes_summary, aes(x = trait, y = factor(community), size = community_size, fill = average_NES)) +
  geom_point(shape = 21, stroke = 0) +
  geom_point(data = nes_summary[nes_summary$p_adj < 0.05, ],
             aes(x = trait, y = factor(community), size = community_size, fill = average_NES),
             shape = 21, stroke = 1, color = "black") +
  scale_size_continuous(range = c(2, 6), breaks = c(10, 30, 50)) +
  scale_fill_gradient2(low = "lightblue", mid = "white", high = "red", midpoint = 2.5, name = "Average NES") +
  labs(title = " ",
       x = "Trait",
       y = "",
       size = "Cluster size",
       fill = "Average NES") +
  theme_minimal() +
  theme(axis.text.y = element_text(colour = 'black'),
        # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        axis.text.x = element_blank()
  ) +
  coord_flip()
print(p1)





df_strain_prop <- fisher_test_df %>%
  count(consensus_cluster, strain, name = "n") %>%
  group_by(consensus_cluster) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

df_tissue_prop <- fisher_test_df %>%
  count(consensus_cluster, tissue, name = "n") %>%
  group_by(consensus_cluster) %>%
  mutate(
    prop = n / sum(n)
  ) %>%
  ungroup()

p_strain <- ggplot(df_strain_prop, aes(x = factor(consensus_cluster), y = prop, fill = strain)) +
  geom_col(width = 0.8) +
  labs(x = "", y = "Proportion", fill = "Strain") +
  theme_classic() + 
  scale_fill_manual(values = strain_col) +
  theme(axis.text.x = element_blank())

p_tissue <- ggplot(df_tissue_prop, aes(x = factor(consensus_cluster), y = prop, fill = tissue)) +
  geom_col(width = 0.8) +
  labs(x = "", y = "Proportion", fill = "Tissue") +
  theme_classic() + 
  scale_fill_manual(values = tissue_col) +
  theme(axis.text.x = element_blank())

p_strain
p_tissue



merge_all = read.table("metascape.clusters.out", sep='\t', header=T, quote = "", check.names = F)
top_list = merge_all %>%
  filter(`Log(q-value)` < log10(0.05)) %>%
  group_by(Description) %>%
  mutate(pathway_count = n_distinct(cluster),
         gene_count = gsub("/.*", "", InTerm_InList)) %>%
  filter(pathway_count<3) %>% ungroup() %>%
  filter(grepl("Summary", GroupID)) %>%
  group_by(cluster) %>% arrange(`Log(q-value)`) %>%
  slice_head(n=5)
merge_use = merge_all[merge_all$Description %in% top_list$Description & 
                        merge_all$`Log(q-value)` < log10(0.05) &
                        grepl("Member", merge_all$GroupID), c("pathway", "cluster", "Log(q-value)")]
sorted_df <- merge_use[order(merge_use$cluster, merge_use$pathway), ]
sorted_df$pathway <- factor(sorted_df$pathway, levels = unique(sorted_df$pathway))

p2 = ggplot(sorted_df, aes(y = pathway, x = cluster, fill = -1*`Log(q-value)`)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low="white", high="purple") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black', size=10)) +
  labs(x="", y="Pathway", fill="-log(q-value)")


design = "A
          B
          C
          D"
combined_plot <- (p_strain + p_tissue + p1 + p2) + 
  plot_layout(design = design, guides = 'collect', 
              heights = c(0.05, 0.05, 0.12, 0.3)) &
  theme(legend.position = 'right', 
        legend.box = 'vertical',
        plot.margin = margin(t = 10, r = 0, b = -10, l = 10, unit = "pt"),
        text = element_text(family = "Arial"))
print(combined_plot)
ggsave("cell_cluster.asso.png", width=600/96, height=800/96, dpi=300)
