
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
### DEG Analysis for Each Trait (Figure 4b and S13a)
################################################################################

# Define function for one-sided z-test
one_sided_z_test <- function(proportion_hit_genes, proportion_background, n_snp_gene_list) {
  z_score <- (proportion_hit_genes - proportion_background) / sqrt(proportion_background * (1 - proportion_background) / n_snp_gene_list)
  p_value <- pnorm(z_score, lower.tail = FALSE)  # one-sided greater test
  return(p_value)
}

# Load necessary data files
gwas_merge <- read.table("gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')
background_500k <- read.table("proximal_merged.500k.txt", header = TRUE, sep = '\t')

proximal_merge <- proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge <- eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge <- e2g_merge[grepl("rs", e2g_merge$SNP),]
background_500k <- background_500k[grepl("rs", background_500k$SNP),]

# Unique SNP-gene pairs for each dataset
snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)

snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))
snp_gene_500k <- data.frame(SNP = background_500k$SNP, gene = background_500k$Gene_ID)
background_df <- unique(rbind(snp_gene_500k, snp_gene_eqtl, snp_gene_e2g))

# Merge SNP-gene with ensembl data
snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)

# Load DEG data and filter by tissue and strain
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]

# Define DEG thresholds
thresholds <- list(
  "p.adj < 0.05" = deg_merged %>% filter(p_val_adj < 0.05),
  "p.adj < 0.05 & |log2(FC)| > 0.25" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25),
  "p.adj < 0.05 & |log2(FC)| > 0.5" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
  "p.adj < 0.05 & |log2(FC)| > 1" = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) >1)
)

# Analysis loop for each trait
results <- data.frame()
for (trait in unique(gwas_merge$trait)) {
  SNPS <- unique(gwas_merge[gwas_merge$trait == trait, ]$SNPS)
  
  snp_gene_df_use <- unique(snp_gene_df[snp_gene_df$SNP %in% SNPS, ])
  snp_gene_list <- unique(snp_gene_df_use$gene)
  
  for (threshold_name in names(thresholds)) {
    deg_filtered <- thresholds[[threshold_name]]
    mapped_gene_list <- unique(c(snp_gene_df_use[!is.na(snp_gene_df_use$gene_name_mouse), ]$gene_name_mouse,
                                 snp_gene_df_use[!is.na(snp_gene_df_use$gene_name_rat), ]$gene_name_rat))
    mapped_gene_list <- setdiff(mapped_gene_list, "")
    deg_list <- unique(deg_filtered$gene_name)
    
    background_gene_list <- unique(c(m2h[m2h$gene_id %in% background_df$gene, ]$gene_name_mouse,
                                     r2h[r2h$gene_id %in% background_df$gene, ]$gene_name_rat))
    proportion_background <- length(intersect(background_gene_list, deg_list)) / length(unique(background_df$gene))
    
    hit_genes <- intersect(mapped_gene_list, deg_list)
    num_hit_genes <- length(hit_genes)
    proportion_hit_genes <- num_hit_genes / length(snp_gene_list)
    num_DEG_snps <- length(unique(snp_gene_df_use[(!is.na(snp_gene_df_use$gene_name_mouse) & snp_gene_df_use$gene_name_mouse %in% deg_list) | 
                                                    (!is.na(snp_gene_df_use$gene_name_rat) & snp_gene_df_use$gene_name_rat %in% deg_list), ]$SNP))
    
    p_value <- one_sided_z_test(proportion_hit_genes, proportion_background, length(snp_gene_list))
    
    results <- rbind(results, data.frame(
      trait = trait,
      SNPS = length(unique(snp_gene_df_use$SNP)),
      SNP_genes = length(snp_gene_list),
      mapped_genes = length(mapped_gene_list),
      threshold = threshold_name,
      DEGs = num_hit_genes,
      DEG_SNPS = num_DEG_snps,
      proportion_DEG = proportion_hit_genes,
      proportion_background = proportion_background,
      proportion_SNP = num_DEG_snps / length(unique(snp_gene_df_use$SNP)),
      p_value = p_value
    ))
    
  }
}

# Adjust p-values and add labels
results$p.adj <- p.adjust(results$p_value, method = "BH")
results$label <- paste(results$DEGs, 
                       "(", round(results$proportion_DEG * 100, 2), "%)", 
                       ifelse(results$p.adj < 0.05, "*", ""), sep = "")

results$threshold <- factor(results$threshold, levels = rev(names(thresholds)))

write.table(results, "trait.deg_prop.z_test.out", sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)




### visualize DEG-hit results
results_use <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/trait.deg_prop.z_test.out", sep = ",", header = TRUE)

# Update threshold labels for better text wrapping in plots
results_use$threshold <- gsub("& ", "&\n", results_use$threshold)
threshold_col <- RColorBrewer::brewer.pal(5, "Blues")[2:5]
names(threshold_col) <- unique(results_use$threshold)

# Reorder traits by total DEG hits for improved visualization order
results_use$trait <- factor(results_use$trait, levels = results_use %>%
                              group_by(trait) %>%
                              summarise(total_hits = sum(proportion_DEG)) %>%
                              arrange(total_hits) %>%
                              pull(trait))

# Add dashed lines indicating background proportion by threshold
vline_data <- results_use %>%
  group_by(threshold) %>%
  summarise(yintercept = unique(porpotion_background))

# Generate bar plot with DEG proportions by trait and threshold (Figure S13a)
ggplot(results_use, aes(x = trait, y = proportion_DEG, fill = threshold, label = label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = threshold_col) +
  geom_text(aes(x = trait, y = proportion_DEG, hjust = 0), position = position_dodge(0.9), size = 3) +
  labs(x = "Trait", y = "DEG Proportion", fill = "DEG Threshold") +
  scale_y_continuous(limits = c(0, 0.85)) +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')
  ) +
  coord_flip() +
  facet_wrap(~threshold, scale = "free_x", nrow = 1) +
  geom_hline(data = vline_data, aes(yintercept = yintercept), linetype = "dashed", color = "grey")




# Filter results for specific threshold (p.adj < 0.05 & |log2(FC)| > 0.25)
results_threshold <- results_use[results_use$threshold == "p.adj < 0.05 &\n|log2(FC)| > 0.25", ]
results_threshold$trait <- factor(results_threshold$trait, levels = results_threshold %>%
                                    group_by(trait) %>%
                                    summarise(total_hits = sum(proportion_SNP)) %>%
                                    arrange(total_hits) %>%
                                    pull(trait))

# Prepare background line data for this threshold subset
vline_data <- results_threshold %>%
  group_by(threshold) %>%
  summarise(yintercept = unique(porpotion_background))

# Generate line plot with SNP-related gene and SNP proportions by trait
p <- ggplot(results_threshold) +
  geom_segment(aes(x = trait, xend = trait, y = proportion_DEG, yend = proportion_SNP), color = "gray", size = 1) +
  geom_point(aes(x = trait, y = proportion_DEG, label = DEGs, color = "SNP-related genes"), size = 3) +
  geom_point(aes(x = trait, y = proportion_SNP, label = DEG_SNPS, color = "SNPs"), size = 3) +
  geom_text(aes(x = trait, y = proportion_DEG - 0.12, label = DEGs), color = "cadetblue4", size = 3) +
  geom_text(aes(x = trait, y = proportion_SNP + 0.12, label = DEG_SNPS), color = "chocolate4", size = 3) +
  scale_color_manual(name = "", values = c("SNP-related genes" = "cadetblue4", "SNPs" = "chocolate4")) +
  scale_y_continuous(limits = c(0.05, 0.85)) +
  labs(x = "Trait", y = "Proportion Hit by DEGs") +
  theme(
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(colour = 'black')
  ) +
  coord_flip()

# Save the plot
print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4b.deg_snp_number.png", width = 460 / 96, height = 178 / 96, dpi = 300)






################################################################################
### Cell Type Enrichment Analysis od SNP-related DEG
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")

gwas_merge <- read.table("gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t')
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)

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
          deg_list <- unique(deg_merged_use[deg_merged_use$p_val_adj < 0.05 & abs(deg_merged_use$avg_log2FC) > 0.25, ]$gene_name)
          
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
write.table(fisher_test_df, "trait.fisher.snp_gene.enrichment_results.w_SNP.out", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)





################################################################################
### Visualize Fisher Test Results Using Dot Plot (Figure 4e)
################################################################################

# Load the Fisher test results
fisher_test_df <- read.table("trait.fisher.snp_gene.enrichment_results.out", header = TRUE, sep = "\t", comment.char = "")

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
### Nominating SNP-Gene Pairs in Trait Associated-Cell Types (Figure 4f and S14a)
################################################################################
snp_gene_df <- read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
gwas_merge <- read.table("gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)

deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC)>0.25)

gwas_merge_use <- gwas_merge[gwas_merge$trait %in% trait, ] %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

generate_dotplot <- function(trait_var, cell_type_var, tissue_var, gwas_data, snp_gene_data, deg_data) {
  
  gwas_filtered <- gwas_data %>%
    filter(trait %in% trait) %>%
    group_by(SNPS) %>%
    summarise(cumulative_significance = -sum(log10(P.VALUE)))
  
  gene_snp_significance <- snp_gene_data %>%
    separate_rows(combined_gene_name, sep = ";") %>%
    inner_join(gwas_filtered, by = c("SNP" = "SNPS")) %>%
    group_by(gene, Gene.name, combined_gene_name) %>%
    summarise(cumulative_significance = sum(cumulative_significance))
  
  deg_filtered <- deg_data %>%
    filter(cell_type == cell_type_var & tissue == tissue_var) %>%
    mutate(weight = -log10(p_val_adj)) %>%
    group_by(gene_name) %>%
    summarise(weighted_log2FC = sum(-avg_log2FC * weight) / sum(weight))
  
  final_data <- gene_snp_significance %>%
    inner_join(deg_filtered, by = c("combined_gene_name" = "gene_name"))
  
  top_genes <- final_data[order(final_data$cumulative_significance, decreasing = T),][1:10,]
  
  p <- ggplot(final_data, aes(x = weighted_log2FC, y = cumulative_significance)) +
    geom_vline(xintercept = 0, color="darkgrey") +
    geom_hline(yintercept = 0, color="darkgrey") +
    geom_point(aes(color = weighted_log2FC), size = 3) +
    geom_point(data = top_genes, color = "black", size = 5, shape = 21) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(x = "Weighted log2FC\n(hypertension vs. normotension)", y = "Cumulative significance of related SNPs",
         title = paste0(tissue_var, " - ", cell_type_var, "\n(", paste(trait_var, collapse = "/"), ")")) +
    theme_minimal() +
    theme(legend.position = "none") +
    geom_text_repel(data = top_genes, aes(label = combined_gene_name), 
                    nudge_y = 0.05, nudge_x = 0.05, 
                    box.padding = 0.35, point.padding = 0.3, 
                    segment.color = 'grey50')
  
  print(p)
}

generate_dotplot(
  trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"),
  cell_type_var = "EC",
  tissue_var = "HYP",
  gwas_data = gwas_merge_use,
  snp_gene_data = snp_gene_df,
  deg_data = deg_filtered
)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4f.hyp_ec.npr3.png", width=400/96, height=329/96, dpi=300)



combinations <- list(
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Astrocyte", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Myelinating OL", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Microglia", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "CM", tissue_var = "LV"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure", "CAD"), cell_type_var = "Fibroblast", tissue_var = "LV"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "EC", tissue_var = "LK"),
  list(trait_var = c("BUN"), cell_type_var = "TAL", tissue_var = "LK"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "VSMC", tissue_var = "MSA")
)

plots <- list()

# Iterate through each combination
for (i in seq_along(combinations)) {
  comb <- combinations[[i]]
  
  plot <- generate_dotplot(
    trait = comb$trait_var,
    cell_type_var = comb$cell_type_var,
    tissue_var = comb$tissue_var,
    gwas_data = gwas_merge_use,
    snp_gene_data = snp_gene_df,
    deg_data = deg_filtered
  )
  
  plots[[paste(comb$cell_type_var, comb$tissue_var, paste(comb$trait_var, collapse = "_"), sep = "_")]] <- plot
}

# Combine all plots into a single layout
combined_plot <- wrap_plots(plots, ncol = 3) 

print(combined_plot)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/figs14a.trait_cell_type.top10.png", width=1400/96, height=1200/96, dpi=300)









