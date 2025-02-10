calculate_enhanced_index <- function(deg_data) {
  # Calculate the proportion of conditions with significant changes (consistency metric)
  consistency <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::summarise(consistency_metric = mean(p_val_adj < 0.05, na.rm = TRUE))
  # Normalize the fold changes across tissues/models using z-score transformation
  deg_data <- deg_data %>%
  group_by(tissue, strain) %>%
  dplyr::mutate(z_score_fc = scale(-avg_log2FC)) %>%
  ungroup()
  # Calculate the standard deviation of z-scored fold changes across tissues/models
  deg_data <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::mutate(variance_fc = ifelse(n() == 1, 0, sd(-avg_log2FC, na.rm = TRUE))) %>%
  ungroup()
  # Join the consistency metric back to the main data
  deg_data <- deg_data %>%
  left_join(consistency, by = "gene_name")
  # Handling p_val_adj = 1 by replacing it with a small number slightly less than 1
  deg_data <- deg_data %>%
  dplyr::mutate(adjusted_p_val = ifelse(p_val_adj == 1, 0.9999, p_val_adj))
  # Combine the normalized fold change, variance, and consistency into an overall index
  index <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::summarise(
  weighted_z_score_fc = sum(-avg_log2FC * -log10(adjusted_p_val), na.rm = TRUE) / sum(-log10(adjusted_p_val), na.rm = TRUE),
  # robustness_penalty = 1 / (1 + variance_fc),
  robustness_penalty = 1 / (1 + variance_fc),
  consistency_metric = first(consistency_metric),
  # index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * robustness_penalty * consistency_metric)
  index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * robustness_penalty)
  ) %>%
  dplyr::mutate(index_value = ifelse(is.nan(index_value), 0, index_value)) %>%
  unique()
  return(index)
}
calculate_enhanced_index_zscore <- function(deg_data) {
  # Calculate the proportion of conditions with significant changes (consistency metric)
  consistency <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::summarise(consistency_metric = mean(p_val_adj < 0.05, na.rm = TRUE),
  consistency_count = sum(p_val_adj < 0.05, na.rm = TRUE))
  # Normalize the fold changes across tissues/models using z-score transformation
  deg_data <- deg_data %>%
  group_by(tissue, strain) %>%
  dplyr::mutate(z_score_fc = scale(-avg_log2FC)) %>%
  ungroup()
  # Calculate the standard deviation of z-scored fold changes across tissues/models
  deg_data <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::mutate(variance_fc = ifelse(n() == 1, 0, sd(z_score_fc, na.rm = TRUE))) %>%
  ungroup()
  # Join the consistency metric back to the main data
  deg_data <- deg_data %>%
  left_join(consistency, by = "gene_name")
  # Handling p_val_adj = 1 by replacing it with a small number slightly less than 1
  deg_data <- deg_data %>%
  dplyr::mutate(adjusted_p_val = ifelse(p_val_adj == 1, 0.9999, p_val_adj))
  # Combine the normalized fold change, variance, and consistency into an overall index
  index <- deg_data %>%
  group_by(gene_name) %>%
  dplyr::summarise(
  weighted_z_score_fc = sum(z_score_fc * -log10(adjusted_p_val), na.rm = TRUE) / sum(-log10(adjusted_p_val), na.rm = TRUE),
  # robustness_penalty = 1 / (1 + variance_fc),
  robustness_penalty = 1 / (1 + variance_fc),
  consistency_metric = first(consistency_metric),
  # index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * robustness_penalty * consistency_metric)
  index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * consistency_metric)
  ) %>%
  dplyr::mutate(index_value = ifelse(is.nan(index_value), 0, index_value)) %>%
  unique()
  return(index)
}


library(dplyr)
library(ggplot2)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))
hg38_mm10 <- rtracklayer::import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg38ToMm10.over.chain")
hg38_rn7 <- rtracklayer::import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg38ToRn7.over.chain")
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID
################################################################################
### gene-level scores
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_hypt <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_ctrl <- deg_merged[deg_merged$strain %in% c("SD", "WKY"), ]
enhanced_index_hypt <- calculate_enhanced_index_zscore(deg_hypt)
enhanced_index_ctrl <- calculate_enhanced_index_zscore(deg_ctrl)
# enhanced_index_hypt <- calculate_enhanced_index(deg_hypt)
# enhanced_index_ctrl <- calculate_enhanced_index(deg_ctrl)
enhanced_index_ctrl[enhanced_index_ctrl$gene_name=="Nrip1",]
enhanced_index_hypt[enhanced_index_hypt$gene_name=="Nrip1",]
enhanced_index_hypt[enhanced_index_hypt$gene_name %in% c("Slc5a3", "Kcne2", "Mrps6"),]
enhanced_index_ctrl[enhanced_index_ctrl$gene_name %in% c("Slc5a3", "Kcne2", "Mrps6"),]
total_tissue_cell_combinations = nrow(unique(deg_hypt[, c("tissue", "cell_type")]))
breadth_score <- deg_hypt %>%
group_by(gene_name) %>%
dplyr::summarise(breadth_score = n_distinct(tissue, cell_type) / total_tissue_cell_combinations)
combined_index <- enhanced_index_hypt %>%
inner_join(enhanced_index_ctrl, by = "gene_name") %>%
inner_join(breadth_score, by = "gene_name") %>%
mutate(combined_score = index_value.x * breadth_score)
combined_index[combined_index$gene_name=="Nrip1",]


gwas_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  dplyr::select(SNPS, trait, P.VALUE)

gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  dplyr::select(SNPS, trait, P.VALUE)

snps_LD = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/SNPS_in_LD.txt", header = T, sep = "\t")
ld_1 = names(table(snps_LD$rsID1)[table(snps_LD$rsID1)==1])
snp_gene_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')
SNPS = SNPS %>%
  mutate(Chromosome = sapply(strsplit(as.character(POS), "_"), `[`, 1),
  Position = as.numeric(sapply(strsplit(as.character(POS), "_"), `[`, 2))) %>%
  filter(rowSums(is.na(.))==0)
gr = GRanges(seqnames = SNPS$Chromosome,
              ranges = IRanges(start = SNPS$Position,
              end = SNPS$Position),
              POS = SNPS$POS)
snp_mm10 <- rtracklayer::liftOver(x = gr, chain = hg38_mm10) %>% unlist()
snp_rn7 <- rtracklayer::liftOver(x = gr, chain = hg38_rn7) %>% unlist()
lifted_mm10 <- data.frame(original_POS = mcols(snp_mm10)$POS,
POS_mm10 = paste0(seqnames(snp_mm10), "_",
start(snp_mm10), "_",
end(snp_mm10)))
lifted_rn7 <- data.frame(original_POS = mcols(snp_rn7)$POS,
POS_mm10 = paste0(seqnames(snp_rn7), "_",
start(snp_rn7), "_",
end(snp_rn7)))
length(unique(lifted_mm10$original_POS))/length(unique(SNPS$SNP))
length(unique(lifted_rn7$original_POS))/length(unique(SNPS$SNP))
# > lifted_rn7[lifted_rn7$original_POS=="chr21_15184047",]
# original_POS                POS_mm10
# 254 chr21_15184047 chr11_15091368_15091368
# Nrip1 in rn7: NC_051346.1 (14895843..14979490, complement)
# 15091368-14979490 = 111878
snp_gene_df_mod = snp_gene_df %>%
separate_rows(combined_gene_name , sep=";") %>%
merge(., SNPS, by="SNP") %>%
merge(., ensembl, by.x="gene", by.y="Gene.stable.ID") %>%
merge(., gwas_merge_use, by.x = "SNP", by.y="SNPS") %>%
merge(., combined_index, by.x = "combined_gene_name", by.y="gene_name") %>%
mutate(TSS = ifelse(Strand == 1, Gene.start..bp., Gene.end..bp.),
Distance_to_TSS = abs(TSS-Position)) %>%
filter(rowSums(is.na(.))==0)
head(snp_gene_df_mod)







