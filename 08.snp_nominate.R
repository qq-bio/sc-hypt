
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

# deg_hypt[deg_hypt$gene_name=="Nrip1" & deg_hypt$p_val_adj<0.05,]

################################################################################
################################################################################
### load LD information
gwas_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
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
                         POS_rn7 = paste0(seqnames(snp_rn7), "_", 
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

# snp_gene_df_mod[snp_gene_df_mod$Gene.name.x=="NRIP1" & snp_gene_df_mod$SNP=="rs1882961",]


snp_gene_df_filter = snp_gene_df_mod %>% 
  filter(trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"))
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()

snp_gene_df_filter = snp_gene_df_filter %>%
  filter(SNP %in% ld_1) 
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()

snp_gene_df_filter = snp_gene_df_filter %>%
  filter(POS %in% lifted_mm10$original_POS &
         POS %in% lifted_rn7$original_POS)
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()

snp_gene_df_filter = snp_gene_df_filter %>%
  group_by(SNP) %>% dplyr::mutate(count_gene = n_distinct(combined_gene_name),
                                  count_trait = n_distinct(trait),
                                  hypt_specific = weighted_z_score_fc.y*weighted_z_score_fc.x<0) %>% 
  ungroup() %>%
  filter(hypt_specific)
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()

snp_gene_df_filter = snp_gene_df_filter %>%
  group_by(SNP) %>% dplyr::mutate(n_gene_per_SNP = n_distinct(gene)) %>% ungroup() %>%
  group_by(gene) %>% dplyr::mutate(n_SNP_per_gene = n_distinct(SNP)) %>% ungroup() %>%
  filter(n_SNP_per_gene==1)
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()

snp_gene_df_filter = snp_gene_df_filter %>%
  filter(Distance_to_TSS>100000)
snp_gene_df_filter %>% dplyr::select(SNP, combined_gene_name) %>% unique() %>% nrow()



snp_gene_df_filter = snp_gene_df_mod %>% 
  filter(SNP %in% ld_1, 
         POS %in% lifted_mm10$original_POS &
         POS %in% lifted_rn7$original_POS,
         trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"),
         Distance_to_TSS>100000) %>%
  group_by(SNP) %>% dplyr::mutate(n_gene_per_SNP = n_distinct(gene)) %>% ungroup() %>%
  group_by(gene) %>% dplyr::mutate(n_SNP_per_gene = n_distinct(SNP)) %>% ungroup() %>%
  filter(n_SNP_per_gene==1) %>%
  group_by(SNP) %>% dplyr::mutate(count_gene = n_distinct(combined_gene_name),
                                  count_trait = n_distinct(trait),
                                  hypt_specific = weighted_z_score_fc.y*weighted_z_score_fc.x<0) %>% ungroup() %>%
  dplyr::select(SNP, POS, Chromosome, Position, trait, P.VALUE, 
                combined_gene_name, weighted_z_score_fc.x, index_value.x, breadth_score,
                weighted_z_score_fc.y, index_value.y, count_gene, count_trait, hypt_specific) %>% unique()

# snp_gene_df_filter = snp_gene_df_filter %>%
#   left_join(lifted_mm10, by = c("POS" = "original_POS")) %>%
#   left_join(lifted_rn7, by = c("POS" = "original_POS"))

# length(unique(snp_gene_df_filter$SNP))
# length(unique(snp_gene_df_filter$gene))

snp_gene_df_filter  = rbind(unique(snp_gene_df_filter[snp_gene_df_filter$combined_gene_name!="Nrip1",]), unique(snp_gene_df_filter[snp_gene_df_filter$combined_gene_name=="Nrip1",]))

n_counts <- snp_gene_df_filter %>% 
  filter(breadth_score > 0.5 & weighted_z_score_fc.x > 0) %>% 
  group_by(trait) %>% 
  dplyr::summarise(n_count = n()) %>% 
  dplyr::mutate(label = paste0("N=", n_count))

snp_gene_df_filter %>% 
  ggplot(., aes(x=breadth_score, y=weighted_z_score_fc.x, 
                color= combined_gene_name=="Nrip1",
                label= combined_gene_name=="Nrip1")) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_point() +
  geom_text_repel(data = . %>% filter(combined_gene_name == "Nrip1"),
                  aes(label = "rs1882961 - Nrip1"),
                  nudge_y = 0.1, nudge_x = 0.4, 
                  size = 4, color = "red") +
  scale_color_manual(values=c("grey", "red")) +
  labs(x = "Expression breadth score", y="Weighted fold change\n(hypertension vs. normotension)") +
  theme(legend.position = "None") +
  facet_wrap(~trait) +
  geom_text(data = n_counts, aes(x = 1, y = 2.5, label = label), size = 4, color = "black")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/snp_gene_nomination.100k.png", width=565/96, height=246/96, dpi=300)


n_counts <- snp_gene_df_filter %>% 
  filter(trait=="systolic_bp" & breadth_score > 0.5 & weighted_z_score_fc.x > 0) %>% 
  summarise(n_count = n()) %>% 
  mutate(label = paste0("N=", n_count))

snp_gene_df_filter %>% 
  filter(trait=="systolic_bp") %>%
  ggplot(., aes(x=breadth_score, y=weighted_z_score_fc.x, 
                color= combined_gene_name=="Nrip1",
                label= combined_gene_name=="Nrip1")) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_point() +
  geom_text_repel(data = . %>% filter(combined_gene_name == "Nrip1"),
                  aes(label = "rs1882961 - Nrip1"),
                  nudge_y = 1, nudge_x = 0.4, 
                  size = 4, color = "red") +
  scale_color_manual(values=c("grey", "red")) +
  labs(x = "Expression breadth score", y="Weighted fold change\n(hypertension vs. normotension)") +
  theme(legend.position = "None") +
  geom_text(data = n_counts, aes(x = 1, y = 2.2, label = label), size = 4, color = "black")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/snp_gene_nomination.sys_bp.png", width=245/96, height=265/96, dpi=300)



################################################################################
snp_gene_df_filter  = unique(snp_gene_df_filter[order(snp_gene_df_filter$SNP=="rs28451064"),])
snp_gene_df_filter_org = snp_gene_df_filter
snp_gene_df_filter = snp_gene_df_filter_org %>% 
  group_by(SNP) %>%
  dplyr::mutate(count_gene = n_distinct(combined_gene_name),
              count_trait = n_distinct(trait))
  

n_counts <- snp_gene_df_filter %>% 
  ungroup() %>%
  filter(weighted_z_score_fc.y * weighted_z_score_fc.x < 0) %>% 
  group_by(trait) %>%
  dplyr::summarise(n_count = n(),
                   n_snp = n_distinct(SNP)) %>% 
  dplyr::mutate(label = paste0(# "Hypertension specific:\n",
    n_snp, " SNPs\n",
    n_count, " SNP-gene pairs"))

snp_id = "rs28451064"
# snp_id = "rs1882961"
p = snp_gene_df_filter %>% 
  # ggplot(., aes(x=index_value.x, y=index_value.y)) + 
  ggplot(., aes(x=weighted_z_score_fc.x, y=weighted_z_score_fc.y)) + 
  geom_rect(aes(xmin = 0, xmax = -Inf, ymin = 0, ymax = Inf), 
            fill = "#FDEDEC", color = NA) + 
  geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = -Inf), 
            fill = "#FDEDEC", color = NA) + 
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(aes(fill= SNP==snp_id,
                 label= SNP==snp_id),
             color= "black", shape=21, alpha=0.8, size=2) +
  geom_text_repel(data = . %>% filter(SNP==snp_id, trait=="diastolic_bp"),
                  aes(label = paste0(snp_id, " - ", combined_gene_name)),
                  nudge_x = c(2, -2.1),nudge_y = c(-0.35, 0.8), 
                  point.padding = 0.5, size = 4, color = "red") +
  geom_text_repel(data = . %>% filter(SNP==snp_id, trait=="pulse_pressure"),
                  aes(label = paste0(snp_id, " - ", combined_gene_name)),
                  nudge_x = c(2, -2),nudge_y = c(-0.3, 0.35), 
                  point.padding = 0.5, size = 4, color = "red") +
  scale_fill_manual(values=c("grey", "red")) +
  labs(x = "Weighted fold change in hypertensive strain\n(treatment vs. control)", 
       y="Weighted fold change in normotensive strain\n(treatment vs. control)") +
  theme(legend.position = "none") +
  facet_wrap(~trait) +
  geom_text(data = n_counts, aes(x = -4, y = 2.2, label = label), hjust = 0, size = 4, color = "black")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig5a.snp_gene_nom.rs064.png", width=930/96, height=403/96, dpi=300)




snp_gene_df_filter %>% filter(breadth_score>0.5 & weighted_z_score_fc.x>0) %>% group_by(trait) %>% summarise(n_count = n())

snp_gene_df_filter %>% filter(lifted_mm10, lifted_rn7) %>%
  ggplot(., aes(x=breadth_score, y=index_value.x, 
                color=combined_gene_name=="Nrip1",
                size= -log10(P.VALUE))) + geom_point()

snp_gene_df_filter %>% filter(n_SNP_per_gene==1, n_gene_per_SNP==1) %>%
  ggplot(., aes(x=weighted_z_score_fc.x, y=weighted_z_score_fc.y, 
                color=combined_gene_name=="Nrip1",
                size= -log10(P.VALUE))) + geom_point()





################################################################################


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


