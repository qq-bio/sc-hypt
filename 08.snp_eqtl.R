library(dplyr)
library(tidyr)
library(biomaRt)

# ensmbl_hg38 = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", header = T, sep = '\t')
r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out", header = T, sep = '\t')
# rsid_loc = read.table("/xdisk/mliang1/qqiu/reference/dbSNP/dbSNP.pos.chr.txt", header = T, sep = '\t')
# ensembl <- useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")

gwas_bp = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas_catalog_blood pressure.tsv", header = T, sep = '\t')
gwas_hypt = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas_catalog_essential hypertension.tsv", header = T, sep = '\t')

eqtl_kidney = read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.processed.txt", header = T, sep = '\t')

gwas_bp_mod = gwas_bp %>% filter(pValue < 5e-8) %>%
  separate_rows(mappedGenes, sep=",") %>%
  separate_rows(efoTraits, sep=",") %>%
  filter(efoTraits %in% c("blood pressure", "diastolic blood pressure", "systolic blood pressure", "hypertension")) %>%
  full_join(r2h, join_by(mappedGenes == Human.gene.name), relationship = "many-to-many") %>%
  mutate(SNP = gsub("-.*", "", riskAllele),
         norm_p = -1*log10(pValue),
         category = "GWAS mapped genes",
         efoTraits = str_to_sentence(efoTraits)) %>%
  filter(Gene.name != "" & !is.na(Gene.name) & !is.na(SNP)) %>%
  dplyr::select(SNP, norm_p, mappedGenes, Gene.name, efoTraits, category) %>%
  unique() %>% as.data.frame()

gwas_bp_beta = gwas_bp %>% filter(pValue < 5e-8) %>%
  separate_rows(mappedGenes, sep=",") %>%
  separate_rows(efoTraits, sep=",") %>%
  filter(efoTraits %in% c("blood pressure", "diastolic blood pressure", "systolic blood pressure", "hypertension")) %>%
  full_join(r2h, join_by(mappedGenes == Human.gene.name), relationship = "many-to-many") %>%
  mutate(SNP = gsub("-.*", "", riskAllele),
         norm_p = -1*log10(pValue),
         category = "GWAS mapped genes") %>%
  filter(Gene.name != "" & !is.na(Gene.name) & !is.na(SNP)) %>%
  dplyr::select(riskAllele, SNP, beta, norm_p, mappedGenes, Gene.name, efoTraits, category) %>%
  unique() %>% as.data.frame()

gwas_hypt_mod = gwas_hypt %>% filter(pValue < 5e-8) %>%
  separate_rows(mappedGenes, sep=",") %>%
  separate_rows(efoTraits, sep=",") %>%
  full_join(r2h, join_by(mappedGenes == Human.gene.name), relationship = "many-to-many") %>%
  mutate(SNP = gsub("-.*", "", riskAllele),
         norm_p = -1*log10(pValue),
         # efoTraits = "essential hypertension",
         category = "GWAS mapped genes",
         efoTraits = str_to_sentence(efoTraits)) %>%
  filter(Gene.name != "" & !is.na(Gene.name) & !is.na(SNP)) %>%
  dplyr::select(SNP, norm_p, mappedGenes, Gene.name, efoTraits, category) %>%
  unique() %>% as.data.frame()


# eqtl_kidney_mod = eqtl_kidney %>%
#   mutate(norm_p = -1*log10(pval_beta),
#          efoTraits = "eQTL",
#          category = "eQTL genes") %>%
#   left_join(r2h, join_by(gene_name == Human.gene.name)) %>%
#   filter(Gene.name != "" & !is.na(Gene.name) & rs_id!= "") %>%
#   dplyr::select(rs_id, norm_p, gene_name, Gene.name, efoTraits, category) %>%
#   unique() %>% as.data.frame()
# colnames(eqtl_kidney_mod) = c("SNP", "norm_p", "mappedGenes", "Gene.name", "efoTraits", "category")
# write.table(eqtl_kidney_mod, "/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.processed.add_rat_symbol.txt", col.names = T, row.names = F, sep = "\t", quote = F)
eqtl_kidney_mod = read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.processed.add_rat_symbol.txt", header = T, sep = "\t")
colnames(eqtl_kidney_mod)[c(1, 3)] = c("SNP", "mappedGenes")

selected_genes = c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2", "Chrm3", # circulatory system process
  "Arhgap24", "Efna5", "Igf1r", "Ptprd", # regulation of plasma membrane bounded cell projection organization (GO)
  "Pkhd1" # negative regulation of intracellular signal transduction (GO)
)

snp_use = rbind(
  gwas_bp_mod[gwas_bp_mod$Gene.name %in% selected_genes, ],
  gwas_hypt_mod[gwas_hypt_mod$Gene.name %in% selected_genes, ]
  #eqtl_kidney_mod[eqtl_kidney_mod$Gene.name %in% selected_genes, ]
)
snp_use = snp_use %>% # filter(norm_p>4) %>%
  group_by(SNP, Gene.name, efoTraits) %>%
  filter(norm_p == max(norm_p)) %>%
  ungroup() %>% 
  mutate(efoTraits = factor(efoTraits, 
                            levels = c("Systolic blood pressure", "Diastolic blood pressure", "Hypertension", "eQTL"))) %>%
  arrange(efoTraits, desc(norm_p)) %>% mutate(SNP=factor(SNP, levels=unique(SNP)))
p2 = ggplot(snp_use, aes(x = efoTraits, y = SNP)) +
  geom_tile(aes(fill = norm_p)) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low = "white", high = "purple") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    legend.title.align = 0.5,
    legend.position = "bottom",
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", fill="-log10(p-value)")

gene_names <- layer_scales(p1)$y$range$range
total_genes <- length(gene_names)
relative_positions <- seq(from = 0.5/total_genes, to = (total_genes-0.5)/total_genes, length.out = total_genes)
gene_positions <- relative_positions; names(gene_positions) = rev(gene_names)

snp_names <- layer_scales(p2)$y$range$range
total_snps <- length(snp_names)
relative_positions <- seq(from = 0.5/total_snps, to = (total_snps-0.5)/total_snps, length.out = total_snps)
snp_positions <- relative_positions; names(snp_positions) = rev(snp_names)

snp_gene_pair = unique(snp_use[,c("SNP", "Gene.name", "category")])
snp_gene_pair$y1 = gene_positions[snp_gene_pair$Gene.name]
snp_gene_pair$y2 = snp_positions[as.character(snp_gene_pair$SNP)]
snp_gene_pair$x1 = 0
snp_gene_pair$x2 = 1

p3 = ggplot(snp_gene_pair, aes(x = x1, y = y1, xend = x2, yend = y2)) +
  geom_segment(color = "black", size = 0.2) +  # Color lines by gene, adjust size as needed
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +  # Apply a theme that removes most elements
  theme(legend.title = element_text(size = 10),
        legend.title.align = 0.5)

p1 <- p1 + theme(plot.margin = unit(c(0,0,0,0), "cm"))
p2 <- p2 + theme(plot.margin = unit(c(0,0,0,-0.3), "cm"))
p3 <- p3 + theme(plot.margin = unit(c(0,0,0,0), "cm"))

combined_plot <- p1 + p3 + p2 + plot_layout(guides = 'collect', widths = c(1, 0.1, 0.08))
print(combined_plot)



snp_use_beta = gwas_bp_beta[gwas_bp_beta$SNP %in% snp_use$SNP, ]
snp_use_beta[order(snp_use_beta$SNP), ]


dat = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/SNP.gnomAD.out", header=T, sep='\t')

sample_size_threshold <- 0.05

snp_significance <- dat %>%
  filter(Category=="race") %>%  
  group_by(SNP, Genetic.Ancestry.Group) %>%
  summarise(
    Allele_Count = sum(Allele.Count),  
    Total = sum(Allele.Number)  
  ) %>%
  filter(Total > mean(Total) * sample_size_threshold) %>%  
  mutate(Mean_Frequency = Allele_Count / Total) %>%
  group_by(SNP) %>%
  mutate(
    SNP_Mean = mean(Mean_Frequency),
    SD_Frequency = sd(Mean_Frequency)
  ) %>%
  ungroup() %>%
  mutate(Significance = case_when(
    (Mean_Frequency - SNP_Mean) > 2*SD_Frequency ~ "High", # 0.05
    (Mean_Frequency - SNP_Mean) < -2*SD_Frequency ~ "Low",
    TRUE ~ "Baseline"
  ))

snp_ids <- snp_significance %>%
  filter(Significance %in% c("High", "Low")) %>%
  distinct(SNP)

filtered_data <- snp_significance %>%
  filter(SNP %in% snp_ids$SNP)

ggplot(filtered_data, aes(x = Genetic.Ancestry.Group, y = Mean_Frequency, fill = Significance)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(SNP ~ ., scales = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = unit(c(0,0,0,1), "cm"),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "right",
        axis.text = element_text(colour = 'black'),
        strip.background = element_rect(colour = "black", fill = NA)) +
  labs(x = "Genetic ancestry group", y = "Mean allele frequency", fill="Relative\nfrequency") +
  scale_fill_manual(values = c("High" = "red", "Low" = "blue", "Baseline" = "grey"))


gwas_bp_beta[gwas_bp_beta$SNP %in% snp_ids$SNP, ]
