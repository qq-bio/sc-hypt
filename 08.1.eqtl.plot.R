library(vcfR)
library(dplyr)
library(tidyr)
library(ggplot2)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))




################################################################################
### load data
consistent_result = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/eqtl.consistent_result.out", sep="\t", header=T)
consistent_result$variant_id = paste(consistent_result$POS, consistent_result$ref_allele, consistent_result$alt_allele, "b38", sep="_")

gwas_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_merged.txt", header = TRUE, sep = '\t', quote = "")

Artery_Aorta = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_artery_aorta.gct", header=T, skip = 2)
Artery_Coronary = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_artery_coronary.gct", header=T, skip = 2)
Artery_Tibial = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_artery_tibial.gct", header=T, skip = 2)
Brain_Hypothalamus = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_brain_hypothalamus.gct", header=T, skip = 2)
Heart_Left_Ventricle = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_heart_left_ventricle.gct", header=T, skip = 2)
Kidney_Cortex = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_kidney_cortex.gct", header=T, skip = 2)

gene_expr_list = list("Artery_Aorta" = Artery_Aorta, "Artery_Coronary" = Artery_Coronary, 
                      "Artery_Tibial" = Artery_Tibial, "Brain_Hypothalamus" = Brain_Hypothalamus, 
                      "Heart_Left_Ventricle" = Heart_Left_Ventricle, "Kidney_Cortex"= Kidney_Cortex)

vcf <- read.vcfR("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/GTEx.eqtl.filtered.vcf.gz")

deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', header=T)
colnames(deg_pseudo)[12] = "strain"
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$strain = factor(deg_pseudo$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "strain", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = c(c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", "Astrocyte", # "Microglia", "Activated microglia", 
                                                             "OPC", "NFO", "Premyelinating OL", "Myelinating OL", "Tanycyte", "Ependymal cell", "Pars tuberalis cell"), c("CM"),
                                                           c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"), c("EC", "VSMC", "E/P transition cell", "Pericyte", "Fibroblast", "Adipocyte"),
                                                           c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils", "NK cells", "NKT", "T cells", "B cells"))
)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$strain = factor(expr_all$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))

################################################################################
### overall plot
consistent_result_plot = consistent_result %>%
  filter(strain %in% c("C57BL/6", "SHR", "SS")) %>%
  mutate(headline = "Gene expression change during hypertension progression",
         expr_trend = ifelse(avg_log2FC>0, "Decrease", "Increase")) %>%
  group_by(trait, SNPS, gene_id_mod) %>%
  arrange(desc(abs(avg_log2FC))) %>%
  slice_head(n=1)

ggplot(consistent_result_plot, aes(x = -avg_log2FC, y = -log10(p_val), color = cell_type)) +
  geom_point() +
  scale_color_manual(values=cell_col) +
  scale_y_continuous(trans='log2') +
  labs(x = "log2FC (hypertension vs. normotension)", y = "-log10(p-value)") +  # No need for fill in labs
  # theme_minimal() +
  theme(
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black', size=10),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  ) +
  facet_nested( ~ headline + expr_trend, scales = "free")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/allele_specific.expr.volcano.color_by_cell.png", width=458/96, height=273/96, dpi=300)








p = ggplot(consistent_result_plot, aes(x = -avg_log2FC, y = -log10(P.VALUE), fill = expr_effect)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = consistent_result_plot, shape = 21, aes(fill = expr_effect), color = "black") +
  scale_fill_manual(values = c("Increase" = "brown3", "Decrease" = "blue3"), 
                    name = "Allele effect of SNP\nto gene expression") +  # Separate name for effect to expression
  scale_color_manual(values = c("Inconsistent" = "grey"), name = "Consistency") +  # Separate legend for consistency
  labs(x = "log2FC of gene in snRNA-seq\n(hypertension vs. normotension)", y = "-log10(p-value) of SNP in GWAS") +  # No need for fill in labs
  # theme_minimal() +
  theme(
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  ) +
  facet_grid(effect_direction ~ ., scales = "free")

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/allele_specific.png", width=460/96, height=350/96, dpi=300)

dim(unique(consistent_result_plot[, c("trait", "SNPS", "gene_id_mod")]))
dim(unique(consistent_result_plot[, c("SNPS", "gene_id_mod")]))
dim(unique(consistent_result_plot[consistent_result_plot$effect_direction=="Adverse effect to trait", c("SNPS", "gene_id_mod")]))
dim(unique(consistent_result_plot[consistent_result_plot$effect_direction=="Protective effect to trait", c("SNPS", "gene_id_mod")]))


### plot for specific snp-gene pair
## HYP based
# tested: 7042 & 6465(bad eqtl), 9903 (bad expr), 7055 (so-so eqtl), 2850 & 9928 & 12550(bad eqtl & expr), 1822 (bad SD)
# 10618, 6213, 10154 (BUN included), 11893

## LK based
# tested: 1565 (not sig in ss, bad eqtl), 5519 (reversed change between shr and wky baseline)
# 8992, 7656

## LV based
# tested: 1980 & 1980 & 8655 & 1580 & 1430(Arhgap42, bad eqtl), 9459 & 5759 (bad eqtl & expr),  3960 (shr only), 6881 & 6940 & 11927 & 5634(bad expr)
# 4866, 3890 & 6954 (angii only), 11449, 7564, 7021, 2762
row_num = 4866
p1 = generate_forest_plot(consistent_result, row_num, gwas_merge)
p2 = generate_boxplot_with_annotation(consistent_result, row_num, vcf, gene_expr_list)
p3 = generate_dot_plot(consistent_result, row_num, deg_pseudo)
combined_plot <- (p1 + p2 + p3) + plot_layout(guides = 'collect', widths = c(0.5, 0.5, 0.5)) &
  theme(legend.position = 'bottom', 
        legend.box = 'vertical')
print(combined_plot)



# rs3768046-G & Tie1
row_num = 7021
p1 = generate_forest_plot(consistent_result, row_num, gwas_merge) + theme(
  axis.text.x = element_text(angle = 0, hjust = 0.5, colour = 'black', size=12))
p2 = generate_boxplot_with_annotation(consistent_result, row_num, vcf, gene_expr_list)
p3 = generate_dot_plot(consistent_result, row_num, deg_pseudo)
design = "AC
          BC"
combined_plot <- (p1 + p2 + p3) + plot_layout(design = design, guides = 'collect', widths = c(0.3, 0.2)) &
  theme(legend.position = 'right', 
        legend.box = 'vertical')
print(combined_plot)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/eqtl.rs3768046_Tie1.png", width=695/96, height=468/96, dpi=300)

consistent_result[row_num, ]


################################################################################
### forest plot for beta value
generate_forest_plot <- function(consistent_result, row_number, data) {
  # Extract relevant row from consistent_result
  selected_row <- consistent_result[row_number, ]
  
  # Extract the SNP ID and other relevant information
  snp_id <- selected_row$SNPS
  
  # Filter data for the given SNP ID
  snp_data <- data %>%
    filter(SNPS == snp_id) %>%
    group_by(trait) %>%
    slice_min(P.VALUE, n = 1) %>%  # Select record with the smallest p-value for each trait
    ungroup()
  
  # Process the confidence intervals and OR/Beta values
  snp_data <- snp_data %>%
    # Extract lower and upper bounds from the CI text
    mutate(
      ci_bounds = str_extract_all(X95..CI..TEXT., "\\[.*?\\]"),  # Extract text between brackets
      ci_bounds = gsub("\\[|\\]", "", ci_bounds),  # Remove brackets
      lower_ci = as.numeric(sapply(ci_bounds, function(x) strsplit(x, "-")[[1]][1])),
      upper_ci = as.numeric(sapply(ci_bounds, function(x) strsplit(x, "-")[[1]][2]))
    ) %>%
    # Calculate beta and handle effect direction
    mutate(
      beta = as.numeric(OR.or.BETA),
      effect_direction = ifelse(grepl("decrease", X95..CI..TEXT., ignore.case = TRUE), -1, 1)
      # beta = beta * effect_direction
      # lower_ci = ifelse(effect_direction == -1, -upper_ci, lower_ci),
      # upper_ci = ifelse(effect_direction == -1, -lower_ci, upper_ci)
    )
  
  if (any(snp_data$effect_direction == -1)) {
    legend_labels <- c("Protective", "Adverse")
  } else {
    legend_labels <- c("Adverse", "Protective")
  }
  
  # Generate the forest plot
  ggplot(snp_data, aes(x = beta, y = reorder(trait, beta))) +
    geom_point(aes(color = as.factor(effect_direction)), size = 3) +
    geom_errorbarh(aes(xmin = lower_ci, xmax = upper_ci, color = as.factor(effect_direction)), height = 0) +
    geom_text(aes(label = format(signif(P.VALUE, 3), scientific = TRUE),
                  x = beta, y = reorder(trait, beta)), 
              vjust = -0.5, size = 3.5) +  # Adjust vjust to position the text above the point
    scale_color_manual(values = c("1" = "red", "-1" = "blue"), labels = legend_labels) +
    labs(
      title = selected_row$STRONGEST.SNP.RISK.ALLELE,
      x = "Effect Size", # "Effect Size\n(OR or Beta)"
      y = "",
      color = "Effect"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size=12),
      axis.text.y = element_text(colour = 'black'),
      title = element_text(size=12)
    ) +
    coord_flip()
}





################################################################################
### pre-processing
# eqtl_merge <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/eqtl_merged.txt", header = TRUE, sep = '\t')
# eqtl_id = unique(eqtl_merge$variant_id)
# write.table(eqtl_id, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/eqtl.variant_id.txt", col.names = F, row.names = F, quote = F, sep = "\t")
# bcftools view -i 'ID=@/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/eqtl.variant_id.txt' /xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_866Indiv.vcf.gz -Oz -o /xdisk/mliang1/qqiu/project/multiomics-hypertension/data/GTEx.eqtl.filtered.vcf.gz

### eqtl plot
generate_boxplot_with_annotation <- function(consistent_result, row_number, vcf, gene_expr_list) {
  # Extract relevant row from consistent_result
  selected_row <- consistent_result[row_number, ]
  
  # Extract variant_id, gene_id, tissue, slope, and p-value from the selected row
  variant_id <- selected_row$variant_id
  gene_id <- selected_row$gene_id_mod
  gene_name <- selected_row$gene_name.x
  tissue <- selected_row$tissue.x
  slope <- selected_row$slope
  p_value <- selected_row$pval_nominal
  
  # Get the gene expression data for the specific tissue
  gene_expr <- gene_expr_list[[tissue]]
  
  # Clean gene expression sample IDs
  colnames(gene_expr) <- gsub("\\.", "-", colnames(gene_expr))
  colnames(gene_expr)[4:ncol(gene_expr)] <- str_extract(colnames(gene_expr)[4:ncol(gene_expr)], "^GTEX-[^-]+")
  gene_expr$Name <- gsub("\\..*", "", gene_expr$Name)
  
  # Extract ref and alt alleles from the variant_id
  ref_allele <- selected_row$ref_allele
  alt_allele <- selected_row$alt_allele
  
  # Extract the genotype data for the specified variant
  genotypes <- extract.gt(vcf, element = "GT")[variant_id, ]
  
  # Process genotype data to convert to ref/ref, ref/alt, alt/alt
  genotypes <- as.data.frame(genotypes) %>%
    rownames_to_column("SampleID") %>%
    mutate(Genotype = case_when(
      genotypes == "0/0" ~ paste(ref_allele, ref_allele, sep = "/"),
      genotypes == "0/1" | genotypes == "1/0" ~ paste(ref_allele, alt_allele, sep = "/"),
      genotypes == "1/1" ~ paste(alt_allele, alt_allele, sep = "/"),
      TRUE ~ NA_character_  # Handle missing or unusual genotypes
    )) %>%
    mutate(Genotype = factor(Genotype, levels = c(
      paste(ref_allele, ref_allele, sep = "/"),   # "0/0"
      paste(ref_allele, alt_allele, sep = "/"),   # "0/1" or "1/0"
      paste(alt_allele, alt_allele, sep = "/")    # "1/1"
    )))
  
  # Extract expression data for the specified gene
  expr <- gene_expr %>%
    filter(Name == gene_id) %>%
    gather(SampleID, Expression, -id, -Name, -Description)
  
  # Merge genotype and expression data
  merged_data <- merge(genotypes, expr, by = "SampleID")
  
  # Generate the boxplot with slope and p-value annotation
  p <- ggplot(merged_data, aes(x = Genotype, y = log2(Expression+1))) +
    geom_boxplot() +
    labs(
      title = tissue,
      subtitle = paste("Slope:", round(slope, 2), "\nP:", 
                    format(signif(p_value, 3), scientific = TRUE)),
      y = paste0(gene_name, "\nlog2(TPM + 1)"),
      x = ""
    ) +
    # theme_minimal() +
    theme(plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
          axis.text.x = element_text(colour = 'black', size=12),
          title = element_text(size=12))  # Style the subtitle
  
  print(p)
}



################################################################################
### dotplot for snRNA-seq
generate_dot_plot <- function(consistent_result, row_number, deg_pseudo) {
  # Extract relevant row from consistent_result
  selected_row <- consistent_result[row_number, ]
  
  # Extract necessary variables from the selected row
  strain <- selected_row$strain
  tissue <- selected_row$tissue.y
  cell_type <- selected_row$cell_type
  gene_name <- selected_row$combined_gene_name  # or `gene_name.x` based on your preference
  
  # Filter the expression data based on the input
  deg_use <- expr_all[
    # expr_all$strain %in% strain & 
      expr_all$tissue %in% tissue &
      expr_all$cell_type %in% cell_type & 
      expr_all$gene_name %in% gene_name,
  ]

  # Generate the dot plot
  p1 <- ggplot(deg_use, aes(x = treatment, y = gene_name)) +
    geom_point(aes(fill = log2(avg_expr + 1), size = pct * 100), shape = 21) +
    scale_y_discrete(limits = rev(levels(deg_use$gene_name))) +
    scale_fill_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major.y = element_blank(),  # No horizontal grid lines
      legend.justification = c(1, 0.5),
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size=12),
      strip.text = element_text(colour = 'black'),
      title = element_text(size=12)
    ) +
    labs(
      x = "",
      y = "",
      fill = "log2(expression+1)",
      size = "Percentage",
      title = paste0(tissue, "-", cell_type)
    ) +
    facet_nested(project + strain ~ ., scales = "free", space = "free") +
    coord_flip()
  
  print(p1)
}

