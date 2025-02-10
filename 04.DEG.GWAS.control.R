# library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
# library(openxlsx)
# library(httr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(patchwork)
library(stringr)
library(igraph)
library(vegan)
library(ggraph)
library(ggsignif)
Sys.setenv(RETICULATE_PYTHON = "/home/u1/qqiu/.conda/envs/stereopy/bin/python")
library(reticulate)
# use_python("/home/u1/qqiu//Python/pyenv_v3.8/bin/python")
library(leiden)
library(lsa)

# options(timeout = 300)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/control_traits")


################################################################################
### load reference files
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

### h2m
r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) = c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
r2h[r2h$gene_name_rat=="",]$gene_name_rat = r2h[r2h$gene_name_rat=="",]$gene_id_rat
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) = c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")
m2h[m2h$gene_name_rat=="",]$gene_name_rat = m2h[m2h$gene_name_rat=="",]$gene_id_rat


gene_gr <- GRanges(
  seqnames = Rle(ensembl$Chromosome.scaffold.name),
  ranges = IRanges(start = ensembl$Gene.start..bp., end = ensembl$Gene.end..bp.),
  strand = Rle(ensembl$Strand),
  gene_id = ensembl$Gene.stable.ID,
  gene_name = ensembl$Gene.name,
  gene_type = ensembl$Gene.type
)



################################################################################
### extract snp information
gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/control_traits/"

gwas_list = c("gwas_catalog_birth_weight.tsv", "gwas_catalog_bone_density.tsv", 
              "gwas_catalog_breast_cancer.tsv", "gwas_catalog_colorectal_cancer.tsv", 
              "gwas_catalog_alzheimer_disease.tsv", "gwas_catalog_parkinson_disease.tsv",
              "gwas_catalog_multiple_sclerosis.tsv", "gwas_catalog_rheumatoid_arthritis.tsv",
              "gwas_catalog_asthma.tsv", "gwas_catalog_psoriasis.tsv", "gwas_catalog_copd.tsv")
outfile = "gwas_catalog_control_traits.snp.txt"

bp_SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

SNPS = c()
for(gwas_file in gwas_list){
  
  output = gsub("tsv", "snp.txt", gwas_file)
  
  dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
  dat = dat[dat$P.VALUE<5e-8, ]
  
  snp_pos = data.frame(SNP=dat$SNPS, POS = paste0("chr", dat$CHR_ID, "_", dat$CHR_POS))
  snp_pos[!(grepl("^rs", snp_pos$SNP)), ]$POS = gsub(":", "_", snp_pos[!(grepl("^rs", snp_pos$SNP)), ]$SNP)
  
  SNPS = rbind(SNPS, unique(snp_pos))
  
}

any(grepl(",", SNPS$SNP))
any(grepl(";", SNPS$SNP))
any(grepl("x", SNPS$SNP))

SNPS = SNPS[! grepl(",", SNPS$SNP), ]
SNPS = SNPS[! grepl(";", SNPS$SNP), ]
SNPS = SNPS[! grepl("x", SNPS$SNP), ]

SNPS = SNPS[! (SNPS$SNP %in% bp_SNPS$SNP), ]

write.table(SNPS, outfile, col.names = T, row.names = F, quote = F, sep = '\t')

#### in shell
# awk 'NR==FNR { pos[$2] = $1; next } FNR==1 { print $0, "POS_dbSNP"; next } { print $0, pos[$1] }' /xdisk/mliang1/qqiu/reference/dbSNP/dbSNP.pos.chr.reformat.txt gwas_catalog_bp_relevant.snp.txt > gwas_catalog_bp_relevant.snp.pos.txt
# mv gwas_catalog_bp_relevant.snp.pos.txt gwas_catalog_bp_relevant.snp.txt

################################################################################
# library(VariantAnnotation)
# 
# vcf <- readVcf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.vcf", "hg38")  # Specify the genome build (e.g., "hg19", "hg38")
# 
# snp_info <- data.frame(
#   SNP = rownames(vcf),
#   POS = paste(seqnames(rowRanges(vcf)), start(rowRanges(vcf)), sep='_')
# )
# 
# snp_info = unique(snp_info)
# rownames(snp_info) = snp_info$SNP
# SNPS$POS_new = snp_info[SNPS$SNP,]
# SNPS[rowSums(is.na(SNPS))>0,]
# all(SNPS[rowSums(is.na(SNPS))==0,]$POS==SNPS[rowSums(is.na(SNPS))==0,]$POS_new.POS)
# # [1] TRUE


################################################################################
### gwas data pre-processing
gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/control_traits/"

gwas_list = c("gwas_catalog_birth_weight.tsv", "gwas_catalog_bone_density.tsv", 
              "gwas_catalog_breast_cancer.tsv", "gwas_catalog_colorectal_cancer.tsv", 
              "gwas_catalog_alzheimer_disease.tsv", "gwas_catalog_parkinson_disease.tsv",
              "gwas_catalog_multiple_sclerosis.tsv", "gwas_catalog_rheumatoid_arthritis.tsv",
              "gwas_catalog_asthma.tsv", "gwas_catalog_psoriasis.tsv", "gwas_catalog_copd.tsv")

gwas_merge = c()
for(gwas_file in gwas_list){

  output = gsub("tsv", "snp.txt", gwas_file)
  trait = gsub("(gwas_catalog_|.tsv)", "", gwas_file)

  dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
  dat = dat[dat$P.VALUE<5e-8, ]
  dat$trait = trait

  gwas_merge = rbind(gwas_merge, unique(dat))

}

gwas_merge = separate_rows(gwas_merge, SNPS, sep = ";")
gwas_merge = gwas_merge[! grepl(" x ", gwas_merge$SNPS), ]
gwas_merge$SNPS = gsub(" ", "", gwas_merge$SNPS)
gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = " - ")
gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = ", ")
gwas_merge = separate_rows(gwas_merge, REPORTED.GENE.S., sep = ", ")
write.table(gwas_merge, "gwas_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')



################################################################################
### snp-gene pairs in 1 mb windows
# SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/control_traits/gwas_catalog_control_traits.snp.txt", header=T, sep='\t')
# 
# range <- 100000
# 
# SNPS <- SNPS %>% filter(grepl("rs", SNP)) %>% 
#   separate(POS, into = c("chr", "pos"), sep = "_", convert = TRUE, remove = FALSE) %>%
#   mutate(pos = as.numeric(pos)) %>%  filter(!is.na(pos) | pos=="")
# snps_gr <- GRanges(seqnames = Rle(SNPS$chr),
#                    ranges = IRanges(start = SNPS$pos - range, 
#                                     end = SNPS$pos + range),
#                    SNP = SNPS$SNP)
# 
# ensembl <- ensembl %>%
#   mutate(Chromosome.scaffold.name = paste0("chr", Chromosome.scaffold.name),
#          TSS = ifelse(Strand == 1, Gene.start..bp., Gene.end..bp.))
# 
# genes_gr <- GRanges(seqnames = Rle(ensembl$Chromosome.scaffold.name),
#                     ranges = IRanges(start = ensembl$TSS,
#                                      end = ensembl$TSS),
#                     gene_id = ensembl$Gene.stable.ID,
#                     gene_name = ensembl$Gene.name)
# 
# overlaps <- findOverlaps(snps_gr, genes_gr)
# snp_gene_pairs <- data.frame(
#   SNP = mcols(snps_gr)$SNP[queryHits(overlaps)],
#   Gene_ID = mcols(genes_gr)$gene_id[subjectHits(overlaps)],
#   Gene_Name = mcols(genes_gr)$gene_name[subjectHits(overlaps)]
# )
# 
# write.table(snp_gene_pairs, "proximal_merged.100k.txt", col.names = T, row.names = F, quote = F, sep = '\t')





################################################################################
### nearest protein coding genes
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/control_traits/gwas_catalog_control_traits.snp.txt", header=T, sep='\t')

SNPS <- SNPS %>% filter(grepl("rs", SNP)) %>% 
  separate(POS, into = c("chr", "pos"), sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(pos = as.numeric(pos)) %>%  filter(!is.na(pos) | pos=="")
snps_gr <- GRanges(seqnames = Rle(SNPS$chr),
                   ranges = IRanges(start = SNPS$pos, 
                                    end = SNPS$pos),
                   SNP = SNPS$SNP)

ensembl <- ensembl %>%
  filter(Gene.type == "protein_coding") %>%
  mutate(Chromosome.scaffold.name = paste0("chr", Chromosome.scaffold.name))

genes_gr <- GRanges(seqnames = Rle(ensembl$Chromosome.scaffold.name),
                    ranges = IRanges(start = ensembl$Gene.start..bp.,
                                     end = ensembl$Gene.end..bp.),
                    gene_id = ensembl$Gene.stable.ID,
                    gene_name = ensembl$Gene.name)

nearest_genes <- distanceToNearest(snps_gr, genes_gr)
snp_gene_pairs <- data.frame(
  SNP = mcols(snps_gr)$SNP[queryHits(nearest_genes)],
  Gene_ID = mcols(genes_gr)$gene_id[subjectHits(nearest_genes)],
  Gene_Name = mcols(genes_gr)$gene_name[subjectHits(nearest_genes)]
)

write.table(snp_gene_pairs, "nearest_protein_coding_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')


snp_gene_pairs = read.table("nearest_protein_coding_merged.txt", header=T, sep='\t')
gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote="")
gwas_merge = separate_rows(gwas_merge, SNP_GENE_IDS, sep = ", ")

snp_gene = data.frame(SNP=c(gwas_merge$SNPS, gwas_merge$SNPS),
                      Gene_Name = c(gwas_merge$MAPPED_GENE, gwas_merge$REPORTED.GENE.S.))
snp_gene = unique(snp_gene[snp_gene$Gene_Name!="", ])
snp_id = data.frame(SNP=c(gwas_merge$SNPS, gwas_merge$SNPS, gwas_merge$SNPS),
                    Gene_ID = c(gwas_merge$UPSTREAM_GENE_ID, gwas_merge$DOWNSTREAM_GENE_ID, gwas_merge$SNP_GENE_IDS))
snp_id = unique(snp_id[snp_id$Gene_ID!="", ])
snp_gene = base::merge(snp_gene, ensembl[, c("Gene.name", "Gene.stable.ID")], by.x="Gene_Name", by.y="Gene.name", all.x=T)
colnames(snp_gene)[3] = "Gene_ID"
snp_id = base::merge(snp_id, ensembl[, c("Gene.name", "Gene.stable.ID")], by.x="Gene_ID", by.y="Gene.stable.ID", all.x=T)
colnames(snp_id)[3] = "Gene_Name"

snp_gene_merge = unique(rbind(snp_gene[c("SNP", "Gene_ID", "Gene_Name")],
                              snp_id[c("SNP", "Gene_ID", "Gene_Name")],
                              snp_gene_pairs[c("SNP", "Gene_ID", "Gene_Name")]))
# snp_gene_merge = unique(rbind(snp_gene[c("SNP", "Gene_ID", "Gene_Name")],
#                               snp_id[c("SNP", "Gene_ID", "Gene_Name")]))
snp_gene_merge = snp_gene_merge[!(is.na(snp_gene_merge$Gene_ID)),]

write.table(snp_gene_merge, "gwas_snp_gene_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')


################################################################################
### eQTL data processing
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/control_traits/gwas_catalog_control_traits.snp.txt", header=T, sep='\t')

eqtl_folder = "/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/"

eqtl_list = c("Artery_Aorta.v8.signif_variant_gene_pairs.txt", 
              "Artery_Coronary.v8.signif_variant_gene_pairs.txt",
              "Artery_Tibial.v8.signif_variant_gene_pairs.txt", 
              "Brain_Hypothalamus.v8.signif_variant_gene_pairs.txt", 
              "Heart_Left_Ventricle.v8.signif_variant_gene_pairs.txt", 
              "Kidney_Cortex.v8.signif_variant_gene_pairs.txt")

eqtl_merged = c()
for(eqtl_file in eqtl_list){
  
  tissue = gsub(".v8.*", "", eqtl_file)
  eqtl = read.table(paste0(eqtl_folder, eqtl_file), header=T, sep='\t', quote = "", fill = T)
  
  eqtl$POS = gsub("_[ATCG].*", "", eqtl$variant_id)
  eqtl$gene_id_mod = gsub("\\..*", "", eqtl$gene_id)
  eqtl$gene_name = ensembl[eqtl$gene_id_mod, ]$Gene.name
  eqtl$tissue = tissue
  
  eqtl = eqtl[eqtl$POS %in% SNPS$POS, ]
  eqtl = base::merge(eqtl, SNPS, by="POS", all.x=T)
  eqtl = base::merge(eqtl, r2h[r2h$gene_id %in% eqtl$gene_id_mod, ], by.x="gene_id_mod", by.y="gene_id", all.x=T)
  eqtl = base::merge(eqtl, m2h[m2h$gene_id %in% eqtl$gene_id_mod, ], by.x="gene_id_mod", by.y="gene_id", all.x=T)
  
  eqtl_merged = rbind(eqtl_merged, eqtl)
  
}

eqtl_merged = eqtl_merged[, !(colnames(eqtl_merged) %in% c("gene_name.y", "gene_name"))]

write.table(eqtl_merged, "eqtl_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')



################################################################################
### epimap data processing
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/control_traits/gwas_catalog_control_traits.snp.txt", header=T, sep='\t')

SNPS <- transform(SNPS, chr = sub("_.*", "", POS), pos = as.numeric(sub(".*_", "", POS)))
SNPS_filter = SNPS[rowSums(is.na(SNPS))==0,]
snps_gr <- GRanges(seqnames = SNPS_filter$chr,
                   ranges = IRanges(start = SNPS_filter$pos, end = SNPS_filter$pos),
                   SNP = SNPS_filter$SNP)

epimap_folder = "/xdisk/mliang1/qqiu/data/EpiMap/"

e2g_list = c("links_by_group.brain.tsv", 
             "links_by_group.endothelial.tsv", 
             "links_by_group.heart.tsv", 
             "links_by_group.kidney.tsv")

chain <- import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg19ToHg38.over.chain")

e2g_merged = c()
for(e2g_file in e2g_list){
  
  tissue = gsub("(links_by_group.|.tsv)", "", e2g_file)
  e2g = read.table(paste0(epimap_folder, e2g_file), header=T, sep='\t', quote = "", fill = T)
  e2g$gene_name = ensembl[e2g$gene, ]$Gene.name
  
  e2g_gr <- GRanges(seqnames = e2g$chr,
                    ranges = IRanges(start = e2g$start, end = e2g$end),
                    gene = e2g$gene,
                    name = e2g$name,
                    score = e2g$score,
                    group = e2g$group,
                    gene_name = e2g$gene_name)
  
  e2g_gr_hg38 <- liftOver(e2g_gr, chain)
  e2g_gr_hg38 <- unlist(e2g_gr_hg38)
  
  overlaps <- findOverlaps(snps_gr, e2g_gr_hg38)
  
  snp_indices <- queryHits(overlaps)
  e2g_indices <- subjectHits(overlaps)
  
  # Combine information from both dataframes
  combined_df <- data.frame(
    SNP = SNPS_filter$SNP[snp_indices],
    SNP_POS = SNPS_filter$POS[snp_indices],
    chr = seqnames(e2g_gr_hg38)[e2g_indices],
    start = start(e2g_gr_hg38)[e2g_indices],
    end = end(e2g_gr_hg38)[e2g_indices],
    gene = mcols(e2g_gr_hg38)$gene[e2g_indices],
    name = mcols(e2g_gr_hg38)$name[e2g_indices],
    score = mcols(e2g_gr_hg38)$score[e2g_indices],
    group = mcols(e2g_gr_hg38)$group[e2g_indices],
    gene_name = mcols(e2g_gr_hg38)$gene_name[e2g_indices]
  )
  
  combined_df = base::merge(combined_df, r2h[r2h$gene_id %in% combined_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
  combined_df = base::merge(combined_df, m2h[m2h$gene_id %in% combined_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
  
  e2g_merged = rbind(e2g_merged, combined_df)
  
}

e2g_merged = e2g_merged[, !(colnames(e2g_merged) %in% c("gene_name.y", "gene_name"))]

write.table(e2g_merged, "e2g_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')







################################################################################
### #deg in each trait
one_sided_z_test <- function(proportion_hit_genes, porpotion_background, n_snp_gene_list) {
  p_hat <- proportion_hit_genes
  p_0 <- porpotion_background
  n <- n_snp_gene_list
  
  z_score <- (p_hat - p_0) / sqrt(p_0 * (1 - p_0) / n)
  p_value <- pnorm(z_score, lower.tail = FALSE)  # one-sided greater test
  
  return(p_value)
}

gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')
# background_500k <- read.table("proximal_merged.500k.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]
# background_500k = background_500k[grepl("rs", background_500k$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)
snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

bp_snp_gene_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
snp_gene_df = snp_gene_df[!(snp_gene_df$gene %in% bp_snp_gene_df$gene), ]

snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")

snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
# snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_name_rat <- snp_gene_df[snp_gene_df$gene_name_rat == "" & !is.na(snp_gene_df$gene_name_rat), ]$gene_id_rat
# snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_name_mouse <- snp_gene_df[snp_gene_df$gene_name_mouse == "" & !is.na(snp_gene_df$gene_name_mouse), ]$gene_id_mouse

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]

deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]

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
  snp_gene_list = unique(snp_gene_df_use$gene)
  
  for (threshold_name in names(thresholds)) {
    deg_filtered <- thresholds[[threshold_name]]
    mapped_gene_list <- unique(c(snp_gene_df_use[!is.na(snp_gene_df_use$gene_name_mouse), ]$gene_name_mouse,
                                 snp_gene_df_use[!is.na(snp_gene_df_use$gene_name_rat), ]$gene_name_rat))
    mapped_gene_list = setdiff(mapped_gene_list, "")
    deg_list <- unique(deg_filtered$gene_name)
    
    hit_genes <- intersect(mapped_gene_list, deg_list)
    num_hit_genes <- length(hit_genes)
    proportion_hit_genes <- num_hit_genes / length(snp_gene_list)
    num_DEG_snps = length(unique(snp_gene_df_use[(! is.na(snp_gene_df_use$gene_name_mouse) & snp_gene_df_use$gene_name_mouse %in% deg_list) | 
                                                   (!is.na(snp_gene_df_use$gene_name_rat) & snp_gene_df_use$gene_name_rat %in% deg_list), ]$SNP))
    
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
                       "(", round(results$proportion_DEG * 100, 2), "%)")
results$threshold <- factor(results$threshold, levels = rev(names(thresholds)))
results = results[results$trait!="parkinson_disease", ]
write.table(results, "trait.deg_prop.z_test.out", sep = ",", col.names = T, row.names = F, quote = F)



################################################################################
results = read.table("trait.deg_prop.z_test.out", header = TRUE, sep = ",", quote = "")
bp_results = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/trait.deg_prop.z_test.out", header = TRUE, sep = ',', quote = "")

results$category = "Random traits"
bp_results$category = "BP-relevant traits"
select_col = c("trait", "threshold", "proportion_DEG", "category")
results_use = rbind(results[select_col], bp_results[select_col])

write.table(results_use, "trait.bp_random.deg_prop.out", sep = ",", col.names = T, row.names = F, quote = F)

threshold_col = RColorBrewer::brewer.pal(5, "Blues")[2:5]
names(threshold_col) = names(thresholds)

calculate_p_values <- function(data) {
  p_values <- data %>%
    group_by(threshold) %>%
    summarise(p_value = wilcox.test(proportion_DEG[category == "BP-relevant traits"],
                                    proportion_DEG[category == "Random traits"],
                                    alternative = "greater")$p.value) %>%
    ungroup()
  return(p_values)
}
p_values <- calculate_p_values(results_use)

ggplot(results_use, aes(x = category, y = proportion_DEG, fill = threshold, label = trait)) +
  geom_point(shape=21) +
  scale_fill_manual(values = threshold_col) +
  labs(x = "", y = "DEG proportion", color = "DEG Threshold") +
  scale_y_continuous(limits = c(0, 0.85)) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')) +
  facet_wrap(~threshold, scale = "free_x", nrow = 1) +
  geom_text(data = p_values, aes(x = 1.5, y = 0.75, label = paste("p =", format(p_value, digits = 3))), inherit.aes = FALSE)






