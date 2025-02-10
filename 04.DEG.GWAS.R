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
Sys.setenv(RETICULATE_PYTHON = "/home/u1/qqiu/.conda/envs/stereopy/bin/python")
library(reticulate)
# use_python("/home/u1/qqiu//Python/pyenv_v3.8/bin/python")
library(leiden)
library(lsa)

# options(timeout = 300)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")


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
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/"

gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv", 
              "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
              "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv", 
              "gwas_catalog_BUN.tsv")
outfile = "gwas_catalog_bp_relevant.snp.txt"

gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv")
outfile = "gwas_catalog_bp.snp.txt"

gwas_list = c("gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
              "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv", 
              "gwas_catalog_BUN.tsv")
outfile = "gwas_catalog_non-bp.snp.txt"

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

if(any(grepl(";", SNPS$SNP))){
  SNPS = SNPS[! grepl(";", SNPS$SNP), ]
  SNPS = rbind(SNPS, 
               data.frame(SNP = c("rs10755578", "rs2048327", "rs7767084", "rs3127599"),
                          POS = c("chr6_160548706", "chr6_160442500", "chr6_160541471", "chr6_160486102")))
}

write.table(SNPS[, 1], outfile, col.names = T, row.names = F, quote = F, sep = '\t')

#### in shell
# awk 'NR==FNR { pos[$2] = $1; next } FNR==1 { print $0, "POS_dbSNP"; next } { print $0, pos[$1] }' /xdisk/mliang1/qqiu/reference/dbSNP/dbSNP.pos.chr.reformat.txt gwas_catalog_bp_relevant.snp.txt > gwas_catalog_bp_relevant.snp.pos.txt
# mv gwas_catalog_bp_relevant.snp.pos.txt gwas_catalog_bp_relevant.snp.txt

################################################################################
library(VariantAnnotation)

vcf <- readVcf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.vcf", "hg38")  # Specify the genome build (e.g., "hg19", "hg38")

snp_info <- data.frame(
  SNP = rownames(vcf),
  POS = paste(seqnames(rowRanges(vcf)), start(rowRanges(vcf)), sep='_')
)

snp_info = unique(snp_info)
rownames(snp_info) = snp_info$SNP
SNPS$POS_new = snp_info[SNPS$SNP,]
SNPS[rowSums(is.na(SNPS))>0,]
all(SNPS[rowSums(is.na(SNPS))==0,]$POS==SNPS[rowSums(is.na(SNPS))==0,]$POS_new.POS)
# [1] TRUE


################################################################################
### gwas data pre-processing
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
# gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/"
# 
# gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv",
#               "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv",
#               "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
#               "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv",
#               "gwas_catalog_BUN.tsv")
# 
# gwas_merge = c()
# for(gwas_file in gwas_list){
# 
#   output = gsub("tsv", "snp.txt", gwas_file)
#   trait = gsub("(gwas_catalog_|.tsv)", "", gwas_file)
# 
#   dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
#   dat = dat[dat$P.VALUE<5e-8, ]
#   dat$trait = trait
# 
#   gwas_merge = rbind(gwas_merge, unique(dat))
# 
# }
# 
# gwas_merge = separate_rows(gwas_merge, SNPS, sep = ";")
# gwas_merge = gwas_merge[! grepl(" x ", gwas_merge$SNPS), ]
# gwas_merge$SNPS = gsub(" ", "", gwas_merge$SNPS)
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = " - ")
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = ", ")
# gwas_merge = separate_rows(gwas_merge, REPORTED.GENE.S., sep = ", ")
# write.table(gwas_merge, "gwas_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')
# 
# 
# gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
#               "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv")
# 
# gwas_merge = c()
# for(gwas_file in gwas_list){
#   
#   output = gsub("tsv", "snp.txt", gwas_file)
#   
#   dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
#   dat = dat[dat$P.VALUE<5e-8, ]
#   
#   gwas_merge = rbind(gwas_merge, unique(dat))
#   
# }
# 
# gwas_merge = separate_rows(gwas_merge, SNPS, sep = ";")
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = " - ")
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = ", ")
# gwas_merge = separate_rows(gwas_merge, REPORTED.GENE.S., sep = ", ")
# write.table(gwas_merge, "bp.gwas_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')
# 
# 
# gwas_list = c("gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
#               "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv", 
#               "gwas_catalog_BUN.tsv")
# 
# gwas_merge = c()
# for(gwas_file in gwas_list){
#   
#   output = gsub("tsv", "snp.txt", gwas_file)
#   
#   dat = read.table(paste0(gwas_folder, gwas_file), header=T, sep='\t', quote = "", fill = T)
#   dat = dat[dat$P.VALUE<5e-8, ]
#   
#   gwas_merge = rbind(gwas_merge, unique(dat))
#   
# }
# 
# gwas_merge = separate_rows(gwas_merge, SNPS, sep = ";")
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = " - ")
# gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = ", ")
# gwas_merge = separate_rows(gwas_merge, REPORTED.GENE.S., sep = ", ")
# write.table(gwas_merge, "non-bp.gwas_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')



################################################################################
### snp-gene pairs in 1 mb windows
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

range <- 100000

SNPS <- SNPS %>% filter(grepl("rs", SNP)) %>% 
  separate(POS, into = c("chr", "pos"), sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(pos = as.numeric(pos)) %>%  filter(!is.na(pos) | pos=="")
snps_gr <- GRanges(seqnames = Rle(SNPS$chr),
                   ranges = IRanges(start = SNPS$pos - range, 
                                    end = SNPS$pos + range),
                   SNP = SNPS$SNP)

ensembl <- ensembl %>%
  mutate(Chromosome.scaffold.name = paste0("chr", Chromosome.scaffold.name),
         TSS = ifelse(Strand == 1, Gene.start..bp., Gene.end..bp.))

genes_gr <- GRanges(seqnames = Rle(ensembl$Chromosome.scaffold.name),
                    ranges = IRanges(start = ensembl$TSS,
                                     end = ensembl$TSS),
                    gene_id = ensembl$Gene.stable.ID,
                    gene_name = ensembl$Gene.name)

overlaps <- findOverlaps(snps_gr, genes_gr)
snp_gene_pairs <- data.frame(
  SNP = mcols(snps_gr)$SNP[queryHits(overlaps)],
  Gene_ID = mcols(genes_gr)$gene_id[subjectHits(overlaps)],
  Gene_Name = mcols(genes_gr)$gene_name[subjectHits(overlaps)]
)

write.table(snp_gene_pairs, "proximal_merged.100k.txt", col.names = T, row.names = F, quote = F, sep = '\t')





################################################################################
### snp-gene pairs in 1 mb windows
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

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
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

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
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

SNPS <- transform(SNPS, chr = sub("_.*", "", POS), pos = as.numeric(sub(".*_", "", POS)))
SNPS_filter = SNPS[rowSums(is.na(SNPS))==0,]
snps_gr <- GRanges(seqnames = SNPS_filter$chr,
                   ranges = IRanges(start = SNPS_filter$pos, end = SNPS_filter$pos),
                   SNP = SNPS_filter$SNP)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
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
### venn plot and enrichment
SNPS = read.table("gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')
bp_SNPS = read.table("gwas_catalog_bp.snp.txt", header=T, sep='\t')
nb_SNPS = read.table("gwas_catalog_non-bp.snp.txt", header=T, sep='\t')

# proximal_merge = read.table("proximal_merged.100k.txt", header=T, sep='\t')
# gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote="")
proximal_merge = read.table("gwas_snp_gene_merged.txt", header=T, sep='\t')
eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

proximal_merge_bp = proximal_merge[proximal_merge$SNP %in% bp_SNPS$x, ]
# gwas_merge_bp = gwas_merge[gwas_merge$SNPS %in% bp_SNPS$x, ]
eqtl_merge_bp = eqtl_merge[eqtl_merge$SNP %in% bp_SNPS$x, ]
e2g_merge_bp = e2g_merge[e2g_merge$SNP %in% bp_SNPS$x, ]

proximal_merge_nb = proximal_merge[proximal_merge$SNP %in% nb_SNPS$x, ]
# gwas_merge_nb = gwas_merge[gwas_merge$SNPS %in% nb_SNPS$x, ]
eqtl_merge_nb = eqtl_merge[eqtl_merge$SNP %in% nb_SNPS$x, ]
e2g_merge_nb = e2g_merge[e2g_merge$SNP %in% nb_SNPS$x, ]

snp_gene_prox = unique(paste(proximal_merge$SNP, proximal_merge$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge$SNP, eqtl_merge$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge$SNP, e2g_merge$gene, sep="-"))
# snp_gene_gwas = unique(c(paste(gwas_merge$SNPS, gwas_merge$MAPPED_GENE, sep="-"), paste(gwas_merge$SNPS, gwas_merge$REPORTED.GENE.S., sep="-")))
# snp_gene_eqtl = unique(paste(eqtl_merge$SNP, eqtl_merge$gene_name.x, sep="-"))
# snp_gene_e2g = unique(paste(e2g_merge$SNP, e2g_merge$gene_name.x, sep="-"))

snp_gene_prox = unique(paste(proximal_merge_bp$SNP, proximal_merge_bp$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge_bp$SNP, eqtl_merge_bp$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge_bp$SNP, e2g_merge_bp$gene, sep="-"))
# snp_gene_gwas = unique(c(paste(gwas_merge_bp$SNPS, gwas_merge_bp$MAPPED_GENE, sep="-"), paste(gwas_merge_bp$SNPS, gwas_merge_bp$REPORTED.GENE.S., sep="-")))
# snp_gene_eqtl = unique(paste(eqtl_merge_bp$SNP, eqtl_merge_bp$gene_name.x, sep="-"))
# snp_gene_e2g = unique(paste(e2g_merge_bp$SNP, e2g_merge_bp$gene_name.x, sep="-"))

# snp_gene_gwas = unique(c(paste(eqtl_merge_nb$SNPS, eqtl_merge_nb$MAPPED_GENE, sep="-"), paste(eqtl_merge_nb$SNPS, eqtl_merge_nb$REPORTED.GENE.S., sep="-")))
# snp_gene_eqtl = unique(paste(eqtl_merge_nb$SNP, eqtl_merge_nb$gene_name.x, sep="-"))
# snp_gene_e2g = unique(paste(e2g_merge_nb$SNP, e2g_merge_nb$gene_name.x, sep="-"))

library(ggVennDiagram)

x = list("Proximal"=snp_gene_prox, 
         "Expressional"=snp_gene_eqtl,
         "Regulatory"=snp_gene_e2g)

ggVennDiagram(x, category.names = "", label="count", label_alpha = 0) +
  scale_fill_gradient("Number of\nSNP-gene\npairs", low="white",high = "red")


snp_gene_prox = unique(proximal_merge$Gene_ID)
snp_gene_eqtl = unique(eqtl_merge$gene_id_mod)
snp_gene_e2g = unique(e2g_merge$gene)

x = list("Proximal"=snp_gene_prox, 
         "Expressional"=snp_gene_eqtl,
         "Regulatory"=snp_gene_e2g)

ggVennDiagram(x, category.names = "", label="both", label_alpha = 0) +
  scale_fill_gradient("Number of genes", low="white",high = "red")


length(unique(c(proximal_merge$SNP, eqtl_merge$SNP, e2g_merge$SNP)))
length(unique(c(proximal_merge$Gene_ID, eqtl_merge$gene_id_mod, e2g_merge$gene)))


snp_gene_prox = unique(paste(proximal_merge_bp$SNP, proximal_merge_bp$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge_bp$SNP, eqtl_merge_bp$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge_bp$SNP, e2g_merge_bp$gene, sep="-"))
# snp_gene_gwas = unique(c(paste(gwas_merge_bp$SNPS, gwas_merge_bp$MAPPED_GENE, sep="-"), paste(gwas_merge_bp$SNPS, gwas_merge_bp$REPORTED.GENE.S., sep="-"))) %>% filter(grepl("-$",.))
# snp_gene_eqtl = unique(paste(eqtl_merge_bp$SNP, eqtl_merge_bp$gene_name.x, sep="-"))
# snp_gene_e2g = unique(paste(e2g_merge_bp$SNP, e2g_merge_bp$gene_name.x, sep="-"))
bp_snp_ev_count = table(c(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_prox = unique(paste(proximal_merge_nb$SNP, proximal_merge_nb$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge_nb$SNP, eqtl_merge_nb$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge_nb$SNP, e2g_merge_nb$gene, sep="-"))
# snp_gene_gwas = unique(c(paste(eqtl_merge_nb$SNPS, eqtl_merge_nb$MAPPED_GENE, sep="-"), paste(eqtl_merge_nb$SNPS, eqtl_merge_nb$REPORTED.GENE.S., sep="-")))
# snp_gene_eqtl = unique(paste(eqtl_merge_nb$SNP, eqtl_merge_nb$gene_name.x, sep="-"))
# snp_gene_e2g = unique(paste(e2g_merge_nb$SNP, e2g_merge_nb$gene_name.x, sep="-"))
nb_snp_list = unique(c(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

table(names(bp_snp_ev_count) %in% nb_snp_list)

# Define the parameters for the Hypergeometric Test
k <- sum(names(bp_snp_ev_count[bp_snp_ev_count>1]) %in% nb_snp_list)
m <- sum(names(bp_snp_ev_count) %in% nb_snp_list)
n <- sum(!(names(bp_snp_ev_count) %in% nb_snp_list))
N <- sum(bp_snp_ev_count>1)

p_value <- phyper(k - 1, m, n, N, lower.tail = FALSE)
print(p_value)
# [1] 0.0002776661

data <- data.frame(
  Evidence = c("Multiple", "Single", "Proximal-only"),
  Count = c(sum(names(bp_snp_ev_count[bp_snp_ev_count > 1]) %in% nb_snp_list), 
            sum(names(bp_snp_ev_count[bp_snp_ev_count == 1]) %in% nb_snp_list))
) %>%
  mutate(Total = c(sum(bp_snp_ev_count > 1), sum(bp_snp_ev_count == 1)),
         Proportion = Count / Total)

ggplot(data, aes(x = Evidence, y = Proportion)) +
  geom_bar(stat = "identity", fill="red") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.01)), vjust = -0.5) +
  annotate("text", x = 1.5, y = max(data$Proportion) + 0.02, 
           label = paste("Hypergeometric test\np-value:", format(p_value, digits = 2)), size = 5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.08)) +
  labs(title = "BP SNP-gene pairs in\nend-organ damage pairs",
       x = "Evidence Type",
       y = "Proportion") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
  )




bp_snp_gene_prox = unique(paste(proximal_merge_bp$SNP, proximal_merge_bp$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge_bp$SNP, eqtl_merge_bp$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge_bp$SNP, e2g_merge_bp$gene, sep="-"))
bp_snp_list = unique(c(bp_snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_prox = unique(paste(proximal_merge_nb$SNP, proximal_merge_nb$Gene_ID, sep="-"))
snp_gene_eqtl = unique(paste(eqtl_merge_nb$SNP, eqtl_merge_nb$gene_id_mod, sep="-"))
snp_gene_e2g = unique(paste(e2g_merge_nb$SNP, e2g_merge_nb$gene, sep="-"))
nb_snp_list = unique(c(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

data <- data.frame(
  Evidence = c("Combined", "Proximal-only"),
  Count = c(sum(bp_snp_list %in% nb_snp_list), 
            sum(bp_snp_gene_prox %in% nb_snp_list))
) %>%
  mutate(Total = c(length(bp_snp_list), length(bp_snp_gene_prox)),
         Proportion = Count / Total)

ggplot(data, aes(x = Evidence, y = Proportion)) +
  geom_bar(stat = "identity", fill="red") +
  geom_text(aes(label = scales::percent(Proportion, accuracy = 0.01)), vjust = -0.5) +
  # annotate("text", x = 1.5, y = max(data$Proportion) + 0.02, 
  #          label = paste("Hypergeometric test\np-value:", format(p_value, digits = 2)), size = 5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.08)) +
  labs(title = "BP SNP-gene pairs in\nend-organ damage pairs",
       x = "Evidence Type",
       y = "Proportion") +
  theme(legend.position = "none",
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")
  )






################################################################################
### examine overlapped/interested SNP-gene pairs
intersect(intersect(snp_gene_prox, snp_gene_eqtl), snp_gene_e2g)

snp_gene_pairs <- list(
  list(snp_id = "rs28451064", gene = c("KCNE2", "SMIM11", "SLC5A3", "RCAN1", "MRPS6")),
  list(snp_id = "rs2782980", gene = c("ADRB1", "AFAP1L2", "NHLRC2", "CASP7")),
  list(snp_id = "rs1882961", gene = c("NRIP1", "SAMSN1", "MRPL23", "SOD1")),
  list(snp_id = "rs1173771", gene = c("NPR3")),
  list(snp_id = "rs217727", gene = c("H19", "MRPL23", "TNNT3", "IGF2", "LSP1", "SYT8", "CD81", "CTSD")),
  list(snp_id = "rs9833313", gene = c("SHOX2"))
)

for (pair in snp_gene_pairs) {
  snp_id <- pair$snp_id
  gene <- pair$gene
  
  gwas_filtered <- gwas_merge[gwas_merge$SNPS == snp_id & (gwas_merge$MAPPED_GENE == gene | gwas_merge$REPORTED.GENE.S. == gene), ]
  eqtl_filtered <- eqtl_merge[eqtl_merge$SNP == snp_id & eqtl_merge$gene_name.x %in% gene, ]
  e2g_filtered <- e2g_merge[e2g_merge$SNP == snp_id & e2g_merge$gene_name.x %in% gene, ]
  
  print(c(snp_id, gene))
  print(unique(gwas_filtered[, c("SNPS", "MAPPED_GENE", "REPORTED.GENE.S.")]))
  if(nrow(eqtl_filtered)>0){
    print(eqtl_filtered)
  }
  if(nrow(e2g_filtered)>0){
    print(e2g_filtered)
  }
  
}


deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]

deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]

snp_gene_pairs <- list(
  list(snp_id = "rs28451064", gene = c("Kcne2", "Smim11", "Slc5a3", "Rcan1", "Mrps6")),
  list(snp_id = "rs2782980", gene = c("Adrb1", "Afap1l2", "Nhlrc2", "Casp7")),
  list(snp_id = "rs1882961", gene = c("Nrip1", "Samsn1", "Mrpl23", "Sod1")),
  list(snp_id = "rs1173771", gene = c("Npr3")),
  list(snp_id = "rs217727", gene = c("H19", "Mrpl23", "Tnnt3", "Igf2", "Lsp1", "Syt8", "Cd81", "Ctsd")),
  list(snp_id = "rs9833313", gene = c("Shox2"))
)

for (pair in snp_gene_pairs) {
  snp_id <- pair$snp_id
  gene <- pair$gene
  
  deg_filtered <- deg_merged[deg_merged$gene_name %in% gene & deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.25, ]
  
  print(c(snp_id, gene))
  if(nrow(deg_filtered)>0){
    print(deg_filtered)
  }
}


################################################################################
### extract snp-trait list for yong's hyplotype analysis 
gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
snp_trait_df = unique(gwas_merge[,c("SNPS", "trait")])
write.table(snp_trait_df, "snp_trait.for_haplotype.out", col.names = T, row.names = F, sep = "\t", quote = F)






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
background_500k <- read.table("proximal_merged.500k.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]
background_500k = background_500k[grepl("rs", background_500k$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)
snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_500k <- data.frame(SNP = background_500k$SNP, gene = background_500k$Gene_ID)
background_df <- unique(rbind(snp_gene_500k, snp_gene_eqtl, snp_gene_e2g))

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
    
    background_gene_list <- unique(c(m2h[m2h$gene_id %in% background_df$gene, ]$gene_name_mouse,
                                     r2h[r2h$gene_id %in% background_df$gene, ]$gene_name_rat))
    porpotion_background = length(intersect(background_gene_list, deg_list)) / length(unique(background_df$gene))
    
    hit_genes <- intersect(mapped_gene_list, deg_list)
    num_hit_genes <- length(hit_genes)
    proportion_hit_genes <- num_hit_genes / length(snp_gene_list)
    num_DEG_snps = length(unique(snp_gene_df_use[(! is.na(snp_gene_df_use$gene_name_mouse) & snp_gene_df_use$gene_name_mouse %in% deg_list) | 
                                                   (!is.na(snp_gene_df_use$gene_name_rat) & snp_gene_df_use$gene_name_rat %in% deg_list), ]$SNP))
    
    p_value <- one_sided_z_test(proportion_hit_genes, porpotion_background, length(snp_gene_list))
    
    results <- rbind(results, data.frame(
      trait = trait,
      SNPS = length(unique(snp_gene_df_use$SNP)),
      SNP_genes = length(snp_gene_list),
      mapped_genes = length(mapped_gene_list),
      threshold = threshold_name,
      DEGs = num_hit_genes,
      DEG_SNPS = num_DEG_snps,
      proportion_DEG = proportion_hit_genes,
      porpotion_background = porpotion_background,
      proportion_SNP = num_DEG_snps/length(unique(snp_gene_df_use$SNP)),
      p_value = p_value
    ))
    
  }
}

results$p.adj = p.adjust(results$p_value, method = "BH")
results$label <- paste(results$DEGs, 
                           "(", round(results$proportion_DEG * 100, 2), "%)", 
                           ifelse(results$p.adj < 0.05, "*", ""), sep = "")

results$threshold <- factor(results$threshold, levels = rev(names(thresholds)))

write.table(results, "trait.deg_prop.z_test.out", sep = ",", col.names = T, row.names = F, quote = F)




results_use = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/trait.deg_prop.z_test.out", sep = ",", header=T)
results_use$threshold = gsub("& ", "&\n", results_use$threshold)
threshold_col = RColorBrewer::brewer.pal(5, "Blues")[2:5]
names(threshold_col) = unique(results_use$threshold)

results_use$trait <- factor(results_use$trait, levels = results_use %>% group_by(trait) %>% summarise(total_hits = sum(proportion_DEG)) %>% arrange(total_hits) %>% pull(trait))
vline_data <- results_use %>%
  group_by(threshold) %>%
  summarise(yintercept = unique(porpotion_background))

ggplot(results_use, aes(x = trait, y = proportion_DEG, fill = threshold, label = label)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_fill_manual(values = threshold_col) +
  geom_text(aes(x = trait, y = proportion_DEG, hjust = 0), position = position_dodge(0.9), size = 3) +
  labs(x = "Trait", y = "DEG proportion", fill = "DEG Threshold") +
  scale_y_continuous(limits = c(0, 0.85)) +
  theme(legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')) +
  coord_flip() +
  facet_wrap(~threshold, scale = "free_x", nrow = 1) +
  geom_hline(data = vline_data, aes(yintercept = yintercept), linetype = "dashed", color="grey")


results_use = results_use[results_use$threshold=="p.adj < 0.05 &\n|log2(FC)| > 0.25",]
results_use$trait <- factor(results_use$trait, levels = results_use %>% group_by(trait) %>% summarise(total_hits = sum(proportion_SNP)) %>% arrange(total_hits) %>% pull(trait))
vline_data <- results_use %>%
  group_by(threshold) %>%
  summarise(yintercept = unique(porpotion_background))

p = ggplot(results_use) +
  geom_segment(aes(x = trait, xend = trait, y = proportion_DEG, yend = proportion_SNP), color = "gray", size = 1) +
  geom_point(aes(x = trait, y = proportion_DEG, label = DEGs, color = "SNP-related genes"), size = 3) +
  geom_point(aes(x = trait, y = proportion_SNP, label = DEG_SNPS, color = "SNPs"), size = 3) +
  geom_text(aes(x = trait, y = proportion_DEG-0.12, label = DEGs), color = "cadetblue4", size = 3) +
  geom_text(aes(x = trait, y = proportion_SNP+0.12, label = DEG_SNPS), color = "chocolate4", size = 3) +
  scale_color_manual(name = "", values = c("SNP-related genes" = "cadetblue4", "SNPs" = "chocolate4")) +
  scale_y_continuous(limits = c(0.05, 0.85)) +
  labs(x = "Trait", y = "Proportion hit by DEGs") +
  theme(
    # legend.position = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(colour = 'black')
  ) +
  coord_flip()

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg_snp_number.png", width=460/96, height=178/96, dpi=300)






################################################################################
### cell type enrichment
gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
proximal_merge = read.table("proximal_merged.100k.txt", header=T, sep='\t', quote = "")
proximal_merge = read.table("gwas_snp_gene_merged.txt", header=T, sep='\t')
eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

fisher_test_df = c()
num_permutations = 1000
for(trait in unique(gwas_merge$trait)){
  
  SNPS = unique(gwas_merge[gwas_merge$trait==trait, ]$SNPS)
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




################################################################################
### visualize fisher results using dotplot
fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]

# permutation_test_df = read.table("trait.permut.out", header = T, sep = "\t", comment.char = "")

test_df = fisher_test_df
# test_df <- test_df %>%
#   mutate(log_p_value = -log10(p.value))
# fit <- lm(log_p_value ~ X.DEG, data = test_df)
# test_df <- test_df %>%
#   mutate(NES = residuals(fit))
test_df$tissue = factor(test_df$tissue, levels = tissue_order)
test_df$strain = factor(test_df$strain, levels = c("C57BL/6", "SS", "SHR"))
test_df$cell_type = factor(test_df$cell_type, levels = cell_order)
test_df$trait = factor(test_df$trait, 
                       levels = c("systolic_bp", "diastolic_bp", "pulse_pressure", "essential_hypertension", "stroke", "CAD", "eGFR", "BUN", "albuminuria"))

test_df$DEG_prop = test_df$X.SNP.DEG/test_df$X.DEG
epsilon <- 1e-6
test_df$p.value <- ifelse(test_df$p.value == 0, epsilon, test_df$p.value)

test_df_use = test_df[test_df$p.adj<0.05 & test_df$X.SNP.DEG>=5, ]
p <- ggplot() +
  geom_point(data = test_df_use,
             aes(x = strain, y = cell_type, size = DEG_prop, fill = NES),
             shape = 21) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits = rev) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    panel.spacing.y = unit(0.2, "lines"),
    panel.spacing.x = unit(0.1, "lines"),
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black')
  ) +
  labs(fill = "NES", size = "Pct. DEG in\nSNP-gene", y = "", x = "") +
  facet_nested(tissue ~ trait, scales = "free", space = "free")
print(p)




test_df_use = test_df[!(is.na(test_df_use$NES)),]
p <- ggplot() +
  geom_point(data = test_df_use[test_df_use$p.adj >= 0.01 & test_df_use$X.SNP.DEG >= 1, ],
             aes(x = strain, y = cell_type, size = X.SNP.DEG/10, fill = NES),
             shape = 21, stroke = 0, color = "black") +
  geom_point(data = test_df_use[test_df_use$p.adj < 0.01 & test_df_use$X.SNP.DEG >= 1, ],
             aes(x = strain, y = cell_type, size = X.SNP.DEG/10, fill = NES),
             shape = 21, stroke = 0.3, color = "black") +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks=c(5, 10, 15), labels = c(50, 100, 150)) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    panel.spacing.y = unit(0.2, "lines"),
    panel.spacing.x = unit(0.1, "lines"),
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size=8),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text.x = element_text(colour = 'black', angle = 90, hjust = 0, margin=margin(l=0))
  ) +
  labs(fill = "NES", size = "Number of\nDE SNP-gene", y = "", x = "Model") +
  facet_nested(tissue ~ trait, scales = "free", space = "free")

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_type.gwas.fisher.png", width=693/96, height=859/96, dpi=300)





################################################################################
### highlight cell type-trait pairs

combine_pvalues <- function(pvalues) {
  chisq_stat <- -2 * sum(log(pvalues))
  df <- 2 * length(pvalues)
  p_value <- pchisq(chisq_stat, df, lower.tail = FALSE)
  return(p_value)
}

fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]

plot_data <- fisher_test_df %>%
  group_by(tissue, cell_type, trait) %>%
  mutate(average_NES = mean(as.numeric(NES), na.rm = TRUE),
         combined_p_value = combine_pvalues(p.value)) %>%
  ungroup() %>%
  filter(p.adj<0.01) %>%
  group_by(tissue, cell_type, trait) %>%
  mutate(significant_count = n_distinct(strain),
         label = paste(tissue, cell_type, trait, sep = "-"),
         Tissue = factor(tissue, levels = names(tissue_col))) %>%
  select(average_NES, combined_p_value, significant_count, label, Tissue, trait) %>%
  distinct()

highlight_points_conv <- plot_data %>%
  filter(significant_count>1)
highlight_points_sig <- plot_data %>%
  group_by(label) %>%
  mutate(sum_p = sum(-log10(combined_p_value))) %>%
  arrange(desc(sum_p))
highlight_points = rbind(highlight_points_conv,
                         highlight_points_sig)

ggplot(plot_data, aes(x = -log10(combined_p_value), y = average_NES, 
                      size = as.factor(significant_count),
                      color = as.factor(significant_count))) +
  geom_point(shape=21) +
  scale_color_manual(values = c("salmon", "red", "darkred")) +
  geom_label_repel(data = highlight_points, 
                   aes(label = label, fill = Tissue), color = 'white', size = 3,
                   segment.color = 'black', 
                   force = 2, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.5,
                   nudge_x = 0.5 * max(-log10(plot_data$combined_p_value)), 
                   nudge_y = 0.2 * max(plot_data$average_NES)) +
  scale_fill_manual(values = tissue_col) +
  labs(x = "-log10(combined p-value)",
       y = "Average NES",
       size = "Number of\nsignificant models",
       color = "Number of\nsignificant models",
       title = "Convergent trait-cell type associations (p.adj<0.01)")


highlight_points <- plot_data %>%
  arrange(combined_p_value) %>%
  group_by(trait) %>%
  slice_head(n = 1)# plot_data %>% filter(average_NES > 20, significant_count>1) 

ggplot(plot_data, aes(x = -log10(combined_p_value), y = average_NES, 
                      size = as.factor(significant_count),
                      color = as.factor(significant_count))) +
  geom_point(shape=21) +
  scale_color_manual(values = c("salmon", "red")) +
  geom_label_repel(data = highlight_points, 
                   aes(label = label, fill = Tissue), color = 'white', size = 3,
                   segment.color = 'black', 
                   force = 2, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.5,
                   nudge_x = 0.5 * max(-log10(plot_data$combined_p_value)), 
                   nudge_y = 0.2 * max(plot_data$average_NES)) +
  scale_fill_manual(values = tissue_col) +
  labs(x = "-log10(combined p-value)",
       y = "Average NES",
       size = "Number of\nsignificant models",
       color = "Number of\nsignificant models",
       title = "Top trait-cell type associations")

# facet_wrap(~ trait)
# width=599&height=421

# combine_pvalues <- function(pvalues) {
#   chisq_stat <- -2 * sum(log(pvalues))
#   df <- 2 * length(pvalues)
#   p_value <- pchisq(chisq_stat, df, lower.tail = FALSE)
#   return(p_value)
# }
# combined_pvalues <- fisher_test_df %>%
#   group_by(tissue, cell_type, trait) %>%
#   summarise(combined_p_value = combine_pvalues(p.value))
# 
# plot_data <- significant_counts %>%
#   inner_join(combined_pvalues, by = c("tissue", "cell_type", "trait"))
# plot_data <- plot_data %>%
#   mutate(label = paste(tissue, cell_type, trait, sep = "-"),
#          Tissue = factor(tissue, levels = names(tissue_col))
#   )
# highlight_points <- plot_data %>% filter(combined_p_value < 0.001)  # Example filter, adjust as needed
# 
# epsilon <- 1e-6
# plot_data$combined_p_value <- ifelse(plot_data$combined_p_value == 0, epsilon, plot_data$combined_p_value)
# 
# ggplot(plot_data, aes(x = significant_count, y = -log10(combined_p_value))) +
#   geom_point() +
#   geom_label_repel(data = subset(plot_data, combined_p_value < 0.001), 
#                    aes(label = label, fill = Tissue), color = 'white', size = 3,
#                    segment.color = 'black') +
#   scale_x_continuous(breaks = c(1, 2, 3), limits = c(1, 3)) +
#   scale_fill_manual(values = tissue_col) +
#   lims(y = c(0, 8)) +
#   labs(x = "Number of significant models",
#        y = "-log10(combined p-value)")

plot_data = plot_data[plot_data$tissue %in% c("HYP", "LK", "MCA"),]
ggplot(plot_data, aes(x = significant_count, y = -log10(combined_p_value))) +
  geom_point() +
  geom_label_repel(data = subset(plot_data, combined_p_value < 0.01), 
                   aes(label = label, fill = Tissue), color = 'white', size = 3,
                   nudge_y = 0.1, segment.color = 'black') +
  scale_x_continuous(breaks = c(1, 2, 3), limits = c(1, 3)) +
  scale_fill_manual(values = tissue_col) +
  lims(y = c(0, 8)) +
  labs(x = "Number of significant models",
       y = "-log10(combined p-value)") +
  facet_nested(. ~ trait, scales = "free", space = "free")











################################################################################
### nominating SNP-gene pairs from cell type-trait association
#### process enrichment result
fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]

#### process snp-gene table
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)
snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")
snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df = snp_gene_df[rowSums(is.na(snp_gene_df))==0,]

snp_gene_df_mod <- snp_gene_df %>%
  mutate(combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(SNP, Gene.name, combined_gene_name) %>%
  unique()  %>%
  mutate(pair = paste0(SNP, " - ", Gene.name))

### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
# deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
# deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', header=T)
colnames(deg_pseudo)[12] = "strain"
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$strain = factor(deg_pseudo$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "strain", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = cell_order)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$strain = factor(expr_all$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))



gwas_merge_use <- gwas_merge[gwas_merge$trait %in% trait, ] %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

deg_snp_df = snp_gene_df_mod %>%
  merge(., gwas_merge_use, by.x="SNP", by.y="SNPS") %>%
  mutate(DEG =ifelse(combined_gene_name %in% deg_filtered$gene_name, "Yes", "No")) %>%
  group_by(Gene.name) %>%
  mutate(n_SNP = n_distinct(SNP),
         n_SNP_trait = n()) 

deg_snp_df_long <- deg_snp_df %>%
  select(Gene.name, DEG, n_SNP, n_SNP_trait) %>%
  unique() %>%
  gather(key = "variable", value = "value", -c(Gene.name, DEG))

plot1 <- ggplot(deg_snp_df_long, aes(x = value, fill = DEG, color = DEG)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
  # geom_density(alpha = 0.2) +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Distribution of SNP and Trait Metrics by DEG Status", x = "Count", y = "Density") +
  theme(legend.position = "top") +
  coord_flip()

plot2 <- deg_snp_df %>%
  select(pair, DEG, P.VALUE) %>%
  unique() %>%
  ggplot(., aes(x = -log10(P.VALUE), fill = DEG, color = DEG)) +
  geom_histogram(aes(y = ..density..), position = "identity", alpha = 0.6, bins = 30) +
  # geom_density(alpha = 0.2) +
  labs(x = "-log10(P.VALUE)", y = "Density") +
  theme(legend.position = "top") +
  coord_flip()

combined_plot <- plot1 + plot2 + plot_layout(nrow = 1, widths = c(0.5, 0.25))
print(combined_plot)




strain = "C57BL/6" #c("C57BL/6", "SS", "SD", "SHR", "WKY")
tissue = "HYP"
cell_type = "EC"
trait = c("systolic_bp", "diastolic_bp", "pulse_pressure")
rm_bg = FALSE

strain = "C57BL/6" #c("C57BL/6", "SS", "SD", "SHR", "WKY")
tissue = "HYP"
cell_type = "Avp+ neuron"
trait = c("systolic_bp", "diastolic_bp", "pulse_pressure")
rm_bg = FALSE

# strain = c("C57BL/6", "SHR")
# tissue = "HYP"
# cell_type = "Astrocyte"
# trait = c("systolic_bp", "diastolic_bp", "pulse_pressure")

strain = c("SS", "SD")
tissue = "MSA"
cell_type = "VSMC"
trait = c("systolic_bp", "diastolic_bp", "pulse_pressure", "CAD", "eGFR")
rm_bg = TRUE

strain = c("C57BL/6", "SS", "SD", "SHR", "WKY")
tissue = "LK"
cell_type = "TAL"
trait = c("BUN")
rm_bg = FALSE


de_snp_gene = strsplit(fisher_test_df[fisher_test_df$p.adj<0.01 & fisher_test_df$strain %in% strain & fisher_test_df$tissue %in% tissue &
                                        fisher_test_df$cell_type %in% cell_type & fisher_test_df$trait %in% trait, ]$gene_list, ", ")
if(rm_bg){
  ctrl_deg = unique(deg_merged[deg_merged$p_val_adj <0.05 & deg_merged$strain %in% c("SD") & 
                                 deg_merged$tissue %in% tissue &
                                 deg_merged$cell_type %in% cell_type, ]$gene_name)
  de_snp_gene <- unique(setdiff(unlist(de_snp_gene), ctrl_deg))
}else{
  de_snp_gene <- unique(unlist(de_snp_gene))
}

snp_trait_df = snp_gene_df_mod %>%
  filter(combined_gene_name %in% de_snp_gene) %>%
  merge(., gwas_merge_use, by.x="SNP", by.y="SNPS") %>%
  group_by(Gene.name) %>% mutate(sum_p = sum(-log10(P.VALUE))) %>% 
  arrange(desc(sum_p))

top10_genes = head(unique(snp_trait_df$combined_gene_name), n=10)
# top5_snps_to_genes = snp_trait_df %>%
#   filter(combined_gene_name %in% top10_genes) %>% 
#   group_by(SNP, Gene.name) %>% summarise(sum_p = sum(-log10(P.VALUE))) %>% 
#   ungroup() %>% group_by(Gene.name) %>% 
#   arrange(desc(sum_p)) %>% slice_head(n = 5)

top5_snps_to_genes = snp_trait_df %>%
  filter(combined_gene_name %in% top10_genes) %>% 
  group_by(SNP, combined_gene_name) %>% 
  summarise(sum_p = sum(-log10(P.VALUE)), .groups = 'drop') %>% 
  ungroup() %>% 
  mutate(combined_gene_name = factor(combined_gene_name, levels = top10_genes)) %>% 
  group_by(combined_gene_name) %>% 
  arrange(combined_gene_name, desc(sum_p)) %>% 
  slice_head(n = 5)

snp_use = unique(top5_snps_to_genes$SNP)


deg_use = expr_all[expr_all$strain %in% strain & expr_all$tissue %in% tissue &
                     expr_all$cell_type %in% cell_type & expr_all$gene_name %in% top10_genes, ]
deg_use$gene_name = factor(deg_use$gene_name, levels = top10_genes)
p1=ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(fill = log2(avg_expr+1), size=pct*100), shape=21) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="Top 10 genes by cumulative\n-log10(GWAS p-value) of associated loci", fill="log2(expression+1)", size="Percentage", title = paste0(tissue, "-", cell_type)) +
  facet_nested( ~ project + strain, scales = "free", space = "free")


snp_df_use = gwas_merge_use %>% 
  filter(SNPS %in% snp_use) %>%
  group_by(SNPS) %>%
  mutate(sum_p = sum(-log10(P.VALUE))) %>% 
  arrange(desc(sum_p))
snp_df_use$SNPS = factor(snp_df_use$SNPS, levels = snp_use)
p2 = ggplot(snp_df_use, aes(x = trait, y = SNPS)) +
  geom_tile(aes(fill = -log10(P.VALUE))) +
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
  labs(x="", y="", fill="-log10(p-value)", title = "Top 5 SNPs")

gene_names <- layer_scales(p1)$y$range$range
total_genes <- length(gene_names)
relative_positions <- seq(from = 0.5/total_genes, to = (total_genes-0.5)/total_genes, length.out = total_genes)
gene_positions <- relative_positions; names(gene_positions) = rev(gene_names)

snp_names <- layer_scales(p2)$y$range$range
total_snps <- length(snp_names)
relative_positions <- seq(from = 0.5/total_snps, to = (total_snps-0.5)/total_snps, length.out = total_snps)
snp_positions <- relative_positions; names(snp_positions) = rev(snp_names)

# snp_gene_pair = top5_snps_to_genes %>% merge(., unique(snp_trait_df[, c("Gene.name", "combined_gene_name")]))
snp_gene_pair = unique(top5_snps_to_genes[, c("SNP", "combined_gene_name")])
snp_gene_pair$y1 = gene_positions[as.character(snp_gene_pair$combined_gene_name)]
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

combined_plot <- (p1 + p3 + p2) + plot_layout(guides = 'collect', widths = c(0.5, 0.2, 0.5)) &
  theme(legend.position = 'bottom', 
        legend.box = 'vertical')
print(combined_plot)
















################################################################################
### nominating SNP-gene pairs from cell type-trait association
#### 
# fisher_test_df = read.table("trait.fisher.cluster.out", header = T, sep = "\t", comment.char = "")
trait = c("systolic_bp", "diastolic_bp", "pulse_pressure")

snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)

deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
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
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/hyp_ec.npr3.png", width=400/96, height=329/96, dpi=300)



combinations <- list(
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "EC", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Astrocyte", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Myelinating OL", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "Microglia", tissue_var = "HYP"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "CM", tissue_var = "LV"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure", "CAD"), cell_type_var = "Fibroblast", tissue_var = "LV"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "EC", tissue_var = "LK"),
  list(trait_var = c("BUN"), cell_type_var = "TAL", tissue_var = "LK"),
  list(trait_var = c("diastolic_bp", "systolic_bp", "pulse_pressure"), cell_type_var = "VSMC", tissue_var = "MSA")
)

# Initialize a list to store the plots
plots <- list()

# Iterate through each combination
for (i in seq_along(combinations)) {
  comb <- combinations[[i]]
  
  # Generate the plot using the generate_dotplot function
  plot <- generate_dotplot(
    trait = comb$trait_var,
    cell_type_var = comb$cell_type_var,
    tissue_var = comb$tissue_var,
    gwas_data = gwas_merge_use,
    snp_gene_data = snp_gene_df,
    deg_data = deg_filtered
  )
  
  # Add the plot to the list, naming it with the combination
  plots[[paste(comb$cell_type_var, comb$tissue_var, paste(comb$trait_var, collapse = "_"), sep = "_")]] <- plot
}

# Combine all plots into a single layout
combined_plot <- wrap_plots(plots, ncol = 3) # Adjust ncol as needed

# Print or save the combined plot
print(combined_plot)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/trait_cell_type.top10.png", width=1400/96, height=1200/96, dpi=300)












################################################################################
### nominating multifunctioal SNPs
#### process enrichment result
fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]

#### process snp-gene table
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)
snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")
snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df = snp_gene_df[rowSums(is.na(snp_gene_df))==0,]

snp_gene_df_mod <- snp_gene_df %>%
  mutate(combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(SNP, Gene.name, combined_gene_name) %>%
  unique()  %>%
  mutate(pair = paste0(SNP, " - ", Gene.name))

### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', header=T)
colnames(deg_pseudo)[12] = "strain"
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$strain = factor(deg_pseudo$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "strain", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = cell_order)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$strain = factor(expr_all$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))


gwas_merge_use <- gwas_merge[gwas_merge$trait %in% trait, ] %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)


snp_deg_df = snp_gene_df_mod %>%
  group_by(SNP) %>%
  mutate(num_gene = n_distinct(Gene.name)) %>%
  ungroup() %>%
  # merge(., gwas_merge_use, by.x="SNP", by.y="SNPS") %>%
  merge(., deg_filtered, by.x="combined_gene_name", by.y="gene_name") %>%
  group_by(SNP) %>%
  mutate(num_deg = n_distinct(Gene.name),
         num_cell = n_distinct(strain, tissue, treatment, cell_type),
         deg_prop = num_deg/num_gene) %>%
  ungroup() %>%
  group_by(num_deg, deg_prop) %>%
  mutate(snp_count = n_distinct(SNP))
  

snp_list = c("rs28451064", "rs2782980", "rs1882961", "rs1173771", "rs217727", "rs9833313")

gene_list = c("KCNE2", "SMIM11", "SLC5A3", "RCAN1", "MRPS6", "ADRB1", "AFAP1L2", "NHLRC2", "CASP7",
              "NRIP1", "SAMSN1", "MRPL23", "SOD1", "NPR3", "H19", "MRPL23", "TNNT3", "IGF2", "LSP1", "SYT8", "CD81", "CTSD", "SHOX2")

snp_deg_df %>%
  select(SNP, deg_prop, num_deg, snp_count) %>%
  unique() %>%
  ggplot(., aes(x = deg_prop, y = num_deg, label = ifelse(SNP %in% snp_list, SNP, ""))) +
  geom_point(aes(size = snp_count), color = "black") +
  geom_text(nudge_y = 0.5, nudge_x = -0.05, color = "red", size = 3) +
  labs(x = "DEG proportion in all associated genes", y = "Number of associated DEG", size = "SNP count") +
  theme(legend.position = "right")











snp_list = c("rs28451064", "rs2782980", "rs1882961", "rs1173771", "rs217727", "rs9833313")
gene_list = c("KCNE2", "SMIM11", "SLC5A3", "RCAN1", "MRPS6", "ADRB1", "AFAP1L2", "NHLRC2", "CASP7",
              "NRIP1", "SAMSN1", "MRPL23", "SOD1", "NPR3", "H19", "MRPL23", "TNNT3", "IGF2", "LSP1", "SYT8", "CD81", "CTSD", "SHOX2")

snp_gene_pairs <- list(
  list(snp_id = "rs28451064", gene = c("KCNE2", "SMIM11", "SLC5A3", "RCAN1", "MRPS6")),
  list(snp_id = "rs2782980", gene = c("ADRB1", "AFAP1L2", "NHLRC2", "CASP7")),
  list(snp_id = "rs1882961", gene = c("NRIP1", "SAMSN1", "MRPL23", "SOD1")),
  list(snp_id = "rs1173771", gene = c("NPR3")),
  list(snp_id = "rs217727", gene = c("H19", "MRPL23", "TNNT3", "IGF2", "LSP1", "SYT8", "CD81", "CTSD")),
  list(snp_id = "rs9833313", gene = c("SHOX2"))
)

for (pair in snp_gene_pairs) {
  snp_id <- pair$snp_id
  gene <- pair$gene
  
  deg_filtered <- snp_gene_df_mod[snp_gene_df_mod$Gene.name %in% gene & snp_gene_df_mod$SNP %in% snp_id, ]
  
  print(c(snp_id, gene))
  if(nrow(deg_filtered)>0){
    print(deg_filtered)
  }
}

snp_gene_df_mod %>%
  filter(DEG == "Yes") %>%
  select(SNP, gene, Gene.name, trait) %>%
  unique() %>%
  group_by(gene, Gene.name) %>%
  mutate(total_traits = n_distinct(trait),
         total_SNPs = n_distinct(SNP)) %>%
  select(Gene.name, total_traits, total_SNPs) %>%
  unique() %>%
  ggplot(., aes(x = total_traits, y = total_SNPs, label = ifelse(Gene.name %in% gene_list, Gene.name, ""))) +
  geom_point(color = "black") +
  geom_text(nudge_y = 0.5, color = "red", size = 3) +
  labs(x = "Number of traits", y = "Number of SNPs") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")


snp_gene_df_mod %>%
  filter(DEG == "Yes") %>%
  select(SNP, gene, Gene.name, trait) %>%
  unique() %>%
  group_by(SNP) %>%
  mutate(total_traits = n_distinct(trait),
         total_genes = n_distinct(gene)) %>%
  select(SNP, total_traits, total_genes) %>%
  unique() %>%
  ggplot(., aes(x = total_traits, y = total_genes, label = ifelse(SNP %in% snp_list, SNP, ""))) +
  geom_point(color = "black") +
  geom_text(nudge_y = 0.5, color = "red", size = 3) +
  labs(x = "Number of traits", y = "Number of SNPs") +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")



for (pair in snp_gene_pairs) {
  snp_id <- pair$snp_id
  gene <- pair$gene
  
  deg_filtered <- trait_snp_count[trait_snp_count$Gene.name %in% gene & trait_snp_count$SNP %in% snp_id, ]
  
  print(c(snp_id, gene))
  if(nrow(deg_filtered)>0){
    print(deg_filtered)
  }
}


trait_counts <- snp_gene_df_mod %>%
  filter(DEG=="Yes") %>%
  select(gene, Gene.name, trait) %>%
  unique() %>%
  group_by(gene, Gene.name) %>%
  summarise(total_traits = n(),
            deg_traits = c())

trait_counts <- trait_counts %>%
  arrange(desc(total_traits)) %>%
  group_by(total_traits) %>%
  mutate(cum_count = n()) %>%
  mutate(cum_count = n())


# Plot the cumulative count of genes with shared traits
ggplot(trait_counts, aes(x = total_traits, y = cum_count)) +
  geom_step(color = "red") +
  geom_area(fill = "gray", alpha = 0.5) +
  labs(x = "Number of shared traits", y = "Number of genes") +
  theme_minimal() +
  theme(legend.position = "none")

# Highlight top 10 genes
top_genes <- trait_counts %>% top_n(10, total_traits)

# Inset plot for top genes
inset_plot <- ggplot(top_genes, aes(x = total_traits, y = cum_count, label = Gene.name)) +
  geom_point(color = "red") +
  geom_text(nudge_y = 2, color = "red", size = 3) +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), legend.position = "none")

# Combine main plot and inset plot
library(cowplot)
main_plot <- ggplot(trait_counts, aes(x = total_traits, y = cum_count)) +
  geom_step(color = "red") +
  geom_area(fill = "gray", alpha = 0.5) +
  labs(x = "Number of shared traits", y = "Number of genes") +
  theme_minimal() +
  theme(legend.position = "none")

final_plot <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(inset_plot, x = 0.5, y = 0.5, width = 0.4, height = 0.4)

print(final_plot)

































################################################################################
# Create a binary matrix for cell types and genes
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(vegan)
library(circlize)

#### process snp-gene table
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID)
snp_gene_eqtl <- data.frame(SNP = eqtl_merge$SNP, gene = eqtl_merge$gene_id_mod)
snp_gene_e2g <- data.frame(SNP = e2g_merge$SNP, gene = e2g_merge$gene)
snp_gene_df <- unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_df$pair <- paste(snp_gene_df$SNP, snp_gene_df$gene, sep = "-")
snp_gene_df <- base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x = "gene", by.y = "Gene.stable.ID", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df <- base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x = "gene", by.y = "gene_id", all.x = TRUE)
snp_gene_df = snp_gene_df[rowSums(is.na(snp_gene_df))==0,]

snp_gene_df_mod <- snp_gene_df %>%
  mutate(combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(SNP, Gene.name, combined_gene_name) %>%
  unique()  %>%
  mutate(pair = paste0(SNP, " - ", Gene.name))

### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
deg_use = deg_filtered[deg_filtered$gene_name %in% snp_gene_df_mod$combined_gene_name, ]

### generate binary table
deg_use <- deg_use %>%
  mutate(condition = paste(project, strain, tissue, control, treatment, cell_type, sep = "_"))
binary_matrix <- deg_use %>%
  select(gene_name, condition) %>%
  mutate(value = 1) %>%
  spread(key = condition, value = value, fill = 0)
row.names(binary_matrix) <- binary_matrix$gene_name
binary_matrix <- binary_matrix %>% select(-gene_name)

jaccard_dist_rows <- vegdist(binary_matrix, method = "jaccard")
jaccard_dist_columns <- vegdist(t(binary_matrix), method = "jaccard")


column_info <- deg_use[match(colnames(binary_matrix), deg_use$condition), ]
ha = HeatmapAnnotation(Strain = column_info$strain, Cell = column_info$cell_type, Tissue = column_info$tissue, 
                       col = list(Strain = species_col, Cell = cell_col, Tissue = tissue_col))
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
heatmap <- Heatmap(binary_matrix, 
                   name = "DEG\npresence",
                   col = col_fun,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = as.dist(jaccard_dist_rows),
                   clustering_distance_columns = as.dist(jaccard_dist_columns),
                   clustering_method_rows = "average",
                   clustering_method_columns = "average",
                   top_annotation = ha,
                   heatmap_legend_param = list(title = "DEG\npresence", legend_direction = "vertical"))

draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")

ordered_row_indices <- row_order(heatmap)
ordered_col_indices <- column_order(heatmap)
ordered_matrix <- binary_matrix[ordered_row_indices, ordered_col_indices]
ordered_row_names <- rownames(ordered_matrix)
ordered_column_names <- colnames(ordered_matrix)
ordered_matrix_with_names <- cbind(ordered_row_names, ordered_matrix)
write.csv(ordered_matrix_with_names, file = "ordered_binary_matrix.csv", row.names = FALSE)
ordered_matrix_df <- as.data.frame(ordered_matrix)
ordered_matrix_df <- cbind(rownames(ordered_matrix), ordered_matrix_df)
colnames(ordered_matrix_df)[1] <- "Genes"
write.csv(ordered_matrix_df, file = "ordered_binary_matrix_with_names.csv", row.names = FALSE)


library(igraph)
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

subset_gene = rownames(ordered_matrix)[545:760]
subset_gene_snp = gwas_merge_use %>%
  merge(., snp_gene_df_mod, by.x="SNPS", by.y="SNP", all.x=T) %>%
  filter(combined_gene_name %in% subset_gene) %>%
  select(trait, Gene.name) %>% unique()
# subset_gene = deg_filtered %>%
#   filter(strain %in% c("SS", "SHR") & tissue=="LK" & cell_type=="TAL") %>%
#   select(gene_name) %>% as.vector()
# subset_snp = gwas_merge[gwas_merge$trait %in% "BUN", ]$SNPS
# subset_gene_snp = snp_gene_df_mod %>%
#   filter(combined_gene_name %in% subset_gene$gene_name &
#            SNP %in% subset_snp) %>%
#   select(SNP, Gene.name) %>% unique()

bipartite_graph <- graph_from_data_frame(subset_gene_snp, directed = FALSE)
V(bipartite_graph)$type <- bipartite_mapping(bipartite_graph)$type
snp_projection <- bipartite_projection(bipartite_graph)$proj1
snp_betweenness <- betweenness(snp_projection, directed = FALSE)
top_snps <- names(sort(snp_betweenness, decreasing = TRUE)[1:10])  # Top 10 SNPs
print(top_snps)
snp_betweenness[snp_list]

par(mfrow = c(1, 2))
plot(bipartite_graph, vertex.label = NA, vertex.size = 3, vertex.color = ifelse(V(bipartite_graph)$type, "lightblue", "orange"),
     main = "Bipartite SNP-Gene Network")
plot(snp_projection, vertex.size = sqrt(snp_betweenness) * 0.5, vertex.color = "orange",
     main = "Projected SNP Network with Betweenness Centrality")
par(mfrow = c(1, 1))  # Reset plotting area


par(mfrow = c(1, 2))
plot(bipartite_graph, vertex.size = 3, vertex.color = ifelse(V(bipartite_graph)$type, "lightblue", "orange"),
     main = "Bipartite SNP-Gene Network")
plot(snp_projection, vertex.size = sqrt(snp_betweenness) * 0.5, vertex.color = "orange",
     main = "Projected SNP Network with Betweenness Centrality")
par(mfrow = c(1, 1))  # Reset plotting area






subset_gene = rownames(ordered_matrix)[545:760]
subset_gene_snp = gwas_merge_use %>%
  merge(., snp_gene_df_mod, by.x="SNPS", by.y="SNP", all.x=T) %>%
  filter(combined_gene_name %in% subset_gene) %>%
  select(trait, Gene.name) %>% unique()
# subset_gene = deg_filtered %>%
#   filter(strain %in% c("SS", "SHR") & tissue=="LK" & cell_type=="TAL") %>%
#   select(gene_name) %>% as.vector()
# subset_snp = gwas_merge[gwas_merge$trait %in% "BUN", ]$SNPS
# subset_gene_snp = snp_gene_df_mod %>%
#   filter(combined_gene_name %in% subset_gene$gene_name &
#            SNP %in% subset_snp) %>%
#   select(SNP, Gene.name) %>% unique()

bipartite_graph <- graph_from_data_frame(subset_gene_snp, directed = FALSE)
V(bipartite_graph)$type <- bipartite_mapping(bipartite_graph)$type
snp_projection <- bipartite_projection(bipartite_graph)$proj1
snp_betweenness <- betweenness(snp_projection, directed = FALSE)
top_snps <- names(sort(snp_betweenness, decreasing = TRUE)[1:10])  # Top 10 SNPs
print(top_snps)
snp_betweenness[snp_list]

par(mfrow = c(1, 2))
plot(bipartite_graph, vertex.label = NA, vertex.size = 3, vertex.color = ifelse(V(bipartite_graph)$type, "lightblue", "orange"),
     main = "Bipartite SNP-Gene Network")
plot(snp_projection, vertex.size = sqrt(snp_betweenness) * 0.5, vertex.color = "orange",
     main = "Projected SNP Network with Betweenness Centrality")
par(mfrow = c(1, 1))  # Reset plotting area


par(mfrow = c(1, 2))
plot(bipartite_graph, vertex.size = 3, vertex.color = ifelse(V(bipartite_graph)$type, "lightblue", "orange"),
     main = "Bipartite SNP-Gene Network")
plot(snp_projection, vertex.size = sqrt(snp_betweenness) * 0.5, vertex.color = "orange",
     main = "Projected SNP Network with Betweenness Centrality")
par(mfrow = c(1, 1))  # Reset plotting area








################################################################################
### cell-type - trait network using NES
library(ggraph)

fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]

fisher_test_df <- fisher_test_df %>%
  mutate(cell_context = paste(strain, treatment, tissue, cell_type, sep = "_"))

# Filter for significant associations
significant_associations <- fisher_test_df %>%
  filter(p.adj < 0.05)

# Create an edge list for the graph
edges <- significant_associations %>%
  select(trait, cell_context, NES)

# Convert to graph
g <- graph_from_data_frame(edges, directed = FALSE)

V(g)$type <- ifelse(V(g)$name %in% edges$trait, "trait", "cell_context")

snp_projection <- bipartite_projection(g)$proj1
snp_betweenness <- betweenness(snp_projection, directed = FALSE)
top_snps <- names(sort(snp_betweenness, decreasing = TRUE)[1:10])  # Top 10 SNPs
print(top_snps)
snp_betweenness[snp_list]

par(mfrow = c(1, 2))
plot(g, vertex.label = NA, vertex.size = 3, vertex.color = ifelse(V(g)$type, "lightblue", "orange"),
     main = "Bipartite SNP-Gene Network")
plot(snp_projection, vertex.size = sqrt(snp_betweenness) * 0.5, vertex.color = "orange",
     main = "Projected SNP Network with Betweenness Centrality")
par(mfrow = c(1, 1))  # Reset plotting area

# Plot the network
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_alpha = NES), show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 5, show.legend = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_alpha(range = c(0.1, 1)) +
  theme_void() +
  labs(title = "Cell Type-Trait Association Network",
       edge_alpha = "NES")



### cell-type community - trait network
snp_gene_df_mod <- snp_gene_df %>%
  mutate(combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(SNP, Gene.name, combined_gene_name) %>%
  unique()  %>%
  mutate(pair = paste0(SNP, " - ", Gene.name))

### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
deg_use = deg_filtered[deg_filtered$gene_name %in% snp_gene_df_mod$combined_gene_name, ]

### generate binary table
deg_use <- deg_use %>%
  mutate(condition = paste(strain, treatment, tissue, cell_type, sep = "_"))
binary_matrix <- deg_use %>%
  select(gene_name, condition) %>%
  mutate(value = 1) %>%
  spread(key = condition, value = value, fill = 0)
row.names(binary_matrix) <- binary_matrix$gene_name
binary_matrix <- binary_matrix %>% select(-gene_name) %>% as.matrix()

# jaccard_dist_rows <- vegdist(binary_matrix, method = "jaccard")
# jaccard_dist_columns <- vegdist(t(binary_matrix), method = "jaccard")
# similarity_matrix <- as.dist(jaccard_dist_columns)
similarity_matrix <- cosine(binary_matrix)
similarity_graph <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "undirected", weighted = TRUE)

### test resolution
resolution_range <- seq(0.1, 2.0, by = 0.1)
num_communities <- numeric(length(resolution_range))
for (i in seq_along(resolution_range)) {
  resolution <- resolution_range[i]
  communities <- leiden(similarity_graph, resolution_parameter = resolution)
  num_communities[i] <- length(unique(communities))
}
plot(resolution_range, num_communities, type = "b", xlab = "Resolution Parameter", ylab = "Number of Communities",
     main = "Effect of Resolution Parameter on Community Detection")

leiden_communities <- leiden(similarity_graph, resolution_parameter = 1)  # Adjust the resolution parameter as needed
names(leiden_communities) <- V(similarity_graph)$name


fisher_test_df <- fisher_test_df %>%
  mutate(cell_context = paste(strain, treatment, tissue, cell_type, sep = "_")) %>%
  filter(cell_context %in% names(leiden_communities)) %>%
  mutate(community = leiden_communities[cell_context])

nes_summary <- fisher_test_df %>%
  group_by(trait) %>%
  mutate(
    avg_NES_others = map2_dbl(community, trait, ~mean(fisher_test_df %>% # purrr
                                                        filter(trait == .y & community != .x) %>%
                                                        pull(average_NES)))
  ) %>%
  ungroup()


plot_data <- fisher_test_df %>%
  group_by(community, trait) %>%
  summarise(
    avg_NES = mean(NES, na.rm = TRUE),
    num_snp_genes = n(),
    .groups = 'drop'
  )
plot_data$community_label <- paste("Community", plot_data$community, sep = " ")

ggplot(plot_data, aes(x = trait, y = community_label, size = num_snp_genes, color = avg_NES)) +
  geom_point() +
  scale_size_continuous(range = c(3, 15)) +
  scale_color_gradient(low = "lightblue", high = "red") +
  theme_minimal() +
  labs(
    title = "Cell Type-Trait Associations Grouped by Community",
    x = "Trait",
    y = "Community",
    size = "Number of SNP-Genes",
    color = "Average NES"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )








community_trait_summary <- fisher_test_df %>%
  group_by(community, trait) %>%
  summarise(mean_NES = mean(NES, na.rm = TRUE),
            total_NES = sum(NES, na.rm = TRUE),
            count = n()) %>%
  ungroup()
community_trait_edges <- community_trait_summary %>%
  filter(mean_NES > 0) %>%  # Adjust the threshold as needed
  select(community, trait, mean_NES, total_NES)
g <- graph_from_data_frame(community_trait_edges, directed = FALSE)
V(g)$type <- ifelse(V(g)$name %in% community_trait_edges$trait, "trait", "community")
ggraph(g, layout = "fr") +
  geom_edge_link(aes(edge_alpha = mean_NES, edge_width = total_NES), show.legend = TRUE) +
  geom_node_point(aes(color = type), size = 5, show.legend = TRUE) +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_alpha(range = c(0.1, 1)) +
  scale_edge_width(range = c(0.5, 5)) +
  theme_void() +
  labs(title = "Community-Trait Association Network",
       edge_alpha = "Mean NES",
       edge_width = "Total NES")















test_list = list(c(strain=c("SS"), tissue=c("MSA"), cell_type=c("VSMC"), trait=c("systolic_bp", "diastolic_bp", "pulse_pressure")),
                 c(strain=c(), tissue=c(), cell_type=c(), trait=c()),
                 c(strain=c(), tissue=c(), cell_type=c(), trait=c()),
                 c(strain=c(), tissue=c(), cell_type=c(), trait=c()),
                 c(strain=c(), tissue=c(), cell_type=c(), trait=c()),
                 c(strain=c(), tissue=c(), cell_type=c(), trait=c()))
for(ti in test_list){
  
  strain = ti$strain
  tissue = ti$tissue
  cell_type = ti$cell_type
  trait = ti$trait
  
  subset_gene = deg_filtered %>%
    filter(strain %in% strain & tissue %in% tissue & cell_type %in% cell_type) %>%
    select(gene_name) %>% as.vector()
  subset_snp = gwas_merge[gwas_merge$trait %in% trait, ]$SNPS
  subset_gene_snp = snp_gene_df_mod %>%
    filter(combined_gene_name %in% subset_gene$gene_name &
             SNP %in% subset_snp) %>%
    select(SNP, Gene.name) %>% unique()
  
  bipartite_graph <- graph_from_data_frame(subset_gene_snp, directed = FALSE)
  V(bipartite_graph)$type <- bipartite_mapping(bipartite_graph)$type
  snp_projection <- bipartite_projection(bipartite_graph)$proj1
  snp_betweenness <- betweenness(snp_projection, directed = FALSE)
  top_snps <- names(sort(snp_betweenness, decreasing = TRUE)[1:10])  # Top 10 SNPs
  # print(top_snps)
  snp_betweenness[snp_list]
  
}

test = unique(deg_filtered[, c("tissue", "cell_type")])
for(ci in 1:nrow(test)){
  for(trait in unique(gwas_merge$trait)){
    tissue = test[ci,]$tissue
    cell_type = test[ci,]$cell_type
    subset_gene = deg_filtered %>%
      filter(tissue %in% tissue & cell_type %in% cell_type) %>%
      select(gene_name) %>% as.vector()
    subset_snp = gwas_merge[gwas_merge$trait %in% trait, ]$SNPS
    subset_gene_snp = snp_gene_df_mod %>%
      filter(combined_gene_name %in% subset_gene$gene_name &
               SNP %in% subset_snp) %>%
      select(SNP, Gene.name) %>% unique()
    
    bipartite_graph <- graph_from_data_frame(subset_gene_snp, directed = FALSE)
    V(bipartite_graph)$type <- bipartite_mapping(bipartite_graph)$type
    snp_projection <- bipartite_projection(bipartite_graph)$proj1
    snp_betweenness <- betweenness(snp_projection, directed = FALSE)
    top_snps <- names(sort(snp_betweenness, decreasing = TRUE)[1:10])  # Top 10 SNPs
    # print(top_snps)
    if(sum(is.na(snp_betweenness[snp_list]))<6){
      print(c(tissue, cell_type, trait))
      print(snp_betweenness[snp_list])
      
    }
    
  }
  
}




Heatmap(ordered_matrix, 
                   name = "DEG\npresence",
                   col = col_fun,
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   top_annotation = ha,
                   heatmap_legend_param = list(title = "DEG\npresence", legend_direction = "vertical"))







################################################################################
### cluster of gene-cell type pattern
gwas_merge = read.table("gwas_merged.txt", header=T, sep='\t', quote = "")
proximal_merge = read.table("proximal_merged.txt", header=T, sep='\t', quote = "")
eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')

SNP = unique(gwas_merge$SNPS)
eqtl_merge_use = eqtl_merge[eqtl_merge$SNP %in% SNPS, ]
e2g_merge_use = e2g_merge[e2g_merge$SNP %in% SNPS, ]
proximal_merge_use = proximal_merge[proximal_merge$SNP %in% SNPS, ]

snp_gene_prox = data.frame(SNP=gwas_merge_use$SNP, gene=gwas_merge_use$Gene_ID)
snp_gene_eqtl = data.frame(SNP=eqtl_merge_use$SNP, gene=eqtl_merge_use$gene_id_mod)
snp_gene_e2g = data.frame(SNP=e2g_merge_use$SNP, gene=e2g_merge_use$gene)
snp_gene_df = unique(rbind(snp_gene_prox, snp_gene_eqtl, snp_gene_e2g))

snp_gene_df = base::merge(snp_gene_df, ensembl[snp_gene_df$gene, c("Gene.stable.ID", "Gene.name")], by.x="gene", by.y="Gene.stable.ID", all.x=T)
snp_gene_df = base::merge(snp_gene_df, r2h[r2h$gene_id %in% snp_gene_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
snp_gene_df = base::merge(snp_gene_df, m2h[m2h$gene_id %in% snp_gene_df$gene, ], by.x="gene", by.y="gene_id", all.x=T)
snp_gene_df[snp_gene_df$gene_name_rat=="" & !(is.na(snp_gene_df$gene_name_rat)), ]$gene_name_rat = snp_gene_df[snp_gene_df$gene_name_rat=="" & !(is.na(snp_gene_df$gene_name_rat)), ]$gene_id_rat
snp_gene_df[snp_gene_df$gene_name_mouse=="" & !(is.na(snp_gene_df$gene_name_mouse)), ]$gene_name_mouse = snp_gene_df[snp_gene_df$gene_name_mouse=="" & !(is.na(snp_gene_df$gene_name_mouse)), ]$gene_id_mouse


deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]



























################################################################################
### load GWAS data
gwas = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas_catalog_blood pressure.tsv", sep = '\t', header=T)

eqtl = read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.txt", header=T)
eqtl = read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Brain_Hypothalamus.v8.signif_variant_gene_pairs.txt", header=T)


### pre-filter gwas info
gwas_use = gwas[gwas$pValue<1e-8, ]
gwas_use$location_mod = paste0("chr", gsub(":", "_", gwas_use$locations))

eqtl$location_mod = gsub("_[ATCG].*", "", eqtl$variant_id)
eqtl$gene_id_mod = gsub("\\..*", "", eqtl$gene_id)
eqtl$gene_name = ensembl[eqtl$gene_id_mod, ]$Gene.name


eqtl_gene = eqtl[eqtl$location_mod %in% gwas_use$location_mod, ]$gene_id_mod
eqtl_gene_r = r2h[r2h$Human.gene.stable.ID %in% eqtl_gene, ]$Gene.name


### load snRNA-seq data
deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep='\t', header=T)
deg_strain[deg_strain$treatment=="SHR-10w",]$treatment = "SHR-WKY"
deg_strain[deg_strain$treatment=="SS-LS",]$treatment = "SS-SD"


deg_strain[deg_strain$p_val_adj<0.05 & abs(deg_strain$avg_log2FC)>0.5 & deg_strain$gene_name %in% eqtl_gene_r, ]


### LK-eqtl
# p_val avg_log2FC pct.1 pct.2 p_val_adj gene_name           cell_type pct.diff        project         strain tissue control treatment control_size treatment_size test_gene
# 4559   2.7e-113      -1.10 0.288 0.450  5.5e-109     Ppm1e   Inhibitory neuron    -0.16    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY         6691           9053      3996
# 16538   1.2e-12      -2.81 0.053 0.321   2.5e-08     Ppm1e Pars tuberalis cell    -0.27    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY          266            106      1646
# 29754  2.3e-102       0.79 0.609 0.352   5.9e-98   Slc17a3                  PT     0.26    Spontaneous    Spontaneous     LK WKY-10w   SHR-WKY         3315           3920      3076
# 104930  1.9e-10       1.69 0.144 0.018   4.9e-06   Slc17a3                  EC     0.13 Salt-sensitive Salt-sensitive     LK   SD-LS     SS-SD          153            491      1426

### HYP-eqtl
# p_val avg_log2FC pct.1 pct.2 p_val_adj gene_name           cell_type pct.diff        project         strain tissue control treatment control_size treatment_size test_gene
# 8636   5.6e-09      -0.59 0.126 0.180   1.1e-04      Rbm6      Myelinating OL   -0.054    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY         2510           2703      1233
# 12107  6.2e-07      -0.62 0.074 0.184   1.3e-02      Rbm6           Microglia   -0.110    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY          585            353      1291
# 13375  3.5e-08      -0.96 0.086 0.179   7.1e-04      Rbm6                 OPC   -0.093    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY          754            716      1549
# 22405  2.7e-10      -0.65 0.242 0.416   5.5e-06     Cadm2            Tanycyte   -0.174    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY          598            681      1619
# 22476  6.9e-07      -0.87 0.082 0.173   1.4e-02      Rbm6            Tanycyte   -0.091    Spontaneous    Spontaneous    HYP WKY-10w   SHR-WKY          598            681      1619
# 27549  3.0e-22       0.95 0.607 0.366   7.6e-18     Cadm2                 DCT    0.241    Spontaneous    Spontaneous     LK WKY-10w   SHR-WKY          631            767      2158
# 40367  2.0e-13       0.55 0.184 0.095   5.2e-09     Cadm2                 TAL    0.089    Spontaneous    Spontaneous     LK WKY-10w   SHR-WKY         1867           1702      2009
# 42049  8.2e-13       0.79 0.255 0.098   2.1e-08     Inpp1                  CD    0.157    Spontaneous    Spontaneous     LK WKY-10w   SHR-WKY          638            664      1910
# 84676  2.8e-12       0.60 0.102 0.021   5.7e-08      Canx           Astrocyte    0.081 Salt-sensitive Salt-sensitive    HYP   SD-LS     SS-SD         2443            749      2839
# 92264  1.6e-07       0.73 0.475 0.217   3.2e-03     Cadm2 Pars tuberalis cell    0.258 Salt-sensitive Salt-sensitive    HYP   SD-LS     SS-SD          236            198      2564
# 104864 1.6e-16       2.00 0.209 0.020   4.0e-12     Cenpp                  EC    0.189 Salt-sensitive Salt-sensitive     LK   SD-LS     SS-SD          153            491      1426
# 118311 1.9e-07       1.57 0.188 0.042   4.9e-03     Cenpp         Macrophages    0.146 Salt-sensitive Salt-sensitive     LK   SD-LS     SS-SD           96            433      1379
# 122127 4.2e-08      -0.56 0.331 0.563   1.1e-03     Cadm2                 DCT   -0.232 Salt-sensitive Salt-sensitive     LK   SD-LS     SS-SD          492            309      1953
# 124051 3.1e-07      -0.63 0.147 0.376   7.9e-03      Rbm6                  CT   -0.229 Salt-sensitive Salt-sensitive     LK   SD-LS     SS-SD          272            226      1962


purified_id = r2h[r2h$Gene.name %in% c("Slc17a3"), ]$Human.gene.stable.ID
apply(eqtl[eqtl$gene_id_mod %in% purified_id & eqtl$location_mod %in% gwas_use$location_mod, c("gene_id", "variant_id")], 1, function(x) paste(x[1], x[2], sep=','))
gwas_use[gwas_use$location_mod %in% eqtl[eqtl$gene_id_mod %in% purified_id, ]$location_mod, ]

# riskAllele pValue pValueAnnotation riskFrequency orValue                 beta          ci mappedGenes                traitName                efoTraits bgTraits  accessionId  locations
# 3819 rs62394271-A  6e-28                -            NR       - 0.3512 unit increase [0.29-0.41]     SLC17A1 Diastolic blood pressure diastolic blood pressure        - GCST90132904 6:25826926
# 7196 rs62394274-T  5e-21                -            NR       - 0.5187 unit decrease [0.41-0.63]     SLC17A3  Systolic blood pressure  systolic blood pressure        - GCST90132903 6:25837627
# pubmedId      author  location_mod
# 3819 35762941 Plotnikov D chr6_25826926
# 7196 35762941 Plotnikov D chr6_25837627


purified_id = r2h[r2h$Gene.name %in% c("Rbm6", "Cadm2"), ]$Human.gene.stable.ID
apply(eqtl[eqtl$gene_id_mod %in% purified_id & eqtl$location_mod %in% gwas_use$location_mod, c("gene_id", "variant_id")], 1, function(x) paste(x[1], x[2], sep=','))
gwas_use[gwas_use$location_mod %in% eqtl[eqtl$gene_id_mod %in% purified_id, ]$location_mod, ]











deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)









# Compile results
results <- data.frame(
  SNP_ID = snp_gr$snpid[snp_ids],
  Gene_ID = gene_gr$gene_id[gene_ids],
  Gene_Name = gene_gr$gene_name[gene_ids],
  Gene_Type = gene_gr$gene_type[gene_ids]
)

SNP_gene = unique(results[, c("Gene_ID", "Gene_Name")])
results_mod = merge(SNP_gene, r2h[], by.x="Gene_ID", by.y="Human.gene.stable.ID", all.x=T)


### enrichment analysis
cell_type_list = unique(deg$cell_type)

de_gene_sum = deg %>% 
  dplyr::group_by(cell_type) %>%
  mutate(
    specificity = case_when(
      p_val_adj.f < 0.05 & !is.na(p_val_adj.f) & avg_log2FC.f > 0.5 ~ "Female-specific",
      p_val_adj.m < 0.05 & !is.na(p_val_adj.m) & avg_log2FC.m > 0.5 ~ "Male-specific",
      p_val_adj.f < 0.05 & p_val_adj.m < 0.05 & abs(avg_log2FC.f) > 0.5 & abs(avg_log2FC.m) > 0.5 & avg_log2FC.f * avg_log2FC.m > 0  & !is.na(p_val_adj.m)  & !is.na(p_val_adj.f) ~ "Overlapping-consistent",
      p_val_adj.f < 0.05 & p_val_adj.m < 0.05 & abs(avg_log2FC.f) > 0.5 & abs(avg_log2FC.m) > 0.5 & avg_log2FC.f * avg_log2FC.m < 0  & !is.na(p_val_adj.m)  & !is.na(p_val_adj.f)  ~ "Overlapping-but-inconsistent",
      TRUE ~ "Not significant"
    )) %>%
  ungroup() %>% as.data.frame()

# de_gene_sum_v2 = de_gene_merged_mod %>% 
#   dplyr::group_by(cell_type) %>%
#   dplyr::summarize(a = sum(p_val_adj.ctrl<0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m>0),
#                    b = sum(p_val_adj.ctrl<0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m<0),
#                    c = sum(p_val_adj.ctrl<0.05 & (p_val_adj.f-0.05)*(p_val_adj.m-0.05)<0),
#                    d = sum(p_val_adj.ctrl>0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m>0),
#                    e = sum(p_val_adj.ctrl>0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m<0),
#                    f = sum(p_val_adj.ctrl>0.05 & (p_val_adj.f-0.05)*(p_val_adj.m-0.05)<0)) %>%
#   ungroup() %>% as.data.frame()

for(i in cell_type_list){
  
  deg_use = de_gene_sum[de_gene_sum$cell_type==i, ]
  
  gwas_genes <- as.character(unique(intersect(results_mod$Gene.name, deg_use$gene_short_name)))
  degs <- as.character(unique(deg_use[deg_use$specificity %in% c("Female-specific", "Male-specific"), ]$gene_short_name))
  degs <- as.character(unique(deg_use[!(deg_use$specificity %in% c("Not significant")), ]$gene_short_name))
  
  all_genes <- as.character(unique(deg_use$gene_short_name))
  
  # Calculate the necessary quantities
  M <- length(all_genes)  # Total genes
  n <- length(gwas_genes)  # Total GWAS genes
  K <- length(intersect(degs, gwas_genes))  # DEGs that are GWAS genes
  N <- M - n  # Non-DEGs
  
  # Perform the hypergeometric test
  p_value <- phyper(K - 1, n, N, length(degs), lower.tail = FALSE)
  
  print(p_value)
  
}








snp_gr <- GRanges(
  seqnames = Rle(ss$chr),
  ranges = IRanges(start = ss$pos, width = 1),  # SNPs are point locations, so width is 1
  snpid = ss$rsid
)

extended_snp_gr <- resize(snp_gr, width = 1000000, fix = "center")
all_overlaps <- findOverlaps(extended_snp_gr, gene_gr)
snp_ids <- queryHits(all_overlaps)
gene_ids <- subjectHits(all_overlaps)












################################################################################
### pathway level GWAS enrichment

### r2h
r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)

### load GWAS data and pre-filter
gwas = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas_catalog_blood pressure.tsv", sep = '\t', header=T)
gwas_bp_mod = gwas %>% filter(pValue < 5e-8) %>%
  separate_rows(mappedGenes, sep=",") %>%
  separate_rows(efoTraits, sep=",") %>%
  filter(efoTraits %in% c("pulse pressure", "diastolic blood pressure", "systolic blood pressure", "hypertension")) %>%
  full_join(r2h, join_by(mappedGenes == Human.gene.name), relationship = "many-to-many") %>%
  mutate(SNP = gsub("-.*", "", riskAllele),
         norm_p = -1*log10(pValue),
         category = "GWAS mapped genes",
         efoTraits = str_to_sentence(efoTraits)) %>%
  filter(Gene.name != "" & !is.na(Gene.name) & !is.na(SNP)) %>%
  dplyr::select(SNP, norm_p, mappedGenes, Gene.name, efoTraits, category) %>%
  unique() %>% as.data.frame()


### load pathway enrichment
merge_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", sep='\t', header = T, quote = "", check.names = F)

traits <- unique(gwas_bp_mod$efoTraits)
for (trait in traits) {
  merge_all[[paste0("gwas_score_", trait)]] <- 0
}

for (i in 1:nrow(merge_all)) {
  pathway <- merge_all[i, ]
  genes <- strsplit(pathway$Symbols, ",")[[1]]
  if (length(genes) <= 10) {
    next
  }
  total_genes <- length(genes)
  
  for (trait in traits) {
    gwas_trait <- gwas_bp_mod[gwas_bp_mod$efoTraits == trait, ]
    mapped_genes <- sum(genes %in% gwas_trait$mappedGenes)
    proportion <- mapped_genes / total_genes
    snp_count <- length(unique(gwas_trait$SNP[gwas_trait$mappedGenes %in% genes]))
    combined_score <- proportion * log1p(snp_count)
    merge_all[i, paste0("gwas_score_", trait)] <- combined_score
  }
}

trait_of_interest <- "Diastolic blood pressure"
p_value <- sum(merge_all[sample(1:nrow(merge_all), 10000), paste0("gwas_score_", trait_of_interest)] >= merge_all[9300, paste0("gwas_score_", trait_of_interest)]) / 10000
p_value
trait_of_interest <- "Systolic blood pressure"
p_value <- sum(merge_all[sample(1:nrow(merge_all), 10000), paste0("gwas_score_", trait_of_interest)] >= merge_all[9300, paste0("gwas_score_", trait_of_interest)]) / 10000
p_value


trait_of_interest <- "Diastolic blood pressure"
p_value <- sum(merge_all[sample(1:nrow(merge_all), 10000), paste0("gwas_score_", trait_of_interest)] >= merge_all[21723, paste0("gwas_score_", trait_of_interest)]) / 10000
p_value
trait_of_interest <- "Systolic blood pressure"
p_value <- sum(merge_all[sample(1:nrow(merge_all), 10000), paste0("gwas_score_", trait_of_interest)] >= merge_all[21723, paste0("gwas_score_", trait_of_interest)]) / 10000
p_value




















