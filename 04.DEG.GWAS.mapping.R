
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")


# SNPS = read.table("gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t') 
gwas_snp = read.table("gwas_merged.txt", header = T, quote = "\"", sep = "\t")


proximal_merge = read.table("gwas_snp_gene_merged.txt", header=T, sep='\t')
eqtl_merge = read.table("eqtl_merged.txt", header=T, sep='\t')
e2g_merge = read.table("e2g_merged.txt", header=T, sep='\t')


length(unique(c(proximal_merge$SNP, eqtl_merge$SNP, e2g_merge$SNP)))
# [1] 10293
length(unique(c(proximal_merge$Gene_ID, eqtl_merge$gene_id_mod, e2g_merge$gene)))
# [1] 10070

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

length(unique(c(proximal_merge$SNP, eqtl_merge$SNP, e2g_merge$SNP)))
# [1] 10253
length(unique(c(proximal_merge$Gene_ID, eqtl_merge$gene_id_mod, e2g_merge$gene)))
# [1] 10041


snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")

length(unique(snp_gene_df$SNP)) 
# [1] 10253
length(unique(snp_gene_df$gene))
# [1] 10041



# > snp_gene_df$SNP[!(snp_gene_df$SNP %in% gwas_snp$SNPS)]
# [1] "rs1333048 x rs1333049"  "rs1333048 x rs1333049"  "rs1165668 x rs1165669"  " rs3127599"             " rs7767084"             " rs2048327"             "rs11924705 x rs6789378" " rs3127599"             " rs2048327"            
# [10] " rs7767084"         

# > snp_gene_df[!(snp_gene_df$SNP %in% gwas_snp$SNPS), ]
# SNP            gene Gene.name combined_gene_name proximal expressional e2g evidence_summary
# 1937  rs1333048 x rs1333049 ENSG00000147889    CDKN2A               <NA>      yes           no  no         proximal
# 1946  rs1333048 x rs1333049 ENSG00000147883    CDKN2B             Cdkn2b      yes           no  no         proximal
# 4978  rs1165668 x rs1165669 ENSG00000166598   HSP90B1            Hsp90b1      yes           no  no         proximal
# 5756              rs3127599 ENSG00000198670       LPA              Lpal2      yes           no  no         proximal
# 5758              rs7767084 ENSG00000198670       LPA              Lpal2      yes           no  no         proximal
# 5762              rs2048327 ENSG00000198670       LPA              Lpal2      yes           no  no         proximal
# 7440 rs11924705 x rs6789378 ENSG00000182447     OTOL1              Otol1      yes           no  no         proximal
# 9674              rs3127599 ENSG00000146477   SLC22A3            Slc22a3      yes           no  no         proximal
# 9677              rs2048327 ENSG00000146477   SLC22A3            Slc22a3      yes           no  no         proximal
# 9684              rs7767084 ENSG00000146477   SLC22A3            Slc22a3      yes           no  no         proximal


unmatched_snp_list = snp_gene_df$SNP[!(snp_gene_df$SNP %in% gwas_snp$SNPS)]



# original way to generate snp-gene table based on proximal info
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

snp_gene_merge = snp_gene_merge[!(is.na(snp_gene_merge$Gene_ID)),]


# > snp_gene_pairs[snp_gene_pairs$SNP %in% unmatched_snp_list, ]
# [1] SNP       Gene_ID   Gene_Name
# <0 rows> (or 0-length row.names)

# > gwas_merge[gwas_merge$SNPS %in% unmatched_snp_list, ]
# # A tibble: 0 × 39
# # ℹ 39 variables: DATE.ADDED.TO.CATALOG <chr>, PUBMEDID <int>, FIRST.AUTHOR <chr>, DATE <chr>, JOURNAL <chr>, LINK <chr>, STUDY <chr>, DISEASE.TRAIT <chr>, INITIAL.SAMPLE.SIZE <chr>, REPLICATION.SAMPLE.SIZE <chr>, REGION <chr>,
# #   CHR_ID <chr>, CHR_POS <chr>, REPORTED.GENE.S. <chr>, MAPPED_GENE <chr>, UPSTREAM_GENE_ID <chr>, DOWNSTREAM_GENE_ID <chr>, SNP_GENE_IDS <chr>, UPSTREAM_GENE_DISTANCE <int>, DOWNSTREAM_GENE_DISTANCE <int>,
# #   STRONGEST.SNP.RISK.ALLELE <chr>, SNPS <chr>, MERGED <int>, SNP_ID_CURRENT <int>, CONTEXT <chr>, INTERGENIC <int>, RISK.ALLELE.FREQUENCY <chr>, P.VALUE <dbl>, PVALUE_MLOG <dbl>, P.VALUE..TEXT. <chr>, OR.or.BETA <dbl>,
# #   X95..CI..TEXT. <chr>, PLATFORM..SNPS.PASSING.QC. <chr>, CNV <chr>, MAPPED_TRAIT <chr>, MAPPED_TRAIT_URI <chr>, STUDY.ACCESSION <chr>, GENOTYPING.TECHNOLOGY <chr>, trait <chr>

# > snp_gene_merge[snp_gene_merge$SNP %in% unmatched_snp_list, ]
# [1] SNP       Gene_ID   Gene_Name
# <0 rows> (or 0-length row.names)






################################################################################
# assumption 1: proximal dataframe was generated before the update of gwas data
# where "x" and " " were removed

gwas_folder = "/xdisk/mliang1/qqiu/data/HT-GWAS/"
gwas_list = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv",
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv",
              "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
              "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv",
              "gwas_catalog_BUN.tsv")

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
# gwas_merge = gwas_merge[! grepl(" x ", gwas_merge$SNPS), ]
# gwas_merge$SNPS = gsub(" ", "", gwas_merge$SNPS)
gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = " - ")
gwas_merge = separate_rows(gwas_merge, MAPPED_GENE, sep = ", ")
gwas_merge = separate_rows(gwas_merge, REPORTED.GENE.S., sep = ", ")
# write.table(gwas_merge, "gwas_merged.txt", col.names = T, row.names = F, quote = F, sep = '\t')


################################################################################
# based on assumption 1, what if we handle "x" and " " in the propper way, how the
# numbers gonna change

snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
snp_gene_df = snp_gene_df[! grepl(" x ", snp_gene_df$SNP), ]
snp_gene_df$SNP = gsub(" ", "", snp_gene_df$SNP)

length(unique(snp_gene_df$SNP)) 
# [1] 10247
length(unique(snp_gene_df$gene))
# [1] 10040


################################################################################
# but the most appropreate way is to re-generate the evi_org table
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
# > length(unique(gwas_merge_use$SNPS))
# [1] 10325
gwas_merge_use <- gwas_merge %>%
  dplyr::filter(grepl("rs", SNPS)) %>%
  dplyr::filter(!(grepl("x", SNPS))) %>%
  mutate(SNPS = gsub(" ", "", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  dplyr::select(SNPS, trait, P.VALUE)

proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
proximal_merge = proximal_merge[! grepl(" x ", proximal_merge$SNP), ]
proximal_merge$SNP = gsub(" ", "", proximal_merge$SNP)

eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]
# > length(unique(c(proximal_merge$SNP, eqtl_merge$SNP, e2g_merge$SNP)))
# [1] 10247
# > length(unique(c(proximal_merge$Gene_ID, eqtl_merge$gene_id_mod, e2g_merge$gene)))
# [1] 10040

library(ggVennDiagram)

snp_gene_prox = unique(proximal_merge$Gene_ID)
snp_gene_eqtl = unique(eqtl_merge$gene_id_mod)
snp_gene_e2g = unique(e2g_merge$gene)

x = list("Proximal"=snp_gene_prox, 
         "Expressional"=snp_gene_eqtl,
         "Regulatory"=snp_gene_e2g)

ggVennDiagram(x, category.names = "", label="both", label_alpha = 0) +
  scale_fill_gradient("Number of genes", low="white",high = "red")




snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID, proximal = "yes")
snp_gene_eqtl <- eqtl_merge %>% dplyr::select(SNP, gene_id_mod, tissue) %>%
  group_by(SNP, gene_id_mod) %>%
  summarise(expressional = paste(unique(tissue), collapse = ", ")) %>%
  rename(gene = gene_id_mod)
snp_gene_e2g <- e2g_merge %>% dplyr::select(SNP, gene, group) %>%
  group_by(SNP, gene) %>%
  summarise(e2g = paste(unique(group), collapse = ", "))
snp_gene_df <- full_join(snp_gene_prox, snp_gene_eqtl, by = c("SNP", "gene")) %>%
  full_join(., snp_gene_e2g, by = c("SNP", "gene"))
snp_gene_df <- snp_gene_df %>%
  rowwise() %>%
  mutate(proximal = ifelse(is.na(proximal), "no", "yes"),
         expressional = ifelse(is.na(expressional), "no", expressional),
         e2g = ifelse(is.na(e2g), "no", e2g),
         evidence_summary = paste(na.omit(c(ifelse(proximal == "yes", "proximal", NA), 
                                            ifelse(expressional != "no", "expressional", NA), 
                                            ifelse(e2g != "no", "e2g", NA))), collapse = ", "))
# > length(unique(snp_gene_df$SNP))
# [1] 10247
# > length(unique(snp_gene_df$gene))
# [1] 10040
r2h_agg <- r2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_rat = paste(unique(gene_name_rat), collapse = ";"))

m2h_agg <- m2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_mouse = paste(unique(gene_name_mouse), collapse = ";"))

snp_gene_df <- snp_gene_df %>%
  mutate(pair = paste(SNP, gene, sep = "-")) %>%
  left_join(ensembl[, c("Gene.stable.ID", "Gene.name")], by = c("gene" = "Gene.stable.ID")) %>%
  left_join(r2h_agg, by = c("gene" = "gene_id")) %>%
  left_join(m2h_agg, by = c("gene" = "gene_id")) %>%
  rowwise() %>%
  mutate(combined_gene_name = case_when(
    is.na(gene_name_mouse) ~ gene_name_rat,
    is.na(gene_name_rat) ~ gene_name_mouse,
    TRUE ~ paste(
      unique(unlist(strsplit(paste(gene_name_mouse, gene_name_rat, sep = ";"), ";"))),
      collapse = ";"
    )
  )) %>%
  dplyr::select(SNP, gene, Gene.name, combined_gene_name, proximal, expressional, e2g, evidence_summary) %>%
  unique()

# > length(unique(snp_gene_df$SNP))
# [1] 10247
# > length(unique(snp_gene_df$gene))
# [1] 10040


snp_gene_df[!(snp_gene_df$SNP %in% gwas_merge$SNPS), ]
# A tibble: 0 × 8
# Rowwise: 
# ℹ 8 variables: SNP <chr>, gene <chr>, Gene.name <chr>, combined_gene_name <chr>, proximal <chr>, expressional <chr>, e2g <chr>, evidence_summary <chr>

gwas_merge = gwas_merge[gwas_merge$SNPS %in% snp_gene_df$SNP, ]

gwas_snp_reform = gwas_merge %>% 
  dplyr::select(SNPS, trait, PUBMEDID, STUDY.ACCESSION) %>%
  group_by(SNPS) %>%
  summarise(
    TRAITS = paste(unique(trait), collapse = ", "),
    PUBMEDID = paste(unique(PUBMEDID), collapse = ", "),
    STUDY.ACCESSION = paste(unique(STUDY.ACCESSION), collapse = ", "),
    .groups = "drop"
  )

gwas_snp_gene = merge(gwas_snp_reform, bp_snp_gene, by.x="SNPS", by.y="SNP")

colnames(gwas_snp_gene) = c("SNP", "Trait", "PubMed_ID", "GWAS_Catalog_study_accession", "Gene_ID_Human", "Gene_symbol_Human",
                            "Gene_symbol_Mouse_Rat", "Proximal_evidence", "Expression_evidence", "Regulatory_evidence", "Evidence_summary")

# > length(unique(gwas_snp_gene$SNP))
# [1] 10247
# > length(unique(gwas_snp_gene$Gene_ID_Human))
# [1] 10040

gwas_snp_gene[is.na(gwas_snp_gene)] <- ""

write.table(gwas_snp_gene, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_snp_gene.summary.out", col.names = T, row.names = F, sep = "\t", quote = F)



length(unique(gwas_snp_gene[grepl("expressional", gwas_snp_gene$Evidence_summary), ]$Gene_ID_Human))
# 3624
length(unique(gwas_snp_gene[grepl("proximal", gwas_snp_gene$Evidence_summary), ]$Gene_ID_Human))
# 6814
length(unique(gwas_snp_gene[grepl("e2g", gwas_snp_gene$Evidence_summary), ]$Gene_ID_Human))
# 2848



