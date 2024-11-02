
library(GenomicRanges)
library(tidyr)
library(dplyr)


################################################################################
### load reference files
### h2m
r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) <- c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
r2h[r2h$gene_name_rat=="",]$gene_name_rat <- r2h[r2h$gene_name_rat=="",]$gene_id_rat

m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) <- c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")
m2h[m2h$gene_name_mouse=="",]$gene_name_mouse <- m2h[m2h$gene_name_mouse=="",]$gene_id_mouse

### grch38
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

gene_gr <- GRanges(
  seqnames = Rle(ensembl$Chromosome.scaffold.name),
  ranges = IRanges(start = ensembl$Gene.start..bp., end = ensembl$Gene.end..bp.),
  strand = Rle(ensembl$Strand),
  gene_id = ensembl$Gene.stable.ID,
  gene_name = ensembl$Gene.name,
  gene_type = ensembl$Gene.type
)




################################################################################
### Extract SNP Information from GWAS Catalog
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
gwas_folder <- "/xdisk/mliang1/qqiu/data/HT-GWAS/"

# Define the GWAS file lists and output files
gwas_file_groups <- list(
  bp_relevant = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
                  "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv", 
                  "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
                  "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv", 
                  "gwas_catalog_BUN.tsv"),
  bp_only = c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv", 
              "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv")
)
output_files <- c("gwas_catalog_bp_relevant.snp.txt", "gwas_catalog_bp.snp.txt")

# Loop through each group of GWAS files
for (group_name in names(gwas_file_groups)) {
  
  gwas_list <- gwas_file_groups[[group_name]]
  outfile <- output_files[match(group_name, names(gwas_file_groups))]
  
  SNPS <- data.frame()
  
  # Process each GWAS file in the group
  for (gwas_file in gwas_list) {
    
    file_path <- paste0(gwas_folder, gwas_file)
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      next
    }
    
    dat <- read.table(file_path, header = TRUE, sep = '\t', quote = "", fill = TRUE)
    dat <- dat[dat$P.VALUE < 5e-8, ]  # Filter significant p-values
    
    # Extract relevant SNP information
    trait <- gsub("gwas_catalog_", "", gsub(".tsv", "", gwas_file))  # Simplify trait name
    snp_pos <- data.frame(
      SNP = dat$SNPS,
      POS = paste0("chr", dat$CHR_ID, "_", dat$CHR_POS),
      P.VALUE = dat$P.VALUE,
      BETA = dat$OR.or.BETA,
      TRAIT = trait,
      TYPE = dat$CONTEXT
    )
    
    # Filter valid SNPs
    snp_pos <- snp_pos[grepl("^rs", snp_pos$SNP) & !grepl("[,|;|x]|NA", snp_pos$SNP) & !is.na(snp_pos$POS),]
    
    # Keep only the lowest P-value for each SNP
    snp_pos_filtered <- snp_pos %>%
      group_by(SNP) %>%
      filter(P.VALUE == min(P.VALUE)) %>%
      ungroup()
    
    # Combine SNPs
    SNPS <- rbind(SNPS, snp_pos_filtered)
  }
  
  write.table(SNPS, outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')
}





################################################################################
### SNP-Gene Information from GWAS Catalog
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
gwas_folder <- "/xdisk/mliang1/qqiu/data/HT-GWAS/"

gwas_list <- c("gwas_catalog_systolic_bp.tsv", "gwas_catalog_diastolic_bp.tsv",
               "gwas_catalog_pulse_pressure.tsv", "gwas_catalog_essential_hypertension.tsv",
               "gwas_catalog_stroke.tsv", "gwas_catalog_CAD.tsv",
               "gwas_catalog_albuminuria.tsv", "gwas_catalog_eGFR.tsv",
               "gwas_catalog_BUN.tsv")

gwas_merge <- data.frame()

for (gwas_file in gwas_list) {
  
  output <- gsub("tsv", "snp.txt", gwas_file)
  trait <- gsub("(gwas_catalog_|.tsv)", "", gwas_file)
  
  dat <- read.table(file.path(gwas_folder, gwas_file), header = TRUE, sep = '\t', quote = "", fill = TRUE)
  dat <- dat[dat$P.VALUE < 5e-8, ]
  dat$trait <- trait
  
  gwas_merge <- rbind(gwas_merge, unique(dat))
}

# Process SNPs and gene columns by separating multiple entries
gwas_merge <- gwas_merge %>%
  separate_rows(SNPS, sep = ";") %>%
  filter(!grepl(" x ", SNPS)) %>%
  mutate(SNPS = gsub(" ", "", SNPS)) %>%
  separate_rows(MAPPED_GENE, sep = " - ") %>%
  separate_rows(MAPPED_GENE, sep = ", ") %>%
  separate_rows(REPORTED.GENE.S., sep = ", ")

write.table(gwas_merge, "gwas_catalog_bp_relevant.snp_gene.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')





################################################################################
### Identify Proximal Genes to SNPs: Including GWAS Catalog Reported Genes and Nearest Protein Coding Genes
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
gwas_folder <- "/xdisk/mliang1/qqiu/data/HT-GWAS/"

SNPS <- read.table("gwas_catalog_bp_relevant.snp.txt", header = TRUE, sep = '\t')

SNPS <- SNPS %>%
  filter(grepl("^rs", SNP)) %>% 
  separate(POS, into = c("chr", "pos"), sep = "_", convert = TRUE, remove = FALSE) %>%
  mutate(pos = as.numeric(pos)) %>%
  filter(!is.na(pos) & pos != "")

snps_gr <- GRanges(
  seqnames = Rle(SNPS$chr),
  ranges = IRanges(start = SNPS$pos, end = SNPS$pos),
  SNP = SNPS$SNP
)

# Filter Ensembl data for protein-coding genes and construct GRanges
ensembl <- ensembl %>%
  filter(Gene.type == "protein_coding") %>%
  mutate(Chromosome.scaffold.name = paste0("chr", Chromosome.scaffold.name))

genes_gr <- GRanges(
  seqnames = Rle(ensembl$Chromosome.scaffold.name),
  ranges = IRanges(start = ensembl$Gene.start..bp., end = ensembl$Gene.end..bp.),
  gene_id = ensembl$Gene.stable.ID,
  gene_name = ensembl$Gene.name
)

# Identify nearest protein-coding genes to SNPs
nearest_genes <- distanceToNearest(snps_gr, genes_gr)
snp_gene_pairs <- data.frame(
  SNP = mcols(snps_gr)$SNP[queryHits(nearest_genes)],
  Gene_ID = mcols(genes_gr)$gene_id[subjectHits(nearest_genes)],
  Gene_Name = mcols(genes_gr)$gene_name[subjectHits(nearest_genes)]
)

# Save the nearest protein-coding gene information
write.table(snp_gene_pairs, "nearest_protein_coding_merged.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')


# Load GWAS SNP-gene data and extract additional information
snp_gene_pairs <- read.table("nearest_protein_coding_merged.txt", header = TRUE, sep = '\t')
gwas_merge <- read.table("gwas_catalog_bp_relevant.snp_gene.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge <- separate_rows(gwas_merge, SNP_GENE_IDS, sep = ", ")

# Combine and clean SNP-Gene associations from multiple sources
snp_gene <- unique(data.frame(
  SNP = c(gwas_merge$SNPS, gwas_merge$SNPS),
  Gene_Name = c(gwas_merge$MAPPED_GENE, gwas_merge$REPORTED.GENE.S.)
)) %>% filter(Gene_Name != "")

snp_id <- unique(data.frame(
  SNP = c(gwas_merge$SNPS, gwas_merge$SNPS, gwas_merge$SNPS),
  Gene_ID = c(gwas_merge$UPSTREAM_GENE_ID, gwas_merge$DOWNSTREAM_GENE_ID, gwas_merge$SNP_GENE_IDS)
)) %>% filter(Gene_ID != "")

# Merge gene names and IDs with Ensembl data for completeness
snp_gene <- merge(snp_gene, ensembl[, c("Gene.name", "Gene.stable.ID")], by.x = "Gene_Name", by.y = "Gene.name", all.x = TRUE)
colnames(snp_gene)[3] <- "Gene_ID"
snp_id <- merge(snp_id, ensembl[, c("Gene.name", "Gene.stable.ID")], by.x = "Gene_ID", by.y = "Gene.stable.ID", all.x = TRUE)
colnames(snp_id)[3] <- "Gene_Name"

# Combine all sources of SNP-Gene associations and filter out incomplete entries
snp_gene_merge <- unique(
  rbind(snp_gene[c("SNP", "Gene_ID", "Gene_Name")],
        snp_id[c("SNP", "Gene_ID", "Gene_Name")],
        snp_gene_pairs[c("SNP", "Gene_ID", "Gene_Name")])
) %>% filter(!is.na(Gene_ID))

write.table(snp_gene_merge, "gwas_snp_gene_merged.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')








################################################################################
### eQTL Data Processing
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")
SNPS <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header = TRUE, sep = '\t')

# Define the folder and list of eQTL files for processing
eqtl_folder <- "/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/"
eqtl_list <- c("Artery_Aorta.v8.signif_variant_gene_pairs.txt", 
               "Artery_Coronary.v8.signif_variant_gene_pairs.txt",
               "Artery_Tibial.v8.signif_variant_gene_pairs.txt", 
               "Brain_Hypothalamus.v8.signif_variant_gene_pairs.txt", 
               "Heart_Left_Ventricle.v8.signif_variant_gene_pairs.txt", 
               "Kidney_Cortex.v8.signif_variant_gene_pairs.txt")

eqtl_merged <- data.frame()

for(eqtl_file in eqtl_list){
  
  tissue <- gsub(".v8.*", "", eqtl_file)
  eqtl <- read.table(paste0(eqtl_folder, eqtl_file), header = TRUE, sep = '\t', quote = "", fill = TRUE)
  
  # Extract position and gene information, add tissue identifier
  eqtl <- eqtl %>%
    mutate(POS = sub("_[ATCG].*", "", variant_id),
           gene_id_mod = sub("\\..*", "", gene_id),
           gene_name = ensembl[gene_id_mod, ]$Gene.name,
           tissue = tissue) %>%
    filter(POS %in% SNPS$POS)
  
  # Merge with SNP data and orthology mappings
  eqtl <- eqtl %>%
    left_join(SNPS, by = "POS") %>%
    left_join(r2h %>% select(gene_id, gene_name_rat), by = c("gene_id_mod" = "gene_id")) %>%
    left_join(m2h %>% select(gene_id, gene_name_mouse), by = c("gene_id_mod" = "gene_id"))
  
  eqtl_merged <- bind_rows(eqtl_merged, eqtl)
}

eqtl_merged <- eqtl_merged %>%
  select(-gene_name.y, -gene_name)

write.table(eqtl_merged, "eqtl_merged.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')





################################################################################
### EpiMap Data Processing
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")
SNPS <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header = TRUE, sep = '\t')

# Process SNP positions and create GenomicRanges object
SNPS <- SNPS %>%
  mutate(chr = sub("_.*", "", POS),
         pos = as.numeric(sub(".*_", "", POS))) %>%
  filter(!is.na(pos))

snps_gr <- GRanges(seqnames = SNPS$chr,
                   ranges = IRanges(start = SNPS$pos, end = SNPS$pos),
                   SNP = SNPS$SNP)

# Set EpiMap folder and list files
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/")
epimap_folder <- "/xdisk/mliang1/qqiu/data/EpiMap/"

e2g_list <- c("links_by_group.brain.tsv", 
              "links_by_group.endothelial.tsv", 
              "links_by_group.heart.tsv", 
              "links_by_group.kidney.tsv")

# Load chain file for liftover
chain <- import.chain("/xdisk/mliang1/qqiu/reference/liftover/hg19ToHg38.over.chain")

e2g_merged <- data.frame()

for(e2g_file in e2g_list) {
  
  tissue <- gsub("(links_by_group.|.tsv)", "", e2g_file)
  e2g <- read.table(paste0(epimap_folder, e2g_file), header = TRUE, sep = '\t', quote = "", fill = TRUE)
  
  e2g <- e2g %>%
    mutate(gene_name = ensembl[e2g$gene, ]$Gene.name)  # Add gene names
  
  # Create GenomicRanges for EpiMap data
  e2g_gr <- GRanges(seqnames = e2g$chr,
                    ranges = IRanges(start = e2g$start, end = e2g$end),
                    gene = e2g$gene,
                    name = e2g$name,
                    score = e2g$score,
                    group = e2g$group,
                    gene_name = e2g$gene_name)
  
  # LiftOver to hg38 and unlist
  e2g_gr_hg38 <- liftOver(e2g_gr, chain) %>% unlist()
  
  # Find overlaps with SNPs
  overlaps <- findOverlaps(snps_gr, e2g_gr_hg38)
  
  # Combine SNP and EpiMap data where overlaps are found
  combined_df <- data.frame(
    SNP = mcols(snps_gr)$SNP[queryHits(overlaps)],
    SNP_POS = SNPS$POS[queryHits(overlaps)],
    chr = as.character(seqnames(e2g_gr_hg38[subjectHits(overlaps)])),
    start = start(e2g_gr_hg38[subjectHits(overlaps)]),
    end = end(e2g_gr_hg38[subjectHits(overlaps)]),
    gene = mcols(e2g_gr_hg38)$gene[subjectHits(overlaps)],
    name = mcols(e2g_gr_hg38)$name[subjectHits(overlaps)],
    score = mcols(e2g_gr_hg38)$score[subjectHits(overlaps)],
    group = mcols(e2g_gr_hg38)$group[subjectHits(overlaps)],
    gene_name = mcols(e2g_gr_hg38)$gene_name[subjectHits(overlaps)]
  )
  
  # Merge with orthology mappings for rat and mouse genes
  combined_df <- combined_df %>%
    left_join(r2h %>% select(gene_id, gene_name_rat), by = c("gene" = "gene_id")) %>%
    left_join(m2h %>% select(gene_id, gene_name_mouse), by = c("gene" = "gene_id"))
  
  e2g_merged <- bind_rows(e2g_merged, combined_df)
}

e2g_merged <- e2g_merged %>% select(-gene_name.y, -gene_name)

write.table(e2g_merged, "e2g_merged.txt", col.names = TRUE, row.names = FALSE, quote = FALSE, sep = '\t')






################################################################################
### Merge SNP-Gene from Multiple Evidence Sources
################################################################################
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

# Filter SNPs with "rs" identifiers
proximal_merge <- proximal_merge[grepl("rs", proximal_merge$SNP), ]
eqtl_merge <- eqtl_merge[grepl("rs", eqtl_merge$SNP), ]
e2g_merge <- e2g_merge[grepl("rs", e2g_merge$SNP), ]

# Create SNP-gene evidence dataframes
snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID, proximal = "yes")
snp_gene_eqtl <- eqtl_merge %>%
  select(SNP, gene_id_mod, tissue) %>%
  group_by(SNP, gene_id_mod) %>%
  summarise(expressional = paste(unique(tissue), collapse = ", "), .groups = "drop") %>%
  rename(gene = gene_id_mod)
snp_gene_e2g <- e2g_merge %>%
  select(SNP, gene, group) %>%
  group_by(SNP, gene) %>%
  summarise(e2g = paste(unique(group), collapse = ", "), .groups = "drop")

# Merge all evidence dataframes
snp_gene_df <- full_join(snp_gene_prox, snp_gene_eqtl, by = c("SNP", "gene")) %>%
  full_join(., snp_gene_e2g, by = c("SNP", "gene")) %>%
  rowwise() %>%
  mutate(
    proximal = ifelse(is.na(proximal), "no", "yes"),
    expressional = ifelse(is.na(expressional), "no", expressional),
    e2g = ifelse(is.na(e2g), "no", e2g),
    evidence_summary = paste(na.omit(c(
      ifelse(proximal == "yes", "proximal", NA), 
      ifelse(expressional != "no", "expressional", NA), 
      ifelse(e2g != "no", "e2g", NA)
    )), collapse = ", ")
  )

# Aggregate rat and mouse gene names for each gene_id
r2h_agg <- r2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_rat = paste(unique(gene_name_rat), collapse = ";"), .groups = "drop")
m2h_agg <- m2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_mouse = paste(unique(gene_name_mouse), collapse = ";"), .groups = "drop")

# Merge with additional information and aggregate gene names
snp_gene_df <- snp_gene_df %>%
  mutate(pair = paste(SNP, gene, sep = "-")) %>%
  left_join(ensembl, by = c("gene" = "Gene.stable.ID")) %>%
  left_join(r2h_agg, by = c("gene" = "gene_id")) %>%
  left_join(m2h_agg, by = c("gene" = "gene_id")) %>%
  rowwise() %>%
  mutate(
    combined_gene_name = coalesce(gene_name_rat, gene_name_mouse, Gene.name)
  ) %>%
  select(SNP, gene, Gene.name, combined_gene_name, proximal, expressional, e2g, evidence_summary) %>%
  unique()

write.table(snp_gene_df, "snp_gene.evi_org.out", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)





################################################################################
### Venn Diagram and Enrichment Analysis (Figure 4a)
################################################################################
library(ggVennDiagram)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")
SNPS <- read.table("gwas_catalog_bp_relevant.snp.txt", header = TRUE, sep = '\t')

# Load SNP-gene relationship data for different categories
proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t')
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

# Filter valid SNP entries with "rs" prefix
proximal_merge <- proximal_merge[grepl("^rs", proximal_merge$SNP),]
eqtl_merge <- eqtl_merge[grepl("^rs", eqtl_merge$SNP),]
e2g_merge <- e2g_merge[grepl("^rs", e2g_merge$SNP),]

# Unique SNP-gene pairs for each category
snp_gene_prox <- unique(na.omit(paste(proximal_merge$SNP, proximal_merge$Gene_ID, sep = "-")))
snp_gene_eqtl <- unique(na.omit(paste(eqtl_merge$SNP, eqtl_merge$gene_id_mod, sep = "-")))
snp_gene_e2g <- unique(na.omit(paste(e2g_merge$SNP, e2g_merge$gene, sep = "-")))

# Prepare list for Venn Diagram
snp_gene_list <- list("Proximal" = snp_gene_prox, 
                      "Expressional" = snp_gene_eqtl,
                      "Regulatory" = snp_gene_e2g)

ggVennDiagram(snp_gene_list, label = "both", label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "red", name = "Number of\nSNP-gene\npairs") +
  theme(legend.position = "right")




