# libsatlas.so.3

library(biomaRt)
library(dplyr)
library(ggplot2)

mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
mouse2human = unique(mouse2human[mouse2human[,4]==1, c("Human.gene.name", "Gene.name")])


DEG = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.multiomics.DEG.out", header=T)

GWAS = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas-association-downloaded_2024-03-13-accessionId_GCST008058.tsv", header=T, sep='\t')
GWAS = tidyr::separate_rows(GWAS, MAPPED_GENE, sep = ", ")
GWAS = tidyr::separate_rows(GWAS, MAPPED_GENE, sep = " - ")
GWAS = as.data.frame(GWAS)
GWAS_columns = c("DISEASE.TRAIT", "CHR_ID", "CHR_POS", "STRONGEST.SNP.RISK.ALLELE", "SNPS", "RISK.ALLELE.FREQUENCY", "OR.or.BETA", "REPORTED.GENE.S.", "MAPPED_GENE")
GWAS = GWAS[, GWAS_columns]

DEG_merged = merge(DEG, mouse2human, by.x="gene_name", by.y="Gene.name", all.x=T)
DEG_merged = merge(DEG_merged, GWAS, by.x="Human.gene.name", by.y="MAPPED_GENE", all.x=T)


SNP = read.table("/xdisk/mliang1/qqiu/data/atSNP/GCST008058.eGRF.atSNP.out", header=T, sep='\t')
SNP = SNP[SNP$pval_rank_bh<0.05, ]
DEG_merged = merge(DEG_merged, SNP, by.x="SNPS", by.y="snpid", all.x=T)



unique(DEG_merged[! is.na(DEG_merged$REPORTED.GENE.S.), ]$gene_name)


DEG_filter = DEG_merged[! is.na(DEG_merged$REPORTED.GENE.S.), ]
DEG_filter[order(DEG_filter$p_val),]

unique(DEG_filter[order(DEG_filter$p_val),]$gene_name)





################################################################################
GWAS = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas-association-downloaded_2024-03-13-accessionId_GCST008058.tsv", header=T, sep='\t')
eqtl = read.table("/xdisk/mliang1/qqiu/data/GTEx/GTEx_Analysis_v8_eQTL/Kidney_Cortex.v8.signif_variant_gene_pairs.txt", header=T, sep='\t')


mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <-  eqtl$gene_id
gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                  values = genes, mart= mart)

GWAS$vairant = paste0("chr", GWAS$CHR_ID, "_", GWAS$CHR_POS)
eqtl$vairant = gsub("_[ATCG]_[ATCG].*", "", eqtl$variant_id)

left_join(eqtl, gene_names, by = c("gene_id"="ensembl_gene_id"))





################################################################################
library(atSNP)
library(qvalue)

GWAS = read.table("/xdisk/mliang1/qqiu/data/HT-GWAS/gwas-association-downloaded_2024-03-14-HP_0012594.tsv", header=T, sep='\t')
GWAS = GWAS[grepl("rs", GWAS$SNPS), ]

data(encode_library)
snp_info1 <- LoadSNPData(snpids = GWAS$SNPS, 
                         genome.lib ="BSgenome.Hsapiens.UCSC.hg38", 
                         snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)

atsnp.scores <- ComputeMotifScore(encode_motif, snp_info1, ncores = 1)
atsnp.scores$snp.tbl

atsnp.result <- ComputePValues(motif.lib = encode_motif, snp.info = snp_info1,
                               motif.scores = atsnp.scores$motif.scores, ncores = 1, testing.mc=TRUE)

head(atsnp.result[order(atsnp.result$pval_rank), c("snpid", "motif", "pval_ref", "pval_snp", "pval_rank")])

qval_rank = qvalue(atsnp.result$pval_rank, pi0=0.1)$qvalues
atsnp.result = cbind(atsnp.result, qval_rank)
atsnp.result$trans_factor = as.character(lapply(strsplit(atsnp.result$motif, "_"), function (x) x[[1]][1]))
atsnp.result$pval_rank_bh = p.adjust(atsnp.result$pval_rank, "BH")

atsnp.result$func = ""
atsnp.result[atsnp.result$pval_ref<=0.05 & atsnp.result$pval_snp>0.05 & atsnp.result$pval_rank<=0.05, ]$func = "loss"
atsnp.result[atsnp.result$pval_ref>0.05 & atsnp.result$pval_snp<=0.05 & atsnp.result$pval_rank<=0.05, ]$func = "gain"

write.table(atsnp.result, "/xdisk/mliang1/qqiu/data/atSNP/HP_0012594.alb.atSNP.out", 
            col.names = T, row.names = F, sep = '\t', quote = F)



atsnp.result = read.table("/xdisk/mliang1/qqiu/data/atSNP/MONDO_0001134.HTN.atSNP.out", header = T, sep = '\t')

de_tf = unique(DEG[DEG$gene_name %in% dep_tf, ]$gene_name)

snp_tf = unique(mouse2human[mouse2human$Human.gene.name %in% snp_tf, ]$Gene.name)

intersect(de_tf, snp_tf)

snpids = atsnp.result[atsnp.result$pval_rank_bh<0.1 & atsnp.result$trans_factor=="BHLHE40",]$snpid

snp_info1 <- LoadSNPData(snpids = snpids, 
                         genome.lib ="BSgenome.Hsapiens.UCSC.hg38", 
                         snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38", half.window.size = 30, default.par = TRUE, mutation = FALSE)

atsnp.scores <- ComputeMotifScore(encode_motif, snp_info1, ncores = 1)

match.subseq_result <- MatchSubsequence(snp.tbl = atsnp.scores$snp.tbl, motif.scores = atsnp.result, 
                                        motif.lib = encode_motif, 
                                        snpids = "rs7763581", motifs = "BHLHE40_disc2", ncores = 1)
match.subseq_result[c("snpid", "motif", "IUPAC", "ref_match_seq", "snp_match_seq")]

match.seq<-dtMotifMatch(atsnp.scores$snp.tbl, atsnp.scores$motif.scores, snpids="rs11249906", motifs="BHLHE40_disc2", motif.lib = encode_motif)
plotMotifMatch(match.seq,  motif.lib = encode_motif)

match.seq<-dtMotifMatch(atsnp.scores$snp.tbl, atsnp.scores$motif.scores, snpids="rs2282143", motifs="BHLHE40_disc2", motif.lib = encode_motif)
plotMotifMatch(match.seq,  motif.lib = encode_motif)

match.seq<-dtMotifMatch(atsnp.scores$snp.tbl, atsnp.scores$motif.scores, snpids="rs2304615", motifs="BHLHE40_disc2", motif.lib = encode_motif)
plotMotifMatch(match.seq,  motif.lib = encode_motif)

match.seq<-dtMotifMatch(atsnp.scores$snp.tbl, atsnp.scores$motif.scores, snpids="rs7763581", motifs="BHLHE40_disc2", motif.lib = encode_motif)
plotMotifMatch(match.seq,  motif.lib = encode_motif)














################################################################################
### https://bookdown.org/content/b298e479-b1ab-49fa-b83d-a57c2b034d49/ranking.html



selected_paths <- unique(DEG_merged[! is.na(DEG_merged$REPORTED.GENE.S.), ]$gene_name)
path_ranking <- dplyr::arrange(DEG_merged, dplyr::desc(p_val))
path_ranking$path_rank <- dplyr::percent_rank(path_ranking[["p_val"]]) * 100
df_sub <- subset(path_ranking, path_ranking$gene_name %in% selected_paths)

ggplot(data = path_ranking, aes(-log10(p_val), path_rank, label = Human.gene.name)) + 
  geom_point(shape = 19, cex = 2, color = ifelse(is.na(path_ranking$REPORTED.GENE.S.), "grey50", "red")) + 
  geom_text_repel() +
  xlab("Qval") + ylab("Pathway rank") + 
  theme(# panel.border = element_rect(fill = NA), 
    axis.text = element_text(size = 10, color="black"),
    axis.title = element_text(size = 12),
    text = element_text(size = 12)) + 
  coord_flip()
