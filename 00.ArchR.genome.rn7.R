ss <- function(x, pattern, slot = 1, ..) {
  sapply(strsplit(x = x, split = pattern, ..), '[', slot)
}
if(FALSE){
  # install the rn7 references
  BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn7")
  BiocManager::install("org.Rn.eg.db")
} 

###########################################
## load in all the required packages ######
suppressPackageStartupMessages(library(ArchR))
library(BSgenome.Rnorvegicus.UCSC.rn7)
library(GenomicFeatures)
library(AnnotationDbi)
library(rtracklayer)

GENOMEDIR='/xdisk/mliang1/qqiu/reference/ArchR/rn7'

###########################################################
###### create the genome annotation files for ArchR #######
# seqlevels(BSgenome.Rnorvegicus.UCSC.rn7) = gsub("chr", "", seqlevels(BSgenome.Rnorvegicus.UCSC.rn7))
genomeAnnotation = createGenomeAnnotation(genome = BSgenome.Rnorvegicus.UCSC.rn7, 
                                          filter = TRUE)
chromSizes = genomeAnnotation$chromSizes
genome(chromSizes) <- "rn7"

# to avoid ends of chromosomes tiled regions
restrict = chromSizes
start(restrict) = start(restrict)+1e6
end(restrict) = end(restrict)-1e6

blacklist = import(file.path(GENOMEDIR,'rn7_liftOver_mm10-blacklist.v2.bed'))
# seqlevels(blacklist) = gsub("chr", "", seqlevels(blacklist))
blacklist = blacklist[!grepl('NW',seqnames(blacklist))]
blacklist <- sort(sortSeqlevels(blacklist), ignore.strand = TRUE)
seqlevels(blacklist) <- seqlevels(chromSizes)
seqlengths(blacklist) = seqlengths(chromSizes)
genome(blacklist) <- "rn7"
blacklist = intersect(blacklist, restrict)

genomeAnnotation$blacklist = blacklist

############################################################
###### load in all the genomes, annotations, for rn7 #######
txdb_sqlite_fn = file.path(GENOMEDIR, 'rn7.sqlite')
if(!file.exists(txdb_sqlite_fn)){
  # load in the gene annotation and save to sqlite 
  txdb = makeTxDbFromGFF(
    "/xdisk/mliang1/qqiu/reference/cellranger-arc/refdata-rat7-2022-sb-r108/Rattus_norvegicus.mRatBN7.2.108.filtered.gtf", 
    organism = 'Rattus norvegicus',
    metadata = c("gene_name")) # for Rattus norvegicus
  seqlevels(txdb)[1:23] <- paste0("chr", seqlevels(txdb)[1:23])
  saveDb(txdb, txdb_sqlite_fn)
}else {
  txdb = loadDb(txdb_sqlite_fn)
}

txdb_rda = file.path(GENOMEDIR, 'rn7.genes_exons_TSS.rda')
if(!file.exists(txdb_rda)){
  # format genes
  library(org.Rn.eg.db)
  genes <- GenomicFeatures::genes(txdb)
  genes = genes[ !duplicated(genes)]
  genes <- keepSeqlevels(genes, 
                         # c(1:20, "X", "Y"),
                         grep('chr',seqlevels(genes), value = T),
                         pruning.mode="coarse")
  
  gtf = import.gff("/xdisk/mliang1/qqiu/reference/cellranger-arc/refdata-rat7-2022-sb-r108/Rattus_norvegicus.mRatBN7.2.108.filtered.gtf")
  
  genes2 = gtf[gtf$type == "gene"]
  genes2 = data.frame(row.names = genes2$gene_id, symbol = genes2$gene_name)
  genes2$symbol[is.na(genes2$symbol)] = rownames(genes2)[is.na(genes2$symbol)]
  
  mcols(genes)$symbol <- genes2[mcols(genes)$gene_id, ]
  genes <- sort(sortSeqlevels(genes), ignore.strand = TRUE)
  seqlevels(genes, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(genes) = end(chromSizes)
  genome(genes) <- "rn7"
  genes = subsetByOverlaps(genes, restrict)
  
  # format exons
  exons <- unlist(GenomicFeatures::exonsBy(txdb, by = "tx"))
  exons <- dropSeqlevels(exons, grep('NW',seqlevels(exons), value = T),
                         pruning.mode="coarse")
  exons$tx_id <- names(exons)
  mcols(exons)$symbol <- genes2[select(txdb, keys = paste0(mcols(exons)$tx_id), 
                                column = "GENEID", keytype = "TXID")[, "GENEID"], ]
  names(exons) <- NULL
  mcols(exons)$exon_id <- NULL
  mcols(exons)$exon_name <- NULL
  mcols(exons)$exon_rank <- NULL
  mcols(exons)$tx_id <- NULL
  exons <- sort(sortSeqlevels(exons), ignore.strand = TRUE)
  exons = exons[!is.na(exons$symbol) & !duplicated(exons) & 
                  exons$symbol %in% genes$symbol]
  seqlevels(exons, pruning.mode="coarse") <- seqlevels(chromSizes)
  seqlengths(exons) = end(chromSizes)
  genome(exons) <- "rn7"
  genes = genes[genes$symbol %in% unique(exons$symbol)]
  # exons = subsetByOverlaps(exons, restrict)
  
  # TSS
  TSS <- resize(genes, 1, "start")
  TSS <- sort(sortSeqlevels(TSS), ignore.strand = TRUE)
  TSS <- dropSeqlevels(TSS, grep('NW',seqlevels(TSS), value = T),
                       pruning.mode="coarse")
  seqlengths(TSS) = end(chromSizes)
  genome(TSS) <- "rn7"
  # TSS = subsetByOverlaps(TSS, restrict)
  
  # save the files
  save(genes, exons, TSS, file = txdb_rda)
} else {
  load(txdb_rda)
}

# make ArchR gene annotation 
geneAnnotation <- createGeneAnnotation(genes = genes, exons = exons, TSS = TSS)

save(genomeAnnotation, geneAnnotation, file = 
       file.path(GENOMEDIR,'rn7.ArchR_annotations.rda'))
