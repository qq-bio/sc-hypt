library(Seurat)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library(dplyr)

set.seed(1)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot")


################################################################################

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 6


sce <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.rds")

keep_cells <- !is.na(slingPseudotime(sce)[, 3])
sce_sub <- sce[, keep_cells]

counts_mat <- assay(sce_sub, "counts")
logcounts_mat <- assay(sce_sub, "logcounts")
seurat_obj <- CreateSeuratObject(counts = counts_mat)
seurat_obj <- SetAssayData(seurat_obj, slot = "data", new.data = logcounts_mat)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
variable_genes <- VariableFeatures(seurat_obj)

counts <- as.matrix(sce_sub@assays@data$counts)
lineage_pseudotime  <- slingPseudotime(sce_sub)[,3]
weight <- slingCurveWeights(sce_sub)[,3]
sce <- fitGAM(counts = counts, 
              pseudotime = lineage_pseudotime,
              cellWeights = weight,
              parallel = T,
              BPPARAM = BPPARAM,
              # conditions = factor(seurat$condition[names(lineage_pseuodtime)]),
              genes = variable_genes)

saveRDS(sce, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.l3.rds")
