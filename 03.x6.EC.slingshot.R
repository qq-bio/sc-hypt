library(Seurat)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library(dplyr)

set.seed(1)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot")


### 1. load data
i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds"
outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/",
                 gsub("rds", "slingshot.rds", basename(i)))

seurat_object <- readRDS(i)
cluster <- "seurat_clusters"
Idents(seurat_object) <- cluster

sce <- as.SingleCellExperiment(seurat_object)

### 2. using slingshot
sce <- slingshot(sce, clusterLabels = seurat_object$seurat_clusters, reducedDim = 'UMAP', allow.breaks = TRUE, approx_points = 5)
saveRDS(sce, outfile)

summary(sce$slingPseudotime_1)

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')

plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


### 3. temporally dynamic genes
sce <- fitGAM(sce)
ATres <- associationTest(sce)

topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
pst.ord <- order(sce$slingPseudotime_1, na.last = NA)
heatdata <- assays(sce)$counts[topgenes, pst.ord]
heatclus <- sce$GMM[pst.ord]

heatmap(log1p(heatdata), Colv = NA,
        ColSideColors = brewer.pal(9,"Set1")[heatclus])


### 3.1 specific genes
plot_differential_expression <- function(feature_id) {
  feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
  cowplot::plot_grid(plotGeneCount(curves, filt_counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                     plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
}

#### genes change with pseudotime
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)

feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
plot_differential_expression(feature_id)

#### genes are different between lineages
different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)


branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)

feature_id <- branch_point_association %>% filter(pvalue < 0.05) %>% arrange(desc(waldStat)) %>% dplyr::slice(1) %>% pull(feature_id)
plot_differential_expression(feature_id)




#### PAGA
library(SCP)

seurat_object <- RunPAGA(srt = seurat_object, assay_X = "RNA", group_by = "seurat_clusters", linear_reduction = NA, nonlinear_reduction = "UMAP")
CellDimPlot(seurat_object, group.by = "seurat_clusters", reduction = "draw_graph_fr")
PAGAPlot(seurat_object, reduction = "UMAP")
CellDimPlot(seurat_object, group.by = "seurat_clusters", reduction = "UMAP", paga = seurat_object@misc$paga)

seurat_object <- RunPAGA(
  srt = seurat_object, group_by = "seurat_clusters", linear_reduction = NA, nonlinear_reduction = "UMAP",
  embedded_with_PAGA = TRUE, infer_pseudotime = TRUE
)
head(seurat_object[[]])
names(seurat_object@reductions)
FeatureDimPlot(seurat_object, features = "dpt_pseudotime", reduction = "PAGAUMAP2D")
PAGAPlot(seurat_object, reduction = "PAGAUMAP2D")
CellDimPlot(seurat_object, group.by = "seurat_clusters", reduction = "PAGAUMAP2D", paga = seurat_object@misc$paga)


