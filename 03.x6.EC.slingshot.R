library(Seurat)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library(dplyr)

set.seed(1)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot")


################################################################################
### 1. load data
i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/",
                 gsub("rds", "slingshot.LV.rds", basename(i)))

seurat_object <- readRDS(i)
cluster <- "seurat_clusters"
Idents(seurat_object) <- cluster

seurat_object <- subset(seurat_object, seurat_clusters %in% c("C1", "C7", "C9", "M0610", "M24", "C18", "M5813"))

sce <- as.SingleCellExperiment(seurat_object)

### 2. using slingshot
sce <- slingshot(sce, clusterLabels = seurat_object$seurat_clusters, reducedDim = 'UMAP', allow.breaks = TRUE)
saveRDS(sce, outfile)



################################################################################
sce <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.rds")

### three branches in total
summary(sce$slingPseudotime_1) # NA in _1 means more than one lineage, weights of pseudotime are accessible via slingCurveWeights
summary(sce$slingPseudotime_2)
summary(sce$slingPseudotime_3)


### use weighted mean, doesnt work well
sce$slingPseudotime_1[is.na(sce$slingPseudotime_1)] <- 0
sce$slingPseudotime_2[is.na(sce$slingPseudotime_2)] <- 0
sce$slingPseudotime_3[is.na(sce$slingPseudotime_3)] <- 0

sce$slingPseudotime = sce$slingPseudotime_1 * slingCurveWeights(sce)[, 1] +
  sce$slingPseudotime_2 * slingCurveWeights(sce)[, 2] +
  sce$slingPseudotime_3 * slingCurveWeights(sce)[, 3]

colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime, breaks=100)]
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')


plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plotcol[is.na(sce$slingPseudotime_1)] <- "gray"
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')
lines(slingCurves(sce)[[1]], lwd=2, col = 'red')



plot_slingshot_lineage <- function(sce, lineage_num, color_palette = 'PiYG', breaks = 100, 
                                   arrow_col = 'red', line_col = 'black', legend_thickness = 0.1) {
  library(fields)
  library(RColorBrewer)
  
  layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1)) 
  
  colors <- colorRampPalette(rev(brewer.pal(11, color_palette)[-6]))(breaks)
  pseudotime_col <- paste0('slingPseudotime_', lineage_num)
  
  plotcol <- colors[cut(sce[[pseudotime_col]], breaks = breaks)]
  plotcol[is.na(sce[[pseudotime_col]])] <- 'gray'
  
  # Plot main panel
  par(mar = c(5, 5, 4, 2))
  plot(reducedDims(sce)$UMAP, col = plotcol, pch = 16, asp = 1,
       xlab = 'UMAP 1', ylab = 'UMAP 2')
  lines(SlingshotDataSet(sce), lwd = 2, type = 'lineages', col = line_col)
  lines(slingCurves(sce)[[lineage_num]], lwd = 3, col = arrow_col)
  
  # Add arrow at the end of the curve
  curve_points <- slingCurves(sce)[[lineage_num]]$s
  end_idx <- nrow(curve_points)
  arrows(
    x0 = curve_points[end_idx - 1, 1], y0 = curve_points[end_idx - 1, 2],
    x1 = curve_points[end_idx, 1], y1 = curve_points[end_idx, 2],
    length = 0.2, angle = 25, col = arrow_col, lwd = 3
  )
  
  # Plot legend in second panel
  par(mar = c(5, 2, 4, 2))
  plot.new()
  fields::image.plot(
    legend.only = TRUE,
    zlim = range(sce[[pseudotime_col]], na.rm = TRUE),
    col = colors,
    horizontal = FALSE,
    legend.args = list(text = 'Pseudotime', side = 3, line = 1),
    smallplot = c(0.2, 0.2 + legend_thickness, 0.3, 0.7)
  )
}

png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.lineage2.png", width = 463/96, height = 374/96, units = "in", res = 300)
plot_slingshot_lineage(sce, 2)
dev.off()

png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.lineage3.png", width = 463/96, height = 374/96, units = "in", res = 300)
plot_slingshot_lineage(sce, 3)
dev.off()




################################################################################
### 3. temporally dynamic genes
sce <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.rds")

gam_l2 <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.l2.rds")
pseudotime_association_l2 <- associationTest(gam_l2)

gam_l3 <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.l3.rds")
pseudotime_association_l3 <- associationTest(gam_l3)

# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]


### 3.1 specific genes
plot_gene_pseudotime <- function(sce, gene, lineage = 1, col = "dodgerblue", remove_outliers = TRUE) {
  pt <- slingPseudotime(sce)[, lineage]
  expr <- log1p(assay(sce, "counts")[gene, ])  # or "logcounts" if you prefer
  valid <- !is.na(pt)
  
  if (remove_outliers) {
    pt_valid <- pt[valid]
    q1 <- quantile(pt_valid, 0.25)
    q3 <- quantile(pt_valid, 0.75)
    iqr <- q3 - q1
    upper <- q3 + 1.5 * iqr
    lower <- q1 - 1.5 * iqr
    keep <- (pt_valid >= lower) & (pt_valid <= upper)
    removed_ratio <- 1 - sum(keep) / length(keep)
    message(sprintf("Outlier removal proportion: %.2f%%", removed_ratio * 100))
    expr <- expr[valid][keep]
    pt <- pt[valid][keep]
  } else {
    expr <- expr[valid]
    pt <- pt[valid]
  }
  
  smoothed <- loess(expr[valid] ~ pt[valid])
  pred <- predict(smoothed, se = TRUE)
  x_sorted <- sort(pt[valid])
  order_idx <- order(pt[valid])
  fit_sorted <- pred$fit[order_idx]
  se_sorted <- pred$se.fit[order_idx]
  
  ylim_range <- range(fit_sorted + 2 * se_sorted, fit_sorted - 2 * se_sorted, na.rm = TRUE)
  
  plot(NA, xlim = range(pt, na.rm = TRUE), ylim = ylim_range,
       xlab = "Pseudotime",
       ylab = "Expression")
  mtext(gene, side = 3, adj = 0, line = 1, cex = 1.2, font = 2)
  
  # Add confidence interval as shaded region
  polygon(c(x_sorted, rev(x_sorted)),
          c(fit_sorted + 1.96 * se_sorted, rev(fit_sorted - 1.96 * se_sorted)),
          col = adjustcolor(col, alpha.f = 0.2), border = NA)
  
  # Add smoothed line
  lines(x_sorted, fit_sorted, col = col, lwd = 2)
}


# keep_cells <- !is.na(slingPseudotime(sce)[, 2])
# sce_sub <- sce[, keep_cells]
# counts <- as.matrix(sce_sub@assays@data$counts)

plot_gene_pseudotime(sce, 'Il1r1', lineage=3)
plot_gene_pseudotime(sce, 'Myo10', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Kcnt2', lineage=2, remove_outliers = TRUE)








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


