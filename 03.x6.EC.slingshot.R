library(Seurat)
library(slingshot)
library(grDevices)
library(RColorBrewer)
library(tradeSeq)
library(dplyr)
library(ggplot2)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

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
### 3. visualize slingshot res
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
                                   arrow_col = 'red', legend_thickness = 0.1, 
                                   remove_outliers = FALSE) {
  library(fields)
  library(RColorBrewer)
  
  pseudotime_col <- paste0('slingPseudotime_', lineage_num)
  pt <- sce[[pseudotime_col]]
  valid <- !is.na(pt)
  
  keep_cells <- rep(TRUE, length(pt))
  if (remove_outliers) {
    pt_valid <- pt[valid]
    q1 <- quantile(pt_valid, 0.25)
    q3 <- quantile(pt_valid, 0.75)
    iqr <- q3 - q1
    upper <- q3 + 1.5 * iqr
    lower <- q1 - 1.5 * iqr
    outlier_mask <- valid & ((pt < lower) | (pt > upper))
    keep_cells[outlier_mask] <- FALSE
    removed_ratio <- sum(outlier_mask) / sum(valid)
    message(sprintf("Outlier removal proportion: %.2f%% (NA values retained)", removed_ratio * 100))
  }
  
  layout(matrix(c(1, 2), ncol = 2), widths = c(4, 1)) 
  
  colors <- colorRampPalette(rev(brewer.pal(11, color_palette)[-6]))(breaks)
  plotcol <- rep("gray", length(pt))
  plotcol[valid & keep_cells] <- colors[cut(pt[valid & keep_cells], breaks = breaks)]
  
  # Plot main panel
  par(mar = c(5, 5, 4, 2))
  plot_coords <- reducedDims(sce)$UMAP[keep_cells, ]
  plot(plot_coords, col = plotcol[keep_cells], pch = 16, asp = 1,
       xlab = 'UMAP 1', ylab = 'UMAP 2')
  
  slingshot_ds <- SlingshotDataSet(sce)
  lines(slingshot_ds, lwd = 2, type = 'lineages', col = "black")
  
  curve_points <- slingCurves(sce)[[lineage_num]]$s
  lines(curve_points, lwd = 3, col = arrow_col)
  
  max_visible_x <- max(plot_coords[, 1], na.rm = TRUE)
  visible_points <- which(curve_points[, 1] <= max_visible_x)
  
  if (length(visible_points) >= 2) {
    end_idx <- tail(visible_points, 1)
    if (end_idx < nrow(curve_points)) {
      arrows(
        x0 = curve_points[end_idx, 1], y0 = curve_points[end_idx, 2],
        x1 = curve_points[end_idx + 1, 1], y1 = curve_points[end_idx + 1, 2],
        length = 0.2, angle = 25, col = arrow_col, lwd = 3
      )
    }
  } else {
    end_idx <- nrow(curve_points)
    arrows(
      x0 = curve_points[end_idx - 1, 1], y0 = curve_points[end_idx - 1, 2],
      x1 = curve_points[end_idx, 1], y1 = curve_points[end_idx, 2],
      length = 0.2, angle = 25, col = arrow_col, lwd = 3
    )
  }
  
  # Plot legend in second panel
  par(mar = c(5, 2, 4, 2))
  plot.new()
  fields::image.plot(
    legend.only = TRUE,
    zlim = range(pt[valid & keep_cells], na.rm = TRUE),
    col = colors,
    horizontal = FALSE,
    legend.args = list(text = 'Pseudotime', side = 3, line = 1),
    smallplot = c(0.2, 0.2 + legend_thickness, 0.3, 0.7)
  )
}

png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.lineage2.png", width = 463/96, height = 374/96, units = "in", res = 300)
plot_slingshot_lineage(sce, 2)
dev.off()

png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.lineage2.rm_outlier.png", width = 463/96, height = 374/96, units = "in", res = 300)
plot_slingshot_lineage(sce, 2, remove_outliers = TRUE)
dev.off()

png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.lineage3.png", width = 463/96, height = 374/96, units = "in", res = 300)
plot_slingshot_lineage(sce, 3)
dev.off()




################################################################################
### 4. temporally dynamic genes (in another code)
sce <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.rds")


### check sample source
tbl = table(colData(sce)$sxtxt, colData(sce)$seurat_clusters)
tbl = tbl[grepl("LV", rownames(tbl)), ]
tbl/rowSums(tbl)*100



gam_l2 <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.l2.rds")
pseudotime_association_l2 <- associationTest(gam_l2)

gam_l3 <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/slingshot/ec.scvi.gene_nb.hvg_1k.refined.merged.slingshot.LV.l3.rds")
pseudotime_association_l3 <- associationTest(gam_l3)

# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]


### 5. visualize specific genes
plot_gene_pseudotime <- function(sce, gene, lineage = 1, col = "dodgerblue", remove_outliers = FALSE) {
  pt <- slingPseudotime(sce)[, lineage]
  expr <- log1p(assay(sce, "counts")[gene, ])
  valid <- !is.na(pt)
  
  if (remove_outliers) {
    pt_valid <- pt[valid]
    q <- quantile(pt_valid, probs = c(0.25, 0.75))
    iqr <- q[2] - q[1]
    keep <- valid & pt >= (q[1] - 1.5 * iqr) & pt <= (q[2] + 1.5 * iqr)
  } else {
    keep <- valid
  }
  
  pt_keep <- pt[keep]
  expr_keep <- expr[keep]
  
  smoothed <- loess(expr_keep ~ pt_keep)
  pred <- predict(smoothed, se = TRUE)
  ord <- order(pt_keep)
  df_plot <- data.frame(
    pt = pt_keep[ord],
    fit = pred$fit[ord],
    upper = pred$fit[ord] + 1.96 * pred$se.fit[ord],
    lower = pred$fit[ord] - 1.96 * pred$se.fit[ord]
  )
  
  ggplot(df_plot, aes(x = pt, y = fit)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = col, alpha = 0.2) +
    geom_line(color = col, size = 1.2) +
    labs(x = "Pseudotime", y = "Expression", title = gene) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0))
}

# width=244&height=179

### some markers
plot_gene_pseudotime(sce, 'Il1r1', lineage=3)
plot_gene_pseudotime(sce, 'Vwf', lineage=3)
plot_gene_pseudotime(sce, 'Kcnt2', lineage=2)
plot_gene_pseudotime(sce, 'St6galnac3', lineage=2)
plot_gene_pseudotime(sce, 'Kcnt2', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'St6galnac3', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Myo10', lineage=2, remove_outliers = TRUE)

### C7: lipoprotein
plot_gene_pseudotime(sce, 'Lpl', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Cdh13', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Cd36', lineage=2, remove_outliers = TRUE)

### C7: interfeson
plot_gene_pseudotime(sce, 'Notch1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Notch1', lineage=2)
plot_gene_pseudotime(sce, 'Notch1', lineage=3)

plot_gene_pseudotime(sce, 'Ifit1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Mx1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Oas2', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Gbp1', lineage=2, remove_outliers = TRUE)

### C7: ECM
plot_gene_pseudotime(sce, 'Col4a1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Lama4', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Sulf1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Hmcn1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Hsp90aa1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Hsp90ab1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Cdc37', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Nrp2', lineage=2, remove_outliers = TRUE)

### C7: blood coagulation
plot_gene_pseudotime(sce, 'F8', lineage=2, remove_outliers = TRUE)

### C7: angiogenesis
plot_gene_pseudotime(sce, 'Cdh13', lineage=2, remove_outliers = TRUE)


### M5813: aerobic respiration
plot_gene_pseudotime(sce, 'Aco2', lineage=3)
plot_gene_pseudotime(sce, 'mt-Atp6', lineage=3)
plot_gene_pseudotime(sce, 'mt-Co1', lineage=3)
plot_gene_pseudotime(sce, 'mt-Co3', lineage=3)
plot_gene_pseudotime(sce, 'mt-Cytb', lineage=3)
plot_gene_pseudotime(sce, 'mt-Nd4', lineage=3)
plot_gene_pseudotime(sce, 'Sdhc', lineage=3)
plot_gene_pseudotime(sce, 'Me3', lineage=3)
plot_gene_pseudotime(sce, 'Atp1a1', lineage=3)
plot_gene_pseudotime(sce, 'Ednra', lineage=3)


plot_gene_pseudotime(sce, 'Ttn', lineage=3)
plot_gene_pseudotime(sce, 'Tpm1', lineage=3)
plot_gene_pseudotime(sce, 'Mapt', lineage=3)

plot_gene_pseudotime(sce, 'Auts2', lineage=3)
plot_gene_pseudotime(sce, 'Ttn', lineage=3)
plot_gene_pseudotime(sce, 'Sema5b', lineage=3)
plot_gene_pseudotime(sce, 'Auts2', lineage=3)
plot_gene_pseudotime(sce, 'Ttn', lineage=3)
plot_gene_pseudotime(sce, 'Sema5b', lineage=3)


### cellchat based
### C7
plot_gene_pseudotime(sce, 'Col4a1', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Col4a2', lineage=2, remove_outliers = TRUE)

### M5813
plot_gene_pseudotime(sce, 'Col4a1', lineage=3)
plot_gene_pseudotime(sce, 'Col4a2', lineage=3)

### M0610
plot_gene_pseudotime(sce, 'Lama4', lineage=2, remove_outliers = TRUE)
plot_gene_pseudotime(sce, 'Plxna4', lineage=2, remove_outliers = TRUE)

plot_gene_pseudotime(sce, 'Lama4', lineage=3)
plot_gene_pseudotime(sce, 'Plxna4', lineage=3)









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


