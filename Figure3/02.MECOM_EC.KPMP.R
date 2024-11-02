library(Seurat)
library(sctransform)
library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)


################################################################################
### MECOM+ EC Identification in KPMP Dataset (Figure 3f)
################################################################################
setwd("/xdisk/mliang1/qqiu/data/KPMP/cluster")

seurat_object_orig <- readRDS("/xdisk/mliang1/qqiu/data/KPMP/cluster/KPMP.sc.query.l3.rds")

DefaultAssay(seurat_object_orig) <- "SCT"

# Subset only endothelial cell types
ec_list <- unique(seurat_object_orig$class.l3)[grepl("^EC", unique(seurat_object_orig$class.l3))]
seurat_object <- subset(seurat_object_orig, class.l3 %in% ec_list)
rm(seurat_object_orig)

# Process data
seurat_object <- FindVariableFeatures(seurat_object) %>%
  ScaleData() %>%
  RunPCA()

ElbowPlot(seurat_object, ndims = 50, reduction = "pca") + labs(x="PCA")

# UMAP and clustering
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:30) %>%
  FindNeighbors(dims = 1:30) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.5, 0.8, 1, 2, 4, 5, 10))

## Select resolution for visualization
assay <- "SCT"
for(resolution in c(0.5)) {
  reso_para <- paste0(assay, '_snn_res.', resolution)
  seurat_object$seurat_clusters <- seurat_object@meta.data[[reso_para]]
  
  p <- DimPlot(seurat_object, reduction = "umap", group.by = reso_para, pt.size = 0.2, raster = FALSE, label = TRUE) + 
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0("Resolution=", resolution), color = "Cluster") +
    theme(legend.position = "right")
  print(p)
}

saveRDS(seurat_object, "KPMP.sc.EC.cluster.rds")



# Load and analyze EC subset for MECOM+ expression visualization
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/metabolism-kidney/cluster/KPMP.sc.EC.cluster.rds")

# Set active identity to the desired clustering resolution
Idents(seurat_object) <- "SCT_snn_res.0.5"
seurat_object$cluster <- seurat_object$SCT_snn_res.0.5

# Find markers for each cluster
markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25) %>%
  dplyr::mutate(pct.diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff), .by_group = TRUE)

write.table(markers, file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC.allmarker.0.25.long.txt", sep = "\t")

# Save markers in wide format
marker_tbl <- markers %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
  group_by(cluster) %>%
  arrange(desc(pct.diff), .by_group = TRUE) %>%
  mutate(id = row_number()) %>%
  select(id, cluster, gene) %>%
  pivot_wider(names_from = cluster, values_from = gene, names_prefix = "Cluster ")

write.table(marker_tbl, file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC.allmarker.0.25.wide.txt", sep = "\t", row.names = FALSE)

# Plot MECOM expression on UMAP
FeaturePlot(seurat_object, "MECOM") + 
  labs(x = "UMAP 1", y = "UMAP 2", title = "MECOM Expression in Human Kidney EC") +
  theme(plot.title = element_text(face = "plain"))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3f.kpmp.ec.mecom.png", width = 338 / 96, height = 285 / 96, dpi = 300)







################################################################################
### Enriched Pathway by Genes Significantly Higher Expressed in MECOM+ EC (Figure 3g)
################################################################################

pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/KPMP.sc.EC_Mecom.metascape/KPMP.sc.EC.metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

# Plot -log10(q-value) for pathways
ggplot(pathway_res, aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 16), oob = scales::squish)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3g.kpmp.mecom_ec.pathway_barplot.png", width = 706 / 96, height = 285 / 96, dpi = 300)





