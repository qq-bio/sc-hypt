library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggh4x)
library(ggtext)
library(patchwork)
library(data.table)
library(tidyr)
library(purrr)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))




################################################################################
ec_order <- c(7, 12, 3, 2, 5, 10, 4, 1, 6, 9, 13, 15, 8, 14, 11)

ec_colors <- c(
  "9" = "#D73027",
  "13" = "#FDD49E",
  
  "7"   = "#1D91C0",
  "12" = "#4575B4",
  
  "4" = "#66C2A5",
  "3" = "#74ADD1",
  "1" =   "#41B6C4",
  "2"    = "#A6DBA0",
  "5"    = "#74C476",
  "10"    = "#31A354",
  "14"   = "#006D2C",
  
  "6"   = "#B8E186",
  "8"   = "#7FC97F",
  
  "11"  = "#B3B3B3",
  "15" = "#9E9AC8",
  "12" = "#807DBA"
)

cluster_names <- c(
  "7" = "C7 (Art-Slc8a1)", "12" = "C12 (Art-Tox)", 
  "3" = "C3 (Arteriolar)", "2" = "C2 (Cap-Nrp1)", "5" = "C5 (Cap-Nrp2)",
  "10" = "C10 (Cap-cycling)", "4" = "C4 (Cap-activated)", 
  "1" = "C1 (Cap-Cd36)", "6" = "C6 (Cap-Fabp4)",
  "9" = "C9 (BBB-Tfrc)", "13" = "C13 (BBB-Abcb1b)",
  "15" = "C15 (AV-intermediate)", "8" = "C8 (Venous)", 
  "14" = "C14 (Endocardial)", "11" = "C11 (Lymphatic)"
)




################################################################################
### Cross-organ EC QC and merge clusters based on marker gene overlap
################################################################################
### further remove potential doublets by checking expression of markers of other major cell types
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.rds")
seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 2, 2.5, 3, 3.5, 4))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)

Idents(seurat_object) = "RNA_snn_res.4"

contamination_list = c("Ptprc", "Col1a1", "Dcn", "Rgs5", "Acta2", "Myh11", 
                       "Myl2", "Mb", "Tnnt2", "Actc1", "Myh6", "Myh7",
                       "Syt1", "Rbfox3", "Slc1a2", "Mbp", "Pdgfra", "Cx3cr1", "Pdgfrb", "Rgs5", "Tgfbr2",
                       "Ptprc", "Col1a1", "Dcn", "Rgs5", "Acta2", "Myh11", 
                       "Lrp2", "Slc34a1", "Slc8a1", "Slc12a1", "Slc12a3", "Aqp1", "Atp6v0d2",
                       "Mecom"
)
contamination_list = unique(contamination_list)
DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T, group.by = "seurat_clusters")
DotPlot(seurat_object, features = contamination_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = c("Pecam1", "Mb", "Lrp2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

seurat_object <- subset(seurat_object, RNA_snn_res.4 %in% c(1:59, 62:65))
seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 2, 2.5, 3, 3.5, 4))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)
Idents(seurat_object) = "RNA_snn_res.2"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.2
DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T)

DotPlot(seurat_object, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(seurat_object, features = contamination_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.rds")




### calculate overlap of top markers 
# functions
# helper to get top-N genes per cluster 
get_top_genes <- function(markers, n = 50) {
  stopifnot(all(c("cluster", "gene") %in% colnames(markers)))
  has_rank <- "rank" %in% colnames(markers)
  
  if (has_rank) {
    markers %>%
      group_by(cluster) %>%
      dplyr::arrange(rank, .by_group = TRUE) %>%
      slice_head(n = n) %>%
      dplyr::summarise(genes = list(unique(gene)), .groups = "drop")
  } else {
    # fall back ordering if you don't have rank
    pv <- if ("p_val_adj" %in% colnames(markers)) "p_val_adj" else "p_val"
    if (!pv %in% colnames(markers)) stop("Need `rank` OR `p_val_adj`/`p_val` to order genes.")
    markers %>%
      group_by(cluster) %>%
      dplyr::arrange(.data[[pv]], desc(avg_log2FC), .by_group = TRUE) %>%
      slice_head(n = n) %>%
      summarise(genes = list(unique(gene)), .groups = "drop")
  }
}

# compute pairwise overlap stats for a given N
pairwise_overlap <- function(top_tbl, n) {
  clus <- as.character(top_tbl$cluster)
  sets <- setNames(top_tbl$genes, clus)
  
  pairs <- combn(clus, 2, simplify = FALSE)
  purrr::map_dfr(pairs, function(p) {
    a <- p[[1]]; b <- p[[2]]
    A <- sets[[a]]; B <- sets[[b]]
    inter <- intersect(A, B)
    union_ <- union(A, B)
    tibble(
      top_n = n,
      cluster1 = a,
      cluster2 = b,
      n1 = length(A),
      n2 = length(B),
      overlap_n = length(inter),
      jaccard = if (length(union_) > 0) length(inter) / length(union_) else NA_real_,
      overlap_genes = paste(inter, collapse = ";")
    )
  })
}


seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.rds")

markers = FindAllMarkers(seurat_object, group.by="seurat_clusters", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(avg_log2FC), .by_group=TRUE) %>%
  dplyr::mutate(rank = row_number()) %>% ungroup()

top_list <- list(
  `50`  = get_top_genes(markers, n = 50),
  `100` = get_top_genes(markers, n = 100)
)

overlap_df <- bind_rows(
  pairwise_overlap(top_list$`50`,  50),
  pairwise_overlap(top_list$`100`, 100)
) %>%
  arrange(top_n, desc(overlap_n), desc(jaccard))



seurat_object$seurat_clusters = as.character(seurat_object$RNA_snn_res.1)
seurat_object@meta.data[seurat_object$RNA_snn_res.2==22, ]$seurat_clusters = "12-1"
seurat_object@meta.data[seurat_object$RNA_snn_res.2==25, ]$seurat_clusters = "12-2"
seurat_object@meta.data[seurat_object$seurat_clusters=="12", ]$seurat_clusters = "10"

seurat_object@meta.data[seurat_object$seurat_clusters %in% c("2", "5"), ]$seurat_clusters <- "M1"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("4", "6", "7"), ]$seurat_clusters <- "M2"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("3", "9"), ]$seurat_clusters <- "M3"

### rename idx
cell_by_size <- table(seurat_object$seurat_clusters)
cell_by_size <- cell_by_size[order(cell_by_size, decreasing = T)]
new_idx <- as.character(1:length(cell_by_size))
names(new_idx) <- names(cell_by_size)

new_idx <- factor(new_idx, levels = new_idx)
seurat_object$new_idx <- as.character(new_idx[as.character(seurat_object$seurat_clusters)])
seurat_object$new_idx <- factor(seurat_object$new_idx, levels = new_idx)

DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T, group.by = "new_idx")
DotPlot(seurat_object, features = marker_list, group.by="new_idx") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cluster_order = c(9, 13, 7, 12, 4, 3, 1, 2, 5, 10, 14, 6, 11, 8, 15)
seurat_object$seurat_clusters <- factor(seurat_object$new_idx, levels=cluster_order)
Idents(seurat_object) = "seurat_clusters"


saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.rds")







################################################################################
### Cross-organ EC cluster overview (Figure 3a-b; Fig. S12a-b)
################################################################################

marker_by_cat <- list(
  # large vessels
  C7  = c("Pcsk5", "Slc8a1", "Sulf1"),
  C12 = c("Tox", "Pgm5", "Zfp521"),
  C14 = c("Npr3", "Ltbp1"),
  
  # capillary
  "Micro-Core" = c("Smad6", "Efnb2", "Nrp1", "Dach1"), 
  
  C3  = c("Epas1"),  
  C2  = c("Mgll", "Myo10"),            
  C5  = c("Nrp2", "Nox4"),             
  C10 = c("Mki67", "Top2a"),           
  C4  = c("Vegfc", "Calcrl", "Thsd7a"),
  C1  = c("Flt1", "Cd36"),             
  C6  = c("Fabp4", "Gpihbp1", "Id1"),  
  
  # BBB
  C9  = c("Tfrc"),
  C13 = c("Abcb1b"),
  
  # remaining
  C15 = c("Pcdh17", "Emcn", "Rapgef4"), # Mixex transition
  C8  = c("Exoc3l2", "Inpp4b"),   # Venous-like
  C11 = c("Flt4", "Prox1")                # Lymphatic
)

i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.rds"
seurat_object <- readRDS(i)

# seurat_object$seurat_clusters = factor(seurat_object$seurat_clusters, levels=as.character(ec_order))
# cluster_names_appd = factor(cluster_names[as.character(seurat_object$seurat_clusters)], levels=cluster_names)
# names(cluster_names_appd) = colnames(seurat_object)
# seurat_object$cluster_name = cluster_names_appd


### umap plot (Figure 3a)
p <- DimPlot(seurat_object, group.by = "new_idx", reduction = "umap", pt.size = 1, label = T, repel = F) + 
  blank_theme + labs(title = "") +
  scale_color_manual(values = ec_colors)

pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.pdf", width = 342 / 96, height = 359 / 96)
print(p)
dev.off()



### plot umap by strain (Fig. S12a)
plot_umap_highlight_by <- function(so, group_key = "tissue",
                                   highlight_col = "#2b6cb0", other_col = "grey85",
                                   pt_size = 1, ncol = 4) {
  levs <- unique(so@meta.data[[group_key]])
  plots <- lapply(levs, function(lv) {
    so$..highlight_tmp <- ifelse(so@meta.data[[group_key]] == lv, "highlight", "other")
    p <- DimPlot(so, reduction = "umap", group.by = "..highlight_tmp",
                 cols = c("other"=other_col, "highlight"=highlight_col), pt.size = pt_size, raster = TRUE) +
      ggplot2::ggtitle(lv) + theme(legend.position = "none")
    so$..highlight_tmp <- NULL
    p
  })
  patchwork::wrap_plots(plots, ncol = ncol)
}

p <- plot_umap_highlight_by(seurat_object, group_key = "strain",
                            highlight_col = "#1f77b4", other_col = "grey90",
                            pt_size = 3, ncol = 5) &
  labs(x = "UMAP 1", y = "UMAP 2")
p
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.per_strain.png", width = 1036 / 96, height = 235 / 96, dpi = 300)



### dot plot (Figure 3b)
marker_list = unique(c(unlist(marker_by_cat)))
seurat_object %>%
  DotPlot(., features = marker_list, group.by = "cluster_name") + 
  # scale_x_discrete() +
  scale_y_discrete(limits=rev) +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    text = element_text(family="Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(size = 10), 
    legend.text  = element_text(size = 10)) +
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    colour = guide_colorbar(title = "Average\nExpression")
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.marker.dotplot.png", width=1017/96, height=384/96, dpi=300)





### tissue enrichment plot (obs/exp and composition)
prop_data_tissue <- seurat_object@meta.data %>%
  group_by(seurat_clusters, tissue) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  left_join(
    seurat_object@meta.data %>%
      dplyr::count(tissue, name = "total_cells") %>%
      dplyr::mutate(expected_prop = total_cells / sum(total_cells)),
    by = "tissue"
  ) %>%
  # Compute log(obs/exp)
  dplyr::mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()

p_tissue <- ggplot(prop_data_tissue[prop_data_tissue$log_obs_exp>0, ], aes(x = tissue, y = seurat_clusters, fill = log_obs_exp)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low = "white", high = "purple", limits = c(0,4.5)) +
  scale_y_discrete(labels = cluster_names) +
  labs(x = "Tissue", y = "", fill = "Log(obs/exp)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
    # axis.text.y = element_blank()
  )

p_tissue_composition <- ggplot(prop_data_tissue, aes(x = tissue, y = seurat_clusters, fill = proportion)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low = "white", high = "darkred", limits = c(0,1)) +
  labs(x = "Tissue", y = "", fill = "Within-cluster\ntissue\ncomposition") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_blank()
  )


### strain composition plot
prop_data_strain <- seurat_object@meta.data %>%
  group_by(seurat_clusters, strain) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup()

p_strain <- ggplot(prop_data_strain, aes(x = strain, y = seurat_clusters, fill = proportion)) +
  geom_tile(color="black") +  
  scale_fill_gradient(low = "white", high = "darkgreen", limits = c(0,1)) +
  labs(x = "Strain", y = "", fill = "Within-cluster\nstrain\ncomposition") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_blank()
  )


### number of degs
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.DEG_all.out")

deg_strain <- unique(deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5, c("gene_name", "cell_type", "strain")])
deg_count <- as.data.frame(table(deg_strain[, c("cell_type", "strain")]))
deg_count$strain <- factor(deg_count$strain, levels = strain_order)
full_df <- expand_grid(strain = strain_order) %>%
  left_join(deg_count, by = c("strain")) %>%
  dplyr::mutate(Freq = replace_na(Freq, 0),
                strain = factor(strain, levels = strain_order),
                cell_type = factor(cell_type, levels = ec_order))

p_deg <- ggplot(full_df, aes(x = Freq, y = cell_type, color = strain)) +
  geom_point(size=3, alpha=0.6) + 
  scale_color_manual(values = strain_col) +
  labs(x = "Number of DEGs", y = "", color = "Strain") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial", color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_blank()
  )


### combine panels (Fig. S12b)
leg_tissue_composition <- get_legend(p_tissue_composition + theme(legend.position = "right"))
leg_tissue <- get_legend(p_tissue + theme(legend.position = "right"))
leg_strain <- get_legend(p_strain + theme(legend.position = "right"))
leg_deg    <- get_legend(p_deg    + theme(legend.position = "right"))

p_tissue_comp_noleg <- p_tissue_composition + theme(legend.position = "none")
p_tissue_noleg <- p_tissue + theme(legend.position = "none")
p_strain_noleg <- p_strain + theme(legend.position = "none")
p_deg_noleg    <- p_deg    + theme(legend.position = "none")

main_panel <- plot_grid(
  p_tissue_comp_noleg,
  p_tissue_noleg,
  p_strain_noleg,
  p_deg_noleg,
  ncol = 3,
  rel_widths = c(1, 1, 1)
)

legend_panel <- plot_grid(
  leg_tissue,
  leg_tissue_composition,
  leg_strain,
  leg_deg,
  ncol = 2
)

p_tissue_noleg + p_tissue_comp_noleg + p_strain_noleg + p_deg_noleg + legend_panel + plot_layout(nrow = 1, widths = c(1, 1, 1, 1, 1.4))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.enrichment_deg.combined.png", width=1100/96, height=400/96, dpi=300)






