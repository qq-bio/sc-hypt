ec_order <- c("C15","C19","C12","C23","C21",
              "C20","C14",
              "M0610","M24","M5813","C1","C7","C9","C18","C22","C16",
              "C3","C11","C17")

ec_colors <- c(
  "C15" = "#D73027",
  "C19" = "#F46D43",
  "C12" = "#FDAE61",
  "C23" = "#FDD49E",
  "C21" = "#FEC44F",
  
  "C20" = "#4575B4",
  "C14" = "#74ADD1",
  
  "M0610" = "#66C2A5",
  "M24"   = "#41B6C4",
  "M5813" = "#1D91C0",
  "C1"    = "#A6DBA0",
  "C7"    = "#74C476",
  "C9"    = "#31A354",
  "C18"   = "#006D2C",
  "C22"   = "#B8E186",
  "C16"   = "#7FC97F",
  
  "C3"  = "#B3B3B3",
  "C11" = "#9E9AC8",
  "C17" = "#807DBA"
)

blank_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)

################################################################################
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
# library(ggsci)
library(ggh4x)
library(ggtext)
library(patchwork)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))




################################################################################


i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
seurat_object <- readRDS(i)

### rename idx
cell_by_size <- table(seurat_object$seurat_clusters)
cell_by_size <- cell_by_size[order(cell_by_size, decreasing = T)]
new_idx <- as.character(1:19)
names(new_idx) <- names(cell_by_size)

new_idx <- factor(new_idx, levels = new_idx[ec_order])
seurat_object$new_idx <- as.character(new_idx[as.character(seurat_object$seurat_clusters)])
seurat_object$new_idx <- factor(seurat_object$new_idx, levels = new_idx[ec_order])

### umap plot
names(ec_colors) <- as.character(new_idx[ec_order])
p <- DimPlot(seurat_object, group.by = "new_idx", reduction = "umap", pt.size = 1, label = T, repel = F) + 
  blank_theme + labs(title = "") +
  scale_color_manual(values = ec_colors)

pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.pdf", width = 342 / 96, height = 359 / 96)
print(p)
dev.off()

# convert -density 300 /xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.pdf /xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.png



### dotplot
marker_list = c("Pecam1", "Egfl7", "Vwf", # EC
                "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                "Adgrl3", "Slc38a3", # BBB EC
                "Plvap", # venous EC
                "Rgcc", # capillary EC
                "Npr3",
                "Ccl21", "Prox1", # lymphatic EC
                "Igfbp5" # kidney capillary EC
                # "Ehd3" # "Emcn", "Sost", "Plat", # glomerular EC
)

DotPlot(seurat_object, features = marker_list, group.by = "new_idx") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


dp <- DotPlot(seurat_object, features = marker_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  geom_vline(xintercept = 0.3, color = "grey30", linewidth = 0.5) + 
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),  
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )

y_levels <- levels(new_idx)
y_df <- data.frame(
  id = y_levels,
  x = -0.5,
  y = y_levels,
  label = y_levels,
  col = ec_colors[y_levels]  # manually assign label color
)

dp +
  geom_tile(
    data = y_df,
    aes(x = -1, y = y),   # place tile at x = -1
    fill = y_df$col,
    width = 0.3,
    height = 0.8,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = y_df,
    aes(x = -0.7, y = y, label = label),  # place text to the right of the tile
    hjust = 0,
    size = 3,
    inherit.aes = FALSE
  ) +
  coord_cartesian(xlim = c(-0.8, length(marker_list)))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.marker.dotplot.png", width=600/96, height=359/96, dpi=300)



### tissue enrichment plot
prop_data <- seurat_object@meta.data %>%
  group_by(new_idx, tissue) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(new_idx) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  # Compute expected proportion based on total cell distribution per tissue
  left_join(
    seurat_object@meta.data %>%
      dplyr::count(tissue, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = "tissue"
  ) %>%
  # Compute log(obs/exp) and adjust by variance
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()

ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = tissue, y = new_idx, fill = log_obs_exp)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "purple") +
  labs(x = "", y = "", fill = "Log(obs/exp)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 12),
    axis.text.y = element_blank()
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.tissue.enrich.png", width=260/96, height=359/96, dpi=300)




### dot plot of metabolism related genes
marker_list = c("Mfsd2a", "Cpt1a", "Pfkfb3")
# marker_list = fatty_trans_genes
seurat_object_bbb <- subset(seurat_object, seurat_clusters %in% c("C14", "C20"))
DotPlot(seurat_object, features = marker_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "None"
  ) + coord_flip()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.meta.dotplot.png", width=500/96, height=131/96, dpi=300)


DotPlot(seurat_object, features = marker_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="", size = "Percent\nExpressed", colour ="Average\nExpression") +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right"
  ) + coord_flip()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.meta.dotplot.legend.png", width=662/96, height=336/96, dpi=300)




### number of degs
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")

deg_count <- as.data.frame(table(deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5, c("cell_type", "strain")]))
deg_count$new_idx <- factor(new_idx[as.character(deg_count$cell_type)], levels = levels(new_idx))
deg_count$strain <- factor(deg_count$strain, levels = strain_order)
new_idx_order <- levels(new_idx)
full_df <- expand_grid(new_idx = levels(new_idx),
                       strain = strain_order) %>%
  left_join(deg_count, by = c("new_idx", "strain")) %>%
  mutate(Freq = replace_na(Freq, 0),
         new_idx = factor(new_idx, levels = new_idx_order),
         strain = factor(strain, levels = strain_order))

ggplot(full_df, aes(x = Freq, y = new_idx, color = strain)) +
  geom_point(size=3, alpha=0.6) + 
  scale_color_manual(values = species_col) +
  labs(x = "Number of DEGs", y = "", color = "Strain") +
  theme_classic() +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 12),
    axis.text.y = element_blank()
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.deg.dotplot.png", width=262/96, height=359/96, dpi=300)


### pathway enrichment result of C14
library(openxlsx)
pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/C14/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

# Plot -log10(q-value) for pathways
ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C14.path.barplot.png", width = 650 / 96, height = 252 / 96, dpi = 300)






### cell communication result of C14
cc_df_diff <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.diff.out", sep = "\t", header = T)

cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type=="ECC14") %>% 
  mutate(EC_type = paste("Cluster", new_idx[gsub("EC", "", EC_type)]))

path_list = cc_df_use[abs(cc_df_use$value)>0.1, ] %>% arrange(desc(abs(value)))
path_list = unique(path_list$pathway_name)[1:15]
cc_df_use %>% 
  group_by(EC_type) %>% arrange(value) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  filter(pathway_name %in% path_list) %>%
  ggplot(aes(x = rank, y = value, colour = Strain)) +
  geom_point() +
  geom_line() +
  xlab("Rank") +
  ylab("Differential communication score\n(treatment vs. control)") +
  scale_y_continuous(
    breaks = c(0, 0.1)  # customize to your range
  ) +
  facet_grid2(EC_type~pathway_name) +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C14.cellchat.dotplot.png", width=600/96, height=165/96, dpi=300)






### MSA result
library(tidyr)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")

deg_long <- deg_merged %>%
  filter(control_size >= 10 & treatment_size >= 10) %>%
  mutate(
    group = paste0(gene_name, "-", cell_type, "-", strain, "-", treatment),
    hypertension = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive"),
    significance = ifelse(p_val_adj < 0.05, "Yes", "No"),
    linewidth = ifelse(p_val_adj < 0.05, 1.2, 0.4),
    new_idx_name = new_idx[cell_type]
  ) %>%
  pivot_longer(
    cols = c(pct.1, pct.2),
    names_to = "condition",
    values_to = "pct_expr"
  ) %>%
  mutate(
    condition = recode(condition, "pct.1" = "Control", "pct.2" = "Treatment"),
    label = ifelse(p_val_adj < 0.05 & condition == "Treatment", paste0("Cluster ", new_idx_name, "\n(", strain, ")"), NA),
    log2FC = -1 * avg_log2FC
  )

target_gene = "Klf2"
target_gene = "Klf4"

plot_data <- deg_long %>% filter(gene_name == target_gene) %>% arrange(desc(significance))
label_data <- plot_data %>% filter(!is.na(label))

ggplot(plot_data, aes(x = condition, y = pct_expr, group = group)) +
  geom_line(aes(color = log2FC, linewidth = significance)) +
  geom_point(size = 2) +
  geom_text(
    data = label_data,
    aes(label = label),
    hjust = 1.1,
    vjust = 1, 
    size = 4
  ) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  scale_linewidth_manual(values = c("Yes" = 1.2, "No" = 0.4)) +
  labs(
    x = NULL,
    y = "Percent Expressed",
    color = "log2FC",
    linewidth = "Significant",
    title = target_gene
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12)
    
  ) +
  facet_grid(~hypertension)

# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.Klf2.dotplot.png", width=410/96, height=260/96, dpi=300)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.Klf4.dotplot.png", width=410/96, height=260/96, dpi=300)




### pathway enrichment result of C19
library(openxlsx)
pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/C19/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

# Plot -log10(q-value) for pathways
ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C19.path.barplot.png", width = 650 / 96, height = 420 / 96, dpi = 300)









### pathway enrichment result of C3
library(openxlsx)
pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/C3/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

# Plot -log10(q-value) for pathways
ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C3.path.barplot.png", width = 750 / 96, height = 213 / 96, dpi = 300)



### cell communication result of C3
cc_df_diff <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.diff.out", sep = "\t", header = T)

cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type=="ECC3") %>% 
  mutate(EC_type = paste("Cluster", new_idx[gsub("EC", "", EC_type)]))

path_list = cc_df_use[abs(cc_df_use$value)>0.1, ] %>% arrange(desc(abs(value)))
path_list = unique(path_list$pathway_name)[1:15]
cc_df_use %>% 
  group_by(EC_type) %>% arrange(value) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  filter(pathway_name %in% path_list) %>%
  ggplot(aes(x = rank, y = value, colour = Strain)) +
  geom_point() +
  geom_line() +
  xlab("Rank") +
  ylab("Differential communication score\n(treatment vs. control)") +
  scale_y_continuous(
    breaks = c(0, 0.2)  # customize to your range
  ) +
  facet_grid2(EC_type~pathway_name) +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C3.cellchat.dotplot.png", width=800/96, height=165/96, dpi=300)




