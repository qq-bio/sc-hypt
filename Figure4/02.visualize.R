library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggh4x)
library(ggtext)
library(patchwork)
library(data.table)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))



################################################################################
### load data
i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
seurat_object <- readRDS(i)


### umap plot
p <- DimPlot(seurat_object, group.by = "new_idx", reduction = "umap", pt.size = 1, label = T, repel = F) + 
  blank_theme + labs(title = "") +
  scale_color_manual(values = ec_colors)

pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.umap.pdf", width = 342 / 96, height = 359 / 96)
print(p)
dev.off()



### dot plot with tile
marker_list <- c("Pecam1", "Egfl7", "Vwf", # EC
                "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                "Adgrl3", "Slc38a3", # BBB EC
                "Plvap", # venous EC
                "Rgcc", # capillary EC
                "Npr3",
                "Ccl21", "Prox1", # lymphatic EC
                "Igfbp5" # kidney capillary EC
)

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
  left_join(
    seurat_object@meta.data %>%
      dplyr::count(tissue, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = "tissue"
  ) %>%
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()

ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = tissue, y = new_idx, fill = log_obs_exp)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low = "white", high = "purple") +
  labs(x = "", y = "", fill = "Log(obs/exp)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 12),
    axis.text.y = element_blank()
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.tissue.enrich.png", width=260/96, height=359/96, dpi=300)





### number of degs
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")

deg_count <- as.data.frame(table(deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5, c("cell_type", "strain")]))
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
  theme_classic()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.deg.dotplot.with_label.png", width=272/96, height=359/96, dpi=300)





### dot plot of selected metabolic genes
marker_list = c("Mfsd2a", "Cpt1a", "Pfkfb3")
seurat_object_bbb <- subset(seurat_object, seurat_clusters %in% c("C14", "C20"))
DotPlot(seurat_object, features = marker_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    legend.position = "None"
  ) + coord_flip()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.meta.dotplot.png", width=500/96, height=131/96, dpi=300)




### lined dot plot of pdgfd, pdgfrb
ligand_expr <- deg_merged[deg_merged$cell_type=="C14" & deg_merged$gene_name=="Pdgfd", ]
receptor_expr <- fread("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out")[cell_type == "Pericyte" & gene_name == "Pdgfrb" & tissue == "HYP"]

ligand_long <- ligand_expr %>%
  filter(control_size >= 10 & treatment_size >= 10,
         strain %in% c("C57BL/6", "SHR", "SS")) %>%
  mutate(
    group = paste0(gene_name, "-", cell_type, "-", strain, "-", treatment),
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
    label = "Pdgfd\n(C10)",
    log2FC = -1 * avg_log2FC
  )

receptor_long <- receptor_expr %>%
  filter(control_size >= 10 & treatment_size >= 10,
         strain %in% c("C57BL/6", "SHR", "SS")) %>%
  mutate(
    group = paste0(gene_name, "-", cell_type, "-", strain, "-", treatment),
    significance = ifelse(p_val_adj < 0.05, "Yes", "No"),
    linewidth = ifelse(p_val_adj < 0.05, 1.2, 0.4)
  ) %>%
  pivot_longer(
    cols = c(pct.1, pct.2),
    names_to = "condition",
    values_to = "pct_expr"
  ) %>%
  mutate(
    condition = recode(condition, "pct.1" = "Control", "pct.2" = "Treatment"),
    label = "Pdgfrb\n(pericyte)",
    log2FC = -1 * avg_log2FC
  )

intersect_col <- intersect(colnames(ligand_long), colnames(receptor_long))
plot_data <- rbind(ligand_long[, intersect_col], receptor_long[, intersect_col])
label_data <- plot_data %>% filter(!is.na(label))

ggplot(plot_data, aes(x = condition, y = pct_expr, group = group)) +
  geom_line(aes(color = log2FC, linewidth = significance)) +
  geom_point(size = 2) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  scale_linewidth_manual(values = c("Yes" = 1.2, "No" = 0.4)) +
  labs(
    x = NULL,
    y = "Percent Expressed",
    color = "log2FC",
    linewidth = "Significant",
    title = "PDGF signaling in\nhypertensive strains"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12)
    
  ) +
  facet_grid(~label)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/PDGF.ligand_receptor.expr.png", width=330/96, height=260/96, dpi=300)





### lined dot plot of Klf2
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

target_gene <- "Klf2"

plot_data <- deg_long %>% filter(gene_name == target_gene) %>% arrange(desc(significance))
label_data <- plot_data %>% filter(!is.na(label))

ggplot(plot_data, aes(x = condition, y = pct_expr, group = group)) +
  geom_line(aes(color = log2FC, linewidth = significance)) +
  geom_point(size = 2) +
  geom_text(
    data = label_data,
    x = -Inf,
    aes(label = label),
    hjust = 0,
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

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.Klf2.dotplot.png", width=410/96, height=260/96, dpi=300)


















