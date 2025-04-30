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
  theme_classic() # +
  # theme(
  #   # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 12),
  #   axis.text.y = element_blank()
  # )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.deg.dotplot.with_label.png", width=272/96, height=359/96, dpi=300)


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




### test Pdgfd, Pdgfrb
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")
ligand_expr <- deg_merged[deg_merged$cell_type=="C14" & deg_merged$gene_name=="Pdgfd", ]
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T)
receptor_expr <- deg_merged[deg_merged$cell_type=="Pericyte" & deg_merged$gene_name=="Pdgfrb" & deg_merged$tissue=="HYP", ]

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




ligand_long <- ligand_expr %>%
  filter(control_size >= 10 & treatment_size >= 10,
         strain %in% c("SD", "WKY")) %>%
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
    label = "Pdgfd (C10)",
    log2FC = -1 * avg_log2FC
  )

receptor_long <- receptor_expr %>%
  filter(control_size >= 10 & treatment_size >= 10,
         strain %in% c("SD", "WKY")) %>%
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
    label = "Pdgfrb (pericyte)",
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
    title = "PDGF signaling\n(normotensive strains)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12)
    
  ) +
  facet_grid(~label)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/PDGF.ligand_receptor.expr.normo.png", width=450/96, height=260/96, dpi=300)









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
# target_gene = "Klf4"

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
# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.Klf4.dotplot.png", width=410/96, height=260/96, dpi=300)




gene_list <- c("Klf6", "Slc38a2", "Nr4a1", "Il15", "Esam", "Clu", "Emp1", "Adcy4", 
               "Smg6", "Zyx", "Ptprj", "Zeb2", "Ubc", "Map3k3", "Hnrnph1", "Fosb", "Ddx5", "Il3ra")
DotPlot(seurat_object, features = gene_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="", size = "Percent\nExpressed", colour ="Average\nExpression") +
  theme(
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    legend.position = "right"
  ) + coord_flip()



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





### tf-regulon and marker genes in C15
library(data.table)
library(stringr)
regulonsfile <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/scenic/rat_scenic_res/rat_10k10k/reg.csv"
regulons <- fread(regulonsfile, header = FALSE)
regulons <- regulons[-c(1:3), ]
colnames(regulons) <- c("TF","MotifID","AUC","NES","MotifSimilarityQvalue",
                        "OrthologousIdentity","Annotation","Context","TargetGenes","RankAtMax")
numeric_cols <- c("AUC", "NES", "MotifSimilarityQvalue", "OrthologousIdentity", "RankAtMax")
regulons[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]
regulons[, TargetGeneList := str_extract_all(TargetGenes, "'(.*?)'")]
regulons[, TargetGeneList := sapply(TargetGeneList, function(x) paste(gsub("'", "", x), collapse = ";"))]

regulons_long <- regulons %>%
  separate_rows(TargetGeneList, sep = ";") %>%
  dplyr::select(Gene = TargetGeneList, TF, MotifID, AUC, NES)

regulons_use <- regulons_long[regulons_long$TF %in% c("Klf2", "Klf4"), ]

# load C15 marker genes
# marker_list <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.cluster_wise.DEG.out", header = T)
# marker_c19 <- marker_list[marker_list$cluster]

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")
deg_c19 <- deg_merged[deg_merged$cell_type=="C19" & deg_merged$p_val_adj<0.05, ]
deg_target <- unique(deg_c19[deg_c19$gene_name %in% regulons_use$Gene, ]$gene_name)


TF_gene_pathway_all <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/scenic/rat_scenic_res/TF_gene_pathway_all.csv")
klf_gene <- TF_gene_pathway_all[TF_gene_pathway_all$TF %in% c("Klf2", "Klf4") & TF_gene_pathway_all$data == "rat_10k10k", ]
gene_list <- setdiff(unique(klf_gene$Gene), c("Klf2", "Klf4"))



deg_long <- deg_merged %>%
  filter(control_size >= 10 & treatment_size >= 10,
         cell_type=="C19" & strain == "WKY",
         gene_name %in% gene_list) %>%
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
    label = gene_name,
    log2FC = -1 * avg_log2FC
  )

plot_data <- deg_long %>% arrange(desc(significance))

ggplot(plot_data, aes(x = condition, y = pct_expr, group = group)) +
  geom_line(aes(color = log2FC, linewidth = significance)) +
  geom_point(size = 2) +
  # ) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0) +
  scale_linewidth_manual(values = c("Yes" = 1.2, "No" = 0.4)) +
  labs(
    x = NULL,
    y = "Percent Expressed",
    color = "log2FC",
    linewidth = "Significant",
    title = "Targets of Klf2/Klf4 in C15 (WKY)"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text.y = element_text(size = 10, angle = 0)
  ) +
  facet_wrap(~ label, scales = "fixed", ncol=3)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/Klf2_4.target.expr_change.png", width = 416 / 96, height = 593 / 96, dpi = 300)










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







### path enrichment results of LK ECs
library(openxlsx)
cell_type_list = c("C7", "M5813", "M0610", "M24")

pathway_res_merged = c()
for(i in cell_type_list){
  
  path_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/", i, "/metascape_result.xlsx")
  pathway_res <- read.xlsx(path_file, sheet = 2)
  # pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
  pathway_res <- pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ]
  pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
  pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")
  pathway_res$cell_type <- paste0("C", new_idx[i])
  pathway_res_merged <- rbind(pathway_res_merged, pathway_res)
  
}

summary_list <- pathway_res_merged %>% group_by(cell_type) %>%
  filter(grepl("Summary", GroupID)) %>% ungroup() %>% dplyr::select(pathway)


pathway_res_use <- pathway_res_merged[pathway_res_merged$pathway %in% summary_list$pathway &
                                        grepl("Member", pathway_res_merged$GroupID), ]

ggplot(pathway_res_use, aes(x = pathway, y = cell_type, fill = -1 * `Log(q-value)`)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "", fill = "-log10(q-value)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black', size = 12)
  ) +
  coord_flip() +
  scale_x_discrete(position = "top")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.path.barplot.png", width = 830 / 96, height = 700 / 96, dpi = 300)


pathway_res_use[pathway_res_use$pathway=="aerobic respiration (GO)", ]







cell_type_list = c("C7", "M5813", "M0610", "M24")

# deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")
pathway_res_use[pathway_res_use$pathway=="aerobic respiration (GO)", ]
deg_merged[deg_merged$cell_type %in% cell_type_list & deg_merged$gene_name %in% c("Mt-co3", "mt-Co3"), ]





### individual 
library(openxlsx)
pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/C7/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme_classic() +
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C7.path.barplot.png", width = 690 / 96, height = 300 / 96, dpi = 300)


pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/M5813/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme_classic() +
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.M5813.path.barplot.png", width = 480 / 96, height = 200 / 96, dpi = 300)


pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/M0610/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme_classic() +
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.M0610.path.barplot.png", width = 670 / 96, height = 300 / 96, dpi = 300)


pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/M24/metascape_result.xlsx", sheet = 2)
pathway_res <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")

ggplot(pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +  
  theme_classic() +
  theme(
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) + 
  scale_y_discrete(position = "right") + 
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.M24.path.barplot.png", width = 480 / 96, height = 200 / 96, dpi = 300)








### cell communication results
cc_df_diff <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.diff.out", sep = "\t", header = T)

lv_ec_list = c("C7", "M5813", "M0610", "M24")
lv_ec_list = c("C7", "M5813", "M0610", "M24", "C22", "C9", "C1", "C18")
lv_ec_list = paste0("EC", lv_ec_list)

cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type %in% lv_ec_list) %>% 
  mutate(EC_type = paste("Cluster", new_idx[gsub("EC", "", EC_type)]))

path_list = cc_df_use[abs(cc_df_use$value)>0.1, ] %>% arrange(desc(abs(value)))
path_list = unique(path_list$pathway_name)[1:15]
cc_df_use %>% 
  group_by(EC_type) %>% arrange(value) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  filter(pathway_name %in% path_list) %>%
  ggplot(aes(x = rank, y = value, colour = Strain)) +
  geom_point(aes(size=factor(abs(value)>0.1))) +
  geom_line() +
  xlab("Rank") +
  ylab("Differential communication score\n(treatment vs. control)") +
  scale_y_continuous(
    breaks = c(-0.1, 0, 0.1),
    labels = c(-0.1, 0, 0.1)
  ) +
  scale_size_manual(values=c(1, 2.5)) +
  facet_grid2(EC_type~pathway_name) +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.cellchat.dotplot.png", width=950/96, height=685/96, dpi=300)



### test Nrg1, Erbb4
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out")
ligand_expr <- deg_merged[deg_merged$cell_type=="C22" & deg_merged$gene_name=="Nrg1", ]
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T)
receptor_expr <- deg_merged[deg_merged$cell_type=="CM" & deg_merged$gene_name=="Erbb4" & deg_merged$tissue=="LV", ]

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
    label = "Nrg1\n(C18)",
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
    label = "Erbb4\n(cardiomyocyte)",
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
    title = "NRG signaling in\nhypertensive strains"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12)
    
  ) +
  facet_grid(~label)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/NRG.ligand_receptor.expr.png", width=450/96, height=260/96, dpi=300)




lv_sub_genes <- c("Kdr", "Epas1", "Eng", "Nrp1", "Nrp2", "Col4a1", "Col4a2", "Il1r1", "Calcrl", "Vwf", "Vegfc", "Smad6", "Efnb2", "Notch1")
lv_ec_list = c("C7", "M5813", "M0610", "M24")

seurat_object_lv_sub <- subset(seurat_object, seurat_clusters %in% c("M0610", "M24", "C7", "M5813"))
seurat_object_lv_sub$new_idx <- paste0("C", seurat_object_lv_sub$new_idx)
DotPlot(seurat_object_lv_sub, features = lv_sub_genes, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )



lv_sub_genes <- c("Nrg1")
lv_ec_list = c("C7", "M5813", "M0610", "M24")
lv_ec_list = c("C7", "M5813", "M0610", "M24", "C22", "C9", "C1", "C18")

seurat_object_lv <- subset(seurat_object, seurat_clusters %in% c("C7", "M5813", "M0610", "M24", "C22", "C9", "C1", "C18"))
seurat_object_lv$new_idx <- paste0("C", seurat_object_lv$new_idx)
DotPlot(seurat_object_lv, features = lv_sub_genes, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
  )
























