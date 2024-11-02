
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(CellChat)
library(ggsci)
library(ggupset)
library(ggh4x)
library(patchwork)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")



################################################################################
### Visualize Enriched Pathway Counts (Figure 3a)
################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = merge_reshape[merge_reshape$tissue=="LK",]

# Define a function to summarize pathway counts based on log q-values for each model
summarize_pathway_counts <- function(data, model_label, comparisons) {
  data %>%
    filter(summary == "Yes") %>%
    group_by(tissue, cell_type) %>%
    summarize(across(all_of(comparisons), ~ sum(. > 1.3), .names = "count_{.col}")) %>%
    pivot_longer(cols = starts_with("count_"), names_to = "Comparison", values_to = "Count") %>%
    mutate(model = model_label) %>%
    filter(Count > 0) %>%
    as.data.frame()
}

# Define comparisons and labels for each model
ss_comparisons <- c("Log.q.value..SS.SD", "Log.q.value..SS.HS.3d", "Log.q.value..SS.HS.21d")
shr_comparisons <- c("Log.q.value..SHR.WKY", "Log.q.value..SHR.26w")
mouse_comparisons <- c("Log.q.value..C57BL.6.AngII.3d", "Log.q.value..C57BL.6.AngII.28d")

# Summarize pathway counts for each model
path_count_ss <- summarize_pathway_counts(merge_reshape, "Salt-sensitive", ss_comparisons)
path_count_shr <- summarize_pathway_counts(merge_reshape, "Spontaneous", shr_comparisons)
path_count_mouse <- summarize_pathway_counts(merge_reshape, "AngII", mouse_comparisons)

# Combine all model data
path_count_merge <- bind_rows(path_count_ss, path_count_shr, path_count_mouse)

# Set factor levels for cell type and tissue
path_count_merge$cell_type <- factor(path_count_merge$cell_type, levels = cell_order)
path_count_merge$tissue <- factor(path_count_merge$tissue, levels = tissue_order)

# Plot pathway count data
ggplot(path_count_merge, aes(y = cell_type, x = Comparison, label = Count)) +
  geom_point(aes(size = Count / 3), alpha = 0.5, color = "red") +
  geom_text() +
  theme_classic() +
  theme(
    axis.text.y = element_text(color = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, color = 'black'),
    legend.text = element_text(color = 'black')
  ) +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(
    breaks = c(1, 3, 6),
    labels = c("3", "9", "18")
  ) +
  labs(x = "", y = "", size = "Number of\nEnriched\nPathways") +
  facet_grid(cols = vars(model), scales = "free", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3a.LK.pathway-count.png", width = 455 / 96, height = 360 / 96, dpi = 300)





################################################################################
### Visualize Gene Expression of Selected Pathways (Figure 3b)
################################################################################
# Load the data and merge with strain-specific DEGs
deg_strain <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep = '\t', header = TRUE)
deg_strain[deg_strain$treatment == "SS-LS", ]$treatment <- "SS-SD"

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- rbind(deg_merged, deg_strain)

# Filter for significant DEGs
deg_merged <- deg_merged[deg_merged$p_val_adj < 0.05 & abs(deg_merged$avg_log2FC) > 0.5, ]
deg_merged$cell_type <- factor(deg_merged$cell_type, levels = cell_order)
deg_merged$tissue <- factor(deg_merged$tissue, levels = tissue_order)
deg_merged$treatment <- factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "SHR-WKY", "10w", "26w", "SS-SD", "LS", "HS 3d", "HS 21d"))
deg_merged$strain <- factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))

# Define selected genes for each pathway category
selected_genes <- c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", "Chrm3", # circulatory system process
  "Arhgap24", "Igf1r", "Pkhd1", # intracellular signal transduction
  "Ank3", "Cacnb4", "Efna5", "Erbb4", "Ptprd", "Sdk1" # cell junction organization
)

# Filter and categorize the data
pathway_gene_use <- selected_genes
deg_use <- deg_merged %>%
  filter(cell_type == "EC", tissue == "LK", gene_name %in% pathway_gene_use) %>%
  mutate(
    project = ifelse(strain == "Salt-sensitive", "Baseline\nDifference", "Treatment vs. Control"),
    strain = ifelse(strain == "Salt-sensitive", "SS vs. SD", as.character(strain)),
    gene_cat = case_when(
      gene_name %in% c("Arhgap24", "Igf1r", "Pkhd1") ~ "Intracellular Signal\nTransduction",
      gene_name %in% c("Ank3", "Cacnb4", "Efna5", "Erbb4", "Ptprd", "Sdk1") ~ "Cell Junction\nOrganization",
      TRUE ~ "Circulatory System\nProcess"
    )
  )

# Set levels for consistent ordering
deg_use$gene_cat <- factor(deg_use$gene_cat, levels = c("Circulatory System\nProcess", "Cell Junction\nOrganization", "Intracellular Signal\nTransduction"))

# Plot gene expression for selected pathways
ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size = -log10(p_val))) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    legend.position = "right",
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.text.y = element_text(angle = 360, colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "", y = "", color = "log2(FC)", size = "-log10(p-value)", 
       title = "LK - EC") +
  facet_nested(gene_cat ~ project + strain, scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3b.LK.pathway-gene.png", width = 637 / 96, height = 477 / 96, dpi = 300)






################################################################################
### Visualize EC Subtype, Adaptive Score, and MECOM+ EC Proportion (Figure 3c and d)
################################################################################
lk_ec_ss <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
lk_ec_shr <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")

lk_ec_ss$sxt <- paste0(lk_ec_ss$strain, " - ", lk_ec_ss$treatment)
lk_ec_shr$sxt <- paste0(lk_ec_shr$strain, " - ", lk_ec_shr$treatment)

# Define cell colors and adaptive score genes
cell_order <- sort(unique(c(lk_ec_ss$subclass_level2, lk_ec_shr$subclass_level2)))
getPalette <- colorRampPalette(brewer.pal(8, "Set2"))
lk_ec_col <- setNames(getPalette(length(cell_order)), cell_order)

ec_adaptive_genes <- c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", # circulatory system process
  "Arhgap24", "Igf1r", "Pkhd1", # intracellular signal transduction
  "Ank3", "Cacnb4", "Efna5", "Erbb4", "Ptprd", "Sdk1" # cell junction organization
)

# Add EC adaptive score and plot UMAP for SS and SHR data
plot_umap_and_score <- function(data, project_name, file_name) {
  data <- AddModuleScore(object = data, features = list(ec_adaptive_genes), name = "EC_Adaptive_Score")
  
  p1 <- DimPlot(data, label = TRUE, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col, repel = TRUE) +
    labs(title = project_name, x = "UMAP 1", y = "UMAP 2") + theme(legend.position = "None")
  
  p2 <- FeaturePlot(data, features = "EC_Adaptive_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") +
    labs(title = "EC adaptive score", x = "UMAP 1", y = "UMAP 2", color = "Score")
  
  p1 + p2
  ggsave(file_name, width = 753 / 96, height = 341 / 96, dpi = 300)
}

plot_umap_and_score(lk_ec_ss, unique(lk_ec_ss$project), "/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3c.LK.ss.EC.umap.png")
plot_umap_and_score(lk_ec_shr, unique(lk_ec_shr$project), "/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3c.LK.shr.EC.umap.png")

# Plot adaptive score distribution
adaptive_score_data <- rbind(lk_ec_ss@meta.data, lk_ec_shr@meta.data) %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SS", "SHR"), "Hypertensive", "Normotensive")) %>%
  filter(subclass_level2 %in% c("EC_Mecom")) %>%
  dplyr::select(sxt, EC_Adaptive_Score1, hypt, strain) %>%
  mutate(sxt = factor(sxt, levels = sxt_order)) %>%
  group_by(sxt) %>%
  mutate(mean_value = mean(EC_Adaptive_Score1))

ggplot(adaptive_score_data, aes(x = sxt, y = EC_Adaptive_Score1, fill = mean_value)) +
  geom_violin(trim = FALSE) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(y = "EC adaptive\nscore", x = "", fill = "Mean") +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.position = "right", text = element_text(size = 12)) +
  facet_nested(~ hypt + strain, scales = "free")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3d.LK.EC.score.vln.png", width = 405 / 96, height = 275 / 96, dpi = 300)

# Plot MECOM+ EC proportion across conditions
mecom_prop_data <- rbind(lk_ec_ss@meta.data, lk_ec_shr@meta.data) %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SS", "SHR"), "Hypertensive", "Normotensive")) %>%
  group_by(seqID2, subclass_level2, hypt, strain, treatment) %>%
  dplyr::summarise(cell_count = n()) %>%
  ungroup() %>%
  group_by(seqID2) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  filter(subclass_level2 == "EC_Mecom") %>%
  group_by(hypt, strain, treatment, subclass_level2) %>%
  dplyr::summarise(mean_proportion = mean(proportion), sd_proportion = sd(proportion))

ggplot(mecom_prop_data, aes(x = treatment, y = mean_proportion, fill = subclass_level2)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = mean_proportion - sd_proportion, ymax = mean_proportion + sd_proportion),
                position = position_dodge(width = 0.8), width = 0.2) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = lk_ec_col) +
  labs(x = "Treatment", y = "Cell proportion") +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.position = "None") +
  facet_nested(subclass_level2 ~ hypt + strain, scales = "free")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig3d.LK.EC.prop-bar.png", width = 351 / 96, height = 240 / 96, dpi = 300)






