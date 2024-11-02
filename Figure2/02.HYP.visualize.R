
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(CellChat)
library(ggh4x)
library(patchwork)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")





################################################################################
### Visualize DEGs by Occurrence Rank (Figure 2a)
################################################################################
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged <- deg_merged[deg_merged$tissue == "HYP" & deg_merged$p_val_adj < 0.05 & abs(deg_merged$avg_log2FC) > 0.25, ]

deg_merged$cell_type <- factor(deg_merged$cell_type, levels = cell_order)
deg_merged$treatment <- factor(deg_merged$treatment, levels = treatment_order)
deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)


selected_genes <- c("Fth1", "Rora", "Apoe", "Ptgds", "Csmd1", "App", "Ncam1", "Hsp90ab1", "Zbtb16", "Lars2")

p <- deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  group_by(gene_name) %>%
  mutate(deg_num = n()) %>% 
  ungroup() %>%
  select(gene_name, deg_num) %>%
  distinct() %>%
  arrange(desc(deg_num)) %>%
  mutate(
    index = row_number(),
    highlight = ifelse(gene_name %in% selected_genes, "Highlighted", "Normal")
  ) %>%
  arrange(desc(highlight)) 

# Plot DEG occurrences
p <- ggplot(plot_data, aes(x = deg_num, y = index, color = highlight)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) + 
  scale_y_continuous(trans = "reverse") +
  ggrepel::geom_label_repel(
    aes(label = ifelse(gene_name %in% selected_genes, gene_name, "")), 
    size = 4, 
    max.overlaps = Inf,
    force = 10,
    box.padding = 1,
    point.padding = 0.5,
    color = "black"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Number of DEG Occurrences", y = "Rank of Genes by DEG Occurrence")

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2a.HYP.deg.dot.png", width=353/96, height=287/96, dpi=300)





################################################################################
### Visualize Hsp90ab1 Expression (Figure 2b)
################################################################################
hyp_m <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds')
hyp_ss <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds')
hyp_shr <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds')

hyp_m$sxt <- factor(paste0(hyp_m$species, " - ", hyp_m$treatment), levels = sxt_order)
hyp_ss$sxt <- factor(paste0(hyp_ss$species, " - ", hyp_ss$treatment), levels = sxt_order)
hyp_shr$sxt <- factor(paste0(hyp_shr$species, " - ", hyp_shr$treatment), levels = sxt_order)

gene <- c("Hsp90ab1")
p1 <- hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment", order = T) & labs(x="", y="") & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"))
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2b.HYP.angii.Hsp90ab1.umap.png", width=498/96, height=155/96, dpi=300)

p2 <- hyp_ss %>% subset(species=="SS") %>% FeaturePlot(., features = gene, split.by = "sxt", order = T) & labs(x="", y="") & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"))
print(p2)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2b.HYP.ss.Hsp90ab1.umap.png", width=498/96, height=155/96, dpi=300)





################################################################################
### Visualize Consensus Changes in snRNA-seq & Bulk RNA-seq (Figure 2c)
################################################################################

bulk <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/bulk/DEG_bulk.csv", header = TRUE, sep = ",")

selected_genes <- c("Snhg11", "Meg3", "Sntg1", "Zbtb16", "Csmd1", "ENSRNOG00000065867",
                    "Cacna1c", "Lars2", "Gpm6b", "Ncam1")

# Plot snRNA-seq differential expression
p1 <- deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_genes) %>%
  mutate(gene_name = factor(gene_name, levels = selected_genes)) %>%
  ggplot(aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size = -log10(p_val_adj))) +
  scale_y_discrete(limits = rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "snRNA-seq", x = "Cell type-wise log2(FC)", y = "", color = "log2(FC)", size = "-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

# Plot bulk RNA-seq differential expression
p2 <- bulk %>%
  filter(tissue == "HYP", strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_genes) %>%
  mutate(gene_name = factor(gene_name, levels = selected_genes),
         strain = factor(strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))) %>%
  ggplot(aes(x = log2FoldChange, y = gene_name)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_point(shape = 21, aes(fill = log2FoldChange, size = -log10(padj))) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical",
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "Bulk RNA-seq", x = "log2(FC) in bulk", y = "", fill = "log2(FC)", size = "-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(fill = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

p1+p2+ plot_layout(widths = c(1,1))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2c.HYP.sn_bulk.consist_expr.png", width=832/96, height=356/96, dpi=300)







################################################################################
### Pathway Enrichment Visualization (Figure 2d)
################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = merge_reshape[merge_reshape$tissue=="HYP",]

sxt_rename <- c(
  "Log.q.value..C57BL.6.AngII.3d" = "C57BL/6 - AngII 3d", 
  "Log.q.value..C57BL.6.AngII.28d" = "C57BL/6 - AngII 28d", 
  "Log.q.value..SS.HS.3d" = "SS - HS 3d", 
  "Log.q.value..SS.HS.21d" = "SS - HS 21d",
  "Log.q.value..SHR.26w" = "SHR - 26w"
)

melted_data <- merge_reshape %>%
  filter(category != "Other", summary == "Yes", model_count > 1) %>%
  pivot_longer(cols = starts_with("Log.q.value"), names_to = "condition", values_to = "log_q_value") %>%
  filter(condition %in% names(sxt_rename)) %>%
  mutate(
    condition = sxt_rename[condition],
    log_q_value = ifelse(log_q_value < 1.3, 0, log_q_value),
    pathway = fct_reorder(pathway, cell_type),
    cell_type = fct_inorder(cell_type),
    condition = factor(condition, levels = c("C57BL/6 - AngII 3d", "C57BL/6 - AngII 28d", "SS - HS 3d", "SS - HS 21d", "SHR - 26w"))
  )


# Load group info and merge with data
group_info <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/HYP.multi_model.pathway.group.out", header = FALSE, sep = '\t')
colnames(group_info) <- c("pathway", "group")
group_info <- unique(group_info)
rownames(group_info) <- group_info$pathway
melted_data$group <- group_info[as.character(melted_data$pathway), "group"]

# Heatmap for pathway enrichment
p_data <- ggplot(melted_data %>% filter(log_q_value > 1.3), aes(x = condition, y = pathway, fill = log_q_value)) +
  geom_tile(color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    plot.margin = margin(0, 0, 0, -2, "cm")
  ) +
  labs(x = "", y = "", fill = "-log(q-val)", title = "Shared Pathways in at Least Two Hypertension Models") +
  facet_grid(cols = vars(cell_type), scales = "free", space = "free")

# Bar plot for pathway groups
p_y <- ggplot(melted_data) +
  geom_col(aes(y = pathway, x = 1, fill = group), width = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(
    axis.text.y = element_text(colour = 'black'),
    legend.position = "right",
    plot.margin = margin(0, -2, 0, 0, "cm")
  ) +
  labs(fill = "Category")

p_y + p_data + plot_layout(widths = c(0.03, 1), guides = "collect")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2d.pathway.heatmap.png", width=1325/96, height=548/96, dpi=300)







################################################################################
### Visualize CellChat Results (Figure 2e and f)
################################################################################
cc_df <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.result.out", 
                    header = TRUE, sep = '\t')

cc_df$model <- factor(cc_df$model, levels = model_order)
cc_df$strain <- factor(cc_df$strain, levels = strain_order)
cc_df$tissue <- factor(cc_df$tissue, levels = tissue_order)
cc_df$treatment <- factor(cc_df$treatment, levels = treatment_order)
cc_df$sxt <- paste0(cc_df$strain, "-", cc_df$treatment)

id_vars <- setdiff(colnames(cc_df), c("sxt", "value"))
cc_df_reshape <- reshape(cc_df, idvar = id_vars, timevar = "sxt", direction = "wide")
cc_df_reshape[is.na(cc_df_reshape)] <- 0

# Calculate differential interaction strength for selected conditions
cc_df_reshape$`diff.C57BL/6-AngII 3d` <- cc_df_reshape$`value.C57BL/6-AngII 3d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.C57BL/6-AngII 28d` <- cc_df_reshape$`value.C57BL/6-AngII 28d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.SS-HS 3d` <- cc_df_reshape$`value.SS-HS 3d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SS-HS 21d` <- cc_df_reshape$`value.SS-HS 21d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SD-HS 3d` <- cc_df_reshape$`value.SD-HS 3d` - cc_df_reshape$`value.SD-LS`
cc_df_reshape$`diff.SHR-26w` <- cc_df_reshape$`value.SHR-26w` - cc_df_reshape$`value.SHR-10w`
cc_df_reshape$`diff.WKY-26w` <- cc_df_reshape$`value.WKY-26w` - cc_df_reshape$`value.WKY-10w`

# Reshape data to long format for plotting
diff_cols <- grep("diff", colnames(cc_df_reshape), value = TRUE)
cc_df_diff <- reshape2::melt(cc_df_reshape[, c(id_vars, diff_cols)], 
                             id.vars = id_vars, measured.vars = diff_cols,
                             variable.name = "sxt")
cc_df_diff$sxt <- gsub("diff.", "", cc_df_diff$sxt)
cc_df_diff$treatment <- factor(sapply(strsplit(cc_df_diff$sxt, "-"), `[`, 2), levels = treatment_order)
cc_df_diff <- cc_df_diff[cc_df_diff$value != 0, ]

# Summarize the differential values for each pathway interaction
cc_df_use <- cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name, interaction_name_2) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# Define selected pathways and interactions for visualization
pathway_selected <- c("NRXN", "NCAM", "NEGR") 
interaction_name_list <- c("Ncam1  - Ncam2", "Negr1  - Negr1", "Nrxn1  - Nlgn1", "Nrxn3  - Nlgn1")
cc_df_use$name_use <- paste0(cc_df_use$pathway_name, ": ", cc_df_use$interaction_name_2)

# Generate heatmap of interaction strengths in selected pathways for HYP tissue
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & 
                   cc_df_use$interaction_name_2 %in% interaction_name_list & 
                   cc_df_use$tissue == "HYP", ]) +
  geom_tile(aes(x = treatment, y = name_use, fill = prob_sum)) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  theme(
    panel.grid.major.y = element_blank(), 
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x = "", y = "", fill = "Differential\nInteraction\nStrength") +
  facet_nested(~ model + strain, scales = "free", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2e.HYP.cellchat.heatmap.png", 
       width = 596/96, height = 241/96, dpi = 300)





e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.rds")

saline3d=e$`mouse.HYP_C57BL.6_Saline 3d`
angii3d=e$`mouse.HYP_C57BL.6_AngII 3d`
angii28d=e$`mouse.HYP_C57BL.6_AngII 28d`
ss_ls=e$rat.ss.HYP_SS_LS
ss_hs3d=e$`rat.ss.HYP_SS_HS 3d`
ss_hs21d=e$`rat.ss.HYP_SS_HS 21d`
sd_ls=e$rat.ss.HYP_SD_LS
sd_hs3d=e$`rat.ss.HYP_SD_HS 3d`
shr_10w=e$rat.sp.HYP_SHR_10w
shr_26w=e$rat.sp.HYP_SHR_26w
wky_10w=e$rat.sp.HYP_WKY_10w
wky_26w=e$rat.sp.HYP_WKY_26w


title = c("saline3d"="Saline 3d", "angii3d"="Ang II 3d")
for( i in c("saline3d", "angii3d") ){
  cellchat = get(i)
  outfile = paste0("fig2f.", i, "cellchat.png")
  
  png(filename = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/", outfile), 
      width = 786/96*300, height = 291/96*300, res = 300)
  
  par(mfrow=c(1, 3))
  
  for(j in c("NCAM", "NEGR", "NRXN") ){
    pathways.show <- j
    
    netVisual_aggregate(cellchat, signaling = pathways.show, 
                        layout = "circle", title.space = 1,
                        remove.isolate = T,
                        weight.scale = F, 
                        edge.weight.max = 1, 
                        edge.width.max = 8 
    )
    
    title(main = paste0(title[i], " - ", j))
  }
  
  dev.off()
  par(mfrow=c(1, 1))
}







################################################################################
### Visualize Avp Expression (Figure 2g)
################################################################################

neuron_list <- c("Inhibitory neurons", "Excitatory neurons", "Avp+ neurons")

hyp_m %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_m$project))
hyp_ss %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_ss$project))
hyp_shr %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_shr$project))


subset_data <- hyp_m %>% subset(subcluster == "Avp+ neurons")
p1 <- VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + theme(legend.position = 'None') + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons\n(AngII)")
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1g.HYP.angii.avp.expr.png", width=356/96, height=307/96, dpi=300)

subset_data <- hyp_ss %>% subset(subcluster=="Glu-15")
p2 <- VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p2)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1g.HYP.ss.avp.expr-1.png", width=356/96, height=270/96, dpi=300)
subset_data <- hyp_ss %>% subset(subcluster=="GABA-10")
p3 <- VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p3)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1g.HYP.ss.avp.expr-2.png", width=356/96, height=270/96, dpi=300)

subset_data <- hyp_shr %>% subset(subcluster=="Glu-15")
p4 <- VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons\n(Spontaneous)") + theme(legend.position = 'None')
print(p4)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1g.HYP.shr.avp.expr.png", width=356/96, height=270/96, dpi=300)







################################################################################
### Correlation of Neuronal Connection Score and AVP Expression in GTEx Data (Figure 2h)
################################################################################
library(GSVA)

# Load and preprocess GTEx data
gtex <- read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_brain_hypothalamus.gct", 
                   header = TRUE, sep = "\t", skip = 2)

# Filter out genes with too many zero values and log-transform data
gtex_use <- gtex[rowSums(gtex == 0) < (ncol(gtex) - 3) * 0.9, ]
gene_name <- gtex_use$Description
gtex_use <- t(gtex_use[, -c(1:3)])
gtex_use <- log2(as.matrix(gtex_use) + 1)
colnames(gtex_use) <- gene_name

# Define gene set for neuronal connection score and calculate GSVA score
gene_sets <- list(com = c("NCAM1", "NCAM2", "FGFR1", "L1CAM", "NEGR1", 
                          "NRXN1", "NRXN2", "NRXN3", "NLGN1", "NLGN2", "NLGN3"))
gsva_scores <- gsva(t(gtex_use), gene_sets, method = "ssgsea")

# Extract AVP expression and calculate correlation with GSVA score
avp_expr <- gtex_use[, "AVP"]
correlation <- cor.test(as.numeric(gsva_scores), as.numeric(avp_expr))

# Prepare data frame for plotting
score_avp <- data.frame(score = as.numeric(gsva_scores),
                        avp = as.numeric(avp_expr))

# Generate scatter plot with regression line and correlation annotation
ggplot(score_avp, aes(x = avp, y = score)) +
  geom_point() + 
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  annotate("text", x = 5, y = 3, label = paste0("R = ", round(correlation$estimate, 2), 
                                                ", p = ", format(correlation$p.value, scientific = TRUE)),
           size = 5) + 
  theme_classic() +
  labs(x = "AVP Expression", y = "Neuronal Connection Score", title = "Human Hypothalamus Bulk RNA-seq")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2h.avp_neu_score.cor.png", 
       width = 348/96, height = 253/96, dpi = 300)






################################################################################
### Visualize SCENIC Results for Ep300 Activity in Avp+ Neurons (Figure 2i)
################################################################################

hyp_m <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds')
hyp_ss <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds')
hyp_shr <- readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds')

# Load AUC data for SCENIC results
mhyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/mhyp_all_regulon_AUC.Ep300.csv", header = TRUE)
sshyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_ssall_regulon_AUC.Ep300.csv", header = TRUE)
shrhyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_shrall_regulon_AUC.Ep300.csv", header = TRUE)


# Prepare Ep300 activity for visualization
# Mouse Data
rownames(mhyp_auc) <- mhyp_auc[,1]
mhyp_auc <- mhyp_auc[,-1]
colnames(mhyp_auc) <- gsub("\\.", "-", colnames(mhyp_auc))
avp_neurons <- colnames(hyp_m)[hyp_m$new.cluster.ids_umap == "Avp+ neurons"]

avp_tf <- data.frame(
  Regulon_activity = as.numeric(mhyp_auc["Ep300(+)", avp_neurons]),
  Regulon = as.numeric(hyp_m@assays$RNA@data["Ep300", avp_neurons]),
  Avp = as.numeric(hyp_m@assays$RNA@data["Avp", avp_neurons]),
  treatment = factor(hyp_m$treatment[avp_neurons], levels = c("Saline 3d", "AngII 3d", "AngII 28d")),
  sxt = paste0("C57BL/6 - ", hyp_m$treatment[avp_neurons])
)

# Plot for mouse Avp+ neurons
ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = TRUE, scale = "width", fill = species_col['C57BL/6']) +
  geom_jitter(size = 0.1, color = "grey") +
  theme_classic() +
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons\n(AngII)") +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.angii.ep300.activity.png", width = 365/96, height = 189/96, dpi = 300)


# Rat Salt-Sensitive Data - Glutamatergic Avp+ Neurons
rownames(sshyp_auc) <- sshyp_auc[,1]
sshyp_auc <- sshyp_auc[,-1]
colnames(sshyp_auc) <- gsub("\\.", "-", colnames(sshyp_auc))
common_columns <- intersect(colnames(sshyp_auc), colnames(hyp_ss))
sshyp_auc <- sshyp_auc[, common_columns]
hyp_ss <- subset(hyp_ss, cells = common_columns)
avp_neurons <- colnames(hyp_ss)[hyp_ss$subcluster == "Glu-15"]

avp_tf <- data.frame(
  Regulon_activity = as.numeric(sshyp_auc["Ep300(+),", avp_neurons]),
  Regulon = as.numeric(hyp_ss@assays$RNA@data["Ep300", avp_neurons]),
  Avp = as.numeric(hyp_ss@assays$RNA@data["Avp", avp_neurons]),
  treatment = factor(hyp_ss$treatment[avp_neurons], levels = c("LS", "HS 3d", "HS 21d")),
  sxt = paste0("SS - ", hyp_ss$treatment[avp_neurons])
)

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = TRUE, scale = "width", fill = species_col['SS']) +
  geom_jitter(size = 0.1, color = "grey") +
  theme_classic() +
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.activity-1.png", width = 365/96, height = 189/96, dpi = 300)


# Rat Salt-Sensitive Data - GABAergic Avp+ Neurons
avp_neurons <- colnames(hyp_ss)[hyp_ss$subcluster == "GABA-10"]

avp_tf <- data.frame(
  Regulon_activity = as.numeric(sshyp_auc["Ep300(+),", avp_neurons]),
  Regulon = as.numeric(hyp_ss@assays$RNA@data["Ep300", avp_neurons]),
  Avp = as.numeric(hyp_ss@assays$RNA@data["Avp", avp_neurons]),
  treatment = factor(hyp_ss$treatment[avp_neurons], levels = c("LS", "HS 3d", "HS 21d")),
  sxt = paste0("SS - ", hyp_ss$treatment[avp_neurons])
)

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = TRUE, scale = "width", fill = species_col['SS']) +
  geom_jitter(size = 0.1, color = "grey") +
  theme_classic() +
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.activity-2.png", width = 365/96, height = 189/96, dpi = 300)




