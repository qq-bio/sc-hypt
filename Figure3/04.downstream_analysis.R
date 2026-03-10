library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggh4x)
library(ggtext)
library(patchwork)
library(data.table)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))




################################################################################
### BBB EC results (Figure 3c-d; Fig. S12e-f)
################################################################################
### dot plot of selected metabolic genes

i <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.rds"
seurat_object <- readRDS(i)

marker_list = c("Mfsd2a", "Cpt1a", "Pfkfb3")
seurat_object_bbb <- subset(seurat_object, seurat_clusters %in% c("C14", "C20"))
DotPlot(seurat_object, features = marker_list, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    legend.position = "None"
  ) + coord_flip()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.metabolite.dotplot.png", width=500/96, height=131/96, dpi=300)


### BBB subtypes comparison
markers = seurat_object %>% subset(seurat_clusters %in% c(9,13)) %>% FindAllMarkers(., group.by="seurat_clusters", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  filter(p_val_adj<0.05) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(rank = row_number()) %>% 
  ungroup()

marker_top = markers %>% group_by(cluster) %>%
  filter(!grepl("ENSRNOG", gene)) %>%
  dplyr::mutate(rank = row_number()) %>% 
  filter(rank<=20) %>%
  select(gene)

marker_list = unique(marker_top$gene)
seurat_object %>% subset( seurat_clusters %in% c(9,13) ) %>%
  DotPlot(., features = marker_list) +
  scale_y_discrete(limits=rev, labels = cluster_names) +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="", size = "Percent\nExpressed", colour ="Average\nExpression") +
  theme(
    text = element_text(family="Arial"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.direction = "vertical", legend.box = "horizontal")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/BBB.EC.marker.dotplot.png", width=1279/96, height=281/96, dpi = 300)



### cell chat results of BBB subtypes
cc_df_diff <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.diff.out", sep = "\t", header = T)

cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  dplyr::mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  dplyr::mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type %in% c("EC9", "EC13")) %>% 
  dplyr::mutate(EC_type = paste("Cluster", gsub("EC", "", EC_type)))

cc_df_use <- cc_df_use %>% 
  group_by(EC_type) %>% arrange(value) %>%
  dplyr::mutate(rank = row_number()) %>% ungroup()

path_list = unique(cc_df_use$pathway_name)
path_list_1 = path_list[1:10]
path_list_2 = path_list[11:20]

p1 <- ggplot(cc_df_use[cc_df_use$pathway_name %in% path_list_1, ], aes(x = rank, y = value, colour = Strain)) +
  geom_point() + geom_line() + theme_classic() +
  xlab("Rank") +
  ylab("Differential communication score\n(treatment vs. control)") +
  scale_y_continuous(breaks = c(0, 0.1)) +
  facet_grid(EC_type~pathway_name) +
  theme(text = element_text(family="Arial")) +
  coord_flip()
p2 <- ggplot(cc_df_use[cc_df_use$pathway_name %in% path_list_2, ], aes(x = rank, y = value, colour = Strain)) +
  geom_point() + geom_line() + theme_classic() +
  xlab("Rank") +
  ylab("Differential communication score\n(treatment vs. control)") +
  scale_y_continuous(breaks = c(0, 0.1)) +
  facet_grid(EC_type~pathway_name) +
  theme(text = element_text(family="Arial")) +
  coord_flip()

p_combined <- arrangeGrob(p1, p2, nrow = 2)



### lined dot plot of pdgfd, pdgfrb
ligand_expr <- deg_merged[deg_merged$cell_type=="9" & deg_merged$gene_name=="Pdgfd", ]
receptor_expr <- fread("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out")[cell_type == "Pericyte" & gene_name == "Pdgfrb" & tissue == "HYP"]

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






################################################################################
### C12 results (Figure 3f-h)
################################################################################
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

target_gene <- "Klf2" # or Klf4

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




### organize scenic results
path <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/scenic/rat_scenic_res/rat_10k10k"

### load regulons file
regulonsfile <- file.path(path, "reg.csv") 
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
  select(Gene = TargetGeneList, TF, MotifID, AUC, NES)

### load AUCell file
aucellfile <- file.path(path,"aucell.csv")
aucell <- read_csv(aucellfile)
aucell <- as.data.frame(aucell)
rownames(aucell) <- aucell$Cell
aucell$Cell <- NULL

overlap_cell <- intersect(colnames(rat), rownames(aucell))
rat <- subset(rat, cell_id %in% overlap_cell)

aucell <- aucell[colnames(rat), ]

### add into seurat

regulon_mat <- t(as.matrix(aucell))
rat[["Regulon"]] <- CreateAssayObject(data = regulon_mat)
DefaultAssay(rat) <- "Regulon"
Idents(rat) <- "seurat_clusters"

markers <- FindAllMarkers(rat, group.by = "seurat_clusters",  only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0 ) 

write.csv(markers, file = paste0(path,"/regulon_markers.csv"))


### Rank plot (Figure 3g)
markers <- read.csv("rat_10k10k/regulon_markers.csv")

markers_c12 <- markers %>%
  filter(cluster == "12", !is.na(p_val), p_val > 0) %>%
  arrange(p_val) %>%
  dplyr::mutate(
    rank = row_number(),
    logp = -log10(p_val),
    point_color = ifelse(gene %in% c("Klf4(+)", "Klf2(+)"), "red", "grey")
  ) %>%
  arrange(gene %in% c("Klf2(+)", "Klf4(+)"))

top10_genes <- markers_c12 %>%
  top_n(-10, p_val) %>%
  mutate(label_color = ifelse(gene %in% c("Klf4(+)", "Klf2(+)"), "red", "black"))


my_theme <- theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    strip.placement = "outside",
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

plot_markers_c12 <- ggplot(markers_c12, aes(x = rank, y = logp)) +
  geom_point(aes(color = point_color), alpha = 1) +
  geom_text_repel(
    data = top10_genes,
    force = 2,   
    force_pull = 0.5,  
    aes(label = gene, color = label_color),
    segment.color = 'grey',
    size = 4,
    show.legend = FALSE,
    max.overlaps = 100
  ) +
  scale_color_identity() +
  labs(
    x = "Rank",
    y = "-log10(p-value)",
    title = "C12 specific TFs"
  ) +
  ylim(0, max(markers_c12$logp, na.rm = TRUE) * 1.1) +
  my_theme


ggsave(paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.C12.regulon.rank.png"), 
       plot = plot_markers_c12, 
       dpi = 300,        
       width = 226/96,  
       height = 363/96)  



### TF-gene network (Figure 3h)
TF_gene_pathway_all <- read.csv("rat_filter/TF_gene_pathway_all.csv")

klf4 <- TF_gene_pathway_all %>% filter(TF == "Klf4", data == "rat_10k10k")
klf2 <- TF_gene_pathway_all %>% filter(TF == "Klf2", data == "rat_10k10k")

combined_regulons <- bind_rows(klf2, klf4) %>%
  select(TF, Gene) %>%
  distinct() %>%
  rename(from = TF, to = Gene)

klf2_klf4_edges <- tibble(
  from = c("Klf2", "Klf4"),
  to   = c("Klf4", "Klf2")
)

all_edges <- bind_rows(combined_regulons, klf2_klf4_edges) %>%
  distinct()  

# Create graph object
graph <- tbl_graph(edges = all_edges, directed = TRUE)

graph <- graph %>%
  mutate(
    type = case_when(
      name %in% c("Klf2", "Klf4") ~ "TF",
      TRUE ~ "Gene"
    )
  )

layout_df <- create_layout(graph, layout = "fr")

layout_df <- layout_df %>%
  mutate(
    x = ifelse(name == "Zyx", x - 0.2, x),  # shift left
    y = ifelse(name == "Zyx", y + 0.2, y),  # shift up
    hjust = ifelse(name %in% c("Klf2", "Klf4"), -0.1, 0.5) 
  )

plot_TF_gene <- ggraph(layout_df) +
  geom_edge_link(color = "gray") +  # Just lines, no arrows
  geom_node_point(aes(color = type), size = 5) +
  geom_node_text(aes(label = name), size = 5, repel = TRUE) +
  scale_color_manual(values = c("TF" = "red3", "Gene" = "skyblue"), name = "Type") +
  theme_void() +
  theme(
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  ) 

ggsave(paste0("/xdisk/mliang1/jingh/multi-omics-RNAseq/SCENIC_hyp/rat_filter/plot_TF_gene.png"), 
       plot = plot_TF_gene, 
       dpi = 300,        
       width = 318/96,  
       height = 363/96)  







################################################################################
### WKY EC and VSMC results (Fig. S12g-h)
################################################################################
msa_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MSA.RNA.anno.L2.rds')
sxt_order = c("SHR-10w", "SHR-26w", "WKY-10w", "WKY-26w")
msa_shr$SxT = factor(msa_shr$SxT, levels = sxt_order)

genes_mech_core <- c("Piezo1", "Pecam1", "Kdr", "Tek", "Trpv4", "Cav1")
genes_adh_junc  <- c("Itgb1", "Itga5", "Cdh5", "Fak")

gene_score_list = list(
  `Primary mechanosensors` = genes_mech_core,
  `Adhesion and juncitonal proteins` = genes_adh_junc
)

msa_shr_score = calc_module_scores_and_summary(msa_shr, gene_score_list, c("Project", "Strain",  "Treatment", "SxT", "Cell_type_L1"))$summary_long

score_merged <- msa_shr_score %>%
  mutate(SxT = factor(SxT, levels = sxt_order),
         Cell_type_L1 = factor(Cell_type_L1, levels = cell_order))

ct_selected = c("EC")

p1 = score_vis(score_merged, score_selected = c("Primary mechanosensors"), project_selected = c("Spontaneous"), cell_type_selected = ct_selected) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 = score_vis(score_merged, score_selected = c("Adhesion and juncitonal proteins"), project_selected = c("Spontaneous"), cell_type_selected = ct_selected)

cowplot::plot_grid(
  p1, p2, 
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1.5)
)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/MSA.shr.mechano_score.dotplot.png", width=800/96, height=370/96, dpi=300)




### Pathway enrichment in VSMC (Fig. S12h)
pathway_res <- read.xlsx("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape.L1/WKY_MSA_10w_26w_VSMC/metascape_result.xlsx", sheet = 2)
outfile <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/WKY_MSA_10w_26w_VSMC.path.barplot.png"

pathway_res_use <- pathway_res[grepl("Summary", pathway_res$GroupID), ]
pathway_res_use$category <- sapply(strsplit(as.character(pathway_res_use$Category), " "), `[`, 1)
pathway_res_use$pathway <- paste0(pathway_res_use$Description, " (", pathway_res_use$category, ")")

ggplot(pathway_res_use[pathway_res_use$`Log(q-value)`< log10(0.05), ], aes(x = -`Log(q-value)`, y = fct_reorder(pathway, -`Log(q-value)`), fill = -`Log(q-value)`)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  labs(x = "-log10(q-value)", y = "", fill = "-log10(q-value)") +
  theme(
    text = element_text(family = "Arial"),
    axis.text = element_text(colour = 'black'),
    legend.position = "left"
  ) +
  labs(title = "WKY MSA VSMC (26w vs. 10w)") +
  scale_y_discrete(position = "right") +
  scale_fill_continuous(low = "white", high = "red", limits = c(0, 2), oob = scales::squish)
ggsave(outfile, width = 800 / 96, height = 370 / 96, dpi = 300)




################################################################################
### C14 results (Fig. S13)
################################################################################
# Fig. S13d
lv_sub_genes <- c("Nrg1")
lv_ec_list = c("11", "6", "1", "4", "10", "5", "2", "3", "14")

seurat_object_lv <- subset(seurat_object, seurat_clusters %in% lv_ec_list)
seurat_object_lv$new_idx <- paste0("C", seurat_object_lv$seurat_clusters)
DotPlot(seurat_object_lv, features = lv_sub_genes, group.by = "new_idx") +
  scale_color_gradient(low = "white", high = "firebrick") +
  labs(y="", x="") +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 12),
    legend.box = "horizontal"
  ) +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/NRG.expr.png", width=770/96, height=180/96, dpi=300)



# Fig. S13a
cc_df_diff <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.diff.out", sep = "\t", header = T)
lv_ec_list = c("EC11", "EC6", "EC1", "EC4", "EC10", "EC5", "EC2", "EC3", "EC14")
cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type %in% lv_ec_list) %>% 
  mutate(EC_type = paste0("C", gsub("EC", "", EC_type)))

path_list = cc_df_use[abs(cc_df_use$value)>0.1, ] %>% arrange(desc(abs(value)))
path_list = unique(path_list$pathway_name)[1:15]
cc_df_use %>% 
  group_by(EC_type) %>% arrange(value) %>%
  mutate(rank = row_number()) %>% ungroup() %>%
  filter(pathway_name %in% path_list) %>%
  ggplot(aes(x = rank, y = value, colour = Strain)) +
  geom_point(aes(size=factor(abs(value)>0.1))) +
  geom_line() +
  labs(x = "Rank", y = "Differential communication score\n(treatment vs. control)",
       size = "|value|>0.1") +
  scale_y_continuous(
    breaks = c(-0.1, 0, 0.1),
    labels = c(-0.1, 0, 0.1)
  ) +
  theme(text = element_text(family = "Arial")) +
  scale_size_manual(values=c(1, 2.5), labels = c("No", "Yes")) +
  facet_grid2(EC_type~pathway_name) +
  coord_flip()
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LV.EC.cellchat.dotplot.png", width=900/96, height=685/96, dpi=300)



# Fig. S13b
cc_df_use <- cc_df_diff %>% 
  filter(grepl("^EC", source) | grepl("^EC", target)) %>%
  filter(!(grepl("^EC", source) & grepl("^EC", target))) %>%
  mutate(EC_type = ifelse(grepl("^EC", source), source, target)) %>% 
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  mutate(Strain = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotensive"))%>% 
  filter(EC_type=="EC14") %>% 
  mutate(EC_type = paste0("C", gsub("EC", "", EC_type))) %>% 
  filter(pathway_name == "NRG")%>% 
  group_by(EC_type) %>% arrange(value) %>%
  mutate(rank = row_number()) %>% ungroup()

label_data <- cc_df_use %>% filter(strain=="SS")

cc_df_use  %>%
  ggplot(aes(x = rank, y = value, colour = Strain)) +
  geom_linerange(aes(x=rank, ymax=value, ymin=0), color = "grey") +
  geom_point(aes(color = strain), size=4) + # size=factor(abs(value)>0.1), 
  geom_text(
    data = label_data,
    aes(x = rank, y = value, label = treatment),
    hjust = 1,
    vjust = 1.5, 
    size = 4.5, color = "black"
  ) +
  labs(x = "", y = "Differential communication score\n(treatment vs. control)",
       title = "NRG signaling") +
  scale_y_continuous(
    breaks = c(-0.1, 0, 0.1),
    labels = c(-0.1, 0, 0.1)
  ) +
  geom_hline(yintercept=0) +
  theme(
    text = element_text(family = "Arial"),
    axis.ticks.y = element_blank(), axis.text.y = element_blank(),
    axis.text.x = element_text(color="black"),
    legend.box = "horizontal") +
  # scale_size_manual(values=c(2, 4), labels = c("No", "Yes")) +
  scale_color_manual(values = strain_col, name = "Strain") +
  # facet_grid2(~Strain) +
  coord_flip()

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/NRG.cellchat.strain_wise.png", width=380/96, height=180/96, dpi=300)



# Fig. S13c
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.DEG_all.out")
ligand_expr <- deg_merged[deg_merged$cell_type=="14" & deg_merged$gene_name=="Nrg1", ]
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = T) 
receptor_expr <- deg_merged[deg_merged$cell_type=="CM" & deg_merged$gene_name=="Erbb4" & deg_merged$tissue=="LV", ]

ligand_long <- ligand_expr %>%
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
    label = "Nrg1\n(C14)",
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
plot_data$label <- factor(plot_data$label, levels = unique(plot_data$label))

label_data <- plot_data %>% filter(condition=="Treatment", strain=="SS")

ggplot(data = plot_data, aes(x = condition, y = pct_expr, group = group)) +
  geom_point(aes(fill = strain),
             size = 5, shape = 21, color = "white") +
  scale_fill_manual(values = strain_col, name = "Strain") +
  
  geom_line(aes(color = log2FC, linewidth = significance)) +
  scale_color_gradient2(low = "blue", mid = "gray90", high = "red", midpoint = 0, name = "log2FC") +
  scale_linewidth_manual(values = c("Yes" = 1.2, "No" = 0.4)) +
  
  geom_text_repel(
    data = label_data,
    aes(x = condition, y = pct_expr, label = treatment),
    hjust = 1,
    vjust = 0.5, 
    size = 4.5, color = "black"
  ) +
  
  labs(
    x = NULL,
    y = "Percent Expressed",
    # color = "log2FC",
    linewidth = "Significant",
    title = "NRG signaling in hypertensive strains"
  ) +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    axis.text.x = element_text(size = 10, color = "black"), #angle = 45, hjust = 1, size = 12, 
    axis.text.y = element_text(size = 10, color = "black"),
    strip.text = element_text(size = 12),
    legend.box = "horizontal"
  ) +
  facet_wrap(~label)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/NRG.ligand_receptor.expr.png", width=770/96, height=260/96, dpi=300)
























