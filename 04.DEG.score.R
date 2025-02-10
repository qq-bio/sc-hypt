library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggsci)
library(ggupset)
library(ggh4x)
library(lemon)

# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DEG_enrichment.R")



################################################################################
### elbow
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)

DEG_df <- deg_merged %>%
  group_by(project, strain, treatment, tissue, cell_type) %>%
  summarize(
    FC_0.1 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 0.1),
    FC_0.25 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25),
    FC_0.5 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5),
    FC_0.75 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 0.75),
    FC_1 = sum(p_val_adj < 0.05 & abs(avg_log2FC) > 1),
    .groups = 'drop'
  ) %>%
  as.data.frame()




################################################################################
### DEG figure
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)

deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]
deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$project %in% c("AngII", "Salt-sensitive")),]

DEG_df = deg_merged %>% 
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>% 
  group_by(project, strain, treatment, tissue, cell_type) %>%
  dplyr::select(project, strain, treatment, tissue, cell_type, control_size, treatment_size, test_gene) %>%
  mutate(
    DEG_num = n(),
    mean_cell_size = (mean(control_size, na.rm = TRUE) + mean(treatment_size, na.rm = TRUE)) / 2,
    min_cell_size = pmin(control_size, treatment_size),
    max_cell_size = pmax(control_size, treatment_size),
    cell_size_diff = abs(control_size - treatment_size),
    combined = factor(interaction(strain, treatment, sep = "-"), 
                      levels = c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d",
                                 "SS-HS 3d", "SS-HS 21d", "SD-HS 3d",
                                 "SHR-26w", "WKY-26w"))
  ) %>%
  unique() %>%
  as.data.frame()

DEG_df$cell_type = factor(DEG_df$cell_type, levels=unique(cell_order))
DEG_df$tissue = factor(DEG_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))

species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
alpha <- setNames(c(0.5, 1, 1, 0.5, 1), unique(DEG_df$treatment))

combined <- unique(DEG_df$combined)
combined_colors <- setNames(species_col[gsub("\\-.*", "", combined)], combined)
combined_alphas <- setNames(alpha[gsub(".*\\-", "", combined)], combined)

ggplot(DEG_df, aes(x = DEG_num, y = cell_type)) +
  geom_segment(aes(yend = cell_type), xend = 0, colour = "grey50") +
  geom_point(aes(colour = combined, alpha = combined, size=log10(max_cell_size)*3)) +
  scale_y_discrete(limits=rev) +
  scale_colour_manual(values = combined_colors) +
  scale_alpha_manual(values = combined_alphas) +
  scale_size_continuous(breaks = c(6, 9, 12), labels = c("100", "1,000", "10,000")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="Number of DEGs", y="", color="Strain x treatment", alpha="Strain x treatment", size="Number of cells") +
  facet_grid2(cols = vars(project), rows = vars(tissue), 
             scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.deg.0.5.png", width=577/96,height=913/96,dpi=300)





library(lemon)

ggplot(DEG_df, aes(x = DEG_num, y = cell_type)) +
  geom_segment(aes(yend = cell_type), xend = 0, colour = "grey50") +
  geom_point(aes(colour = combined, alpha = combined, size=log10(max_cell_size)*3)) +
  scale_y_discrete(limits=rev) +
  scale_colour_manual(values = combined_colors) +
  scale_alpha_manual(values = combined_alphas) +
  scale_size_continuous(breaks = c(6, 9, 12), labels = c("100", "1,000", "10,000")) +
  scale_x_continuous(breaks =seq(0, 500, 100)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="Number of DEGs", y="", color="Strain x treatment", alpha="Strain x treatment", size="Number of cells") +
  facet_grid2(cols = vars(project), rows = vars(tissue), 
              scales = "free", space = "free_y") +
  coord_capped_cart(xlim = c(0, 500))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.deg.0.5.rescale.png", width=577/96,height=913/96,dpi=300)



################################################################################
### cell type wise
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", sep=' ', header=T)
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$tissue = factor(deg_merged$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_merged$species = factor(deg_merged$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

pathway_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", header = T, sep = '\t', quote = "")

mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
mouse2human = unique(mouse2human[, c("Human.gene.name", "Gene.name")])

pathway_gene = pathway_all %>% filter(str_detect(GroupID, "Member")) %>%
  separate_rows(Symbols, sep = ",") %>% distinct() %>%
  left_join(mouse2human, by = c("Symbols" = "Human.gene.name"), relationship = "many-to-many") %>%
  separate_rows(Gene.name, sep = ",") %>% distinct() %>%
  group_by(pathway) %>% dplyr::summarise(Symbols_h = list(unique(Symbols)), Symbols_m = list(unique(Gene.name)))

pathway_gene <- setNames(pathway_gene$Symbols_m, pathway_gene$pathway)

# deg_use = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>1 & deg_merged$cell_type=="EC", ]
# deg_use = deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$gene_name %in% c("Heg1", "Pdgfb", "Slc40a1", "Flt1", "Pecam1", "Sema3c", "Ndk1", "Eng", "Prom1", "Tek",
#                                                                               "Clic5", "Epas1", "Robo4", "Bsg", "Id1"), ]
# 
# deg_use = deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$gene_name %in% c("Smad3", "Pparg", "Lmna", "Foxo1", "Errfi1", "Atp2b4", 
#                                                                                                           "Mef2c", "Flt1", "Pparg", "Tsc2", "Col4a3",
#                                                                                                           "Sox17", "Acvr1", "Smad4", "Nrg1", "Notch1",
#                                                                                                           "Rims2", "Erc2", "Erc1", "Rimbp2", "Rims1"), ]
# 
# deg_use = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>1 & deg_merged$cell_type=="EC" & 
#                        deg_merged$gene_name %in% c("Lamb2", "Ntrk3", "Tnr", 
#                                                    "Flt1", "Pparg", "Tsc2", "Col4a3",
#                                                    "Pparg", "Lmna", "Foxo1", "Errfi1"), ]

pathway_id = "Protein-protein interactions at synapses (Reactome)"; cell_type="Inhibitory neurons"
pathway_id = "Protein-protein interactions at synapses (Reactome)"; cell_type="OPC"
pathway_id = "glutamate receptor signaling pathway (GO)"; cell_type="OPC"
pathway_id = "glutamate receptor signaling pathway (GO)"; cell_type="Microglia"
pathway_id = "NABA CORE MATRISOME (Canonical)"; cell_type="FIB"
pathway_id = "blood vessel development (GO)"; cell_type="FIB"
pathway_id = "blood vessel morphogenesis (GO)"; cell_type="FIB"
pathway_id = "receptor clustering (GO)"; cell_type="OPC"
pathway_gene_use = pathway_gene[pathway_id][[1]]
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       deg_merged$cell_type==cell_type & 
                       deg_merged$tissue=="HYP" & 
                       deg_merged$gene_name %in% pathway_gene_use, ]


ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  # geom_segment(aes(yend = cell_type), xend = 0, colour = "grey50") +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  # scale_size_continuous(breaks = c(6, 9, 12), labels = c("100", "1,000", "10,000")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", title = paste0(pathway_id, "\n", cell_type)) +
  facet_nested(~tissue + species, scales = "free", space = "free")





for(i in c("mean_cell_size", "min_cell_size", "max_cell_size", "test_gene")){
  p = ggplot(DEG_df, aes(x = DEG_num, y = get(i))) +
    geom_point(aes(color = species)) + 
    geom_smooth(method = lm, se = FALSE) +
    labs(y=i) +
    scale_color_manual(values = species_col)
  print(p)
}







deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header = T, row.names = 1)

enrichment_analysis(deg_df)







################################################################################
#### gene wise
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$tissue = factor(deg_merged$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_merged$species = factor(deg_merged$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.pseudo.merged_all.out", sep='\t', header=T)
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$species = factor(deg_pseudo$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "species", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "species", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "species", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = cell_order)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$species = factor(expr_all$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

mouse2human[mouse2human$Human.gene.name=="MRPS6", ]
gene_name = c("Rras", "Npr1", "Ppara", "Serpinc1", "Malat1")
gene_name = c("Npr3", "Nrip1", "Samsn1", "Mrpl23", "Sod1", "Smim11", "Smim34", "Slc5a3", "Rcan1", "Mrps6")
gene_name = "Npr1"
deg_use = deg_merged[  deg_merged$p_val<0.05 &
                       deg_merged$gene_name %in% gene_name, ]
p1=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
  geom_point(aes(fill = -avg_log2FC, size=-log10(p_val),
                 shape = p_val_adj<0.05)) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  scale_shape_manual(values = c(23, 21)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="log2(FC)", size="-log10(p-value)", shape = "adjusted.p < 0.05", title = paste(gene_name, "(DEG-wilcox)")) +
  facet_nested(tissue ~ project + species, scales = "free", space = "free")

deg_use = deg_pseudo[  deg_pseudo$p_val<0.05 &
                         deg_pseudo$gene_name %in% gene_name, ]
p2=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
  geom_point(aes(fill = -avg_log2FC, size=-log10(p_val),
                 shape = p_val_adj<0.05)) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                       high = "red", space = "Lab" ) +
  scale_shape_manual(values = c(23, 21)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="log2(FC)", size="-log10(p-value)", shape = "adjusted.p < 0.05", title = paste(gene_name, "(DEG-pseudo+deseq2)")) +
  facet_nested(tissue ~ project + species, scales = "free", space = "free")

deg_use = expr_all[expr_all$gene_name %in% gene_name, ]
p3=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
  geom_point(aes(color = log2(avg_expr+1), size=pct*100)) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient(low = "white", high = "red") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", color="log2(expression+1)", size="Percentage", title = paste(gene_name, "(overall expression)")) +
  facet_nested(tissue ~ project + species, scales = "free", space = "free")

# library(cowplot)
pp = plot_grid(plotlist=list(p1, p2), ncol=1, align='v')
plot_grid(plotlist=list(pp, p3), ncol=2, align='h', )

p1+p2+p3+plot_layout(design="
                     AC
                     BC")

gene_name_list = c("Rras", "Npr1", "Ppara", "Serpinc1", "Malat1", "Npr3", "Nrip1", "Samsn1", "Mrpl23", "Sod1", "Smim11", "Smim34", "Slc5a3", "Rcan1", "Mrps6")
for(gene_name in gene_name_list){
  deg_use = expr_all[expr_all$gene_name %in% gene_name, ]
  p3=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
    geom_point(aes(color = log2(avg_expr+1), size=pct*100)) +
    scale_y_discrete(limits=rev) +
    scale_color_gradient(low = "white", high = "red") +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),   # No horizontal grid lines
      legend.justification = c(1, 0.5),
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      strip.text = element_text(colour = 'black')
    ) +
    labs(x="", y="", color="log2(expression+1)", size="Percentage", title = paste(gene_name, "(overall expression)")) +
    facet_nested(tissue ~ project + species, scales = "free", space = "free")
  print(p3)
}


tmp1 = deg_use[,-c(4, 13)]; tmp2 = deg_use[,-c(3, 12)]; names(tmp2) = names(tmp1)
deg_use2 = rbind(tmp1, tmp2)
ggplot(deg_use2, aes(x = control, y = cell_type)) +
  geom_point(aes(fill = -avg_log2FC, size=pct.1,
                 shape = p_val_adj<0.05)) +
  scale_y_discrete(limits=rev) +
  # scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
  #                      high = "red", space = "Lab" ) +
  scale_shape_manual(values = c(23, 21)) +
  # scale_size_continuous(breaks = c(6, 9, 12), labels = c("100", "1,000", "10,000")) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="log2(FC)", size="-log10(p-value)", shape = "adjusted.p < 0.05", title = gene_name) +
  facet_nested(tissue + cell_type ~ project + species, scales = "free", space = "free")



################################################################################
### upset plot (https://blog.csdn.net/weixin_45822007/article/details/121413844)


SS_HYP_EC = deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$species=="C57BL/6" & deg_merged$treatment=="AngII 3d", ]$gene_name
deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$species=="C57BL/6" & deg_merged$treatment=="AngII 28d", ]$gene_name

deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$species=="C57BL/6" & deg_merged$treatment=="AngII 28d", ]$gene_name
deg_merged[deg_merged$p_val_adj<0.05 & deg_merged$cell_type=="EC" & deg_merged$species=="C57BL/6" & deg_merged$treatment=="AngII 28d", ]$gene_name


input = list(ang_c15=unique(ang_marker$Cluster.15), ang_c19=unique(ang_marker$Cluster.19),
             sp_c13=unique(sp_marker$Cluster.13), sp_c14=unique(sp_marker$Cluster.14))

upset(fromList(input), order.by = "freq")























de_gene_sum_v2 = de_gene_merged_mod %>% 
  dplyr::group_by(cell_type) %>%
  dplyr::summarize(a = sum(p_val_adj.ctrl<0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m>0),
                   b = sum(p_val_adj.ctrl<0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m<0),
                   c = sum(p_val_adj.ctrl<0.05 & (p_val_adj.f-0.05)*(p_val_adj.m-0.05)<0),
                   d = sum(p_val_adj.ctrl>0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m>0),
                   e = sum(p_val_adj.ctrl>0.05 & p_val_adj.f<0.05 & p_val_adj.m<0.05 & avg_log2FC.f*avg_log2FC.m<0),
                   f = sum(p_val_adj.ctrl>0.05 & (p_val_adj.f-0.05)*(p_val_adj.m-0.05)<0)) %>%
  ungroup() %>% as.data.frame()

de_gene_sum_v2$cell_type = factor(de_gene_sum_v2$cell_type, levels = c("CM-1", "CM-2", "CM-3", "CM-4", 
                                                                       "EC-1", "EC-2", "EC-3", "EC-4", 
                                                                       "FB-1", "FB-2", "FB-3", "PC", 
                                                                       "IMM", "MACRO-1", "MACRO-2", "DC", 
                                                                       "NEU", "UNKNOWN"
))
# de_gene_sum_v2$sex_biased_gene = de_gene_sum_v2$a + de_gene_sum_v2$b + de_gene_sum_v2$c
# de_gene_sum_v2$sex_ind_gene = de_gene_sum_v2$d + de_gene_sum_v2$e + de_gene_sum_v2$f

de_gene_sum_v2_mlt = melt(de_gene_sum_v2[,1:7])
ggplot(de_gene_sum_v2_mlt, aes(y = cell_type, x = variable, label = value)) +
  geom_point(aes(size=log10(value)*3), alpha = .5, color="red") + geom_text() +
  theme(axis.text.y = element_text(size=12, colour = 'black'),
        axis.text.x = element_text(size=12, colour = 'black'),
        legend.text = element_text(size=12, colour = 'black')) +
  # scale_size_continuous(breaks = c(1, 2, 3)*3) +
  scale_size_continuous(breaks = c(1, 2, 3)*3,
                        labels = c("10", "100", "1000")) +
  labs(x="Expression pattern", y="", size="Number of genes")















cor.test(DEG_df$DEG_num, DEG_df$mean_cell_size)
cor.test(DEG_df$DEG_num, pmin(DEG_df$control_size, DEG_df$treatment_size))
cor.test(DEG_df$DEG_num, pmax(DEG_df$control_size, DEG_df$treatment_size))
cor.test(DEG_df$DEG_num, DEG_df$test_gene)


### TAS
model <- lm(DEG_num ~ mean_cell_size + test_gene, data = DEG_df)
summary(model)
DEG_df$tas_scores <- predict(model, newdata = DEG_df)
print(DEG_df$tas_scores)






