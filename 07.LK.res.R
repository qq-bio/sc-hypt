dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

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

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))









################################################################################
# https://stackoverflow.com/questions/71986699/plotting-heatmap-with-triangular-split-tiles-of-more-than-one-categorical-variab
# https://stackoverflow.com/questions/71147521/split-overlapping-tiles-by-facet-in-geom-tile


################################################################################
### preprocess cellchat results
cc_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.result.out", header = T, sep = '\t')

cc_df$model = factor(cc_df$model, levels=c("AngII", "Salt-sensitive", "Spontaneous"))
cc_df$strain = factor(cc_df$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
cc_df$tissue = factor(cc_df$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
cc_df$treatment = factor(cc_df$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w"))
cc_df$sxt = paste0(cc_df$strain, "-", cc_df$treatment)

id_vars = setdiff(colnames(cc_df), c("sxt", "value"))
cc_df_reshape = reshape(cc_df, idvar = id_vars, timevar = "sxt", direction = "wide")
cc_df_reshape[is.na(cc_df_reshape)] <- 0

cc_df_reshape$`diff.C57BL/6-AngII 3d` = cc_df_reshape$`value.C57BL/6-AngII 3d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.C57BL/6-AngII 28d` = cc_df_reshape$`value.C57BL/6-AngII 28d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.SS-HS 3d` = cc_df_reshape$`value.SS-HS 3d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SS-HS 21d` = cc_df_reshape$`value.SS-HS 21d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SD-HS 3d` = cc_df_reshape$`value.SD-HS 3d` - cc_df_reshape$`value.SD-LS`
cc_df_reshape$`diff.SHR-26w` = cc_df_reshape$`value.SHR-26w` - cc_df_reshape$`value.SHR-10w`
cc_df_reshape$`diff.WKY-26w` = cc_df_reshape$`value.WKY-26w` - cc_df_reshape$`value.WKY-10w`

diff_cols = colnames(cc_df_reshape)[grepl("diff", colnames(cc_df_reshape))]
dis_cols = c("treatment", colnames(cc_df_reshape)[grepl("value.", colnames(cc_df_reshape))])
id_vars = setdiff(colnames(cc_df_reshape), c(dis_cols, diff_cols))
cc_df_diff = reshape2::melt(cc_df_reshape[, c(id_vars, diff_cols)], id.vars = id_vars, measured.vars=diff_cols,
                            variable.name = "sxt")
cc_df_diff$sxt = gsub("diff.", "", cc_df_diff$sxt)
cc_df_diff$treatment = as.character(lapply(strsplit(cc_df_diff$sxt, "-"), function(x) x[2]))
cc_df_diff$treatment = factor(cc_df_diff$treatment, c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w"))
cc_df_diff = cc_df_diff[cc_df_diff$value!=0, ]



e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.rds")
mouse_lk_list = list(saline3d=e$`mouse.LK_C57BL.6_Saline 3d`, angii3d=e$`mouse.LK_C57BL.6_AngII 3d`, angii28d=e$`mouse.LK_C57BL.6_AngII 28d`)
mouse_lk = mergeCellChat(mouse_lk_list, add.names = names(mouse_lk_list))

ss_lk_list = list(ss_ls=e$rat.ss.LK_SS_LS, ss_hs3d=e$`rat.ss.LK_SS_HS 3d`, ss_hs21d=e$`rat.ss.LK_SS_HS 21d`, sd_ls=e$rat.ss.LK_SD_LS, sd_hs3d=e$`rat.ss.LK_SD_HS 3d`)
ss_lk = mergeCellChat(ss_lk_list, add.names = names(ss_lk_list))

sp_lk_list = list(shr_10w=e$rat.sp.LK_SHR_10w, shr_26w=e$rat.sp.LK_SHR_26w, wky_10w=e$rat.sp.LK_WKY_10w, wky_26w=e$rat.sp.HYP_WKY_26w)
sp_lk = mergeCellChat(sp_lk_list, add.names = names(sp_lk_list))



################################################################################
### all changes
### all changes
cell_order = c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", # "Microglia", "Activated microglia", 
    "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("Art EC", "Cap EC", "A Cap EC", "Prolif Cap EC", "Ven Cap EC", "Ven EC", "Ang EC", "CR EC", "EC_Mecom", "Glom EC", "Lymph EC", "Imm EC", "IFN EC", 
    "EC_Il1r1", "EC_Prdm16", "EC_Lama2", "EC_Syt1", "EC_Zeb2", "EC-NM", "EC_Acta2", "EC_Lgr5", "EC_Adgrg6", "EC_Nrp1", "EC_Plin1"),
  c("VSMC", "E/P transition cell", "Pericyte", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
    "NK cells", "NKT", "T cells", "B cells")
)

cc_df_diff$source = factor(cc_df_diff$source, cell_order)
cc_df_diff$target = factor(cc_df_diff$target, cell_order)
cc_df_use = cc_df_diff %>% filter(!(is.na(source) | is.na(target))) %>%
  group_by(sxt, model, strain, tissue, treatment, source, target) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

ggplot(cc_df_use[cc_df_use$tissue=="LK", ]) +
  geom_tile(mapping = aes(x = target, y = source, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested( ~ strain + treatment , scales = "free", space = "free")




cc_df_diff$source = factor(cc_df_diff$source, cell_order)
cc_df_diff$target = factor(cc_df_diff$target, cell_order)

for(cell_selected in c("A Cap EC")){
  for(i in c("source", "target")){
    cc_df_use = cc_df_diff %>% filter(!(is.na(source) | is.na(target))) %>%
      group_by(sxt, model, strain, tissue, treatment, !!sym(i), pathway_name) %>%
      dplyr::summarise(prob_sum = sum(value)) %>%
      as.data.frame()
    
    title = paste0(cell_selected, " (", i, ")")
    p = ggplot(cc_df_use[cc_df_use[,i] %in% cell_selected & cc_df_use$tissue=="LK", ]) +
      geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
      # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
      scale_y_discrete(limits=rev) +
      scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
      theme(
        panel.grid.major.y = element_blank(),   # No horizontal grid lines
        # legend.position = c(1, 0.55),           # Put legend inside plot area
        legend.justification = c(1, 0.5),
        axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        # text = element_text(colour = 'black'),
        strip.text = element_text(colour = 'black')
      ) +
      labs(x="", y="", fill="Differential\ninteraction\nstrength", title=title) +
      facet_nested( ~ model + strain, scales = "free", space = "free")
    print(p)
  }
}



cc_df_use = cc_df_diff %>% filter(!(is.na(source) | is.na(target))) %>%
  group_by(sxt, model, strain, tissue, treatment, source, target, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
pathway_selected = c("PTPRM") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & cc_df_use$tissue=="LK", ]) +
  geom_tile(mapping = aes(x = target, y = source, fill=prob_sum)) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested( ~ strain + treatment , scales = "free", space = "free")

# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.cellchat.heatmap.png", width=494/96, height=227/96, dpi=300)

pathway_selected = c("APP") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & cc_df_use$tissue=="LK", ]) +
  geom_tile(mapping = aes(x = target, y = source, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested( ~ strain + treatment, scales = "free", space = "free")



cell_selected = c("Prolif Cap EC") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$target %in% cell_selected & cc_df_use$tissue=="LK", ]) +
  geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")


cc_df_use = cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name, interaction_name_2) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
pathway_selected = c("NRXN", "NCAM", "NEGR") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
interaction_name_list = c("Ncam1  - Ncam2", "Negr1  - Negr1", "Nrxn1  - Nlgn1", "Nrxn3  - Nlgn1")
cc_df_use$name_use = paste0(cc_df_use$pathway_name, ": ", cc_df_use$interaction_name_2)
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & 
                   cc_df_use$interaction_name_2 %in% interaction_name_list & 
                   cc_df_use$tissue=="HYP", ]) +
  geom_tile(mapping = aes(x = treatment, y = name_use, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested( ~ model + strain , scales = "free", space = "free")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.cellchat.heatmap.png", width=596/96, height=241/96, dpi=300)


gg1 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,2), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(C57BL/6 - AngII 3d vs Saline 3d)")
gg2 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,3), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(C57BL/6 - AngII 28d vs Saline 3d)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/mouse.HYP.cellchat.all.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(1,2), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(SHR - 26w vs 10w)")
gg2 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(3,4), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(WKY - 26w vs 10w)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/sp.HYP.cellchat.NCAM.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,2), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(SS - HS 3d vs LS)")
gg2 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,3), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(SS - HS 21d vs LS)")
gg3 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(4,5), cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength\n(SD - HS 3d vs LS)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/ss.HYP.cellchat.NCAM.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2 + gg3)
dev.off()




gg1 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,2), signaling = "NCAM", cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(C57BL/6 - AngII 3d vs Saline 3d)")
gg2 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,3), signaling = "NCAM",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(C57BL/6 - AngII 28d vs Saline 3d)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/mouse.HYP.cellchat.NCAM.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(1,2), signaling = "NCAM", cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(SHR - 26w vs 10w)")
gg2 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(3,4), signaling = "NCAM",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(WKY - 26w vs 10w)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/sp.HYP.cellchat.NCAM.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()


gg1 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,2), signaling = "NCAM", cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(SS - HS 3d vs LS)")
gg2 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,3), signaling = "NCAM",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(SS - HS 21d vs LS)")
gg3 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(4,5), signaling = "NCAM",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NCAM\n(SD - HS 3d vs LS)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/ss.HYP.cellchat.NCAM.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2 + gg3)
dev.off()



gg1 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,2), signaling = "NRXN", cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NRXN\n(C57BL/6 - AngII 3d vs Saline 3d)")
gg2 <- netVisual_heatmap(mouse_hyp, measure = "weight", comparison = c(1,3), signaling = "NRXN",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NRXN\n(C57BL/6 - AngII 28d vs Saline 3d)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/mouse.HYP.cellchat.NRXN.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()

gg1 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(1,2), signaling = "NRXN", cluster.rows = T, cluster.cols = F,
                         title.name = "Differential interaction strength of NRXN\n(SHR - 26w vs 10w)")
gg2 <- netVisual_heatmap(sp_hyp, measure = "weight", comparison = c(3,4), signaling = "NRXN",  cluster.rows = T, cluster.cols = F,
                         title.name = "Differential interaction strength of NRXN\n(WKY - 26w vs 10w)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/sp.HYP.cellchat.NRXN.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2)
dev.off()


gg1 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,2), signaling = "NRXN", cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NRXN\n(SS - HS 3d vs LS)")
gg2 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(1,3), signaling = "NRXN",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NRXN\n(SS - HS 21d vs LS)")
gg3 <- netVisual_heatmap(ss_hyp, measure = "weight", comparison = c(4,5), signaling = "NRXN",  cluster.rows = T, cluster.cols = T,
                         title.name = "Differential interaction strength of NRXN\n(SD - HS 3d vs LS)")
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/ss.HYP.cellchat.NRXN.png", width=786*3, height=448*3, res=300)
print(gg1 + gg2 + gg3)
dev.off()




################################################################################
### prepocess deg results
deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep='\t', header=T)
deg_strain[deg_strain$treatment=="SHR-10w",]$treatment = "SHR-WKY"
deg_strain[deg_strain$treatment=="SS-LS",]$treatment = "SS-SD"

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = rbind(deg_merged, deg_strain)

deg_merged = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5, ]
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$tissue = factor(deg_merged$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "SHR-WKY", "10w", "26w", "SS-SD", "LS", "HS 3d", "HS 21d"))
deg_merged$strain = factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))


pathway_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", header = T, sep = '\t', quote = "")

mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
mouse2human = unique(mouse2human[, c("Human.gene.name", "Gene.name")])

pathway_gene = pathway_all %>% filter(str_detect(GroupID, "Member")) %>%
  separate_rows(Symbols, sep = ",") %>% distinct() %>%
  left_join(mouse2human, by = c("Symbols" = "Human.gene.name"), relationship = "many-to-many") %>%
  separate_rows(Gene.name, sep = ",") %>% distinct() %>%
  group_by(pathway) %>% dplyr::summarise(Symbols_h = list(unique(Symbols)), Symbols_m = list(unique(Gene.name)))

pathway_gene <- setNames(pathway_gene$Symbols_m, pathway_gene$pathway)


################################################################################
pathway_id = "cell junction organization (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "regulation of plasma membrane bounded cell projection organization (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "negative regulation of intracellular signal transduction (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "circulatory system process (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "response to hormone (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "cell junction organization (GO)"; cell_type="PT"; tissue="LK"
pathway_id = "Chemical carcinogenesis - reactive oxygen species (KEGG)"; cell_type="PT"; tissue="LK"
pathway_id = "intracellular chemical homeostasis (GO)"; cell_type="PT"; tissue="LK"
pathway_id = "enzyme-linked receptor protein signaling pathway (GO)"; cell_type="PT"; tissue="LK"
pathway_id = "Diseases of metabolism (Reactome)"; cell_type="PT"; tissue="LK"
pathway_id = "Chemical carcinogenesis - reactive oxygen species (KEGG)"; cell_type="TAL"; tissue="LK"
pathway_id = "Electron transport chain OXPHOS system in mitochondria (WikiPathways)"; cell_type="TAL"; tissue="LK"
pathway_id = "cellular response to lipid (GO)"; cell_type="TAL"; tissue="LK"
pathway_id = "cellular response to hormone stimulus (GO)"; cell_type="TAL"; tissue="LK"

pathway_gene_use = pathway_gene[pathway_id][[1]]
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.5 &
                       deg_merged$cell_type==cell_type & 
                       deg_merged$tissue==tissue & 
                       deg_merged$gene_name %in% pathway_gene_use, ]

ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0(tissue, " - ", cell_type, "\n", pathway_id)) +
  facet_nested(~tissue + strain, scales = "free", space = "free")


selected_genes = c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", "Chrm3", # circulatory system process
  "Arhgap24", "Igf1r", "Pkhd1", # negative regulation of intracellular signal transduction (GO)
  "Ank3","Cacnb4","Efna5","Erbb4","Ptprd","Sdk1" # cell junction organization
)

pathway_gene_use = selected_genes; cell_type="EC"; tissue="LK"
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.5 &
                       deg_merged$cell_type==cell_type & 
                       deg_merged$tissue==tissue & 
                       deg_merged$gene_name %in% pathway_gene_use, ]
deg_use$project = "Treatment vs. control"
deg_use[deg_use$strain=="Salt-sensitive",]$project = "Baseline\ndifference"
levels(deg_use$strain) = c(levels(deg_use$strain), "SS vs. SD")
deg_use[deg_use$strain=="Salt-sensitive",]$strain = "SS vs. SD"
deg_use$gene_cat = "Circulatory\nsystem process"
deg_use[deg_use$gene_name %in% c("Arhgap24", "Igf1r", "Pkhd1"), ]$gene_cat = "Intracellular signal\ntransduction"
deg_use[deg_use$gene_name %in% c("Ank3","Cacnb4","Efna5","Erbb4","Ptprd","Sdk1"), ]$gene_cat = "Cell junction\norganization"
deg_use$gene_cat = factor(deg_use$gene_cat, levels = c("Circulatory\nsystem process", "Cell junction\norganization", "Intracellular signal\ntransduction"))
p1 = ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  # scale_y_discrete(limits=rev) +
  # scale_y_discrete(limits = rev(levels(fct_reorder(deg_use$gene_name, deg_use$gene_cat)))) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    legend.direction = "vertical",
    legend.title = element_text(size = 10),
    # legend.title.align = 0.5,
    legend.position = "right",
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.text.y = element_text(angle = 360, colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0("\n", tissue, " - ", cell_type)) +
  facet_nested(gene_cat ~ project + strain, scales = "free", space = "free_y")
p1
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.pathway-gene.png", width=637/96, height=477/96, dpi=300)

# for i in Sdk1 Ptprd Erbb4 Efna5 Cacdb4 Ank3 Slit2 Slc4a4 2 Ptprj Nox4 Nedd4l Lrp2 Col4a3 Chrm3 Atp1b1 Pkhd1 Igf1r Arhgap24; do grep $i DEG_bulk.csv | grep SS | grep LK; done



deg_pseudo = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', header=T)
colnames(deg_pseudo)[12] = "strain"
deg_pseudo$cell_type = factor(deg_pseudo$cell_type, levels = cell_order)
deg_pseudo$tissue = factor(deg_pseudo$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_pseudo$treatment = factor(deg_pseudo$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_pseudo$strain = factor(deg_pseudo$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

expr_all = unique(rbind(data.frame(unname(deg_pseudo[, c("pct.1", "avg_expr.1", "gene_name", "cell_type", "project", "strain", "tissue", "control")])),
                        data.frame(unname(deg_pseudo[, c("pct.2", "avg_expr.2", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")]))))
names(expr_all) = c("pct", "avg_expr", "gene_name", "cell_type", "project", "strain", "tissue", "treatment")
expr_all$cell_type = factor(expr_all$cell_type, levels = cell_order)
expr_all$tissue = factor(expr_all$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
expr_all$treatment = factor(expr_all$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
expr_all$strain = factor(expr_all$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))

# deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.pseudo.DEG.all.out", sep='\t', header=T)
# deg_strain$cell_type = factor(deg_strain$cell_type, levels = cell_order)
# deg_strain$tissue = factor(deg_strain$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
# deg_strain$treatment = factor(deg_strain$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
# deg_strain$strain = factor(deg_strain$strain, levels = c("Salt-sensitive", "Spontaneous"))
# levels(deg_strain$strain) = c("SS", "SHR")
# deg_strain$project = "Baseline"
# 
# deg_pseudo = rbind(deg_pseudo, deg_strain)
# deg_pseudo$project = factor(deg_pseudo$project, levels = c("AngII", "Salt-sensitive", "Spontaneous", "Baseline"))
# 

selected_genes = c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2", "Chrm3", # circulatory system process
  "Arhgap24", "Efna5", "Igf1r", "Ptprd", # regulation of plasma membrane bounded cell projection organization (GO)
  "Pkhd1" # negative regulation of intracellular signal transduction (GO)
)
selected_genes = c("Ank3","Cacnb4","Efna5","Erbb4","Gpc4","Igf1r","Pkhd1","Tns1","Sdk1","Eng","Flt1","Magt1","Prkch","Ehd4","Mast4","Ptprb","Ldb2","Pbx1")

cell_type="EC"; tissue="LK"
deg_use = expr_all[expr_all$gene_name %in% selected_genes &
                     expr_all$cell_type==cell_type & 
                     expr_all$tissue==tissue, ]
ggplot(deg_use, aes(x = paste(strain, "-", treatment), y = log2(avg_expr+1))) +
  geom_line(aes(group=gene_name)) + 
  geom_point(aes(color = pct, size=pct*100)) +
  # scale_y_discrete(limits=rev) +
  scale_color_gradient(low = "white", high = "red") +
  # theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", color="log2(expression+1)", size="Percentage", title = paste(gene_name, "(overall expression)")) +
  facet_nested( ~ project, scales = "free", space = "free")

################################################################################

deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/EC.strain_wise.DEG.all.out", sep='\t', header=T)
deg_strain[deg_strain$treatment=="SHR-10w",]$treatment = "SHR-WKY"
deg_strain[deg_strain$treatment=="SS-LS",]$treatment = "SS-SD"

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/EC.DEG.all.out", sep='\t', header=T)
deg_merged = rbind(deg_merged, deg_strain)

deg_merged$tissue = factor(deg_merged$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "SHR-WKY", "10w", "26w", "SS-SD", "LS", "HS 3d", "HS 21d"))
deg_merged$strain = factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY", "Salt-sensitive", "Spontaneous"))

deg_merged = deg_merged[!(deg_merged$cell_type %in% c("PT Cell", "UnID Cell", "Imm Cell")),]


selected_genes = c("Ptprm", "Nrp1", "Sbk2", "Lrrc30", "Ctnnb1", "Ctnnd1", "Ptpn6", "Ptprq", "Ptpn3", "Pts", "Ptpra")
selected_genes = c(
  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2", "Chrm3", # circulatory system process
  "Arhgap24", "Efna5", "Igf1r", "Ptprd", # regulation of plasma membrane bounded cell projection organization (GO)
  "Pkhd1" # negative regulation of intracellular signal transduction (GO)
)

pathway_gene_use = selected_genes; cell_type="Ang EC"; tissue="LK"
deg_use = deg_merged[deg_merged$p_val<0.05 &
                       # abs(deg_merged$avg_log2FC)>0.5 &
                       # deg_merged$cell_type==cell_type & 
                       deg_merged$tissue==tissue & 
                       deg_merged$gene_name %in% pathway_gene_use, ]
deg_use$project = "Treatment vs. control"
# deg_use[deg_use$strain=="Salt-sensitive",]$project = "Baseline difference"
p1 = ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    # legend.title.align = 0.5,
    legend.position = "bottom",
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0("\n", tissue, " - ", cell_type)) +
  facet_nested(cell_type ~ project + strain, scales = "free", space = "free_y")

p1




ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds")
sp = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds")

marker_list = c("Efna5", "Lgf1r", "Nedd4l", "Ptprd", "Ptprj", "Nox4", "Atp1b1", "Chrm3")
DotPlot(ss, features = marker_list, group.by = "RNA_snn_res.0.5")


  
pathway_id = "EC-status"; tissue="LK"; 
ec_list=c("Art EC", "Cap EC", "A Cap EC", "Prolif Cap EC", "Ven Cap EC", "Ven EC", "Ang EC", "CR EC", "EC_Mecom", "Glom EC", "Lymph EC", "Imm EC", "IFN EC", 
            "EC_Il1r1", "EC_Prdm16", "EC_Lama2", "EC_Syt1", "EC_Zeb2", "EC-NM", "EC_Acta2", "EC_Lgr5", "EC_Adgrg6", "EC_Nrp1", "EC_Plin1")

permeability_genes <- c("Cdh5", "Cldn5", "Ocln", "Pecam1")  # Permeability
inflammatory_genes <- c("Vcam1", "Icam1", "Sele")  # Inflammatory Status
thrombogenic_genes <- c("F3", "Vwf")  # Thrombogenic Status
angiogenesis_genes <- c("Vegfa", "Fgf1", "Fgf2", "Angpt1")  # Angiogenesis
antioxidant_genes <- c("Nos3", "Sod1", "Sod2", "Sod3")  # Antioxidant Status
senescence_genes <- c("Trp53", "Cdkn1a")  # Senescence/Agingmarker_list = c("Efna5", "Lgf1r", "Nedd4l", "Ptprd", "Ptprj", "Nox4", "Atp1b1", "Chrm3")
marker_list <- c(permeability_genes, inflammatory_genes, thrombogenic_genes, 
                 angiogenesis_genes, antioxidant_genes, senescence_genes)
ss_ec = subset(ss, subclass_level2 %in% ec_list)
DotPlot(ss_ec, features = marker_list, group.by = "subclass_level2")
  





pathway_gene_use = c()
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.5 &
                       deg_merged$cell_type %in% cell_type & 
                       deg_merged$tissue==tissue & 
                       deg_merged$gene_name %in% pathway_gene_use, ]

ggplot(deg_use, aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0(pathway_id, "\n", tissue, " - ", cell_type)) +
  facet_nested(~tissue + strain, scales = "free", space = "free")







################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = merge_reshape[merge_reshape$tissue=="LK",]

path_count_ss = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`SS-LS vs. SD-LS` = sum(Log.q.value..SS.SD>1.3),
                   `SS-HS 3d vs. SS-LS` = sum(Log.q.value..SS.HS.3d>1.3),
                   `SS-HS 21d vs. SS-LS` = sum(Log.q.value..SS.HS.21d>1.3)
  ) %>% melt() %>% 
  mutate(model = "Salt-sensitive") %>% filter(value>0) %>% as.data.frame()

path_count_shr = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`SHR-10w vs. WKY-10w` = sum(Log.q.value..SHR.WKY>1.3),
                   `SHR-26w vs. SHR-10w` = sum(Log.q.value..SHR.26w>1.3)
  ) %>% melt() %>% 
  mutate(model = "Spontaneous") %>% filter(value>0) %>% as.data.frame()

path_count_mouse = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`C57BL/6-AngII 3d vs. Saline 3d` = sum(Log.q.value..C57BL.6.AngII.3d>1.3),
                   `C57BL/6-AngII 28d vs. Saline 3d` = sum(Log.q.value..C57BL.6.AngII.28d>1.3)
  ) %>% melt() %>% 
  mutate(model = "AngII") %>% filter(value>0) %>% as.data.frame()

path_count_merge = rbind(path_count_ss, path_count_shr, path_count_mouse)
path_count_merge$cell_type = factor(path_count_merge$cell_type, levels = cell_order)
path_count_merge$tissue = factor(path_count_merge$tissue, levels = tissue_order)

ggplot(path_count_merge, aes(y = cell_type, x = variable, label = value)) +
  geom_point(aes(size=value/3), alpha = .5, color="red") + geom_text() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  scale_y_discrete(limits=rev) +
  scale_size_continuous(breaks = c(1, 3, 6),
                        labels = c("3", "9", "18")) +
  labs(x="", y="", size="Number of\nenriched\npathways") +
  facet_grid(cols = vars(model), 
             scales = "free", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.pathway-count.png", width=455/96, height=360/96, dpi=300)


length(intersect(merge_reshape[merge_reshape$Log.q.value..SS.SD>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                 merge_reshape[merge_reshape$Log.q.value..SS.HS.3d>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway))
# 7
# [1] "circulatory system process (GO)"                                         "regulation of plasma membrane bounded cell projection organization (GO)"
# [3] "response to hormone (GO)"                                                "negative regulation of intracellular signal transduction (GO)"          
# [5] "Proximal tubule transport (WikiPathways)"                                "sodium ion transport (GO)"                                              
# [7] "regulation of anatomical structure size (GO)"        

length(intersect(merge_reshape[merge_reshape$Log.q.value..SHR.WKY>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="PT" & merge_reshape$summary=="Yes", ]$pathway,
                 merge_reshape[merge_reshape$Log.q.value..SHR.26w>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="PT" & merge_reshape$summary=="Yes", ]$pathway))
# 10
# [1] "electron transport chain (GO)"                            "purine-containing compound metabolic process (GO)"       
# [3] "cellular homeostasis (GO)"                                "ribonucleotide metabolic process (GO)"                   
# [5] "Chemical carcinogenesis - reactive oxygen species (KEGG)" "cell junction assembly (GO)"                             
# [7] "intracellular chemical homeostasis (GO)"                  "Cardiac muscle contraction (KEGG)"                       
# [9] "cell junction organization (GO)"                          "enzyme-linked receptor protein signaling pathway (GO)"   

length(intersect(merge_reshape[merge_reshape$Log.q.value..SS.SD>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                 intersect(merge_reshape[merge_reshape$Log.q.value..SS.HS.3d>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                           merge_reshape[merge_reshape$Log.q.value..SHR.26w>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway)))
# [1] "circulatory system process (GO)"                                         "regulation of plasma membrane bounded cell projection organization (GO)"





################################################################################
### EC subtype

################################################################################
lk_ec_mouse = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds")
lk_ec_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
lk_ec_shr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")

lk_ec_mouse = subset(lk_ec_mouse, subclass_level2 %in% setdiff(unique(lk_ec_mouse$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_ss = subset(lk_ec_ss, subclass_level2 %in% setdiff(unique(lk_ec_ss$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_shr = subset(lk_ec_shr, subclass_level2 %in% setdiff(unique(lk_ec_shr$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

lk_ec_mouse$sxt = paste0(lk_ec_mouse$strain, " - ", lk_ec_mouse$treatment)
lk_ec_ss$sxt = paste0(lk_ec_ss$strain, " - ", lk_ec_ss$treatment)
lk_ec_shr$sxt = paste0(lk_ec_shr$strain, " - ", lk_ec_shr$treatment)

cell_order = sort(unique(c(lk_ec_mouse$subclass_level2, lk_ec_ss$subclass_level2, lk_ec_shr$subclass_level2)))
getPalette = colorRampPalette(brewer.pal(8, "Set2"))
lk_ec_col = getPalette(length(cell_order))
names(lk_ec_col) = cell_order

# marker_list = c("Icam1", "Vcam1", "Pecam1", "Vegfa",
#                 "Aqp1", "Scin", "Sox17", 
#                 "Unc5c", "Sema3a", "Sgcd",
#                 "Kdr", 
#                 "Rbpms", "Plscr2.1", "Nf1",
#                 "Mecom",
#                 "Plat", "Ehd3", "Nrp1", "Smad6", "Nostrin",
#                 "Mmrn1", "Tbx1", "Pkhd1l1", "Prox1",
#                 "Mt-nd4l", "Mt-nd5", "Mt-nd1", "Mt-co2")
# DotPlot(lk_ec_mouse, features = marker_list, group.by = "subclass_level2") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = unique(lk_ec_mouse$project))
# DotPlot(lk_ec_ss, features = marker_list, group.by = "subclass_level2") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = unique(lk_ec_ss$project))
# DotPlot(lk_ec_shr, features = marker_list, group.by = "subclass_level2") + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="", title = unique(lk_ec_shr$project))
# 
# DimPlot(lk_ec_mouse, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col) + labs(title = unique(lk_ec_mouse$project))
# DimPlot(lk_ec_ss, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col) + labs(title = unique(lk_ec_ss$project))
# DimPlot(lk_ec_shr, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col) + labs(title = unique(lk_ec_shr$project))
# 
# DimPlot(lk_ec_ss, label = T, group.by = "strain", reduction = "wnn.umap.harmony")
# DimPlot(lk_ec_shr, label = T, group.by = "strain", reduction = "wnn.umap.harmony")
# DimPlot(lk_ec_ss, label = T, group.by = "sxt", reduction = "wnn.umap.harmony")
# DimPlot(lk_ec_shr, label = T, group.by = "sxt", reduction = "wnn.umap.harmony")

ec_damage_genes <- c(  "Atp1b1", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Ptprj", "Slc16a12", "Slc4a4", "Slit2", # circulatory system process
                       "Arhgap24", "Igf1r", "Pkhd1", # negative regulation of intracellular signal transduction (GO)
                       "Ank3","Cacnb4","Efna5","Erbb4","Ptprd","Sdk1" # cell junction organization
)

lk_ec_mouse <- AddModuleScore(
  object = lk_ec_mouse,
  features = list(ec_damage_genes),
  name = "EC_Adaptive_Score"
)
p1 = DimPlot(lk_ec_mouse, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col, repel = T) + labs(title = unique(lk_ec_mouse$project), x="UMAP 1", y="UMAP 2") + theme(legend.position = "None")
p2 = FeaturePlot(lk_ec_mouse, features = "EC_Adaptive_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(title = "EC adaptive score", x="UMAP 1", y="UMAP 2", color="Score")
p1 + p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.mouse.EC.umap.png", width=753/96, height=341/96, dpi=300)
# RidgePlot(lk_ec_mouse, features = "EC_Deficiency_Score1", group.by = "subclass_level2")

lk_ec_ss <- AddModuleScore(
  object = lk_ec_ss,
  features = list(ec_damage_genes),
  name = "EC_Adaptive_Score"
)
p1 = DimPlot(lk_ec_ss, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col, repel = T) + labs(title = unique(lk_ec_ss$project), x="UMAP 1", y="UMAP 2") + theme(legend.position = "None")
p2 = FeaturePlot(lk_ec_ss, features = "EC_Adaptive_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(title = "EC adaptive score", x="UMAP 1", y="UMAP 2", color="Score")
p1 + p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.ss.EC.umap.png", width=753/96, height=341/96, dpi=300)
# RidgePlot(lk_ec_ss, features = "EC_Deficiency_Score1", group.by = "subclass_level2")


lk_ec_shr <- AddModuleScore(
  object = lk_ec_shr,
  features = list(ec_damage_genes),
  name = "EC_Adaptive_Score"
)
p1 = DimPlot(lk_ec_shr, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony", cols = lk_ec_col, repel = T) + labs(title = unique(lk_ec_shr$project), x="UMAP 1", y="UMAP 2") + theme(legend.position = "None")
p2 = FeaturePlot(lk_ec_shr, features = "EC_Adaptive_Score1", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(title = "EC adaptive score", x="UMAP 1", y="UMAP 2", color="Score")
p1 + p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.shr.EC.umap.png", width=753/96, height=341/96, dpi=300)
# RidgePlot(lk_ec_shr, features = "EC_Deficiency_Score1", group.by = "subclass_level2")


rbind(lk_ec_ss@meta.data, lk_ec_shr@meta.data) %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SS", "SHR"), "Hypertensive", "Normotensive")) %>%
  filter(subclass_level2 %in% c("EC_Mecom")) %>%
  dplyr::select(sxt, EC_Adaptive_Score1, hypt, strain) %>%
  mutate(sxt = factor(sxt, levels = sxt_order)) %>% 
  group_by(sxt) %>%
  mutate(mean_value = mean(EC_Adaptive_Score1)) %>%
  ggplot(aes(x = sxt, y = EC_Adaptive_Score1, fill = mean_value)) +
  geom_violin(trim = FALSE) +
  # geom_jitter(color = "black", size = 0.5, width = 0.2) +  # Add points on top of the violins
  scale_fill_gradient(low = "white", high = "red") +
  labs(y = "EC adaptive\nscore", x = "", fill = "Mean") +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.position = "right", text = element_text(size = 12)) +
  facet_nested( ~ hypt + strain, scales = "free") 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.EC.score.vln.png", width=405/96, height=275/96, dpi=300)


# lk_ec_ss@meta.data %>%
#   group_by(sxt, subclass_level2) %>%               # Group by group and cell type
#   dplyr::summarise(cell_count = n()) %>%              # Count the number of cells in each group-cell type combination
#   dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each group
#   ggplot(aes(x = sxt, y = proportion, fill = subclass_level2)) +  # Plot the results
#   geom_bar(stat = "identity", position = "fill") +  # Create a stacked bar plot with proportions
#   scale_y_continuous(labels = scales::percent) +    # Format the y-axis as percentages
#   labs(x = "Group", y = "Proportion of Cell Types", fill = "Cell Type",
#        title = "Proportion of Cell Types in Different Groups") +  # Add labels and title
#   theme_minimal()

subclass_order <- c("EC_Mecom", "Ang EC")
subclass_order <- c("EC_Mecom")
rbind(lk_ec_ss@meta.data, lk_ec_shr@meta.data) %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SS", "SHR"), "Hypertensive", "Normotensive")) %>%
  group_by(seqID2, subclass_level2, hypt, strain, treatment) %>%               # Group by group and cell type
  dplyr::summarise(cell_count = n()) %>%              # Count the number of cells in each group-cell type combination
  ungroup() %>% group_by(seqID2) %>% 
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%  # Calculate the proportion of each cell type within each group
  filter(subclass_level2 %in% subclass_order) %>%
  dplyr::mutate(subclass_level2 = factor(subclass_level2, levels = subclass_order)) %>%
  group_by(hypt, strain, treatment, subclass_level2) %>%  # Group for mean and SD calculation
  dplyr::summarise(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion)
  ) %>% 
  ggplot(aes(x = treatment, y = mean_proportion, fill = subclass_level2)) +  # Plot the results
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Create a bar plot with means
  geom_errorbar(aes(ymin = mean_proportion - sd_proportion, ymax = mean_proportion + sd_proportion),
                position = position_dodge(width = 0.8), width = 0.2) +  # Add error bars
  scale_y_continuous(labels = scales::percent) +  # Format the y-axis as percentages
  scale_fill_manual(values = lk_ec_col) +
  labs(x = "Treatment", y = "Cell proportion") +  # Add labels and title
  theme(    axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
            legend.position = "None"
  ) +
  facet_nested(subclass_level2 ~ hypt + strain, scales = "free") 

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/LK.EC.prop-bar.png", width=351/96, height=240/96, dpi=300)






lk_ec_mouse = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds")
lk_ec_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
lk_ec_shr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")

lk_ec_mouse = subset(lk_ec_mouse, subclass_level2 %in% setdiff(unique(lk_ec_mouse$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_ss = subset(lk_ec_ss, subclass_level2 %in% setdiff(unique(lk_ec_ss$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_shr = subset(lk_ec_shr, subclass_level2 %in% setdiff(unique(lk_ec_shr$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

seurat_object = lk_ec_mouse
Idents(seurat_object) = "subclass_level2"
markers = FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE)
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/mouse.LK.EC.allmarker.0.25.long.txt", sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.25,] %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                               v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/mouse.LK.EC.allmarker.0.25.wide.txt", sep="\t", row.names = F)


seurat_object = lk_ec_ss
Idents(seurat_object) = "subclass_level2"
markers = FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE)
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/rat.ss.LK.EC.allmarker.0.25.long.txt", sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.25,] %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                               v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/rat.ss.LK.EC.allmarker.0.25.wide.txt", sep="\t", row.names = F)


seurat_object = lk_ec_shr
Idents(seurat_object) = "subclass_level2"
markers = FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE)
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/rat.sp.LK.EC.allmarker.0.25.long.txt", sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.25,] %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                               v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/rat.sp.LK.EC.allmarker.0.25.wide.txt", sep="\t", row.names = F)







FeaturePlot(lk_ec_ss, features = "Mecom", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(x="UMAP 1", y="UMAP 2")
FeaturePlot(lk_ec_shr, features = "Mecom", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(x="UMAP 1", y="UMAP 2")

FeaturePlot(lk_ec_ss, features = "Ptprj", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(x="UMAP 1", y="UMAP 2")
FeaturePlot(lk_ec_shr, features = "Ptprj", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + labs(x="UMAP 1", y="UMAP 2")

FeaturePlot(lk_ec_ss, features = "Vegfa", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") +
FeaturePlot(lk_ec_shr, features = "Vegfa", cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") & labs(x="UMAP 1", y="UMAP 2")

gene ="Cgnl1"
FeaturePlot(lk_ec_ss, features = gene, cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") + 
  FeaturePlot(lk_ec_shr, features = gene, cols = c("lightblue", "red"), pt.size = 0.5, reduction = "wnn.umap.harmony") & labs(x="UMAP 1", y="UMAP 2")


################################################################################
### trajectory analysis
library(monocle3)
source("/xdisk/mliang1/qqiu/project/multiomics-CKD/src/function/seurat-wrappers.monocle3.R")
source("/xdisk/mliang1/qqiu/project/multiomics-CKD/src/function/seurat-wrappers.internal.R")


lk_ec_mouse = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds")
lk_ec_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
lk_ec_shr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")

lk_ec_mouse = subset(lk_ec_mouse, subclass_level2 %in% setdiff(unique(lk_ec_mouse$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_ss = subset(lk_ec_ss, subclass_level2 %in% setdiff(unique(lk_ec_ss$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_shr = subset(lk_ec_shr, subclass_level2 %in% setdiff(unique(lk_ec_shr$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

lk_ec_mouse$sxt = paste0(lk_ec_mouse$strain, " - ", lk_ec_mouse$treatment)
lk_ec_ss$sxt = paste0(lk_ec_ss$strain, " - ", lk_ec_ss$treatment)
lk_ec_shr$sxt = paste0(lk_ec_shr$strain, " - ", lk_ec_shr$treatment)

lk_ec_mouse@reductions[["UMAP"]] <- lk_ec_mouse@reductions[["wnn.umap.harmony"]]
lk_ec_ss@reductions[["UMAP"]] <- lk_ec_ss@reductions[["wnn.umap.harmony"]]
lk_ec_shr@reductions[["UMAP"]] <- lk_ec_shr@reductions[["wnn.umap.harmony"]]


seurat_object = lk_ec_shr
cds <- as.cell_data_set(seurat_object)
cds <- estimate_size_factors(cds)
rowData(cds) = data.frame(gene_short_name=rownames(cds))
cds = cluster_cells(cds, resolution=4e-3)
# manually assign pre-defined clsuters
cds@clusters[["UMAP"]]$clusters <- as.factor(seurat_object$subclass_level2)
plot_cells(cds, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "subclass_level2", show_trajectory_graph = FALSE)
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
plot_cells(cds,
           color_cells_by = "subclass_level2",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) %in% c("Cap EC-1")]))
plot_cells(cds,
           cell_size = 1,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60",
           trajectory_graph_segment_size = 1) +
  labs(colour="Pseudotime")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/rat.sp.LK.EC.monocle3.png", width=366/96, height=245/96, dpi=300)

gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
saveRDS(cds, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.sp.LK.EC.monocle3.rds")

track_genes_pg = graph_test(cds, neighbor_graph = "principal_graph", cores=1) #kkn
track_genes_pg = track_genes_pg[order(track_genes_pg$morans_I, decreasing = T),]
write.table(track_genes_pg, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.sp.LK.EC.trajectory.principal_graph.genes.out", sep='\t', quote=F, col.names=T, row.names=F)

track_genes_knn = graph_test(cds, neighbor_graph = "knn", cores=1) #kkn
track_genes_knn = track_genes_knn[order(track_genes_knn$morans_I, decreasing = T),]
write.table(track_genes_knn, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.sp.LK.EC.trajectory.knn.genes.out", sep='\t', quote=F, col.names=T, row.names=F)




seurat_object = lk_ec_ss
cds <- as.cell_data_set(seurat_object)
cds <- estimate_size_factors(cds)
rowData(cds) = data.frame(gene_short_name=rownames(cds))
cds = cluster_cells(cds, resolution=4e-3)
# manually assign pre-defined clsuters
cds@clusters[["UMAP"]]$clusters <- as.factor(seurat_object$subclass_level2)
plot_cells(cds, show_trajectory_graph = FALSE)
plot_cells(cds, color_cells_by = "subclass_level2", show_trajectory_graph = FALSE)
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
plot_cells(cds,
           color_cells_by = "subclass_level2",
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=FALSE)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) %in% c("Cap EC-1")]))
plot_cells(cds,
           cell_size = 1,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60",
           trajectory_graph_segment_size = 1) +
  labs(colour="Pseudotime")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/rat.ss.LK.EC.monocle3.png", width=366/96, height=245/96, dpi=300)

gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
saveRDS(cds, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.ss.LK.EC.monocle3.rds")

track_genes_pg = graph_test(cds, neighbor_graph = "principal_graph", cores=1) #kkn
track_genes_pg = track_genes_pg[order(track_genes_pg$morans_I, decreasing = T),]
write.table(track_genes_pg, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.ss.LK.EC.trajectory.principal_graph.genes.out", sep='\t', quote=F, col.names=T, row.names=F)

track_genes_knn = graph_test(cds, neighbor_graph = "knn", cores=1) #kkn
track_genes_knn = track_genes_knn[order(track_genes_knn$morans_I, decreasing = T),]
write.table(track_genes_knn, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/monocle3/rat.ss.LK.EC.trajectory.knn.genes.out", sep='\t', quote=F, col.names=T, row.names=F)




# plot_genes_in_pseudotime(cds,
#                          color_cells_by="subclass_level2",
#                          min_expr=0.5)
# 




