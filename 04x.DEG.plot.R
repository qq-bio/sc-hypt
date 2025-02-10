dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggupset)
library(ggh4x)
library(patchwork)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))


mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
mouse2human = unique(mouse2human[, c("Human.gene.name", "Gene.name")])




################################################################################
#### gene wise
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[!(deg_merged$tissue=="MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$tissue = factor(deg_merged$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_merged$strain = factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep='\t', header=T)
deg_strain$cell_type = factor(deg_strain$cell_type, levels = cell_order)
deg_strain$tissue = factor(deg_strain$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_strain$treatment = factor(deg_strain$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_strain$strain = factor(deg_strain$strain, levels = c("Salt-sensitive", "Spontaneous"))
levels(deg_strain$strain) = c("SS", "SHR")
deg_strain$project = "Baseline"

deg_merged = rbind(deg_merged, deg_strain)
deg_merged$project = factor(deg_merged$project, levels = c("AngII", "Salt-sensitive", "Spontaneous", "Baseline"))

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

deg_strain = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.pseudo.DEG.all.out", sep='\t', header=T)
deg_strain$cell_type = factor(deg_strain$cell_type, levels = cell_order)
deg_strain$tissue = factor(deg_strain$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
deg_strain$treatment = factor(deg_strain$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_strain$strain = factor(deg_strain$strain, levels = c("Salt-sensitive", "Spontaneous"))
levels(deg_strain$strain) = c("SS", "SHR")
deg_strain$project = "Baseline"

deg_pseudo = rbind(deg_pseudo, deg_strain)
deg_pseudo$project = factor(deg_pseudo$project, levels = c("AngII", "Salt-sensitive", "Spontaneous", "Baseline"))

# mouse2human[mouse2human$Human.gene.name=="MRPS6", ]
gene_name_list = c("Rras", "Npr1", "Ppara", "Serpinc1", "Malat1", "Npr3", "Nrip1", 
                   "Samsn1", "Mrpl23", "Sod1", "Smim11", "Smim11a", "Smim34", "Slc5a3", "Rcan1", "Mrps6", "ENSRNOG00000059381")
gene_name_list = c("Adrb1", "Afap1l2", "Nhlrc2", "Casp7")
gene_name_list = c("Adra1a", "Arhgap42", "Ctsz", "Lnpep", "Ndst2", "Tpm1") # overlap with helen & yong's list
gene_name_list = c("Asic2", "Atp1a2", "Bmpr2", "Cacna1c", "Nav2", "Nedd4l", "Pde4d", "Ppara", "Smad3", "Stk39", "Tpm1", "Dnajc13", "Fgfr2", "Foxn3", "Gpc6", "Acta2", "Scaf11") # multi strain or tissue
gene_name_list = c("Bin1", "Fgfr2", "Grb10", "Fubp1", "Pcnx", "Ube3c", "Trim33", "Gtf2ird1", 
                   "Clip2", "Foxn3", "Mef2d", "Actn4", "Scaf11", "Ckb", "Slc39a10", "Abcc8", 
                   "Amz1", "Slc15a2", "Notch4", "Klhl23", "Tmem51", "Ankh", "Myl12a", "Casq2")
gene_name_list = c("Shox2", "H19", "Mrpl23", "Tnnt3", "Igf2", "Lsp1", "Syt8", "Cd81", "Ctsd")
gene_name_list = c("Kcne2", "Slc5a3")
gene_name_list = c("Fdft1")
for(gene_name in gene_name_list){
  
  if(gene_name %in% expr_all$gene_name){
    
    deg_use = deg_merged[  ! is.na(deg_merged$p_val) &
                             deg_merged$p_val<0.05 &
                             deg_merged$gene_name %in% gene_name, ]
    if(nrow(deg_use)>0){
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
        facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    }else{
      p1=NULL
    }
    
    deg_use = deg_pseudo[  ! is.na(deg_pseudo$p_val) &
                           deg_pseudo$p_val<0.05 &
                           deg_pseudo$gene_name %in% gene_name, ]
    if(nrow(deg_use)>0){
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
        facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    }else{
      p2=NULL
    }
    
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
      facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    
    p = p3+p1+p2+plot_layout(design="
                     AB
                     AC", widths=c(1.5,1), heights=c(2,1))
    print(p)
  }
  
}


# for(gene_name in gene_name_list){
#   deg_use = expr_all[expr_all$gene_name %in% gene_name, ]
#   p3=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
#     geom_point(aes(color = log2(avg_expr+1), size=pct*100)) +
#     scale_y_discrete(limits=rev) +
#     scale_color_gradient(low = "white", high = "red") +
#     theme_bw() +
#     theme(
#       panel.grid.major.y = element_blank(),   # No horizontal grid lines
#       legend.justification = c(1, 0.5),
#       axis.text.y = element_text(colour = 'black'),
#       axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
#       strip.text = element_text(colour = 'black')
#     ) +
#     labs(x="", y="", color="log2(expression+1)", size="Percentage", title = paste(gene_name, "(overall expression)")) +
#     facet_nested(tissue ~ project + species, scales = "free", space = "free")
#   print(p3)
# }






################################################################################
### cell type wise
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

pathway_id = "Protein-protein interactions at synapses (Reactome)"; cell_type="Inhibitory neurons"
pathway_id = "Protein-protein interactions at synapses (Reactome)"; cell_type="OPC"; tissue="HYP"
pathway_id = "glutamate receptor signaling pathway (GO)"; cell_type="OPC"
pathway_id = "glutamate receptor signaling pathway (GO)"; cell_type="Microglia"
pathway_id = "NABA CORE MATRISOME (Canonical)"; cell_type="FIB"
pathway_id = "blood vessel development (GO)"; cell_type="FIB"
pathway_id = "blood vessel morphogenesis (GO)"; cell_type="FIB"
pathway_id = "receptor clustering (GO)"; cell_type="OPC"

pathway_id = "regulation of cell projection organization (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "enzyme-linked receptor protein signaling pathway (GO)"; cell_type="TAL"; tissue="LK"
pathway_id = "Cardiac muscle contraction (KEGG)"; cell_type="CM"; tissue="LV"
pathway_id = "Diabetic cardiomyopathy (KEGG)"; cell_type="CM"; tissue="LV"
pathway_id = "regulation of postsynaptic membrane potential (GO)" ; cell_type="OPC"; tissue="HYP"
pathway_id = "regulation of system process (GO)" ; cell_type="Astrocyte"; tissue="HYP"
pathway_id = "glutamate receptor signaling pathway (GO)"; cell_type="Microglia"; tissue="HYP"

### pathways enriched in three models
pathway_id = "cell junction organization (GO)"; cell_type="Microglia"; tissue="HYP"
pathway_id = "Smooth Muscle Contraction (Reactome)"; cell_type="VSMC"; tissue="MCA"
pathway_id = "enzyme-linked receptor protein signaling pathway (GO)"; cell_type="VSMC"; tissue="MCA"
pathway_id = "regulation of sodium ion transport (GO)"; cell_type="EC"; tissue="LK"
pathway_id = "NABA CORE MATRISOME (Canonical)"; cell_type="Fibroblast"; tissue="LV"
pathway_id = "regulation of anatomical structure size (GO)"; cell_type="EC"; tissue="LK"

# circulatory system process (GO)
# regulation of plasma membrane bounded cell projection organization (GO)
# response to hormone (GO)
# negative regulation of intracellular signal transduction (GO)
# regulation of anatomical structure size (GO)

pathway_id = "positive regulation of protein phosphorylation (GO)"; cell_type="Astrocyte"; tissue="HYP"
# positive regulation of phosphorus metabolic process (GO)
# Physico chemical features and toxicity associated pathways (WikiPathways)
# positive regulation of protein phosphorylation (GO)

pathway_id = "cell junction organization (GO)"; cell_type="Microglia"; tissue="HYP"
pathway_id = "cell junction organization (GO)"; cell_type="Myelinating OL"; tissue="HYP"


pathway_gene_use = pathway_gene[pathway_id][[1]]
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.5 &
                       deg_merged$cell_type==cell_type & 
                       deg_merged$tissue==tissue & 
                       deg_merged$gene_name %in% pathway_gene_use &
                       deg_merged$strain %in% strain_order, ]

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
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0(pathway_id, "\n", tissue, " - ", cell_type)) +
  facet_nested(~tissue + strain, scales = "free", space = "free")






pathway_gene_use = c("Nrxn1", "Nlgn1", "Nrxn2", "Nrxn3", "Nlgn2"); tissue="HYP"; pathway_id="NRXN"
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.25 &
                       deg_merged$tissue==tissue & 
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
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
       title = paste0(pathway_id, "\n", tissue, " - ", cell_type)) +
  facet_nested(cell_type~tissue + strain, scales = "free", space = "free")




merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", sep='\t', header=T, quote="")
merge_reshape_filter = merge_reshape[merge_reshape$model_count>1 |
                                     ((merge_reshape$Log.q.value..SS.SD > 1.3 |
                                         merge_reshape$Log.q.value..SHR.WKY > 1.3) &
                                        merge_reshape$summary=="Yes"), c("pathway", "cell_type", "tissue")]

pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/metascape.shared.gene_expr.pdf", width=10, height=10)

for(i in 1:nrow(merge_reshape_filter)){
  
  cell_type = merge_reshape_filter$cell_type[i]
  tissue = merge_reshape_filter$tissue[i]
  pathway_id = merge_reshape_filter$pathway[i]
  
  pathway_gene_use = setdiff(pathway_gene[pathway_id][[1]], NA)
  deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                         abs(deg_merged$avg_log2FC)>0.5 &
                         deg_merged$cell_type==cell_type & 
                         deg_merged$tissue==tissue & 
                         deg_merged$gene_name %in% pathway_gene_use, ]
  if(nrow(deg_use)>0){
    
    p = ggplot(deg_use, aes(x = treatment, y = gene_name)) +
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
      labs(x="", y="", color="log2(FC)", size="-log10(p-value)", 
           title = paste0(pathway_id, "\n", tissue, " - ", cell_type)) +
      facet_nested(~tissue + strain, scales = "free", space = "free")
    
    print(p)
    
  }
  
}

dev.off()






deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header = T, row.names = 1)

enrichment_analysis(deg_df)







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








################################################################################
#### vlnplot: pairwise comparison between three groups of specific genes
library(ggpubr)
vp_case1 <- function(seurat_object, cell_type, gene, file_name, test_sign, y_max){
  seurat_object_tmp = subset(seurat_object, subclass_level1 %in% cell_type)
  plot_case1 <- function(signature){
    VlnPlot(seurat_object_tmp, features = signature,
            pt.size = 1, 
            group.by = "sxt", # cols = type_col,
            y.max = y_max, # add the y-axis maximum value - otherwise p-value hidden
    ) + labs(x="", title=paste0(signature, "-", cell_type)) + 
      theme(legend.position = "None",
            axis.text.x = element_text(angle = 0, hjust = 0.5)) +
      stat_compare_means(comparisons = test_sign, label = "p.format", step.increase = 0.2, method = "wilcox")
  }
  purrr::map(gene, plot_case1) %>% cowplot::plot_grid(plotlist = .)
  # file_name <- paste0(file_name, "_r.png")
  # ggsave(file_name, width = 14, height = 8)
}

ss.LK = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds")
sp.LK = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds")
ss.HYP = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds")
sp.HYP = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds")

seurat_object = ss.HYP
seurat_object$sxt = paste(seurat_object$strain, seurat_object$treatment, sep = "-")
seurat_object$sxt = factor(seurat_object$sxt, levels = c("SD-LS", "SD-HS 3d", "SS-LS", "SS-HS 3d", "SS-HS 21d"))
comparisons <- list(c("SD-LS", "SS-LS"), c("SS-LS", "SS-HS 3d"), c("SS-LS", "SS-HS 21d"))

seurat_object = sp.HYP
seurat_object$sxt = paste(seurat_object$strain, seurat_object$treatment, sep = "-")
seurat_object$sxt = factor(seurat_object$sxt, levels = c("WKY-10w", "WKY-26w", "SHR-10w", "SHR-26w"))
comparisons <- list(c("SHR-10w", "WKY-10w"), c("SHR-10w", "SHR-26w"))



gene_list = c("Atp1b1", "Nedd4l", "Nox4", "Ptprj", "Chrm3", "Slc4a4")
gene_list = c("Igf1r", "Pkhd1", "Chrm3")
gene_list = c("Pbx1","Adgrl4","Ldb2","AABR07066700.1","Meis2","Chrm3","Inpp4b","Rbms3","Tek","ENSRNOG00000000940","Ebf1","Ptprm","Rapgef4","Plpp1","Prex2","Ptprb","Ptprg","Emcn","Ehd4","Fgf12","Unc5c","Exoc3l2","Maml2","Tcf4","Dysf","ENSRNOG00000070284","Tgfbr2","Akt3","Cped1","Arhgap31","Fam117b","Pecam1")

gene_list = c("Slc17a3")
gene_list = c("Cadm2", "Rbm6")

DotPlot(seurat_object, features = gene_list, group.by = "subclass_level1")

cell_type = c("Astrocyte")
# pdf("../DEG/CDC42.vln.pdf", height=9, width = 6)
vp_case1(seurat_object, cell_type = cell_type, gene = gene_list, test_sign = comparisons, y_max = 8)
# dev.off()












################################################################################
seurat_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds")
seurat_object$sxt = paste(seurat_object$strain, seurat_object$treatment, sep = "-")

gene_list = c(rownames(seurat_object)[grep("Il6|Tnf|Ccl2|Cxcl10|Gfap|Nos2|Ptgs2|Nfkb1|Rela|Mapk", rownames(seurat_object))])
DotPlot(seurat_object, features = gene_list, group.by = "subclass_level1", split.by = "strain") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'))




