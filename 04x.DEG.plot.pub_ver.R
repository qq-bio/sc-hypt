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


gene_name_list = c("Npr3")
gene_name_list = c("Kcne2", "Smim11", "Smim11a", "Smim34", "Slc5a3", "Rcan1", "Mrps6", "ENSRNOG00000059381")
for(gene_name in gene_name_list){
  
  if(gene_name %in% expr_all$gene_name){
    
    # deg_use = deg_merged[  deg_merged$p_val<0.05 &
    #                          deg_merged$gene_name %in% gene_name, ]
    # if(nrow(deg_use)>0){
    #   p2=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
    #     geom_point(aes(fill = -avg_log2FC, size=-log10(p_val),
    #                    shape = p_val_adj<0.05)) +
    #     scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
    #                          high = "red", space = "Lab" ) +
    #     scale_shape_manual(values = c(23, 21)) +
    #     theme_bw() +
    #     theme(
    #       panel.grid.major.y = element_blank(),   # No horizontal grid lines
    #       legend.justification = c(1, 0.5),
    #       axis.text.y = element_text(colour = 'black'),
    #       axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    #       strip.text = element_text(colour = 'black')
    #     ) +
    #     labs(x="", y="", fill="log2(FC)", size="-log10(p-value)", shape = "adjusted.p < 0.05", title = paste(gene_name, "(DEG-wilcox)")) +
    #     facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    # }else{
    #   p2=NULL
    # }

    deg_use = expr_all[expr_all$gene_name %in% gene_name, ]
    p1=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
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
    # 
    # p = p1+p2+plot_layout(design="
    #                  AB", widths=c(1.5,1))
    print(p1)
  }
  
}


gene_name_list = c("Kcne2", "Smim11", "Smim11a", "Smim34", "Slc5a3", "Rcan1", "Mrps6", "ENSRNOG00000059381")
for(gene_name in gene_name_list){
  
  if(gene_name %in% expr_all$gene_name){
    
    deg_use = deg_merged[  deg_merged$p_val_adj<0.05 &
                             deg_merged$gene_name %in% gene_name, ]
    if(nrow(deg_use)>0){
      p2=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
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
        guides(shape = "none") + 
        labs(x="", y="", fill="log2(FC)", size="-log10(p-value)", title = paste(gene_name, "(DEG-wilcox)")) +
        facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    }else{
      p2=NULL
    }
    
    # deg_use = expr_all[expr_all$gene_name %in% gene_name, ]
    # p1=ggplot(deg_use, aes(x = treatment, y = cell_type)) +
    #   geom_point(aes(color = log2(avg_expr+1), size=pct*100)) +
    #   scale_y_discrete(limits=rev) +
    #   scale_color_gradient(low = "white", high = "red") +
    #   theme_bw() +
    #   theme(
    #     panel.grid.major.y = element_blank(),   # No horizontal grid lines
    #     legend.justification = c(1, 0.5),
    #     axis.text.y = element_text(colour = 'black'),
    #     axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    #     strip.text = element_text(colour = 'black')
    #   ) +
    #   labs(x="", y="", color="log2(expression+1)", size="Percentage", title = paste(gene_name, "(overall expression)")) +
    #   facet_nested(tissue ~ project + strain, scales = "free", space = "free")
    # 
    # # p = p1+p2+plot_layout(design="
    # #                  AB", widths=c(1.5,1))
    print(p2)
  }
  
}
