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
library(CSCORE)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))



################################################################################
# https://stackoverflow.com/questions/71986699/plotting-heatmap-with-triangular-split-tiles-of-more-than-one-categorical-variab
# https://stackoverflow.com/questions/71147521/split-overlapping-tiles-by-facet-in-geom-tile


sxt_order = c("C57BL/6 - Saline 3d", "C57BL/6 - AngII 3d", "C57BL/6 - AngII 28d", 
              "SS - LS", "SS - HS 3d", "SS - HS 21d", "SD - LS", "SD - HS 3d",
              "SHR - 10w", "SHR - 26w", "WKY - 10w", "WKY - 26w")

sxt_col = c(rep(species_col['C57BL/6'], 3), rep(species_col['SS'], 3), rep(species_col['SD'], 2),
            rep(species_col['SHR'], 2), rep(species_col['WKY'], 2))
names(sxt_col) = sxt_order
################################################################################
### Avp expression
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.final.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.final.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.final.rds')

hyp_m$sxt = factor(paste0(hyp_m$species, " - ", hyp_m$treatment), levels = sxt_order)
hyp_ss$sxt = factor(paste0(hyp_ss$species, " - ", hyp_ss$treatment), levels = sxt_order)
hyp_shr$sxt = factor(paste0(hyp_shr$species, " - ", hyp_shr$treatment), levels = sxt_order)

neuron_list = c("Inhibitory neurons", "Excitatory neurons", "Avp+ neurons")
hyp_m %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_m$project))
hyp_ss %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_ss$project))
hyp_shr %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_shr$project))



subset_data = hyp_m %>% subset(subcluster == "Avp+ neurons")
p1 = VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + theme(legend.position = 'None') + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons\n(AngII)")
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.angii.avp.expr.png", width=356/96, height=307/96, dpi=300)

subset_data = hyp_ss %>% subset(subcluster=="Glu-15")
p2 = VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p2)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.ss.avp.expr-1.png", width=356/96, height=270/96, dpi=300)
subset_data = hyp_ss %>% subset(subcluster=="GABA-10")
p3 = VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p3)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.ss.avp.expr-2.png", width=356/96, height=270/96, dpi=300)

subset_data = hyp_shr %>% subset(subcluster=="Glu-15")
p4 = VlnPlot(subset_data, features = "Avp", group.by = "sxt", cols = sxt_col[unique(subset_data$sxt)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons\n(Spontaneous)") + theme(legend.position = 'None')
print(p4)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.shr.avp.expr.png", width=356/96, height=270/96, dpi=300)



### subcluster analysis
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.astrocyte.subcluster.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.astrocyte.subcluster.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.astrocyte.subcluster.rds')

hyp_m$sxt = factor(paste0(hyp_m$species, " - ", hyp_m$treatment), levels = sxt_order)
hyp_ss$sxt = factor(paste0(hyp_ss$species, " - ", hyp_ss$treatment), levels = sxt_order)
hyp_shr$sxt = factor(paste0(hyp_shr$species, " - ", hyp_shr$treatment), levels = sxt_order)

gene = c("Nlgn1", "Nrxn1", "Nrxn3", "Ncam1", "Ncam2", "Negr1")
gene = c("Gfap", "Gphn", "Asic2", "Adgrl3", "Pitpnc1", "Atp1a2", "Msi2", "Cables1", "Itih3") # angii
gene = c("Gfap", "Gphn", "Cdk8", "Scd2", "Gpm6a", "Map2", "Celf2", "Grm3", "Gria2", "Sash1") # up-regulated in angii at astrocyte level
gene = c("Cntnap2", "Camk1d", "Ldb2", "Kcnn2", "Tmem117", "Map2", "Glis3", "Alk", "Phactr1", "Asic2") # ss
gene = c("Gfap", "Ptprd", "Apoe", "Grik2", "Cst3", "Dlg2", "Fgf14", "Robo1", "Csmd3", "Grid2", "Kdelr3") # up-regulated in ss at astrocyte level
gene = c("Cntnap2", "Camk1d", "Ldb2", "Kcnn2", "Temem117", "Map2", "Glis3", "Alk", "Phactr1", "Asic2") # sp
hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment")
hyp_ss %>% FeaturePlot(., features = gene, split.by = "sxt")
hyp_shr %>% FeaturePlot(., features = gene, split.by = "sxt")





hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.microglia.subcluster.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.microglia.subcluster.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.microglia.subcluster.rds')

hyp_m$sxt = factor(paste0(hyp_m$species, " - ", hyp_m$treatment), levels = sxt_order)
hyp_ss$sxt = factor(paste0(hyp_ss$species, " - ", hyp_ss$treatment), levels = sxt_order)
hyp_shr$sxt = factor(paste0(hyp_shr$species, " - ", hyp_shr$treatment), levels = sxt_order)

gene = c("Marchf1", "Apba1", "Ptgds", "Cst3", "Cdk14", "B2m", "Apoe", "C1qb", "Cadm1")
hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment")
hyp_ss %>% FeaturePlot(., features = gene, split.by = "sxt")
hyp_shr %>% FeaturePlot(., features = gene, split.by = "sxt")





activated_microglia_marker = c("Mrc1", "Marchf1", "Apba1", "Ptgds", "Cst3", "Cdk14", "B2m", "Apoe", "C1qb", "Cadm1", "Cd74")
core_microglial_markers <- c("C1qa", "C1qb", "C1qc", "B2m", "Trem2", "Itgam", "Aif1", "Cd68", "Spp1", "Tnf", "Il1b", "Cx3cr1", "Cd86", "H2-Aa", "H2-Ab1")
inflammatory_genes <- c("Ptgs2", "Nos2", "Nlrp3", "Csf1r", "Fgl2")
metabolic_stress_response <- c("Hmox1", "Sod2", "Nfe2l2", "Hsp90ab1")
chemokines_cytokines <- c("Ccl2", "Ccl3", "Cxcl10", "Ccl5")
phagocytosis_synaptic_pruning <- c("Axl", "Mertk", "Gpnmb", "Laptm5", "Ctss")
proliferation_cell_cycle <- c("Cd14", "Csf1", "Csf2", "Cd44")
from_paper = c("Cd14", "Cd86", "Fcgr2a", "Fcgr3a", "Havcr2", "Il18", "Cd163", "Apoe", "Gnas", "Csf1r", "Cx3cr1", "Tmem119")
mg=c("Microglia", "Activated microglia")

gene_list = c(activated_microglia_marker,
              core_microglial_markers, inflammatory_genes, metabolic_stress_response,
              chemokines_cytokines, phagocytosis_synaptic_pruning, proliferation_cell_cycle,
              from_paper)
gene_list = c("C1qa", "C1qb", "C1qc", "B2m")
p = deg_merged %>%
  filter(gene_name %in% gene_list,
         cell_type %in% mg) %>%
  mutate(gene_name = factor(gene_name, levels=rev(unique(gene_list))),
         p_val = ifelse(p_val==0, 1.6e-303, p_val)) %>%
  ggplot(., aes(x = treatment, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.justification = c(1, 0.5),
    legend.direction = "horizontal",
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x="", y="", color="log2(FC)", size="-log10(p-value)") +
  facet_nested(~ cell_type + strain, scales = "free", space = "fixed")

print(p)


### expr profile across samples
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

expr_all = expr_all[expr_all$tissue=="HYP", ]

pathway_gene_use = c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3")
deg_use = expr_all[expr_all$gene_name %in% pathway_gene_use &
                   expr_all$cell_type %in% mg, ]
ggplot(deg_use, aes(x = treatment, y = cell_type)) +
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
  labs(x="", y="", color="log2(expression+1)", size="Percentage", title = "Overall expression") +
  facet_nested(gene_name + tissue ~ project + strain, scales = "free", space = "free")



################################################################################
### top DEGs in each cell type
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[deg_merged$tissue=="HYP", ]

deg_merged = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.25, ]
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_merged$strain = factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))

cell_type=c("Astrocyte", "Microglia", "Myelinating OL", "OPC")
cell_type=c("Inhibitory neuron", "Excitatory neuron", "Activated microglia")
for(ci in cell_type){
  
  for(pi in c("AngII", "Salt-sensitive", "Spontaneous")){
    
    top_genes <- deg_merged %>%
      filter(cell_type == ci, project == pi, 
             strain %in%  c("C57BL/6", "SS", "SHR")) %>%
      group_by(gene_name) %>%
      slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%  # Keep the row with the highest log2FC per gene
      arrange(desc(avg_log2FC)) %>%  # Sort by log2FC in descending order
      ungroup() %>%
      slice(c(1:10, (n() - 9):n()))
    
    p = deg_merged %>%
      filter(gene_name %in% top_genes$gene_name,
             cell_type == ci, project==pi) %>%
      mutate(gene_name = factor(gene_name, levels=rev(unique(top_genes$gene_name))),
             p_val = ifelse(p_val==0, 1.6e-303, p_val)) %>%
      ggplot(., aes(x = treatment, y = gene_name)) +
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
      labs(x="", y="", color="log2(FC)", size="-log10(p-value)") +
      facet_nested(~ cell_type + strain, scales = "free", space = "free")
    
    print(p)
    
  }
  
}


### consensus + bulk
bulk = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/bulk/DEG_bulk.csv", header = T, sep = ",")
bulk[bulk$project=="Spontaneous ", ]$project = "Spontaneous"
deg_bulk <- deg_merged %>%
  group_by(gene_name, strain, treatment) %>%
  mutate(
    deg_num = n(),  # Count the total number of DEGs in each group
    up_num = sum(avg_log2FC < 0),  # Count how many genes are upregulated (log2FC > 0)
    down_num = sum(avg_log2FC > 0),  # Count how many genes are downregulated (log2FC < 0)
    fc_weighted = sum(-avg_log2FC * (control_size + treatment_size))
  ) %>%
  ungroup() %>%  # Ungroup to apply row_number across the entire dataset
  select(gene_name, strain, treatment, deg_num, up_num, down_num, fc_weighted) %>%  # Ensure "treatment" is included
  distinct() %>%
  arrange(desc(deg_num)) %>%
  mutate(index = row_number()) %>%
  left_join(bulk[bulk$tissue=="HYP", ], by = c("gene_name", "strain", "treatment"))

deg_bulk %>% 
  filter(! is.na(deg_bulk$padj), 
         deg_bulk$padj<0.05,
         strain %in% c("C57BL/6", "SS", "SHR")) %>% 
  # filter(fc_weighted*log2FoldChange>0) %>% 
  filter(((down_num>up_num) & log2FoldChange<0) |
         ((down_num<up_num) & log2FoldChange>0)) %>%
  View()


### bulk-validated candidates
gene = c("Snhg11", "Meg3", "Csmd1", "Cacna1c", "Csmd2", "Nav1", "Pkd1", "Ptprd") # neuron markers, down-regulated in Angii -> decrease prop of neuron
gene = c("Zbtb16", "Lars2", "Gpm6b", "Sept7")
hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())
# gene = c("Mapt", "Kcnc3", "Adgrb1", "Pip5k1c", "Bin1")
# hyp_ss %>% FeaturePlot(., features = gene, split.by = "sxt") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())
gene = c("ENSRNOG00000065867")
hyp_shr %>% FeaturePlot(., features = gene, split.by = "sxt") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())


### merged bulk-list
selected_gene = c("Snhg11", "Meg3", "Sntg1", "Zbtb16", "Csmd1", "ENSRNOG00000065867",
                  "Cacna1c", "Lars2", "Gpm6b", "Ncam1"#, "Csmd2", "Nav1", "Sept7"
                  )
p1 = deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_gene) %>%
  mutate(gene_name = factor(gene_name, levels=selected_gene)) %>%
  ggplot(., aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val_adj))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical", 
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "snRNA-seq", 
       x="Cell type-wise log2(FC)", y="", 
       color="log2(FC)", size="-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1)) 

p2 = bulk %>%
  filter(tissue=="HYP", strain %in% c("C57BL/6", "SS", "SHR"),
         # padj<0.05,
         gene_name %in% selected_gene) %>%
  mutate(gene_name = factor(gene_name, levels=selected_gene),
         strain = factor(strain, levels=c("C57BL/6", "SS", "SD", "SHR", "WKY"))) %>%
  ggplot(., aes(x = log2FoldChange, y = gene_name)) +
  geom_vline(xintercept = 0, colour="black") +
  geom_point(shape = 21, aes(fill = log2FoldChange, size=-log10(padj))) +
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical", 
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "Bulk RNA-seq", 
       x="log2(FC) in bulk", y="", 
       fill="log2(FC)", size="-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(color = guide_legend(nrow = 1), 
         size = guide_legend(nrow = 1)) 

# width=678&height=463
p1+p2+ plot_layout(widths = c(1,1))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.sn_bulk.consist_expr.png", width=832/96, height=356/96, dpi=300)





selected_gene = c("Fth1", "Rora", "Apoe", "Ptgds", "Csmd1", "App", "Ncam1",
  "Hsp90ab1", "Zbtb16", "Lars2")
p = deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  group_by(gene_name) %>%
  mutate(deg_num = n()) %>%
  ungroup() %>%
  dplyr::select(gene_name, deg_num) %>%
  distinct() %>%
  arrange(desc(deg_num)) %>%
  mutate(
    index = row_number(),
    highlight = ifelse(gene_name %in% selected_gene, "Highlighted", "Normal")
  ) %>%
  # Reorder to make highlighted genes appear last (on top in the plot)
  arrange(desc(highlight)) %>%
  ggplot(aes(x = deg_num, y = index, color = highlight)) +
  geom_point(size = 3) +  # Plot points
  scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) +  # Color the highlighted genes red
  scale_y_continuous(trans = "reverse") +  # Reverse the index ordering
  ggrepel::geom_label_repel(
    aes(label = ifelse(gene_name %in% selected_gene, gene_name, "")), 
    size = 4, 
    max.overlaps = Inf,  # Allow infinite overlaps
    force = 10,  # Stronger force to repel labels
    box.padding = 1,  # Increase box padding to give more room
    point.padding = 0.5,  # Increase point padding
    color = "black"
  ) +  # Add non-repel text annotations for selected genes
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.position = "none",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Number of DEG occurrence", y = "Rank of genes by DEG occurrence")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2a.HYP.deg.dot.png", width=353/96, height=287/96, dpi=300)


deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_gene) %>%
  mutate(gene_name = factor(gene_name, levels=selected_gene)) %>%
  ggplot(., aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val_adj))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.justification = c(1, 0.5),
    legend.position = "bottom",
    # legend.direction = "vertical",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "", 
       x="Cell type-wise log2(FC)", y="", 
       color="log2(FC)", size="-log10(p-adj)") +
  facet_nested(~ strain, scales = "free")



# top_30_genes <- deg_merged %>%
#   filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
#   group_by(gene_name) %>%
#   mutate(deg_num = n()) %>%
#   ungroup() %>%  # Ungroup to apply row_number across the entire dataset
#   select(gene_name, deg_num) %>%
#   distinct() %>%
#   arrange(desc(deg_num)) %>%
#   mutate(index = row_number()) %>%
#   slice(1:100)

deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% top_30_genes) %>%
  mutate(gene_name = factor(gene_name, levels=top_30_genes$gene_name)) %>%
  ggplot(., aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val_adj))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.justification = c(1, 0.5),
    legend.position = "bottom",
    legend.direction = "vertical",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "snRNA-seq", 
       x="Cell type-wise log2(FC)", y="", 
       color="log2(FC)", size="-log10(p-adj)") +
  facet_nested(~ strain, scales = "free")


gene = c("Gphn", "Hsp90ab1", "Psap") # pan-cell type in angii
gene = c("Ntm", "Nkain2", "Robo1", "Rora", "Cntnap2", "Kcnd2")
hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())
hyp_ss %>% FeaturePlot(., features = gene, split.by = "sxt") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())
hyp_shr %>% FeaturePlot(., features = gene, split.by = "sxt") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())

gene = c("Gphn", "Hsp90ab1", "Psap")
gene = c("Hsp90ab1")
# width=577&height=523
p1=hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment") & labs(x="", y="") & theme(axis.text = element_blank(), axis.ticks = element_blank())

gene = c("Hsp90ab1")
p1=hyp_m %>% FeaturePlot(., features = gene, split.by = "treatment", order = T) & labs(x="", y="") & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"))
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2b.HYP.angii.Hsp90ab1.umap.png", width=498/96, height=155/96, dpi=300)

p2=hyp_ss %>% subset(species=="SS") %>% FeaturePlot(., features = gene, split.by = "sxt", order = T) & labs(x="", y="") & 
  theme(axis.text = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.1, 0, 0, 0), "cm"))
print(p2)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2b.HYP.ss.Hsp90ab1.umap.png", width=498/96, height=155/96, dpi=300)



deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  group_by(gene_name) %>%
  mutate(deg_num = n(unique(strain))) %>%
  ungroup() %>%  # Ungroup to apply row_number across the entire dataset
  select(gene_name, deg_num) %>%
  distinct() %>%  # Use distinct to get unique gene names
  arrange(desc(deg_num)) %>%
  mutate(index = row_number()) %>%  # Now row_number() will apply correctly across all rows
  ggplot(aes(x = deg_num, y = index)) +
  geom_point() +
  scale_y_continuous(trans = "reverse") +  # Reverse the index ordering
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Number of\nDEG occurrence", y = "Rank of genes by DEG occurrence")


top_30_genes <- deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  group_by(gene_name) %>%
  mutate(deg_num = n()) %>%
  ungroup() %>%  # Ungroup to apply row_number across the entire dataset
  select(gene_name, deg_num) %>%
  distinct() %>%
  arrange(desc(deg_num)) %>%
  mutate(index = row_number()) %>%
  slice(1:30)

deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% top_30_genes$gene_name) %>%
  mutate(gene_name = factor(gene_name, levels=top_30_genes$gene_name)) %>%
  ggplot(., aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size=-log10(p_val))) +
  scale_y_discrete(limits=rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" ) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    legend.justification = c(1, 0.5),
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "Top30 genes with\nmost DEG occurrence", 
       x="Cell type-wise log2(FC)\nacross models", y="", 
       color="log2(FC)", size="-log10(p-value)")


### expr profile across samples
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

expr_all = expr_all[expr_all$tissue=="HYP", ]

pathway_gene_use = c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3")
deg_use = expr_all[expr_all$gene_name %in% pathway_gene_use, ]
ggplot(deg_use, aes(x = treatment, y = cell_type)) +
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
  labs(x="", y="", color="log2(expression+1)", size="Percentage", title = "Overall expression") +
  facet_nested(gene_name + tissue ~ project + strain, scales = "free", space = "free")













################################################################################
### pathway enrichment result
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = merge_reshape[merge_reshape$tissue=="HYP",]

path_count = merge_reshape %>% filter(category != "Other", summary == "Yes") %>% 
  group_by(tissue, cell_type, category) %>%
  dplyr::summarise(count = n()) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)
path_count$category = factor(path_count$category, levels = c("Only AngII", "Only SS", "Only SHR", "AngII & SS",
                                                             "AngII & SHR", "SS & SHR", "AngII & SS & SHR"))
path_count$cell_type_mod = factor(paste0(path_count$cell_type, " •"), levels = paste0(cell_order, " ●"))
cell_col_mod = cell_col
names(cell_col_mod) = paste0(names(cell_col), " ●")
ggplot(path_count, aes(y = cell_type, x = category, label = count)) +
  geom_rect(aes(xmin = 4 - 0.5, xmax = 7 + 0.5, ymin = 0, ymax = Inf), 
            fill = "#FDEDEC", color = NA) + 
  geom_point(aes(size = count / 3), alpha = 0.5, color = "red") +
  geom_text() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  scale_y_discrete(limits = rev) +
  scale_size_continuous(breaks = c(1, 3, 6),
                        labels = c("3", "9", "18")) +
  labs(x = "", y = "", size = "Number of\npathways") +
  facet_grid(rows = vars(tissue), scales = "free", space = "free")



sxt_rename = c("Log.q.value..C57BL.6.AngII.3d"="C57BL/6 - AngII 3d", 
               "Log.q.value..C57BL.6.AngII.28d"="C57BL/6 - AngII 28d", 
               "Log.q.value..SS.HS.3d"="SS - HS 3d", 
               "Log.q.value..SS.HS.21d"="SS - HS 21d",
               "Log.q.value..SHR.26w"="SHR - 26w")

melted_data <- merge_reshape %>% filter(category != "Other", summary == "Yes") %>% 
  filter(model_count>1) %>%
  melt(id.vars = c('pathway', 'cell_type'),
                    measure.vars = grep("^Log.q.value", names(merge_reshape), value = TRUE),
                    variable.name = 'condition', 
                    value.name = 'log_q_value') %>% 
  filter(condition %in% names(sxt_rename)) %>% 
  mutate(condition = sxt_rename[as.character(condition)],
         log_q_value = ifelse(log_q_value<1.3, 0, log_q_value))
melted_data$pathway <- fct_reorder(melted_data$pathway, melted_data$cell_type)
melted_data$cell_type <- fct_inorder(melted_data$cell_type)  # Ensure cell type is in order
melted_data$condition <- factor(melted_data$condition, levels=c("C57BL/6 - AngII 3d", "C57BL/6 - AngII 28d", 
                                                               "SS - HS 3d", "SS - HS 21d", "SHR - 26w"))

group_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/HYP.multi_model.pathway.group.out", header=F, sep='\t')
group_info = unique(group_info)
rownames(group_info) = group_info[,1]
melted_data$group = group_info[as.character(melted_data$pathway), 2]

p_data <- ggplot(melted_data[melted_data$log_q_value>1.3,], aes(x = condition, y = pathway, fill = log_q_value)) +
  geom_tile(shape=21, color="black") +
  scale_y_discrete(limits=rev) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        # legend.text = element_text(colour = 'black'),
        plot.margin = margin(0, 0, 0, -2, "cm")) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "", fill = "-log(q-val)", title="Shared pathways in at least two hypertension models") +
  facet_grid(cols = vars(cell_type), scales = "free", space = "free")

# Second plot: vertical annotation for pathway groups
p_y <- ggplot(melted_data) +
  geom_col(aes(y = pathway, x = 1, fill = group), width = 0.8) +  # Constrain x to 1 and limit the width
  scale_x_continuous(limits = c(0, 1)) +  # Limit x-axis to 1
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
  scale_y_discrete(limits=rev) +
  # theme_minimal() + 
  labs(fill="Category") +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        axis.text.y = element_text(colour = 'black'),
        axis.title.y = element_blank(),
        axis.ticks.length.y = unit(0, "pt"),
        panel.background = element_blank(), 
        plot.background = element_blank(),
        legend.position = "right",
        plot.margin = margin(0, -2, 0, 0, "cm"))

# Combine the two plots
p_y + p_data + plot_layout(widths = c(0.03, 1), 
                                           guides = "collect")  # Adjust widths if needed
# width=1325&height=548
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2d.pathway.heatmap.png", width=1325/96, height=548/96, dpi=300)



################################################################################
### gene expr change of selected pathways
pathway_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", header = T, sep = '\t', quote = "")

mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
rat2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")
ani2human = unique(rbind(mouse2human[, c("Human.gene.name", "Gene.name")],
                   rat2human[, c("Human.gene.name", "Gene.name")]))
ani2human = ani2human[rowSums(ani2human=="")==0,]

pathway_gene = pathway_all %>% filter(str_detect(GroupID, "Member")) %>%
  separate_rows(Symbols, sep = ",") %>% distinct() %>%
  left_join(ani2human, by = c("Symbols" = "Human.gene.name"), relationship = "many-to-many") %>%
  separate_rows(Gene.name, sep = ",") %>% distinct() %>%
  group_by(pathway) %>% dplyr::summarise(Symbols_h = list(unique(Symbols)), Symbols_m = list(unique(Gene.name)))

pathway_gene <- setNames(pathway_gene$Symbols_m, pathway_gene$pathway)


deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged = deg_merged[deg_merged$tissue=="HYP", ]

deg_merged = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.25, ]
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$treatment = factor(deg_merged$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
deg_merged$strain = factor(deg_merged$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))


pathway_id = "Neuronal System (Reactome)"; cell_type=c("Microglia", "Myelinating OL")
pathway_id = "Protein-protein interactions at synapses (Reactome)"; cell_type=c("Microglia", "Myelinating OL")
pathway_id = "cell junction assembly (GO)"; cell_type=c("Microglia", "Myelinating OL")
pathway_id = "cell junction organization (GO)"; cell_type=c("Microglia", "Myelinating OL")
pathway_id = "regulation of system process (GO)"; cell_type=c("OPC", "Myelinating OL")

pathway_id = "dendrite development (GO)"; cell_type=c("Astrocyte")
pathway_id = "import into cell (GO)"; cell_type=c("Astrocyte")
pathway_id = "monoatomic cation homeostasis (GO)"; cell_type=c("Astrocyte")
pathway_id = "neuromuscular process (GO)"; cell_type=c("Astrocyte")

pathway_gene_use = pathway_gene[pathway_id][[1]]
pathway_gene_use = c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3")
deg_use = deg_merged[deg_merged$p_val_adj<0.05 & 
                       abs(deg_merged$avg_log2FC)>0.25 &
                       deg_merged$cell_type %in% cell_type &
                       deg_merged$gene_name %in% pathway_gene_use, ]

deg_use = deg_use %>% group_by(gene_name, cell_type) %>%
  mutate(strain_count = n_distinct(strain)) %>%
  filter(strain_count >= 2)

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
       title = pathway_id) +
  facet_nested(~cell_type + strain, scales = "free", space = "free")












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
# mouse_hyp_list = list(saline3d=e$`mouse.HYP_C57BL.6_Saline 3d`, angii3d=e$`mouse.HYP_C57BL.6_AngII 3d`, angii28d=e$`mouse.HYP_C57BL.6_AngII 28d`)
# mouse_hyp = mergeCellChat(mouse_hyp_list, add.names = names(mouse_hyp_list))
# 
# ss_hyp_list = list(ss_ls=e$rat.ss.HYP_SS_LS, ss_hs3d=e$`rat.ss.HYP_SS_HS 3d`, ss_hs21d=e$`rat.ss.HYP_SS_HS 21d`, sd_ls=e$rat.ss.HYP_SD_LS, sd_hs3d=e$`rat.ss.HYP_SD_HS 3d`)
# ss_hyp = mergeCellChat(ss_hyp_list, add.names = names(ss_hyp_list))
# 
# sp_hyp_list = list(shr_10w=e$rat.sp.HYP_SHR_10w, shr_26w=e$rat.sp.HYP_SHR_26w, wky_10w=e$rat.sp.HYP_WKY_10w, wky_26w=e$rat.sp.HYP_WKY_26w)
# sp_hyp = mergeCellChat(sp_hyp_list, add.names = names(sp_hyp_list))


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
  outfile = paste0("fig2d.", i, "cellchat.png")
  
  png(filename = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/", outfile), 
      width = 786/96*300, height = 291/96*300, res = 300)  # Set resolution and size
  
  par(mfrow=c(1, 3))  # 1 row, 3 columns
  
  for(j in c("NCAM", "NEGR", "NRXN") ){
    pathways.show <- j
    
    # Plot each signaling pathway
    netVisual_aggregate(cellchat, signaling = pathways.show, 
                        # signaling.name = NULL,
                        layout = "circle", title.space = 1,
                        remove.isolate = T,
                        weight.scale = F,  # Scale edge weights for comparability
                        edge.weight.max = 1,  # Set maximum edge weight
                        edge.width.max = 8  # Keep consistent maximum edge width
    )
    
    # Add a title to each plot
    title(main = paste0(title[i], " - ", j))
  }
  
  # ggsave(paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/",outfile), width=1248/96, height=548/96, dpi=300)
  dev.off()
  par(mfrow=c(1, 1))
}




for( i in c("saline3d", "angii3d", "angii28d", "ss_ls", "ss_hs3d", "ss_hs21d", "sd_ls", "sd_hs3d", "shr_10w", "shr_26w", "wky_10w", "wky_26w") ){
  cellchat = get(i)
  
  # Set up the layout to display three plots in one row
  par(mfrow=c(1, 3))  # 1 row, 3 columns
  
  for(j in c("NCAM", "NEGR", "NRXN") ){
    pathways.show <- j
    
    # Plot each signaling pathway
    netVisual_aggregate(cellchat, signaling = pathways.show, 
                        layout = "circle", title.space = 1,
                        weight.scale = F,  # Scale edge weights for comparability
                        edge.weight.max = 1,  # Set maximum edge weight
                        edge.width.max = 8  # Keep consistent maximum edge width
    )
    
    # Add a title to each plot
    title(main = paste0(i, " - ", j))
  }
  
  # Reset layout to 1 plot per page after each set of 3
  par(mfrow=c(1, 1))
}


################################################################################
### all changes
cc_df_use = cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
pathway_selected = c("NRXN", "NCAM", "NEGR", "PTN") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & cc_df_use$tissue=="HYP", ]) +
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

# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.cellchat.heatmap.png", width=494/96, height=227/96, dpi=300)


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
library(scCorr)
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.final.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.final.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.final.rds')

seurat_object = hyp_m
seurat_use_ = subset(seurat_object, subcluster=="Avp+ neuron", species %in% c("C57BL/6", "SS", "SHR"))

seurat_object = hyp_ss
seurat_use = subset(seurat_object, subcluster %in% c("Glu-15", "Gaba-10"), species %in% c("C57BL/6", "SS", "SHR"))

seurat_object = hyp_shr
seurat_use = subset(seurat_object, subcluster %in% c("Glu-15", "Gaba-5"), species %in% c("C57BL/6", "SS", "SHR"))

# VlnPlot(seurat_use, features = c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3", "Avp"), group.by = "treatment")

seurat_dat = as.matrix(seurat_use@assays$RNA@data)
seurat_dat = seurat_dat[rowSums(seurat_dat==0)<ncol(seurat_dat)*0.9, ]
seurat_dat = t(seurat_dat)

cor_merged = c()
# 
for( i in colnames(seurat_dat) ){
  cor.tmp = cor.test(seurat_dat[, "Avp"], seurat_dat[,i])
  cor_merged = rbind(cor_merged, c(i, cor.tmp$estimate, cor.tmp$p.value))
}

cor_merged = as.data.frame(cor_merged)
cor_merged$cor = as.numeric(cor_merged$cor)
cor_merged$V3 = as.numeric(cor_merged$V3)
cor_merged[order(cor_merged$V3),]

cor_merged[cor_merged$V1 %in% c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3" ), ]

for( i in c("Ncam1", "Ncam2", "Negr1", "Nrxn1", "Nlgn1", "Nrxn3" ) ){
  p = ggplot(as.data.frame(seurat_dat), aes(x = Avp, y = Nlgn1)) +
    geom_point()
  print(p)
}

c_gene = c_list(seurat_dat)

# V1                  cor                   V3
# 575  Nlgn1    0.192670990360377 4.45287918035042e-08
# 741  Negr1   -0.148701036871298 2.58889163747576e-05
# 2009 Ncam1 -0.00777541335673668    0.826841857776606
# 2796 Nrxn3   -0.212168820987683 1.56086872892559e-09
# 3420 Ncam2   0.0134994324953393    0.704090088643402
# 3625 Nrxn1 -0.00259318525502355    0.941841296567928



gtex = read.table("/xdisk/mliang1/qqiu/data/GTEx/bulk_RNA-seq/gene_tpm_2017-06-05_v8_brain_hypothalamus.gct", header=T, sep="\t", skip = 2)

# gtex_use = gtex[gtex$Description %in% c("NCAM1", "NCAM2", "NEGR1", "NRXN1", "NLGN1", "NRXN3", "AVP"), -c(1, 2)]
gtex_use = gtex[rowSums(gtex==0)<(ncol(gtex)-3)*0.9, ]
gene_name = gtex_use$Description
# gtex_avp = gtex_use[gtex_use$Description=="AVP", -c(1:3)]
gtex_use = t(gtex_use[, -c(1:3)])
gtex_use = log2(as.matrix(gtex_use) + 1)

cor_merged = c()
for( i in 1:ncol(gtex_use) ){
  cor.tmp = cor.test(gtex_use[, gene_name=="AVP"], gtex_use[,i])
  cor_merged = rbind(cor_merged, c(i, cor.tmp$estimate, cor.tmp$p.value))
  
}
cor_merged = as.data.frame(cor_merged)
cor_merged$gene_name = gene_name
cor_merged$order = order(-1*abs(cor_merged$cor))
cor_merged$fdr = p.adjust(cor_merged$V3, method = "BH")
cor_merged[cor_merged$gene_name %in% c("NCAM1", "NCAM2", "NEGR1", "NRXN1", "NLGN1", "NRXN3"), ]
cor_merged[cor_merged$gene_name %in% c("EP300"), ]
cor_merged[cor_merged$gene_name %in% c("EBF3"), ]

cor_merged = c()
for( i in c("NCAM1", "NCAM2", "NEGR1", "NRXN1", "NLGN1", "NRXN3") ){
  p = ggplot(gtex_use, aes(x = AVP, y = i)) +
    geom_point()
  print(p)
  cor.tmp = cor.test(gtex_use[, "AVP"], gtex_use[,i])
  cor_merged = rbind(cor_merged, c(i, cor.tmp$estimate, cor.tmp$p.value))
  
}


library(GSVA)
library(ggpubr)

colnames(gtex_use) = gene_name
# gene_sets = list(com=c("NCAM1", "NCAM2", "NEGR1", "NRXN1", "NLGN1", "NRXN3"))
gene_sets = list(com=c("NCAM1", "NCAM2", "FGFR1", "L1CAM", "NEGR1", 
                       "NRXN1", "NRXN2", "NRXN3",
                       "NLGN1", "NLGN2", "NLGN3"))
gsva_scores <- gsva(t(gtex_use), gene_sets, method = "ssgsea")

# Correlate with AVP expression
avp_expr <- gtex_use[, "AVP"]  # AVP expression values
correlation <- cor.test(as.numeric(gsva_scores), as.numeric(avp_expr))  # Correlation between AVP and gene set scores

score_avp = data.frame(score = as.numeric(gsva_scores),
                       avp = as.numeric(avp_expr))

ggplot(score_avp, aes(x = avp, y = score)) +
  geom_point() +  # Plot the points
  geom_smooth(method = "lm", color = "red", se = FALSE) +  # Add linear regression line
  annotate("text", x = 5, y = 3, label = paste0("R = ", round(correlation$estimate, 2), 
                                                ", p = ", format(correlation$p.value, scientific = TRUE)),
           size = 5) +  # Manually annotate the Pearson correlation coefficient and p-value
  theme_classic() +
  labs(x = "AVP expression", y = "Neuronal connection score", title = "Human hypothalamus bulk RNA-seq")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2h.avp_neu_score.cor.png", width=348/96, height=253/96, dpi=300)







hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.final.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.final.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.final.rds')

mhyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/mhyp_all_regulon_AUC.csv", header = TRUE)
sshyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_ssall_regulon_AUC.csv", header = TRUE)
shrhyp_auc <- read.csv("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_shrall_regulon_AUC.csv", header = TRUE)

rownames(mhyp_auc) <- mhyp_auc[,1]
mhyp_auc <- mhyp_auc[,-1]
colnames(mhyp_auc) <- gsub("\\.", "-", colnames(mhyp_auc))
all(colnames(mhyp_auc)==colnames(hyp_m))
avp_neurons = colnames(hyp_m)[hyp_m$new.cluster.ids_umap=="Avp+ neurons"]
avp_tf = data.frame(Regulon_activity = as.numeric(mhyp_auc["Ep300(+)", avp_neurons]),
                    Regulon = as.numeric(hyp_m@assays$RNA@data["Ep300", avp_neurons]),
                    Avp = as.numeric(hyp_m@assays$RNA@data["Avp", avp_neurons]),
                    treatment = factor(hyp_m$treatment[avp_neurons], 
                                       levels = c("Saline 3d", "AngII 3d", "AngII 28d")),
                    sxt = paste0("C57BL/6 - ", hyp_m$treatment[avp_neurons]))

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = T, scale = "width", fill = species_col['C57BL/6']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons\n(AngII)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.angii.ep300.activity.png", width=365/96, height=189/96, dpi=300)

ggplot(avp_tf, aes(x = treatment, y = Regulon)) +
  geom_violin(trim = T, scale = "width", fill = species_col['C57BL/6']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300\nexpression", x = "", title = "Avp+ neurons\n(AngII)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.angii.ep300.expr.png", width=365/96, height=189/96, dpi=300)



rownames(sshyp_auc) <- sshyp_auc[,1]
sshyp_auc <- sshyp_auc[,-1]
colnames(sshyp_auc) <- gsub("\\.", "-", colnames(sshyp_auc))
all(colnames(sshyp_auc)==colnames(hyp_ss))
common_columns <- intersect(colnames(sshyp_auc), colnames(hyp_ss))
rhyp_auc_filtered <- sshyp_auc[, common_columns]
hyp_rsp_filtered <- subset(hyp_ss, cells = common_columns)
sshyp_auc <- rhyp_auc_filtered
hyp_ss <- hyp_rsp_filtered
all(colnames(sshyp_auc)==colnames(hyp_ss))
avp_neurons = colnames(hyp_ss)[hyp_ss$subcluster=="Glu-15"]
avp_tf = data.frame(Regulon_activity = as.numeric(sshyp_auc["Ep300(+),", avp_neurons]),
                    Regulon = as.numeric(hyp_ss@assays$RNA@data["Ep300", avp_neurons]),
                    Avp = as.numeric(hyp_ss@assays$RNA@data["Avp", avp_neurons]),
                    treatment = factor(hyp_ss$treatment[avp_neurons], 
                                       levels = c("LS","HS 3d","HS 21d")),
                    sxt = paste0("SS - ", hyp_ss$treatment[avp_neurons]))

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SS']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.activity-1.png", width=365/96, height=189/96, dpi=300)

ggplot(avp_tf, aes(x = treatment, y = Regulon)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SS']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300\nexpression", x = "", title =  "Avp+ neurons, glutamatergic\n(Salt-sensitive)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.expr-1.png", width=365/96, height=189/96, dpi=300)



avp_neurons = colnames(hyp_ss)[hyp_ss$subcluster=="GABA-10"]
avp_tf = data.frame(Regulon_activity = as.numeric(sshyp_auc["Ep300(+),", avp_neurons]),
                    Regulon = as.numeric(hyp_ss@assays$RNA@data["Ep300", avp_neurons]),
                    Avp = as.numeric(hyp_ss@assays$RNA@data["Avp", avp_neurons]),
                    treatment = factor(hyp_ss$treatment[avp_neurons], 
                                       levels = c("LS","HS 3d","HS 21d")),
                    sxt = paste0("SS - ", hyp_ss$treatment[avp_neurons]))

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SS']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.activity-2.png", width=365/96, height=189/96, dpi=300)

ggplot(avp_tf, aes(x = treatment, y = Regulon)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SS']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300\nexpression", x = "", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.ss.ep300.expr-2.png", width=365/96, height=189/96, dpi=300)




rownames(shrhyp_auc) <- shrhyp_auc[,1]
shrhyp_auc <- shrhyp_auc[,-1]
colnames(shrhyp_auc) <- gsub("\\.", "-", colnames(shrhyp_auc))
all(colnames(shrhyp_auc)==colnames(hyp_shr))
common_columns <- intersect(colnames(shrhyp_auc), colnames(hyp_shr))
rhyp_auc_filtered <- shrhyp_auc[, common_columns]
hyp_rsp_filtered <- subset(hyp_shr, cells = common_columns)
shrhyp_auc <- rhyp_auc_filtered
hyp_shr <- hyp_rsp_filtered
all(colnames(shrhyp_auc)==colnames(hyp_shr))
avp_neurons = colnames(hyp_shr)[hyp_shr$subcluster=="Glu-15"]
avp_tf = data.frame(Regulon_activity = as.numeric(shrhyp_auc["Ep300(+),", avp_neurons]),
                    Regulon = as.numeric(hyp_shr@assays$RNA@data["Ep300", avp_neurons]),
                    Avp = as.numeric(hyp_shr@assays$RNA@data["Avp", avp_neurons]),
                    treatment = factor(hyp_shr$treatment[avp_neurons], 
                                       levels = c("10w","26w")),
                    sxt = paste0("SHR - ", hyp_shr$treatment[avp_neurons]))

ggplot(avp_tf, aes(x = treatment, y = Regulon_activity)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SHR']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300(+)\nactivity", x = "", title = "Avp+ neurons\n(Spontaneous)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.shr.ep300.activity.png", width=365/96, height=189/96, dpi=300)

ggplot(avp_tf, aes(x = treatment, y = Regulon)) +
  geom_violin(trim = T, scale = "width", fill = species_col['SHR']) +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Ep300\nexpression", x = "", title = "Avp+ neurons\n(Spontaneous)") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2i.shr.ep300.expr.png", width=365/96, height=189/96, dpi=300)


# cor.test(avp_tf$Zfp595_activity, avp_tf$Avp)
# cor.test(avp_tf$Zfp595, avp_tf$Avp)
# ggplot(avp_tf, aes(x = Avp, y = Zfp595_activity)) +
#   geom_density_2d(bandwidth = 0.2, color = "black") +
#   geom_point(color = "steelblue") +
#   labs(x="Avp expression", y="Zfp595(+) activity") +
#   theme(axis.text = element_text(colour = "black"),
#         legend.position = "none") +
#   facet_wrap(vars(treatment))
# 
# ggplot(avp_tf, aes(x = Avp, y = Zfp595)) +
#   geom_point(color = "steelblue") +
#   labs(x="Avp expression", y="Zfp595 expression") +
#   theme(axis.text = element_text(colour = "black"),
#         legend.position = "none") +
#   facet_wrap(vars(treatment))
# 

ggplot(avp_tf, aes(x = treatment, y = Zfp595_activity)) +
  geom_violin(trim = T, scale = "width") +    # Create violin plot
  geom_jitter(size = 0.1, color="grey") +  # Add scatter points (jittered)
  theme_classic() +    # Classic theme for a clean look
  labs(y = "Zfp595(+)\nactivity", x = "", title = "Avp+ neurons") +  # Axis labels
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold")) 








