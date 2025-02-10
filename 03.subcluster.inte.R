source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DE_enrichment.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/00.initial_setting.R")

library(Seurat)
# library(SeuratDisk)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)
# library(gridExtra)
library(patchwork)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster")





################################################################################
###

plot = function(seurat_object, reso = "RNA_snn_res.0.3"){
  
  if(DefaultDimReduc(seurat_object)=="pca"){
    UMAP_plot = DimPlot(seurat_object, label = T, group.by = "subclass_level2", reduction = "wnn.umap.harmony") + theme(legend.position = "none")
  }else{
    UMAP_plot = DimPlot(seurat_object, label = T, group.by = "subclass_level2") + theme(legend.position = "none")
  }
  
  meta_table = seurat_object@meta.data
  barplot = stacked_barplot(meta_table, 
                            cluster = "subclass_level2", 
                            tissue = NULL,
                            treatment = "treatment")
  barplot_combined <- wrap_plots(barplot, nrow = 1, widths = c(2,2,1))
  
  markers %>%
    group_by(cluster) %>%
    dplyr::filter(p_val_adj<0.05, avg_log2FC>0.1) %>%
    slice_head(n = 10) %>%
    ungroup() -> top10
  heatmap = DoHeatmap(seurat_object, features = top10$gene, size=3) + NoLegend()
  
  heatmap+UMAP_plot+barplot_combined+
    plot_layout(
      # widths = c(1, 1, 5),
      design = "
    12
    13
    "
    )
  
}

stacked_barplot = function(meta_table, 
                           cluster = "subclass_level2", 
                           species = "strain",
                           project = NULL,
                           tissue = "tissue",
                           treatment = "sxt"){
  
  dat = meta_table[, c(cluster, species, project, tissue, treatment)]
  pct = dat %>% add_count(get(cluster)) %>% 
    group_by_all() %>% 
    dplyr::summarise(Count = n()) %>% ungroup() %>% 
    mutate(percentage = Count/n) %>% as.data.frame()
  colnames(pct)[1] = c("Cluster")
  
  # cell_col = cell_color(meta_table, cluster)
  
  # pdf("stackedbarplot.pdf", width=4, height = 4)
  fill_list = c(tissue, treatment)
  p_list = vector("list", length(fill_list)) 
  for(i in 1:length(fill_list)){
    p_list[[i]] = local({
      i <- i
      fill_variable <- fill_list[i]
      
      transitions <- dat %>%
        group_by(strain, !!sym(fill_variable)) %>%
        summarise(count = n(), .groups = 'drop') %>%
        left_join(dat %>% 
                    group_by(strain) %>%
                    summarise(total = n(), .groups = 'drop'), 
                  by = "strain") %>%
        mutate(percent = count / total) %>%
        arrange(strain, !!sym(fill_variable)) %>%
        group_by(strain) %>%
        mutate(cumsum = cumsum(percent))
      
      ggplot(data = pct, aes(y = percentage, x = Cluster, fill = get(fill_variable))) + 
        geom_bar(stat="identity", 
                 position = position_fill(reverse = TRUE)) +
        geom_hline(data = transitions, aes(yintercept = cumsum), color = "black", linetype = "dashed") +
        labs(title="", fill='',
             x='', y='Cell Proportion') +
        theme(legend.position = "bottom",
              plot.title = element_text(size=15, colour = 'black'),
              axis.text.y = element_text(size=10, colour = 'black'),
              axis.text.x = element_text(size=10, colour = 'black'),
              axis.title = element_text(size=15, colour = 'black'),
              strip.text = element_text(size=10)) +
        scale_x_discrete(limits = rev(levels(pct$Cluster))) + 
        coord_flip() +
        # scale_fill_manual(values = cell_col) + 
        facet_nested(. ~ strain, scales = "free_y")
    })
  }
  
  transitions <- dat %>%
    group_by(strain) %>%
    summarize(total = n()/nrow(dat), .groups = 'drop') %>%
    mutate(cumulative = cumsum(total)) %>%
    dplyr::select(strain, cumulative)
  
  p_list[[length(p_list)+1]] =       
    ggplot(data = pct, aes(y = percentage, x = Cluster, fill = strain)) + 
    geom_bar(stat="identity", 
             position = position_fill(reverse = TRUE)) +
    geom_hline(data = transitions, aes(yintercept = cumulative), color = "black", linetype = "dashed") +
    labs(title="", fill='',
         x='', y='Cell Proportion') +
    theme(legend.position = "bottom",
          plot.title = element_text(size=15, colour = 'black'),
          axis.text.y = element_text(size=10, colour = 'black'),
          axis.text.x = element_text(size=10, colour = 'black'),
          axis.title = element_text(size=15, colour = 'black'),
          strip.text = element_text(size=10)) +
    scale_x_discrete(limits = rev(levels(pct$Cluster))) + 
    coord_flip()
  
  return(p_list)
}


################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Astrocyte"); outfile = "mouse.HYP.astrocyte."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.3"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Astrocyte"); outfile = "rat.ss.HYP.astrocyte."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Astrocyte"); outfile = "rat.sp.HYP.astrocyte."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Myelinating OL"); outfile = "mouse.HYP.myelinating_ol."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.3"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Myelinating OL"); outfile = "rat.ss.HYP.myelinating_ol."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Microglia"); outfile = "mouse.HYP.microglia."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.5"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Microglia", "Activated microglia"); outfile = "rat.ss.HYP.microglia."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.5"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Microglia"); outfile = "rat.sp.HYP.microglia."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()




################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Microglia"); outfile = "mouse.MCA.microglia."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Microglia"); outfile = "rat.ss.MCA.microglia."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.6"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Macrophages"); outfile = "mouse.MCA.macro."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.8"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Macrophages"); outfile = "rat.ss.MCA.macro."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()




################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("CM"); outfile = "rat.sp.LV.CM."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.6"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()





################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("PT"); outfile = "mouse.LK.PT."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
DefaultAssay(seurat_object_use) = "RNA"
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.8"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("Macrophages"); outfile = "rat.ss.MCA.macro."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()




################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("CM"); outfile = "rat.sp.LV.CM."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.6"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()




################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("PT"); outfile = "mouse.LK.PT."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony_mo(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.1"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("PT"); outfile = "rat.ss.LK.PT."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("PT"); outfile = "rat.sp.LK.PT."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()







################################################################################
### subcluster for EC

################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("EC"); outfile = "mouse.LK.EC."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony_mo(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.1"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()


EC_ref = read.table("/xdisk/mliang1/qqiu/reference/single_cell_marker_ref/Cell-EC-2020.txt", header = T, sep = '\t')
EC_list = unique(EC_ref$cell.type)
kidney_EC_list = EC_list[grepl("kidney", EC_list)]

ref_list <- EC_ref %>%
  group_by(cell.type) %>%
  summarise(genes = list(official.gene.symbol), .groups = 'drop') %>%
  deframe()

ref_list_ues = ref_list[names(ref_list) %in% kidney_EC_list]
seurat_object_use = AddModuleScore(seurat_object_use, features=ref_list_ues)

################################################################################
mouse_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds"
ss_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds"
sp_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds"
mouse_so = readRDS(mouse_file)
ss_so = readRDS(ss_file)
sp_so = readRDS(sp_file)

cell_type = c("EC"); outfile = "LK.EC."
mouse_so = subset(mouse_so, subclass_level1 %in% cell_type)
ss_so = subset(ss_so, subclass_level1 %in% cell_type)
sp_so = subset(sp_so, subclass_level1 %in% cell_type)

mouse_so = FindVariableFeatures(mouse_so)
ss_so = FindVariableFeatures(ss_so)
sp_so = FindVariableFeatures(sp_so)

total_gene = unique(c(VariableFeatures(mouse_so),
                  VariableFeatures(ss_so),
                  VariableFeatures(sp_so)))

mouse_so = subset(mouse_so, features=total_gene)
ss_so = ss_so[overlap_gene, ]
sp_so = sp_so[overlap_gene, ]

seurat_object_use = merge(mouse_so, list(ss_so, sp_so))

seurat_object_use = QC_harmony_mo(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.1"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("EC"); outfile = "rat.sp.LK.EC."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony_mo(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "wsnn_res.")
reso = "wsnn_res.0.4"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("EC"); outfile = "mouse.LV.EC."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.1"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
# pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
# dev.off()



################################################################################
input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds"
seurat_object = readRDS(input_file)

cell_type = c("EC"); outfile = "rat.ss.LV.EC."
seurat_object_use = subset(seurat_object, subclass_level1 %in% cell_type)
seurat_object_use = QC_harmony(seurat_object_use, nPC = 20)
seurat_meta = seurat_object_use@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")
reso = "RNA_snn_res.0.5"

seurat_object_use$sxt = paste0(seurat_object_use$strain, "-", seurat_object_use$treatment)
seurat_object_use$subclass_level2 = paste0("C", seurat_object_use@meta.data[, reso])
Idents(seurat_object_use) = "subclass_level2"
saveRDS(seurat_object_use, paste0(outfile, "subcluster.rds"))

markers = DE_analysis(seurat_object_use, "subclass_level2", outfile)
pdf(paste0(outfile, "subcluster.pdf"), width = 10, height = 10)
plot(seurat_object_use, reso)
dev.off()

















































