
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
library(miloR)

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")






################################################################################
milo_mod = function(seurat_object, cell_type){
  
  meta_table = seurat_object@meta.data
  
  species_list = unique(meta_table$strain)
  
  for(si in species_list){
    
    seurat_object_use = subset(seurat_object, strain == si)
    
    if(si=="C57BL/6"){si = "mouse"}
    
    outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/",
                     si, ".", cell_type, ".miloR.rds")
    
    reduced.dim = "HARMONY"
    
    # construct KNN graph -> representative neighbourhoods
    sce = as.SingleCellExperiment(seurat_object_use, assay="RNA")
    milo_object = Milo(sce)
    milo_object = buildGraph(milo_object, k = 30, d = 30, reduced.dim = reduced.dim)
    milo_object = makeNhoods(milo_object, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = reduced.dim,
                             refinement_scheme="graph")
    plotNhoodSizeHist(milo_object) # average neighbourhood size should be over 5 x N_samples/ 50-100
    
    milo_object <- countCells(milo_object, meta.data = seurat_object_use@meta.data, sample="orig.ident")
    
    milo_design <- seurat_object_use@meta.data[,c("orig.ident", "treatment")]
    milo_design <- distinct(milo_design)
    rownames(milo_design) <- milo_design$orig.ident
    
    milo_design$orig.ident = as.factor(milo_design$orig.ident)
    milo_design$treatment = as.factor(milo_design$treatment)
    
    # computing neighbourhood connectivity
    milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = reduced.dim)
    saveRDS(milo_object, outfile)
    
  }
}


milo_da_merged = function(input_file, cell_type, outfile, coldata_col="subclass_level2"){
  
  da_merge = c()
  for( i in input_file ){
    
    reduced.dim = "HARMONY"
    
    milo_object = readRDS(i)
    
    milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
    milo_design = data.frame(colData(milo_object))[,c("treatment", "strain", "orig.ident")]
    milo_design = distinct(milo_design)
    rownames(milo_design) = milo_design$orig.ident
    milo_design$treatment.new = as.factor(gsub(" ", "", milo_design$treatment))
    milo_design$strain = as.factor(milo_design$strain)
    
    species_list = unique(milo_design$strain)
    treatment_list = intersect(levels(colData(milo_object)$treatment), unique(milo_design$treatment))
    treatment.new_list = gsub(" ", "", treatment_list)
    
    for( si in species_list ){
      
      for( j in 1:length(treatment.new_list[-1]) ){
        
        control.new = treatment.new_list[1]
        treatment.new = treatment.new_list[j+1]
        control = treatment_list[1]
        treatment = treatment_list[j+1]
        
        model.contrasts = paste0("treatment.new", control.new, " - ", "treatment.new", treatment.new)
        
        da_tmp = try(testNhoods(milo_object, design = ~ 0 + treatment.new, design.df = milo_design,
                                fdr.weighting="graph-overlap", reduced.dim = reduced.dim,
                                model.contrasts = model.contrasts), silent = T)
        
        if(class(da_tmp) != "try-error"){
          
          da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = coldata_col)
          da_tmp$subcluster = da_tmp[, coldata_col]
          da_tmp$cell_type = cell_type
          da_tmp$species = si
          da_tmp$control = control
          da_tmp$treatment = treatment
          da_merge = rbind(da_merge, da_tmp)
          
        }else{
          print(c(i, treatment))
        }
        
      }
      
    }
    
  }
  
  write.table(da_merge, outfile, sep='\t', quote=F, col.names = T, row.names = F)
  
}


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
### subcluster analysis
hyp_marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", "Avp",
         "Slc1a2", "Cx3cr1", "P2ry12", "Tgfbr1", "Cspg4", "Pdgfra","Mbp", "St18", 
         "Col23a1", "Tmem212", "Flt1", "Pecam1", "Ebf1",
         "Atp13a5", "Ptgds")

### astrocyte
hyp_astro_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.astrocyte.subcluster.rds')
hyp_astro_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.astrocyte.subcluster.rds')
hyp_astro_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.astrocyte.subcluster.rds')

hyp_astro_m$sxt = factor(paste0(hyp_astro_m$species, " - ", hyp_astro_m$treatment), levels = sxt_order)
hyp_astro_ss$sxt = factor(paste0(hyp_astro_ss$species, " - ", hyp_astro_ss$treatment), levels = sxt_order)
hyp_astro_shr$sxt = factor(paste0(hyp_astro_shr$species, " - ", hyp_astro_shr$treatment), levels = sxt_order)

hyp_astro_m %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_astro_m %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hyp_astro_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_astro_ss %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hyp_astro_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_astro_shr %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


hyp_astro_m = subset(hyp_astro_m, RNA_snn_res.1 %in% c(0:5, 8:9))
hyp_astro_ss = subset(hyp_astro_ss, RNA_snn_res.1 %in% c(0:8, 11))
hyp_astro_shr = subset(hyp_astro_shr, RNA_snn_res.1 %in% c(0:6))

hyp_astro_m = QC_harmony(hyp_astro_m)
hyp_astro_ss = QC_harmony(hyp_astro_ss)
hyp_astro_shr = QC_harmony(hyp_astro_shr)

hyp_astro_m %>% DimPlot(., group.by = "RNA_snn_res.0.8", label = T)
hyp_astro_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_astro_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)

hyp_astro_m_marker = FindAllMarkers(hyp_astro_m, group.by = "RNA_snn_res.0.8")
hyp_astro_ss_marker = FindAllMarkers(hyp_astro_ss, group.by = "RNA_snn_res.1")
hyp_astro_shr_marker = FindAllMarkers(hyp_astro_shr, group.by = "RNA_snn_res.1")

e <- new.env()
assign("astrocyte_m_marker", hyp_astro_m_marker, envir = e)
assign("astrocyte_ss_marker", hyp_astro_ss_marker, envir = e)
assign("astrocyte_shr_marker", hyp_astro_shr_marker, envir = e)
saveRDS(e, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.astrocyte.marker.rds")


# e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.astrocyte.marker.rds")
# hyp_astro_m_marker = e$astrocyte_m_marker
# hyp_astro_ss_marker = e$astrocyte_ss_marker
# hyp_astro_shr_marker = e$astrocyte_shr_marker

hyp_astro_m_marker = hyp_astro_m_marker[hyp_astro_m_marker$p_val_adj<0.05 & hyp_astro_m_marker$avg_log2FC>0.25, ]
hyp_astro_ss_marker = hyp_astro_ss_marker[hyp_astro_ss_marker$p_val_adj<0.05 & hyp_astro_ss_marker$avg_log2FC>0.25, ]
hyp_astro_shr_marker = hyp_astro_shr_marker[hyp_astro_shr_marker$p_val_adj<0.05 & hyp_astro_shr_marker$avg_log2FC>0.25, ]


hyp_astro_m_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_astro_m, return.seurat = TRUE, group.by = "RNA_snn_res.0.8"), ., draw.lines = FALSE) & coord_flip()

hyp_astro_ss_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_astro_ss, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) & coord_flip()

hyp_astro_shr_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_astro_shr, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) & coord_flip()










### microglia
hyp_micro_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.microglia.subcluster.rds')
hyp_micro_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.microglia.subcluster.rds')
hyp_micro_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.microglia.subcluster.rds')

hyp_micro_m$sxt = factor(paste0(hyp_micro_m$species, " - ", hyp_micro_m$treatment), levels = sxt_order)
hyp_micro_ss$sxt = factor(paste0(hyp_micro_ss$species, " - ", hyp_micro_ss$treatment), levels = sxt_order)
hyp_micro_shr$sxt = factor(paste0(hyp_micro_shr$species, " - ", hyp_micro_shr$treatment), levels = sxt_order)

hyp_micro_m %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_micro_m %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hyp_micro_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_micro_ss %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hyp_micro_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_micro_shr %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


hyp_micro_m = subset(hyp_micro_m, RNA_snn_res.1 %in% c(0:2, 5))
hyp_micro_ss = subset(hyp_micro_ss, RNA_snn_res.1 %in% c(0:6, 8:10))
hyp_micro_shr = subset(hyp_micro_shr, RNA_snn_res.1 %in% c(0:5, 7))

hyp_micro_m = QC_harmony(hyp_micro_m)
hyp_micro_ss = QC_harmony(hyp_micro_ss)
hyp_micro_shr = QC_harmony(hyp_micro_shr)

hyp_micro_m %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_micro_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_micro_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)

hyp_micro_m_marker = FindAllMarkers(hyp_micro_m, group.by = "RNA_snn_res.1")
hyp_micro_ss_marker = FindAllMarkers(hyp_micro_ss, group.by = "RNA_snn_res.1")
hyp_micro_shr_marker = FindAllMarkers(hyp_micro_shr, group.by = "RNA_snn_res.1")

e <- new.env()
assign("microglia_m_marker", hyp_micro_m_marker, envir = e)
assign("microglia_ss_marker", hyp_micro_ss_marker, envir = e)
assign("microglia_shr_marker", hyp_micro_shr_marker, envir = e)
saveRDS(e, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.microglia.marker.rds")

# e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.microglia.marker.rds")
# hyp_micro_m_marker = e$microglia_m_marker
# hyp_micro_ss_marker = e$microglia_ss_marker
# hyp_micro_shr_marker = e$microglia_shr_marker

hyp_micro_m_marker = hyp_micro_m_marker[hyp_micro_m_marker$p_val_adj<0.05 & hyp_micro_m_marker$avg_log2FC>0.25, ]
hyp_micro_ss_marker = hyp_micro_ss_marker[hyp_micro_ss_marker$p_val_adj<0.05 & hyp_micro_ss_marker$avg_log2FC>0.25, ]
hyp_micro_shr_marker = hyp_micro_shr_marker[hyp_micro_shr_marker$p_val_adj<0.05 & hyp_micro_shr_marker$avg_log2FC>0.25, ]


hyp_micro_m_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_micro_m, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) & coord_flip()

hyp_micro_ss_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_micro_ss, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) & coord_flip()

hyp_micro_shr_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_micro_shr, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) & coord_flip()









### myelinating_ol
hyp_ol_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.myelinating_ol.subcluster.rds')
hyp_ol_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.myelinating_ol.subcluster.rds')
# hyp_ol_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.myelinating_ol.subcluster.rds')

hyp_ol_m$sxt = factor(paste0(hyp_ol_m$species, " - ", hyp_ol_m$treatment), levels = sxt_order)
hyp_ol_ss$sxt = factor(paste0(hyp_ol_ss$species, " - ", hyp_ol_ss$treatment), levels = sxt_order)
# hyp_ol_shr$sxt = factor(paste0(hyp_ol_shr$species, " - ", hyp_ol_shr$treatment), levels = sxt_order)

hyp_ol_m %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_ol_m %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
hyp_ol_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
hyp_ol_ss %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# hyp_ol_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
# hyp_ol_shr %>% DotPlot(., features = hyp_marker_list, group.by = "RNA_snn_res.1") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


hyp_ol_m = subset(hyp_ol_m, RNA_snn_res.1 %in% c(0:7, 9:10, 12))
# hyp_ol_ss = subset(hyp_ol_ss, RNA_snn_res.1 %in% c(0:6, 8:11))
# hyp_ol_shr = subset(hyp_ol_shr, RNA_snn_res.1 %in% c(0:7))

hyp_ol_m = QC_harmony(hyp_ol_m)
# hyp_ol_ss = QC_harmony(hyp_ol_ss)
# hyp_ol_shr = QC_harmony(hyp_ol_shr)

hyp_ol_m %>% DimPlot(., group.by = "RNA_snn_res.0.8", label = T)
# hyp_ol_ss %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)
# hyp_ol_shr %>% DimPlot(., group.by = "RNA_snn_res.1", label = T)

hyp_ol_m_marker = FindAllMarkers(hyp_ol_m, group.by = "RNA_snn_res.0.8")
hyp_ol_ss_marker = FindAllMarkers(hyp_ol_ss, group.by = "RNA_snn_res.1")
# hyp_ol_shr_marker = FindAllMarkers(hyp_ol_shr, group.by = "RNA_snn_res.1")

e <- new.env()
assign("olglia_m_marker", hyp_ol_m_marker, envir = e)
assign("olglia_ss_marker", hyp_ol_ss_marker, envir = e)
# assign("olglia_shr_marker", hyp_ol_shr_marker, envir = e)
saveRDS(e, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.myelinating_ol.marker.rds")

# e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/HYP.myelinating_ol.marker.rds")
# hyp_ol_m_marker = e$olglia_m_marker
# hyp_ol_ss_marker = e$olglia_ss_marker
# hyp_ol_shr_marker = e$olglia_shr_marker

hyp_ol_m_marker = hyp_ol_m_marker[hyp_ol_m_marker$p_val_adj<0.05 & hyp_ol_m_marker$avg_log2FC>0.25, ]
hyp_ol_ss_marker = hyp_ol_ss_marker[hyp_ol_ss_marker$p_val_adj<0.05 & hyp_ol_ss_marker$avg_log2FC>0.25, ]
# hyp_ol_shr_marker = hyp_ol_shr_marker[hyp_ol_shr_marker$p_val_adj<0.05 & hyp_ol_shr_marker$avg_log2FC>0.25, ]


hyp_ol_m_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_ol_m, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE) + facet_wrap()

hyp_ol_ss_marker %>% 
  mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
  arrange(desc(pct.diff)) %>%
  slice(seq_len(10)) %>% ungroup() %>%
  pull(gene) %>%
  DoHeatmap(AverageExpression(hyp_ol_ss, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE)

# hyp_ol_shr_marker %>% 
#   mutate(pct.diff = pct.1 - pct.2) %>% group_by(cluster) %>%
#   arrange(desc(pct.diff)) %>%
#   slice(seq_len(10)) %>% ungroup() %>%
#   pull(gene) %>%
#   DoHeatmap(AverageExpression(hyp_ol_shr, return.seurat = TRUE, group.by = "RNA_snn_res.1"), ., draw.lines = FALSE)










