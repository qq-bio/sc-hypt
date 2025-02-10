dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(sctransform)
library(harmony)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(R.utils)
library(ggsci)
library(presto)
library(clustree)
library(ROGUE)
library(UpSetR)
library(miloR)
library(ggh4x)
library(ggbeeswarm)
library(ggtext)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/00.initial_setting.R")



################################################################################
### strain wise
# 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
# )
# #
# #
# cluster = "subclass_level1"
# for(i in input_file){
# 
#   tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
# 
#   seurat_object = readRDS(i)
#   DefaultAssay(seurat_object) = "RNA"
# 
#   seurat_object@active.ident = seurat_object@meta.data[, cluster]
#   Idents(seurat_object) = cluster
# 
#   seurat_object$strain = gsub(" ", "", seurat_object$strain)
#   meta_table = seurat_object@meta.data
# 
#   si = unique(meta_table$project)
#   treatment_list = unique(meta_table$treatment)
# 
#   for(ti in treatment_list){
# 
#     seurat_object_use = subset(seurat_object, treatment == ti)
#     strain_list = unique(seurat_object_use$strain)
#     
#     if(length(strain_list)>1){
#       
#       outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/",
#                        si, ".", ti, ".", tissue, ".miloR.rds")
#       
#       if(tissue=="LK"){
#         reduced.dim = "HARMONY.RNA"
#       }else{
#         reduced.dim = "HARMONY"
#       }
#       
#       # construct KNN graph -> representative neighbourhoods
#       sce = as.SingleCellExperiment(seurat_object_use, assay="RNA")
#       milo_object = Milo(sce)
#       milo_object = buildGraph(milo_object, k = 30, d = 30, reduced.dim = reduced.dim)
#       milo_object = makeNhoods(milo_object, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = reduced.dim,
#                                refinement_scheme="graph")
#       plotNhoodSizeHist(milo_object) # average neighbourhood size should be over 5 x N_samples/ 50-100
#       
#       milo_object <- countCells(milo_object, meta.data = seurat_object_use@meta.data, sample="orig.ident")
#       
#       milo_design <- seurat_object_use@meta.data[,c("orig.ident", "treatment")]
#       milo_design <- distinct(milo_design)
#       rownames(milo_design) <- milo_design$orig.ident
#       
#       milo_design$orig.ident = as.factor(milo_design$orig.ident)
#       milo_design$treatment = as.factor(milo_design$treatment)
#       
#       # computing neighbourhood connectivity
#       milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = reduced.dim)
#       saveRDS(milo_object, outfile)
#       
#       
#     }
# 
#   }
# 
# }









################################################################################
input_file = c(

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.LS.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.LS.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.LS.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.LS.MSA.miloR.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.10w.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.10w.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.10w.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.HS\ 3d.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.HS\ 3d.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.HS\ 3d.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Salt-sensitive.HS\ 3d.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.26w.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.26w.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/Spontaneous.26w.MSA.miloR.rds"
  
)


################################################################################
library(ggtext)
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR")


da_merge = c()
for( i in input_file ){

  tissue = sub(".*\\.([^.]*)\\.(miloR).*", "\\1", basename(i))
  reduction = ifelse(tissue=="LK", "HARMONY.RNA", "HARMONY")

  milo_object = readRDS(i)

  milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
  milo_design = data.frame(colData(milo_object))[,c("treatment", "strain", "orig.ident")]
  milo_design = distinct(milo_design)
  rownames(milo_design) = milo_design$orig.ident
  milo_design$treatment.new = as.factor(gsub(" ", "", milo_design$treatment))
  milo_design$strain = as.factor(milo_design$strain)

  strain_list = unique(milo_design$strain)
  pi = unique(colData(milo_object)$project)
  ti = unique(colData(milo_object)$treatment)
  
  responsive_strain = intersect(strain_list, c("SS", "SHR"))
  control_strain = intersect(strain_list, c("SD", "WKY"))
  
  model.contrasts = paste0("strain", control_strain, " - ", "strain", responsive_strain)

  da_tmp = try(testNhoods(milo_object, design = ~ 0 + strain, design.df = milo_design,
                          fdr.weighting="graph-overlap", reduced.dim = reduction,
                          model.contrasts = model.contrasts), silent = T)

  if(class(da_tmp) != "try-error"){

    da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = "subclass_level1")
    da_tmp$SpatialFDR.org = da_tmp$SpatialFDR
    da_tmp$SpatialFDR = p.adjust(da_tmp$PValue, "fdr")
    da_tmp$subcluster = factor(da_tmp$subclass_level1, levels = cell_order)
    da_tmp$species = pi
    da_tmp$tissue = tissue
    da_tmp$control = paste0(control_strain, " - ", ti)
    da_tmp$treatment = paste0(responsive_strain, " - ", ti)
    da_merge = rbind(da_merge, da_tmp)

  }else{
    print(c(i, treatment))
  }

}

write.table(da_merge, "milo.strain_wise.da_result.out", sep='\t', quote=F, col.names = T, row.names = F)
  



















################################################################################
 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
# )
# # 
# # 
# cluster = "subclass_level1"
# for(i in input_file){
# 
#   # deg_merged = c()
# 
#   tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
# 
#   seurat_object = readRDS(i)
#   DefaultAssay(seurat_object) = "RNA"
# 
#   if(grepl("HYP", i)){
#     seurat_object$subclass_level2 = seurat_object$subclass_level1
#   }
# 
#   seurat_object@active.ident = seurat_object@meta.data[, cluster]
#   Idents(seurat_object) = cluster
#   
#   seurat_object$strain = gsub(" ", "", seurat_object$strain)
#   meta_table = seurat_object@meta.data
# 
#   # project_list = unique(meta_table$project)
#   # for(pi in project_list){
#   #
#     species_list = unique(meta_table$strain)
#     # species_list = intersect(strains, species_list)
# 
#     for(si in species_list){
# 
#       seurat_object_use = subset(seurat_object, strain == si)
# 
#       if(si=="C57BL/6"){si = "mouse"}
# 
#       outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/",
#                        si, ".", tissue, ".miloR.rds")
# 
#       if(tissue=="LK"){
#         reduced.dim = "HARMONY.RNA"
#       }else{
#         reduced.dim = "HARMONY"
#       }
# 
#       # construct KNN graph -> representative neighbourhoods
#       sce = as.SingleCellExperiment(seurat_object_use, assay="RNA")
#       milo_object = Milo(sce)
#       milo_object = buildGraph(milo_object, k = 30, d = 30, reduced.dim = reduced.dim)
#       milo_object = makeNhoods(milo_object, prop = 0.2, k = 30, d=30, refined = TRUE, reduced_dims = reduced.dim,
#                                refinement_scheme="graph")
#       plotNhoodSizeHist(milo_object) # average neighbourhood size should be over 5 x N_samples/ 50-100
# 
#       milo_object <- countCells(milo_object, meta.data = seurat_object_use@meta.data, sample="orig.ident")
#       # head(nhoodCounts(milo_object))
# 
#       milo_design <- seurat_object_use@meta.data[,c("orig.ident", "treatment")]
#       # Sequencing_batch = c(1, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 2)
#       # names(Sequencing_batch) = names(table(milo_design$orig.ident))
#       # milo_design$Sequencing_batch = as.factor(Sequencing_batch[milo_design$orig.ident])
#       milo_design <- distinct(milo_design)
#       rownames(milo_design) <- milo_design$orig.ident
# 
#       milo_design$orig.ident = as.factor(milo_design$orig.ident)
#       milo_design$treatment = as.factor(milo_design$treatment)
# 
#       # computing neighbourhood connectivity
#       milo_object <- calcNhoodDistance(milo_object, d=30, reduced.dim = reduced.dim)
#       saveRDS(milo_object, outfile)
# 
#     }
# 
#   # }
# 
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ################################################################################ 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.MSA.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.MSA.miloR.rds",
#   # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MSA.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.HYP.miloR.rds",
#   # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.MSA.miloR.rds"
#   # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.MCA.miloR.rds"
# )
# 
# 
# 
# 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LV.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.ss.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.ss.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.ss.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.ss.MSA.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.ss.MCA.miloR.rds",
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.sp.HYP.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.sp.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.sp.LK.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.sp.MSA.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/rat.sp.MCA.miloR.rds"
# )



# ################################################################################
# library(ggtext)
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR")
# 
# 
# da_merge = c()
# for( i in input_file ){
# 
#   # tissue = strsplit(basename(i), "\\.")[[1]][2]
#   tissue = sub(".*\\.([^.]*)\\.(miloR).*", "\\1", basename(i))
#   reduction = ifelse(tissue=="LK", "HARMONY.RNA", "HARMONY")
# 
#   milo_object = readRDS(i)
# 
#   milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
#   milo_design = data.frame(colData(milo_object))[,c("treatment", "strain", "orig.ident")]
#   milo_design = distinct(milo_design)
#   rownames(milo_design) = milo_design$orig.ident
#   milo_design$treatment.new = as.factor(gsub(" ", "", milo_design$treatment))
#   milo_design$strain = as.factor(milo_design$strain)
# 
#   species_list = unique(milo_design$strain)
#   treatment_list = intersect(levels(colData(milo_object)$treatment), unique(milo_design$treatment))
#   treatment.new_list = gsub(" ", "", treatment_list)
# 
#   for( si in species_list ){
# 
#     for( j in 1:length(treatment.new_list[-1]) ){
# 
#       control.new = treatment.new_list[1]
#       treatment.new = treatment.new_list[j+1]
#       control = treatment_list[1]
#       treatment = treatment_list[j+1]
# 
#       model.contrasts = paste0("treatment.new", control.new, " - ", "treatment.new", treatment.new)
# 
#       da_tmp = try(testNhoods(milo_object, design = ~ 0 + treatment.new, design.df = milo_design,
#                           fdr.weighting="graph-overlap", reduced.dim = reduction,
#                           model.contrasts = model.contrasts), silent = T)
# 
#       # if(length(species_list)==1){
#         # model.contrasts = paste0("treatment.new", control.new, " - ", "treatment.new", treatment.new)
#         # 
#         # da_tmp = try(testNhoods(milo_object, design = ~ 0 + treatment.new, design.df = milo_design,
#         #                         fdr.weighting="graph-overlap", reduced.dim = reduction,
#         #                         model.contrasts = model.contrasts), silent = T)
# 
#       # }else{
#       #   model.contrasts = paste0("species", si, ":treatment.new", control.new, " - ", "species", si, ":treatment.new", treatment.new)
#       # 
#       #   da_tmp = try(testNhoods(milo_object, design = ~ 0 + species*treatment.new, design.df = milo_design,
#       #                           fdr.weighting="graph-overlap", reduced.dim = reduction,
#       #                           model.contrasts = model.contrasts), silent = T)
#       # 
#       # }
# 
#       if(class(da_tmp) != "try-error"){
# 
#         # p = plotDAbeeswarm(da_tmp, group.by = "subcluster", alpha=0.1) + scale_x_discrete(limits=rev) +
#         #   labs(
#         #     x = "", y = "Log Fold Change",
#         #     title = paste0(control, " vs. ", treatment, " (", species, ")<br>
#         #             <span style='font-size:11pt'>
#         #             <span style='color:#832424;'>Enriched in ", treatment, "<br></span>
#         #             <span style='color:#3A3A98;'>Depleted in ", treatment, "</span>
#         #             </span>")
#         #   ) +
#         #   theme(plot.title = element_markdown(lineheight = 1.1),
#         #         legend.text = element_markdown(size = 11),
#         #         text=element_text(size=11, color="black"),
#         #         axis.text.y = element_text(size=11, colour = 'black'),
#         #         axis.text.x = element_text(size=11, colour = 'black'),
#         #         axis.title = element_text(size=11, colour = 'black'))
#         # print(p)
# 
#         da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = "subclass_level1")
#         da_tmp$subcluster = factor(da_tmp$subclass_level1, levels = cell_order)
#         da_tmp$species = si
#         da_tmp$tissue = tissue
#         da_tmp$control = control
#         da_tmp$treatment = treatment
#         da_merge = rbind(da_merge, da_tmp)
# 
#       }else{
#         print(c(i, treatment))
#       }
# 
#     }
# 
#   }
# 
# }
# 
# write.table(da_merge, "milo.da_result.out", sep='\t', quote=F, col.names = T, row.names = F)
# 
# 
# 
# 
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/")
da_merge = read.table("milo.da_result.out", sep='\t', header=T)
da_merge = da_merge[!(da_merge$tissue=="MCA" & da_merge$species %in% c("C57BL/6", "SS")), ]
da_merge$subcluster = factor(da_merge$subcluster, levels = cell_order)
da_merge$species = factor(da_merge$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
da_merge$tissue = factor(da_merge$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
da_merge$treatment = factor(da_merge$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))
da_merge$SpatialFDR.org = da_merge$SpatialFDR
da_merge$SpatialFDR = p.adjust(da_merge$PValue, "fdr")
da_merge_plot <- da_merge %>%
  mutate(order = case_when(
    SpatialFDR > 0.05 ~ 1,
    SpatialFDR <= 0.05 ~ 2
  ))
da_merge_plot <- da_merge_plot %>% arrange(order)
# plotDAbeeswarm(da_merge_plot, group.by = "subcluster", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "Strain - treatment (vs. control)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in treatment<br></span>
#     <span style='color:#3A3A98;'>Depleted in treatment</span>
#     </span>"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11, color="black"),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=9, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank(),
#         strip.background = element_rect(colour = "black", fill = NA)) +
#   facet_nested(tissue ~ species + treatment , scales = "free", space = "free_x") +
#   geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
#   geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
#   stat_summary(fun = median, geom = "point", color = "black", size = 2) +
#   scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
#   scale_color_gradient2( na.value = "lightgrey" )


group.by = "subcluster"; alpha=0.1
da.res = mutate(da_merge_plot, group_by = da_merge_plot[, group.by])
da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
  arrange(group_by) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
  ggplot(aes(group_by, -logFC, color = logFC_color)) + guides(color = "none") +
  xlab(group.by) + ylab("Log Fold Change") + geom_quasirandom(alpha = 1) +
  coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
  labs(
    x = "", y = "Log Fold Change",
    title = "Strain - treatment (vs. control)<br>
    <span style='font-size:11pt'>
    <span style='color:#832424;'>Enriched in treatment<br></span>
    <span style='color:#3A3A98;'>Depleted in treatment</span>
    </span>"
  ) +
  theme_bw() +
  theme(plot.title = element_markdown(lineheight = 1.1),
        legend.text = element_markdown(size = 11),
        text=element_text(size=11, color="black"),
        axis.text.y = element_text(size=11, colour = 'black'),
        axis.text.x = element_text(size=9, colour = 'black'),
        axis.title = element_text(size=11, colour = 'black'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(colour = "black", fill = NA)) +
  facet_nested(tissue ~ species + treatment , scales = "free", space = "free") +
  geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
  scale_color_gradient2( na.value = "lightgrey" )
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/milo.png", width=731/96, height=957/96, dpi=300)







setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/")
da_merge = read.table("milo.strain_wise.da_result.out", sep='\t', header=T)
da_merge$subcluster = factor(da_merge$subcluster, levels = cell_order)
da_merge$species = factor(da_merge$species, levels = c("Salt-sensitive", "Spontaneous"))
da_merge$tissue = factor(da_merge$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
da_merge$treatment = factor(da_merge$treatment, levels = c("SS - LS", "SS - HS 3d", "SHR - 10w", "SHR - 26w"))
# da_merge$SpatialFDR.org = da_merge$SpatialFDR
# da_merge$SpatialFDR = p.adjust(da_merge$PValue, "fdr")
da_merge_plot <- da_merge %>%
  mutate(order = case_when(
    SpatialFDR > 0.05 ~ 1,
    SpatialFDR <= 0.05 ~ 2
  ))
da_merge_plot <- da_merge_plot %>% arrange(order)

group.by = "subcluster"; alpha=0.1
da.res = mutate(da_merge_plot, group_by = da_merge_plot[, group.by])
da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
  arrange(group_by) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
  ggplot(aes(group.by, -logFC, color = logFC_color)) + guides(color = "none") +
  xlab(group.by) + ylab("Log Fold Change") + geom_quasirandom(alpha = 1) +
  coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
  labs(
    x = "", y = "Log Fold Change",
    title = "Strain - treatment (vs. control)<br>
    <span style='font-size:11pt'>
    <span style='color:#832424;'>Enriched in hypertensive strain<br></span>
    <span style='color:#3A3A98;'>Depleted in hypertensive strain</span>
    </span>"
  ) +
  theme_bw() +
  theme(plot.title = element_markdown(lineheight = 1.1),
        legend.text = element_markdown(size = 11),
        text=element_text(size=11, color="black"),
        axis.text.y = element_text(size=11, colour = 'black'),
        axis.text.x = element_text(size=9, colour = 'black'),
        axis.title = element_text(size=11, colour = 'black'),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        strip.background = element_rect(colour = "black", fill = NA)) +
  facet_nested(tissue ~ species + treatment , scales = "free", space = "free") +
  geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
  scale_color_gradient2( na.value = "lightgrey" )
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/strain_wise.milo.png", width=456/96, height=793/96, dpi=300)








# 
# 
# 
# library(ggbeeswarm)
# da_merge_plot$logFC = da_merge_plot$logFC * -1
# da_merge_plot$logFC_color = ifelse(da_merge_plot$SpatialFDR<0.05, ifelse(da_merge_plot$logFC>0, "up", "down"), NA)
# ggplot(da_merge_plot, aes(subcluster, logFC, color = logFC_color, alpha=0.6)) + geom_quasirandom() +
#   guides(color = "none") +
#   scale_color_manual(values = c(
#     "up" = "#832424",
#     "down" = "#3A3A98"
#   ),
#   na.value = "lightgrey") +
#   ylab("Log Fold Change") +
#   coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "Strain - treatment (vs. control)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in treatment<br></span>
#     <span style='color:#3A3A98;'>Depleted in treatment</span>
#     </span>"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11, color="black"),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=9, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank()) +
#   facet_nested(tissue ~ species + treatment , scales = "free", space = "free_x") +
#   geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
#   geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
#   stat_summary(fun = median, geom = "point", color = "black", size = 2) +
#   scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# da_merge_plot$logFC = da_merge_plot$logFC * -1
# da_merge_plot$logFC_color = ifelse(da_merge_plot$SpatialFDR<0.05, ifelse(da_merge_plot$logFC>0, "up", "down"), NA)
# ggplot(da_merge_plot, aes(subcluster, logFC, color = logFC_color, alpha=0.6)) + geom_quasirandom() +
#   guides(color = "none") +
#   scale_color_manual(values = c(
#     "up" = "#832424",
#     "down" = "#3A3A98"
#   ),
#   na.value = "lightgrey") +
#   ylab("Log Fold Change") +
#   coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "Strain - treatment (vs. control)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in treatment<br></span>
#     <span style='color:#3A3A98;'>Depleted in treatment</span>
#     </span>"
#   ) +
#   theme_bw() +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11, color="black"),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=9, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'),
#         panel.grid.minor.x = element_blank(),
#         panel.grid.major.x = element_blank()) +
#   facet_nested(tissue ~ species + treatment , scales = "free", space = "free") +
#   geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
#   geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
#   stat_summary(fun = median, geom = "point", color = "black", size = 2) +
#   scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# milo_object = readRDS(input_file[1])
# milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
# milo_design = data.frame(colData(milo_object))[,c("treatment", "species", "orig.ident")]
# milo_design = distinct(milo_design)
# rownames(milo_design) = milo_design$orig.ident
# milo_design$treatment = as.factor(gsub(" ", "", milo_design$treatment))
# milo_design$species = as.factor(milo_design$species)
# 
# da_3d = testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                    fdr.weighting="graph-overlap", reduced.dim = "HARMONY",
#                    model.contrasts = c("treatmentSaline3d - treatmentAngII3d"))
# da_28d = testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                     fdr.weighting="graph-overlap", reduced.dim = "HARMONY",
#                     model.contrasts = c("treatmentSaline3d - treatmentAngII28d"))
# 
# da_3d <- annotateNhoods(milo_object, da_3d, coldata_col = "subclass_level1")
# da_3d$subcluster = factor(da_3d$subclass_level1, levels = cell_order)
# plotDAbeeswarm(da_3d, group.by = "subcluster", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "AngII 3d vs. Saline 3d (C57BL/6)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in AngII 3d<br></span>
#     <span style='color:#3A3A98;'>Depleted in AngII 3d</span>
#     </span>"
#   ) +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11, color="black"),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=11, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'))
# 
# 
# da_28d <- annotateNhoods(milo_object, da_28d, coldata_col = "subclass_level1")
# da_28d$subcluster = factor(da_28d$subclass_level1, levels=cell_order)
# plotDAbeeswarm(da_28d, group.by = "subcluster", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "AngII 28d vs. Saline 3d (C57BL/6)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in AngII 28d<br></span>
#     <span style='color:#3A3A98;'>Depleted in AngII 28d</span>
#     </span>"
#   ) +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=11, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'))
# 
# 
# 
# # DA in rat-ss
# ss_subcluster_order = c(paste0("GABA-", 1:18), paste0("Glu-", 1:17), paste0("Astro-", 1:6),
#                         levels(seurat_ss$new.cluster.ids_umap)[5:19])
# 
# milo_object = readRDS("rat.ss.HYP.miloR.p_0.2.k_30.HARMONY.rds")
# milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
# milo_design = data.frame(colData(milo_object))[,c("treatment", "species", "orig.ident")]
# milo_design = distinct(milo_design)
# rownames(milo_design) = milo_design$orig.ident
# milo_design$treatment = as.factor(gsub(" ", "", milo_design$treatment))
# milo_design$species = as.factor(milo_design$species)
# 
# da_3d = testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                    fdr.weighting="graph-overlap", reduced.dim = reduction,
#                    model.contrasts = c("treatmentLS - treatmentHS3d"))
# da_21d = testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                     fdr.weighting="graph-overlap", reduced.dim = reduction,
#                     model.contrasts = c("treatmentLS - treatmentHS21d"))
# 
# da_3d <- annotateNhoods(milo_object, da_3d, coldata_col = "subcluster")
# da_3d$subcluster = factor(da_3d$subcluster, levels=ss_subcluster_order)
# plotDAbeeswarm(da_3d, group.by = "subcluster", alpha=0.5) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "HS 3d vs. LS (SS)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in HS 3d<br></span>
#     <span style='color:#3A3A98;'>Depleted in HS 3d</span>
#     </span>"
#   ) +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=11, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'))
# 
# 
# da_21d <- annotateNhoods(milo_object, da_21d, coldata_col = "subcluster")
# da_21d$subcluster = factor(da_21d$subcluster, levels=ss_subcluster_order)
# plotDAbeeswarm(da_21d, group.by = "subcluster", alpha=0.5) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "HS 21d vs. LS (SS)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in HS 21d<br></span>
#     <span style='color:#3A3A98;'>Depleted in HS 21d</span>
#     </span>"
#   ) +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11),
#         axis.text.y = element_text(size=11, colour = 'black'),
#         axis.text.x = element_text(size=11, colour = 'black'),
#         axis.title = element_text(size=11, colour = 'black'))
# 
# 
# 
# 
# milo_object = readRDS("rat.sd.HYP.miloR.p_0.2.k_30.HARMONY.rds")
# milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
# milo_design = data.frame(colData(milo_object))[,c("treatment", "species", "orig.ident")]
# milo_design = distinct(milo_design)
# rownames(milo_design) = milo_design$orig.ident
# milo_design$treatment = as.factor(gsub(" ", "", milo_design$treatment))
# milo_design$species = as.factor(milo_design$species)
# 
# da_3d = testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                    fdr.weighting="graph-overlap", reduced.dim = reduction,
#                    model.contrasts = c("treatmentLS - treatmentHS3d"))
# 
# da_3d <- annotateNhoods(milo_object, da_3d, coldata_col = "subcluster")
# da_3d$subcluster = factor(da_3d$subcluster, levels=ss_subcluster_order)
# plotDAbeeswarm(da_3d, group.by = "subcluster", alpha=0.9997) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "HS 3d vs. LS (SD)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in HS 3d<br></span>
#     <span style='color:#3A3A98;'>Depleted in HS 3d</span>
#     </span>"
#   ) +
#   theme(plot.title = element_markdown(lineheight = 1.1),
#         legend.text = element_markdown(size = 11),
#         text=element_text(size=11))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # # testing
# i = input_file[1]
# milo_object = readRDS(i)
# angii_3d <- testNhoods(milo_object, design = ~ 0 + treatment, design.df = milo_design,
#                        fdr.weighting="graph-overlap", reduced.dim = "HARMONY",
#                        model.contrasts = c("Saline 3d - AngII 3d"))
# angii_28d <- testNhoods(milo_object, design = ~ 0 + Type, design.df = milo_design,
#                         fdr.weighting="graph-overlap", reduced.dim = "HARMONY",
#                         model.contrasts = c("Saline 3d - AngII 28d"))
# angii_3d %>%
#   arrange(SpatialFDR) %>%
#   head()
# angii_28d %>%
#   arrange(SpatialFDR) %>%
#   head()
# 
# # inspecting testing results
# ggplot(da_results, aes(PValue)) + geom_histogram(bins=50)
# ggplot(da_results, aes(logFC, -log10(SpatialFDR))) +
#   geom_point() +
#   geom_hline(yintercept = 1)
# 
# milo_object <- buildNhoodGraph(milo_object)
# 
# 
# ## assign umap value?
# ## Plot neighbourhood graph
# ### set higher alpha value when no significant result
# library(ggraph)
# milo_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-CKD/miloR/multiomics-CKD.14sample.miloR.k_45.harmony.rds")
# milo_object = readRDS("/xdisk/mliang1/qqiu/project/multiomics-CKD/miloR/multiomics-CKD.14sample.miloR.k_45.pca.rds")
# reduction = "PCA"
# dn_da <- testNhoods(milo_object, design = ~ 0 + Type, design.df = milo_design,
#                     fdr.weighting="graph-overlap", reduced.dim = reduction,
#                     model.contrasts = c("TypeNC - TypeDN"))
# hn_da <- testNhoods(milo_object, design = ~ 0 + Type, design.df = milo_design,
#                     fdr.weighting="graph-overlap", reduced.dim = reduction,
#                     model.contrasts = c("TypeNC - TypeHN"))
# dn_hn <- testNhoods(milo_object, design = ~ 0 + Type, design.df = milo_design,
#                     fdr.weighting="graph-overlap", reduced.dim = reduction,
#                     model.contrasts = c("TypeDN - TypeHN"))
# 
# milo_object <- buildNhoodGraph(milo_object)
# 
# 
# plotNhoodGraphDA_mod = function(milo_object, da, annotateText, title){
#   annotations <- data.frame(
#     xpos = c(-Inf, -Inf), ypos =  c(-Inf, -Inf),
#     hjustvar = c(0, 0), vjustvar = c(-3, -1),
#     annotateText = c("", ""))
#   annotations$annotateText = annotateText
#   plotNhoodGraphDA(milo_object, da, layout="UMAP", alpha=1.0) + labs(title=title) +
#     scale_edge_width(range = c(0, 0)) + scale_fill_gradient2(breaks=c(-5, 0, 5)) +
#     geom_text(data=annotations, aes(x=xpos,y=ypos,hjust=hjustvar, color=annotateText,
#                                     vjust=vjustvar,label=annotateText), size=6) +
#     scale_color_manual(values = c("#3A3A98", "#832424")) +
#     guides(fill=guide_colorbar(title="logFC"), edge_width=FALSE, size=FALSE, color=FALSE)
# }
# # plotNhoodGraphDA_mod(milo_object, dn_da, c("Enriched in DN", "Depleted in DN"), "DN vs NC")
# # plotNhoodGraphDA_mod(milo_object, hn_da, c("Enriched in HN", "Depleted in HN"), "HN vs NC")
# plotNhoodGraphDA_mod(milo_object, dn_hn[order(dn_hn$PValue, decreasing = T),],
#                      c("Enriched in HN", "Enriched in DN"), "HN vs DN")
# 
# 
# ## visualize da result
# dn_da <- annotateNhoods(milo_object, dn_da, coldata_col = "cluster_layer2")
# hn_da <- annotateNhoods(milo_object, hn_da, coldata_col = "cluster_layer2")
# dn_hn <- annotateNhoods(milo_object, dn_hn, coldata_col = "cluster_layer2")
# 
# dn_da$FDR = p.adjust(dn_da$PValue, "fdr")
# hn_da$FDR = p.adjust(hn_da$PValue, "fdr")
# dn_hn$FDR = p.adjust(dn_hn$PValue, "fdr")
# 
# # da_results$celltype <- ifelse(da_results$cluster_layer2_fraction < 0.7, "Mixed", da_results$cluster_layer2)
# 
# dn_da$cluster_layer2 = factor(dn_da$cluster_layer2, levels=levels(seurat_object$cluster_layer2))
# dn_da$SpatialFDR = p.adjust(dn_da$PValue, "fdr")
# plotDAbeeswarm(dn_da, group.by = "cluster_layer2", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(x="", y="Log Fold Change\nHigher in DN <- -> Higher in NC")  +
#   theme(axis.text.y=element_text(size=12, color="black"),
#         axis.text.x=element_text(size=12, color="black"),
#         axis.title.x=element_text(size=12, color="black")) +
#   geom_hline(yintercept=c(0), colour = "black")
# hn_da$cluster_layer2 = factor(hn_da$cluster_layer2, levels=levels(seurat_object$cluster_layer2))
# hn_da$SpatialFDR = p.adjust(hn_da$PValue, "fdr")
# plotDAbeeswarm(hn_da, group.by = "cluster_layer2", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(x="", y="Log Fold Change\nHigher in HN <- -> Higher in NC")  +
#   theme(axis.text.y=element_text(size=12, color="black"),
#         axis.text.x=element_text(size=12, color="black"),
#         axis.title.x=element_text(size=12, color="black")) +
#   geom_hline(yintercept=c(0), colour = "black")
# dn_hn$cluster_layer2 = factor(dn_hn$cluster_layer2, levels=levels(seurat_object$cluster_layer2))
# dn_hn$SpatialFDR = p.adjust(dn_hn$PValue, "fdr")
# plotDAbeeswarm(dn_hn, group.by = "cluster_layer2", alpha=0.1) + scale_x_discrete(limits=rev) +
#   labs(x="", y="Log Fold Change\nHigher in HN <- -> Higher in DN") +
#   theme(axis.text.y=element_text(size=12, color="black"),
#         axis.text.x=element_text(size=12, color="black"),
#         axis.title.x=element_text(size=12, color="black")) +
#   geom_hline(yintercept=c(0), colour = "black")
# 








