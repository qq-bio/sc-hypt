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


# 
# ###############################################################################
# ## 
# seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
# 
# cluster = "seurat_clusters"
# 
# seurat_object@active.ident = seurat_object@meta.data[, cluster]
# Idents(seurat_object) = cluster
# 
# meta_table = seurat_object@meta.data
# species_list <- unique(meta_table$strain)
# 
# # Loop through each species
# for (si in species_list) {
#   
#   # Subset Seurat object by strain
#   seurat_object_use <- subset(seurat_object, strain == si)
#   
#   # Adjust strain label for mouse
#   if (si == "C57BL/6") { si <- "mouse" }
#   
#   # Set output file path
#   outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/",
#                     si, ".miloR.rds")
#   
#   # Set dimensionality reduction method based on tissue type
#   reduced.dim <- "SCVI"
#   
#   # Construct KNN graph and identify representative neighborhoods
#   sce <- as.SingleCellExperiment(seurat_object_use, assay = "RNA")
#   milo_object <- Milo(sce)
#   milo_object <- buildGraph(milo_object, k = 15, d = 5, reduced.dim = reduced.dim)
#   milo_object <- makeNhoods(milo_object, prop = 0.2, k = 15, d = 5, refined = TRUE, 
#                             reduced_dims = reduced.dim, refinement_scheme = "graph")
#   
#   # Visualize neighborhood size histogram
#   plotNhoodSizeHist(milo_object) # Ideal average size: 5 x N_samples / 50-100
#   
#   # Count cells within each neighborhood and create design matrix
#   milo_object <- countCells(milo_object, meta.data = seurat_object_use@meta.data, sample = "orig.ident")
#   
#   # Create and format design matrix for Milo analysis
#   milo_design <- distinct(seurat_object_use@meta.data[, c("orig.ident", "treatment")])
#   rownames(milo_design) <- milo_design$orig.ident
#   milo_design$orig.ident <- as.factor(milo_design$orig.ident)
#   milo_design$treatment <- as.factor(milo_design$treatment)
#   
#   # Calculate neighborhood connectivity
#   milo_object <- calcNhoodDistance(milo_object, d = 5, reduced.dim = reduced.dim)
#   
#   # Save the Milo object for downstream analysis
#   saveRDS(milo_object, outfile)
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
# 
# 
# ################################################################################
# input_file = c(
# 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/mouse.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/SS.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/SD.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/SHR.miloR.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/WKY.miloR.rds"
#   
# )
# 
# 
# ################################################################################
# ### pairwise comparison
# library(ggtext)
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR")
# 
# 
# da_merge = c()
# for( i in input_file ){
# 
#   # tissue = strsplit(basename(i), "\\.")[[1]][2]
#   # tissue = sub(".*\\.([^.]*)\\.(miloR).*", "\\1", basename(i))
#   reduction = "SCVI"
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
#   treatment_list = intersect(treatment_order, unique(milo_design$treatment))
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
#                           fdr.weighting="k-distance", reduced.dim = reduction,
#                           model.contrasts = model.contrasts), silent = T)
#       
#       if(class(da_tmp) != "try-error"){
# 
#         da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = "seurat_clusters")
#         da_tmp$subcluster = da_tmp$seurat_clusters
#         da_tmp$species = si
#         # da_tmp$tissue = tissue
#         da_tmp$control = control
#         da_tmp$treatment = treatment
#         da_merge = rbind(da_merge, da_tmp)
# 
#       }else{
#         print(i)
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


################################################################################
### strain-wise comparison

seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/EC.all.miloR.rds")

cluster = "seurat_clusters"

seurat_object@active.ident = seurat_object@meta.data[, cluster]
Idents(seurat_object) = cluster

meta_table = seurat_object@meta.data
species_list <- unique(meta_table$strain)

# Set dimensionality reduction method based on tissue type
reduced.dim <- "SCVI"

# Construct KNN graph and identify representative neighborhoods
sce <- as.SingleCellExperiment(seurat_object, assay = "RNA")
milo_object <- Milo(sce)
milo_object <- buildGraph(milo_object, k = 30, d = 10, reduced.dim = reduced.dim)
milo_object <- makeNhoods(milo_object, prop = 0.1, k = 30, d = 10, refined = TRUE, 
                          reduced_dims = reduced.dim, refinement_scheme = "graph")

# Visualize neighborhood size histogram
plotNhoodSizeHist(milo_object) # Ideal average size: 5 x N_samples / 50-100

# Count cells within each neighborhood and create design matrix
milo_object <- countCells(milo_object, meta.data = seurat_object@meta.data, sample = "orig.ident")

# Create and format design matrix for Milo analysis
milo_design <- distinct(seurat_object@meta.data[, c("orig.ident", "treatment")])
rownames(milo_design) <- milo_design$orig.ident
milo_design$orig.ident <- as.factor(milo_design$orig.ident)
milo_design$treatment <- as.factor(milo_design$treatment)

# Calculate neighborhood connectivity
milo_object <- calcNhoodDistance(milo_object, d = 10, reduced.dim = reduced.dim)

# Save the Milo object for downstream analysis
saveRDS(milo_object, outfile)





input_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/EC.all.miloR.rds"


reduction = "SCVI"

milo_object = readRDS(input_file)
milo_object$sxt = gsub("[/ ]", "", paste0(milo_object$strain, "_", milo_object$treatment))

milo_object = countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample="orig.ident")
milo_design = data.frame(colData(milo_object))[,c("treatment", "strain", "orig.ident", "sxt")]
milo_design = distinct(milo_design)
rownames(milo_design) = milo_design$orig.ident
milo_design$treatment.new = as.factor(gsub(" ", "", milo_design$treatment))
milo_design$strain = as.factor(milo_design$strain)
milo_design$sxt = as.factor(milo_design$sxt)

comp_list = list(c("C57BL6_Saline3d", "C57BL6_AngII3d"),
                 c("C57BL6_Saline3d", "C57BL6_AngII28d"),
                 c("SD_LS", "SD_HS3d"),
                 c("SS_LS", "SS_HS3d"),
                 c("SS_LS", "SS_HS21d"),
                 c("WKY_10w", "WKY_26w"),
                 c("SHR_10w", "SHR_26w")
)

da_merge = c()
for( i in comp_list ){

  control = i[1]
  treatment = i[2]

  model.contrasts = paste0("sxt", control, " - ", "sxt", treatment)

  da_tmp = try(testNhoods(milo_object, design = ~ 0 + sxt, design.df = milo_design,
                          fdr.weighting="k-distance", reduced.dim = reduction,
                          model.contrasts = model.contrasts), silent = T)

  if(class(da_tmp) != "try-error"){

    da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = "seurat_clusters")
    da_tmp$SpatialFDR.org = da_tmp$SpatialFDR
    da_tmp$SpatialFDR = p.adjust(da_tmp$PValue, "fdr")
    da_tmp$subcluster = da_tmp$seurat_clusters
    # da_tmp$species = pi
    # da_tmp$tissue = tissue
    da_tmp$control = control
    da_tmp$treatment = treatment
    da_merge = rbind(da_merge, da_tmp)

  }else{
    print(i)
  }

}


write.table(da_merge, "EC.all.milo.da_result.out", sep='\t', quote=F, col.names = T, row.names = F)









################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/miloR/")
da_merge = read.table("EC.all.milo.da_result.out", sep='\t', header=T)
da_merge$species = gsub("_.*", "", da_merge$control)
da_merge[da_merge$species=="C57BL6", ]$species = "C57BL/6"
da_merge$species = factor(da_merge$species, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
da_merge$treatment = gsub(".*_", "", da_merge$treatment)
da_merge$treatment = factor(da_merge$treatment, levels = c("Saline3d", "AngII3d", "AngII28d", "10w", "26w", "LS", "HS3d", "HS21d"))
da_merge$SpatialFDR.org = da_merge$SpatialFDR
da_merge$SpatialFDR = p.adjust(da_merge$PValue, "fdr")
da_merge_plot <- da_merge %>%
  mutate(order = case_when(
    SpatialFDR > 0.05 ~ 1,
    SpatialFDR <= 0.05 ~ 2
  ))
da_merge_plot <- da_merge_plot %>% arrange(order)

group.by = "seurat_clusters"; alpha=0.1
da.res = mutate(da_merge_plot, group_by = da_merge_plot[, group.by])
da.res %>% mutate(is_signif = ifelse(SpatialFDR.org < alpha, 1, 0)) %>%
  mutate(logFC_color = ifelse(is_signif == 1, as.numeric(logFC), NA)) %>%
  arrange(group_by) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
  ggplot(aes(as.character(group_by), -logFC, color = logFC_color)) + guides(color = "none") +
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
  facet_nested( ~ species + treatment , scales = "free", space = "free") +
  geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
  geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
  scale_color_gradient2( na.value = "lightgrey" )
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.milo.png", width=731/96, height=957/96, dpi=300)


# 
# 
# 
# 
# 
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/")
# da_merge = read.table("milo.strain_wise.da_result.out", sep='\t', header=T)
# da_merge$subcluster = factor(da_merge$subcluster, levels = cell_order)
# da_merge$species = factor(da_merge$species, levels = c("Salt-sensitive", "Spontaneous"))
# da_merge$tissue = factor(da_merge$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
# da_merge$treatment = factor(da_merge$treatment, levels = c("SS - LS", "SS - HS 3d", "SHR - 10w", "SHR - 26w"))
# # da_merge$SpatialFDR.org = da_merge$SpatialFDR
# # da_merge$SpatialFDR = p.adjust(da_merge$PValue, "fdr")
# da_merge_plot <- da_merge %>%
#   mutate(order = case_when(
#     SpatialFDR > 0.05 ~ 1,
#     SpatialFDR <= 0.05 ~ 2
#   ))
# da_merge_plot <- da_merge_plot %>% arrange(order)
# 
# group.by = "subcluster"; alpha=0.1
# da.res = mutate(da_merge_plot, group_by = da_merge_plot[, group.by])
# da.res %>% mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
#   mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
#   arrange(group_by) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
#   ggplot(aes(group_by, -logFC, color = logFC_color)) + guides(color = "none") +
#   xlab(group.by) + ylab("Log Fold Change") + geom_quasirandom(alpha = 1) +
#   coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
#   labs(
#     x = "", y = "Log Fold Change",
#     title = "Strain - treatment (vs. control)<br>
#     <span style='font-size:11pt'>
#     <span style='color:#832424;'>Enriched in hypertensive strain<br></span>
#     <span style='color:#3A3A98;'>Depleted in hypertensive strain</span>
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
#   facet_nested(tissue ~ species + treatment , scales = "free", space = "free") +
#   geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
#   geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
#   stat_summary(fun = median, geom = "point", color = "black", size = 2) +
#   scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
#   scale_color_gradient2( na.value = "lightgrey" )
# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/strain_wise.milo.png", width=456/96, height=793/96, dpi=300)
# 
# 








