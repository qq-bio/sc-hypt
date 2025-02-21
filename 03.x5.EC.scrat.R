.libPaths("/xdisk/mliang1/qqiu/R/library_Rv4.2_puma9/", include.site = F)

library(scrat)
# library(Seurat)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scrat/")


# ### downsample
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
#   # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
#   # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds", 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds", 
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
# )
# 
# excluded_cell_types <- c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")
# 
# so_list <- list()
# vg_list <- c()
# for( i in input_file){
#   
#   name_tmp <- gsub(".EC.anno.rds", "", basename(i))
#   so_tmp <- readRDS(i)
#   so_tmp <- subset(so_tmp, subclass_level2 %in% setdiff(unique(so_tmp$subclass_level2), excluded_cell_types))
#   
#   so_tmp$species <- ifelse(grepl("mouse", basename(i)), "mouse", "rat")
#   so_tmp$assays <- ifelse(grepl("LK", basename(i)), "snMultiome", "snRNA")
#   
#   so_tmp <- FindVariableFeatures(so_tmp, selection.method = "vst")
#   
#   so_list[[name_tmp]] <- so_tmp
#   vg_list <- c(vg_list, head(VariableFeatures(so_tmp), 2000))
#   
# }
# 
# 
# gene_lists <- lapply(so_list, function(obj) rownames(obj))
# common_genes <- Reduce(intersect, gene_lists)
# 
# vg_list <- intersect(unique(vg_list), common_genes)
# 
# 
# seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
# 
# set.seed(42)
# 
# frac <- 0.5
# subset_cells <- sample(colnames(seurat_object), size = floor(frac * ncol(seurat_object)))
# seurat_subset <- subset(seurat_object, cells = subset_cells, features = vg_list)
# saveRDS(seurat_subset, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.subset.hvg.0.5.rds")
# 
# frac <- 0.25
# subset_cells <- sample(colnames(seurat_object), size = floor(frac * ncol(seurat_object)))
# seurat_subset <- subset(seurat_object, cells = subset_cells, features = vg_list)
# saveRDS(seurat_subset, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.subset.hvg.0.25.rds")


# seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.subset.hvg.0.25.rds")

env <- scrat.new(list(dataset.name="Cross-organ EC"))
# recommend to use expression values in logarithmic scale
env$indata <- seurat_object@assays$RNA@data
env$group.labels <- seurat_object$seurat_clusters

scrat.run(env)
saveRDS(env, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/scrat/ec.scvi.subset.hvg.0.25.scrat.sbatch.rds")

