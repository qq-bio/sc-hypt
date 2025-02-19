.libPaths("/xdisk/mliang1/qqiu/R/library_Rv4.2_puma9/", include.site = F)

library(Seurat)
library(Signac)
library(sceasy)
library(reticulate)

reticulate::use_condaenv("cell2location", conda = "/opt/ohpc/pub/apps/anaconda/2022.05/bin/conda", required = TRUE)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

torch <- import("torch")
pd <- import("pandas")

### Load and merge datasets of EC
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

excluded_cell_types <- c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")

so_list <- list()
vg_list <- c()
for( i in input_file){
  
  name_tmp <- gsub(".EC.anno.rds", "", basename(i))
  so_tmp <- readRDS(i)
  so_tmp <- subset(so_tmp, subclass_level2 %in% setdiff(unique(so_tmp$subclass_level2), excluded_cell_types))
  
  so_tmp$species <- ifelse(grepl("mouse", basename(i)), "mouse", "rat")
  so_tmp$assays <- ifelse(grepl("LK", basename(i)), "snMultiome", "snRNA")
  
  so_tmp <- FindVariableFeatures(so_tmp, selection.method = "vst")

  so_list[[name_tmp]] <- so_tmp
  vg_list <- c(vg_list, head(VariableFeatures(so_tmp), 2000))
  
}


gene_lists <- lapply(so_list, function(obj) rownames(obj))
common_genes <- Reduce(intersect, gene_lists)

vg_list <- intersect(unique(vg_list), common_genes)

seurat_object <- merge(so_list[[1]], y = so_list[-1])
seurat_object$sxtxt <- paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)

seurat_object_vg <- subset(seurat_object, features = vg_list)
print(seurat_object_vg)


### Convert Seurat object to Anndata
adata <- convertFormat(seurat_object_vg, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)


### Training
scvi$model$SCVI$setup_anndata(adata, batch_key="strain",
                              categorical_covariate_keys=list("assays"), 
                              continuous_covariate_keys=list("nCount_RNA", "nFeature_RNA", "percent.mt"))

# model = scvi$model$SCVI(adata, n_layers =2)
# model <- scvi$model$SCVI(adata) # ec.scvi.rds
# model <- scvi$model$SCVI(adata, n_layers=as.integer(2), n_latent= as.integer(30), gene_likelihood="nb") # ec.scvi.rec.rds
model <- scvi$model$SCVI(adata, gene_likelihood="nb") # ec.scvi.gene_nb.rds

model$train(use_gpu = TRUE) # model$train(max_epochs = as.integer(400))

latent <- model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_object)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_2k.rds")
saveRDS(model, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_2k.model.rds")




so_list <- list()
vg_list <- c()
for( i in input_file){
  
  name_tmp <- gsub(".EC.anno.rds", "", basename(i))
  so_tmp <- readRDS(i)
  so_tmp <- subset(so_tmp, subclass_level2 %in% setdiff(unique(so_tmp$subclass_level2), excluded_cell_types))
  
  so_tmp$species <- ifelse(grepl("mouse", basename(i)), "mouse", "rat")
  so_tmp$assays <- ifelse(grepl("LK", basename(i)), "snMultiome", "snRNA")
  
  so_tmp <- FindVariableFeatures(so_tmp, selection.method = "vst")
  
  so_list[[name_tmp]] <- so_tmp
  vg_list <- c(vg_list, head(VariableFeatures(so_tmp), 3000))
  
}


gene_lists <- lapply(so_list, function(obj) rownames(obj))
common_genes <- Reduce(intersect, gene_lists)

vg_list <- intersect(unique(vg_list), common_genes)

seurat_object <- merge(so_list[[1]], y = so_list[-1])
seurat_object$sxtxt <- paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)

seurat_object_vg <- subset(seurat_object, features = vg_list)
print(seurat_object_vg)


### Convert Seurat object to Anndata
adata <- convertFormat(seurat_object_vg, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)


### Training
scvi$model$SCVI$setup_anndata(adata, batch_key="strain",
                              categorical_covariate_keys=list("assays"), 
                              continuous_covariate_keys=list("nCount_RNA", "nFeature_RNA", "percent.mt"))

# model = scvi$model$SCVI(adata, n_layers =2)
# model <- scvi$model$SCVI(adata) # ec.scvi.rds
# model <- scvi$model$SCVI(adata, n_layers=as.integer(2), n_latent= as.integer(30), gene_likelihood="nb") # ec.scvi.rec.rds
model <- scvi$model$SCVI(adata, gene_likelihood="nb") # ec.scvi.gene_nb.rds

model$train(use_gpu = TRUE) # model$train(max_epochs = as.integer(400))

latent <- model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_object)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_3k.rds")
saveRDS(model, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_3k.model.rds")




# 
# seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.rds")
# seurat_object_nb <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.rds")
# seurat_object_rec <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.rec.rds")
# 
# # ###Find clusters, then run UMAP, and visualize
# # seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
# # seurat_object <- FindClusters(seurat_object, resolution =1)
# 
# seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)
# DimPlot(seurat_object, group.by = "species", reduction = "umap", pt.size = 3)
# DimPlot(seurat_object, split.by = "tissue", reduction = "umap", pt.size = 3)
# FeaturePlot(seurat_object, features = c("Pecam1", "Egfl7", "Vwf", "Sema3g", "Flt1", "Plvap", "Pdpn", "Prox1", "Slc38a3"))
# 
# 
# seurat_object_nb <- RunUMAP(seurat_object_nb, dims = 1:10, reduction = "scvi", n.components = 2)
# DimPlot(seurat_object_nb, group.by = "species", reduction = "umap", pt.size = 3)
# DimPlot(seurat_object_nb, split.by = "tissue", reduction = "umap", pt.size = 3)
# FeaturePlot(seurat_object_nb, features = c("Pecam1", "Egfl7", "Vwf", # EC
#                                            "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
#                                            "Flt1",
#                                            "Plvap", # venous EC
#                                            "Rgcc", # capillary EC
#                                            "Ccl21", "Prox1", # lymphatic EC
#                                            "Igfbp5", # EC capillary
#                                            "Ehd3", # glomerular EC
#                                            "Slc38a3" # BBB EC
#                                            ))

# negative marker list: "Ackr1", "Gpihbp1", "Fcn3", "Ca4", "Lyve1", "Pdpn", "Sost",

DimPlot(seurat_object, group.by = "subclass_level2", split.by = "tissue", reduction = "umap", pt.size = 3)
DimPlot(seurat_object, group.by = "tissue", split.by = "subclass_level2", reduction = "umap", pt.size = 3)

FeaturePlot(seurat_object, features = c("Vwf", "Pecam1", "Egfl7", "Flt1", "Sema3g", "Postn", "Plvap", "Pdpn", "Prox1", "Slc38a3"))

FeaturePlot(seurat_object, features = c("Tagln", "Acta2", "Pdgfra", "Pdgfrb", "Notch3", "Myh11", "Kcnj8")) # mural cell
FeaturePlot(seurat_object, features = c("Smoc1", "Inhba", "Npr3")) # endocardial_ec

FeaturePlot(seurat_object, features = c("Slco1a2", "Slc38a3", "Slc7a8"))
FeaturePlot(seurat_object, features = c("Cdh5", "Sema3g", "Rgcc", "Ackr1", "Ccl21"))

FeaturePlot(seurat_object, features = c("Mecom"))



