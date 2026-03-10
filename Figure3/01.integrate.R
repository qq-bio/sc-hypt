library(Seurat)
library(Signac)
library(sceasy)
library(reticulate)


sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

torch <- import("torch")
pd <- import("pandas")



################################################################################
### Extract EC across datasets
################################################################################

LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
print(ncol(LK_EC))
saveRDS(LK_EC, "mouse.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
print(ncol(LV_EC))
saveRDS(LV_EC, "mouse.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
print(ncol(MCA_EC))
saveRDS(MCA_EC, "mouse.MCA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC", "E/P transition cell"))
print(ncol(HYP_EC))
saveRDS(HYP_EC, "mouse.HYP.EC.cluster.rds")



LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds"
MSA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
print(ncol(LK_EC))
saveRDS(LK_EC, "rat.ss.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
print(ncol(LV_EC))
saveRDS(LV_EC, "rat.ss.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
print(ncol(MCA_EC))
saveRDS(MCA_EC, "rat.ss.MCA.EC.cluster.rds")

MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, subclass_level1 %in% c("EC"))
print(ncol(MSA_EC))
saveRDS(MSA_EC, "rat.ss.MSA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC", "E/P transition cell"))
print(ncol(HYP_EC))
saveRDS(HYP_EC, "rat.ss.HYP.EC.cluster.rds")




LK_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds"
MCA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
MSA_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds"
LV_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds"
HYP_infile="/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds"

LK_object = readRDS(LK_infile); DefaultAssay(LK_object) = "RNA"
MCA_object = readRDS(MCA_infile)
MSA_object = readRDS(MSA_infile)
LV_object = readRDS(LV_infile)
HYP_object = readRDS(HYP_infile)

LK_object@active.ident <- factor(LK_object@active.ident)
names(LK_object@active.ident) <- colnames(LK_object)
LK_EC = subset(LK_object, subclass_level1 %in% c("EC"))
print(ncol(LK_EC))
saveRDS(LK_EC, "rat.sp.LK.EC.cluster.rds")

LV_object@active.ident <- factor(LV_object@active.ident)
names(LV_object@active.ident) <- colnames(LV_object)
LV_EC = subset(LV_object, subclass_level1 %in% c("EC"))
print(ncol(LV_EC))
saveRDS(LV_EC, "rat.sp.LV.EC.cluster.rds")

MCA_object@active.ident <- factor(MCA_object@active.ident)
names(MCA_object@active.ident) <- colnames(MCA_object)
MCA_EC = subset(MCA_object, subclass_level1 %in% c("EC"))
print(ncol(MCA_EC))
saveRDS(MCA_EC, "rat.sp.MCA.EC.cluster.rds")

MSA_object@active.ident <- factor(MSA_object@active.ident)
names(MSA_object@active.ident) <- colnames(MSA_object)
MSA_EC = subset(MSA_object, subclass_level1 %in% c("EC"))
print(ncol(MSA_EC))
saveRDS(MSA_EC, "rat.sp.MSA.EC.cluster.rds")

HYP_object@active.ident <- factor(HYP_object@active.ident)
names(HYP_object@active.ident) <- colnames(HYP_object)
HYP_EC = subset(HYP_object, subclass_level1 %in% c("EC", "E/P transition cell"))
print(ncol(HYP_EC))
saveRDS(HYP_EC, "rat.sp.HYP.EC.cluster.rds")





################################################################################
### Load and merge datasets of EC
################################################################################
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.cluster.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.cluster.rds"
)


anno_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/EC.anno.rm_db.txt", header=T, sep='\t')

so_list <- list()
vg_list <- c()
for( i in input_file){
  
  dataset = gsub(".cluster.rds", "", basename(i))
  EC_object = readRDS(i)
  
  anno_info_use = anno_info[anno_info$dataset == dataset, ]
  reso = unique(anno_info_use$resolution)
  cluster_use = anno_info_use[anno_info_use$note!="remove", ]$cluster
  
  Idents(EC_object) = reso
  EC_object$seurat_clusters = EC_object@meta.data[, reso]
  
  if(length(unique(EC_object$seurat_clusters))==nrow(anno_info_use)){
    
    EC_object <- subset(EC_object, seurat_clusters %in% cluster_use)
    
    EC_object$species <- ifelse(grepl("mouse", basename(i)), "mouse", "rat")
    EC_object$assays <- ifelse(grepl("LK", basename(i)), "snMultiome", "snRNA")
    
    EC_object <- FindVariableFeatures(EC_object, selection.method = "vst")
    
    so_list[[dataset]] <- EC_object
    vg_list <- c(vg_list, head(VariableFeatures(EC_object), 2000))
    
  }else{
    sprintf("Clusters don't match: %s", i)
  }
  
}


gene_lists <- lapply(so_list, function(obj) rownames(obj))
common_genes <- Reduce(intersect, gene_lists)

vg_list <- intersect(unique(vg_list), common_genes)

seurat_object <- merge(so_list[[1]], y = so_list[-1])
seurat_object$sxtxt <- paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)

seurat_object <- JoinLayers(
  object = seurat_object
)

seurat_object_vg <- subset(seurat_object, features = vg_list)
print(seurat_object_vg)




################################################################################
### SCVI
################################################################################
### Convert Seurat object to Anndata
assay_name <- DefaultAssay(seurat_object_vg)
seurat_object_vg[[assay_name]] <- as(object = seurat_object_vg[[assay_name]], Class = "Assay")
print(class(seurat_object_vg[[assay_name]])) 
adata <- sceasy::convertFormat(seurat_object_vg, from="seurat", to="anndata", 
                               main_layer="counts", drop_single_values=FALSE)
print(adata)


### Training
scvi$model$SCVI$setup_anndata(adata, batch_key="orig.ident",
                              categorical_covariate_keys=list(list("strain","tissue")))

model <- scvi$model$SCVI(adata, gene_likelihood="nb")

model$train(use_gpu = TRUE) # model$train(max_epochs = as.integer(400))

latent <- model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_object)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.rds")
