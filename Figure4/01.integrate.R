library(Seurat)
library(Signac)
library(sceasy)
library(reticulate)


sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)

torch <- import("torch")
pd <- import("pandas")




################################################################################
### load and merge datasets of EC
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
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
  vg_list <- c(vg_list, head(VariableFeatures(so_tmp), 1000))
  
}


gene_lists <- lapply(so_list, function(obj) rownames(obj))
common_genes <- Reduce(intersect, gene_lists)

vg_list <- intersect(unique(vg_list), common_genes)

seurat_object <- merge(so_list[[1]], y = so_list[-1])
seurat_object$sxtxt <- paste0(seurat_object$strain, "-", seurat_object$treatment, "-", seurat_object$tissue)

seurat_object_vg <- subset(seurat_object, features = vg_list)
print(seurat_object_vg)


### convert Seurat object to Anndata
adata <- convertFormat(seurat_object_vg, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
print(adata)


### training
scvi$model$SCVI$setup_anndata(adata, batch_key="strain",
                              categorical_covariate_keys=list("assays"), 
                              continuous_covariate_keys=list("nCount_RNA", "nFeature_RNA", "percent.mt"))

model <- scvi$model$SCVI(adata, gene_likelihood="nb")

model$train(use_gpu = TRUE)

latent <- model$get_latent_representation()

latent <- as.matrix(latent)
rownames(latent) = colnames(seurat_object)
seurat_object[["scvi"]] <- CreateDimReducObject(embeddings = latent, key = "scvi_", assay = DefaultAssay(seurat_object))


seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 2))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)

Idents(seurat_object) = "RNA_snn_res.1.5"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.1.5


marker_list= c("Nav3",'Ablim3',"Ccdc85a",'Ncald',"Nrp1",'Kitlg',"Dach1","Arhgap18",'Ank3',
               "Hdac9","Nrp2","Lamc1",'Ano4',
               "Nav2",'Diaph3',"Mki67","Rad51b","Sdk1",
               'Fhod3',"Unc5c","Il1r1","Hmcn1","Vwf",
               "Myo10","Lrrc3b","Btnl9",'Kcnt2',
               "Nebl",'St6galnac3',
               "Prdm16","Ptprj","Sulf1","Auts2l1",'Eln',"Myof","Pcdh7","Alcam",
               "Slco1a4","Lef1","Gpcpd1",'Ccdc141',
               "Zfp521","Nuak1","Mctp1",
               'Pkhd1l1','Cgnl1',"Cdh11",
               "Meis2","Pbx1","Ldb2",
               "Rbms3","Inpp4b","Fmnl2","Chrm3","Malat1",
               'Pkhd1',"Erbb4",'Ca12')

DotPlot(seurat_object, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


seurat_object$seurat_clusters <- paste0("C", seurat_object$seurat_clusters)
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C5", "C8", "C13"), ]$seurat_clusters <- "M5813"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C2", "C4"), ]$seurat_clusters <- "M24"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C0", "C6", "C10"), ]$seurat_clusters <- "M0610"


### rename idx
cell_by_size <- table(seurat_object$seurat_clusters)
cell_by_size <- cell_by_size[order(cell_by_size, decreasing = T)]
new_idx <- as.character(1:19)
names(new_idx) <- names(cell_by_size)

new_idx <- factor(new_idx, levels = new_idx[ec_order])
seurat_object$new_idx <- as.character(new_idx[as.character(seurat_object$seurat_clusters)])
seurat_object$new_idx <- factor(seurat_object$new_idx, levels = new_idx[ec_order])



saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.merged.rds")

