library(CellChat)
library(Seurat)



setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/")


input <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.rds"
output <- "cross_organ.EC.cellchat.rds"
ec_seurat_object <- readRDS(input)
ec_seurat_object$seurat_clusters <- paste0("EC", ec_seurat_object$seurat_clusters)
ec_meta_tbl <- ec_seurat_object@meta.data

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

CellChatDB = CellChatDB.mouse
future::plan("multicore", workers = 4) 


e <- new.env()
for(i in input_file){
  
  tissue = gsub("\\.(RNA|multiomics)\\.anno\\.v2\\.rds", "", basename(i))
  
  seurat_object = readRDS(i)
  
  seurat_object$subclass_level2.1 = seurat_object$subclass_level2
  EC_list = rownames(seurat_object@meta.data[seurat_object$subclass_level1=="EC", ])
  discard_list = setdiff(EC_list, ec_meta_tbl$cell_id)
  remain_list = setdiff(seurat_object$cell_id, discard_list)
  
  
  seurat_object = subset(seurat_object, cell_id %in% remain_list)
  seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$subclass_level2.1 = ec_meta_tbl[seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$cell_id, ]$seurat_clusters
  
  
  Idents(seurat_object) <- "subclass_level2.1"
  seurat_object@active.ident = factor(seurat_object@meta.data[, "subclass_level2.1"])
  Idents(seurat_object) = "subclass_level2.1"
  
  strain_list = unique(seurat_object$strain)
  for(si in strain_list){
    
    treatment_list = unique(seurat_object@meta.data[seurat_object$strain==si, ]$treatment)
    
    for(ti in treatment_list){
      
      seurat_object_tmp = subset(seurat_object, strain==si & treatment==ti)
      cellchat = createCellChat(seurat_object_tmp, group.by = "subclass_level2.1")
      cellchat@DB = CellChatDB
      cellchat = subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      
      output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
      
      with(e, {
        assign(output_name, cellchat)
      })
      
    }
    
  }
  
}

saveRDS(e, output)











input <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds"
output <- "cross_organ.EC.refined.cellchat.rds"
ec_seurat_object <- readRDS(input)
ec_seurat_object$seurat_clusters <- paste0("EC", ec_seurat_object$seurat_clusters)
ec_meta_tbl <- ec_seurat_object@meta.data

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

CellChatDB = CellChatDB.mouse
future::plan("multicore", workers = 4) 


e <- new.env()
for(i in input_file){
  
  tissue = gsub("\\.(RNA|multiomics)\\.anno\\.v2\\.rds", "", basename(i))
  
  seurat_object = readRDS(i)
  
  seurat_object$subclass_level2.1 = seurat_object$subclass_level2
  EC_list = rownames(seurat_object@meta.data[seurat_object$subclass_level1=="EC", ])
  discard_list = setdiff(EC_list, ec_meta_tbl$cell_id)
  remain_list = setdiff(seurat_object$cell_id, discard_list)
  
  
  seurat_object = subset(seurat_object, cell_id %in% remain_list)
  seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$subclass_level2.1 = ec_meta_tbl[seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$cell_id, ]$seurat_clusters
  
  
  Idents(seurat_object) <- "subclass_level2.1"
  seurat_object@active.ident = factor(seurat_object@meta.data[, "subclass_level2.1"])
  Idents(seurat_object) = "subclass_level2.1"
  
  strain_list = unique(seurat_object$strain)
  for(si in strain_list){
    
    treatment_list = unique(seurat_object@meta.data[seurat_object$strain==si, ]$treatment)
    
    for(ti in treatment_list){
      
      seurat_object_tmp = subset(seurat_object, strain==si & treatment==ti)
      cellchat = createCellChat(seurat_object_tmp, group.by = "subclass_level2.1")
      cellchat@DB = CellChatDB
      cellchat = subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      
      output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
      
      with(e, {
        assign(output_name, cellchat)
      })
      
    }
    
  }
  
}

saveRDS(e, output)



