library(CellChat)
library(Seurat)



setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/")


safe_computeCentrality <- function(cellchat) {
  tryCatch({
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    return(cellchat)  # Return the updated object
  }, 
  warning = function(w) {
    return(cellchat)  # Return the original object if there's a warning
  }, 
  error = function(e) {
    return(cellchat)  # Return the original object if there's an error
  })
}

# 
# input <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
# output <- "cross_organ.EC.cellchat.rds"
# ec_seurat_object <- readRDS(input)
# ec_seurat_object$seurat_clusters <- paste0("EC", ec_seurat_object$seurat_clusters)
# ec_meta_tbl <- ec_seurat_object@meta.data
# 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
#   
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
# )
# 
# CellChatDB = CellChatDB.mouse
# #future::plan("multicore", workers = 4) 
# 
# 
# e <- new.env()
# for(i in input_file){
#   
#   tissue = gsub("\\.(RNA|multiomics)\\.anno\\.v2\\.rds", "", basename(i))
#   
#   seurat_object = readRDS(i)
#   
#   seurat_object$subclass_level2.1 = seurat_object$subclass_level2
#   EC_list = rownames(seurat_object@meta.data[seurat_object$subclass_level1=="EC", ])
#   discard_list = setdiff(EC_list, ec_meta_tbl$cell_id)
#   remain_list = setdiff(seurat_object$cell_id, discard_list)
#   
#   
#   seurat_object = subset(seurat_object, cell_id %in% remain_list)
#   seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$subclass_level2.1 = ec_meta_tbl[seurat_object@meta.data[seurat_object$cell_id %in% EC_list, ]$cell_id, ]$seurat_clusters
#   
#   
#   Idents(seurat_object) <- "subclass_level2.1"
#   seurat_object@active.ident = factor(seurat_object@meta.data[, "subclass_level2.1"])
#   Idents(seurat_object) = "subclass_level2.1"
#   
#   strain_list = unique(seurat_object$strain)
#   for(si in strain_list){
#     
#     treatment_list = unique(seurat_object@meta.data[seurat_object$strain==si, ]$treatment)
#     
#     for(ti in treatment_list){
#       
#       seurat_object_tmp = subset(seurat_object, strain==si & treatment==ti)
#       cellchat = createCellChat(seurat_object_tmp, group.by = "subclass_level2.1")
#       cellchat@DB = CellChatDB
#       cellchat = subsetData(cellchat)
#       cellchat <- identifyOverExpressedGenes(cellchat)
#       cellchat <- identifyOverExpressedInteractions(cellchat)
#       cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
#       cellchat <- filterCommunication(cellchat, min.cells = 10)
#       cellchat <- computeCommunProbPathway(cellchat)
#       cellchat <- aggregateNet(cellchat)
#       
#       # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
#       cellchat <- safe_computeCentrality(cellchat)
#       
#       output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
#       
#       with(e, {
#         assign(output_name, cellchat)
#       })
#       
#     }
#     
#   }
#   
# }
# 
# saveRDS(e, output)
# 
# 
# 
# 
# 






input <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/ec.scvi.gene_nb.hvg_1k.refined.merged.rds"
output <- "cross_organ.EC.refined.merged.cellchat.rds"
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
#future::plan("multicore", workers = 4) 


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
      
      # cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      cellchat <- safe_computeCentrality(cellchat)
      
      output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
      
      with(e, {
        assign(output_name, cellchat)
      })
      
    }
    
  }
  
}

saveRDS(e, output)





################################################################################



model <- c("AngII", "Salt-sensitive", "Spontaneous")
names(model) <- c("mouse", "rat.ss", "rat.sp")

cellchat_list <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.rds")
cc_df <- data.frame()

for (output_name in names(cellchat_list)) {
  
  cellchat <- cellchat_list[[output_name]]
  
  # Extract tissue, strain, and treatment information
  parts <- strsplit(output_name, "_")[[1]]
  tissue <- parts[1]
  strain <- gsub("\\.", "/", parts[2])
  treatment <- parts[3]
  
  # Extract ligand-receptor interactions
  LR_list <- rownames(cellchat@LR$LRsig)
  
  for (lr in LR_list) {
    prob_tmp <- cellchat@net$prob[, , lr]
    
    lr_df_tmp <- reshape2::melt(prob_tmp, value.name = "value")
    colnames(lr_df_tmp)[1:2] <- c("source", "target")
    lr_df_tmp$LR_pair <- lr
    lr_df_tmp <- lr_df_tmp[lr_df_tmp$value > 0, ]
    
    if (nrow(lr_df_tmp) > 0) {
      # Add metadata and append to results
      lr_df_tmp$model <- model[gsub("(.*)\\.([A-Z]+).*", "\\1", tissue, perl = TRUE)]
      lr_df_tmp$tissue <- gsub("(.*)\\.([A-Z]+).*", "\\2", tissue, perl = TRUE)
      lr_df_tmp$strain <- strain
      lr_df_tmp$treatment <- treatment
      lr_df_tmp <- cbind(lr_df_tmp, cellchat@LR$LRsig[lr, , drop = FALSE])
      cc_df <- rbind(cc_df, lr_df_tmp)
    }
  }
}


write.table(cc_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/cellchat/cross_organ.EC.refined.merged.cellchat.result.out", col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)



cc_df = cc_df[order(cc_df$value, decreasing = T), ]

cc_df[cc_df$source=="ECC14" & cc_df$tissue=="MCA", ]


cc_df[cc_df$source=="ECM0610" & cc_df$tissue=="LV", ]
cc_df[cc_df$source=="ECM0610" & cc_df$tissue=="LV" & cc_df$LR_pair %in% c("APP_CD74"), ]














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










