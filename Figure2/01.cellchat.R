
library(CellChat)
library(Seurat)


################################################################################
### CellChat Analysis for Multiomics Hypertension Project
################################################################################

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/")

input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

# Set up CellChat database and parallel plan
CellChatDB <- CellChatDB.mouse
future::plan("multicore", workers = 4)

# Store cellchat objects in a list
cellchat_list <- list()

for (file_path in input_file) {
  
  tissue <- gsub("\\.(RNA|multiomics)\\.anno\\.v2\\.rds", "", basename(file_path))
  seurat_object <- readRDS(file_path)
  
  # Set cluster identities
  Idents(seurat_object) <- "subclass_level2"
  
  # Iterate through strains and treatments
  for (strain in unique(seurat_object$strain)) {
    treatment_list <- unique(seurat_object@meta.data[seurat_object$strain == strain, "treatment"])
    
    for (treatment in treatment_list) {
      
      # Subset data and create CellChat object
      seurat_object_tmp <- subset(seurat_object, strain == strain & treatment == treatment)
      if (ncol(seurat_object_tmp) == 0) next
      
      cellchat <- createCellChat(seurat_object_tmp, group.by = "subclass_level2")
      cellchat@DB <- CellChatDB
      
      # Perform CellChat pipeline steps
      cellchat <- cellchat %>%
        subsetData() %>%
        identifyOverExpressedGenes() %>%
        identifyOverExpressedInteractions() %>%
        computeCommunProb(raw.use = TRUE) %>%
        filterCommunication(min.cells = 10) %>%
        computeCommunProbPathway() %>%
        aggregateNet() %>%
        netAnalysis_computeCentrality(slot.name = "netP")
      
      output_name <- paste0(tissue, "_", gsub("/", ".", strain), "_", treatment)
      cellchat_list[[output_name]] <- cellchat
      
    }
  }
}

saveRDS(cellchat_list, "cellchat.rds")




################################################################################
### Process CellChat Results and Export for Visualization
################################################################################

model <- c("AngII", "Salt-sensitive", "Spontaneous")
names(model) <- c("mouse", "rat.ss", "rat.sp")

cellchat_list <- readRDS("cellchat.rds")
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


write.table(cc_df, "cellchat.result.out", col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)


