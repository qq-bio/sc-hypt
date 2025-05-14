library(Seurat)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/panglao_anno.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/cluster_sample_sum.R")



QC_harmony = function(seurat_object, 
                      nPC = 20, 
                      reso = c(seq(0.1, 1, 0.1), 1, 1.5, 2)){
  
  
  black_list = c(rownames(seurat_object)[grep("^Rp[sl][[:digit:]]", rownames(seurat_object))],
                 rownames(seurat_object)[grep("^Ig", rownames(seurat_object))])
  rm_black_list = setdiff(rownames(seurat_object), black_list)
  seurat_object = subset(seurat_object, nFeature_RNA>200 & 
                           nCount_RNA<20000 & 
                           percent.mt<20 & 
                           doublet_pc.2=="Singlet",
                         features=rm_black_list)
  
  ### cell cycle score & to keep the signals separating non-cycling and cycling cells; reference: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
  s.genes = c(str_to_title(cc.genes$s.genes), "Cenpu")
  g2m.genes = c(str_to_title(cc.genes$g2m.genes), "Pimreg", "Jpt1")
  seurat_object = CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat_object$CC.Difference <- seurat_object$S.Score - seurat_object$G2M.Score
  
  
  ### basic normalization
  seurat_object = NormalizeData(seurat_object)
  seurat_object = FindVariableFeatures(seurat_object)
  seurat_object = ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident")
  
  seurat_object = RunUMAP(seurat_object, reduction = "harmony", dims=1:nPC)
  seurat_object = FindNeighbors(seurat_object, reduction = "harmony", dims=1:nPC)
  seurat_object = FindClusters(seurat_object, resolution = reso)

  return(seurat_object)
  
}




QC_harmony_mo = function(seurat_object, 
                      nPC = 20, 
                      reso = c(seq(0.1, 1, 0.1), 1, 1.5, 2)){

  DefaultAssay(seurat_object) = "RNA"
  seurat_object = NormalizeData(seurat_object)
  seurat_object = FindVariableFeatures(seurat_object)
  seurat_object = ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "pca",
                             reduction.save = "harmony.rna",assay.use = "RNA", project.dim = F)
  
  ## process the ATAC data
  DefaultAssay(seurat_object) = "ATAC"
  seurat_object@assays$ATAC@counts
  seurat_object= RunTFIDF(seurat_object)
  seurat_object = FindTopFeatures(seurat_object, min.cutoff='q0')
  seurat_object = RunSVD(seurat_object)
  seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "lsi",
                             reduction.save = "harmony.atac", assay.use = "ATAC", project.dim = F)
  
  ### integration & cluster
  seurat_object = FindMultiModalNeighbors(seurat_object, reduction.list = list("harmony.rna", "harmony.atac"),
                                          dims.list = list(1:nPC, 2:nPC))
  seurat_object = RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", reduction.key = "wnnUMAPHarmony_")
  seurat_object = FindClusters(seurat_object, graph.name = "wsnn", algorithm = 3,
                               resolution = reso)
  
  DefaultAssay(seurat_object) = "RNA"
  return(seurat_object)
  
}







DE_analysis = function(seurat_object,
                       cluster = "seurat_clusters",
                       outfile){
  
  # seurat_object@active.ident = cluster
  Idents(seurat_object) = cluster
  
  # cluster_sample_tbl = table(seurat_object$seqID2, seurat_object@meta.data[, cluster])
  # cluster_sample_sum(cluster_sample_tbl, paste0(outfile, ".cluster_sample_sum.xlsx"))

  markers = FindAllMarkers(seurat_object, group.by=cluster, only.pos = TRUE, min.pct = 0.25)
  markers$pct.diff = markers$pct.1 - markers$pct.2
  markers = markers %>% group_by(cluster) %>%
      dplyr::arrange(desc(pct.diff), .by_group=TRUE)
  write.table(markers,file=paste0(outfile, ".allmarker.0.25.long.txt"), sep="\t")
  
  marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
    dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
    dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
    dplyr::select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                          v.names="gene", sep=" ", direction = "wide")
  colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
  write.table(marker_tbl,file=paste0(outfile, ".allmarker.0.25.wide.txt"), sep="\t", row.names = F)
  
  # panglao_anno(marker_tbl, paste0(outfile, ".allmarker.0.25.wide.xlsx"))
  
  return(markers)
  
}



sample_DE_analysis = function(seurat_object,
                       cluster = "seurat_clusters",
                       comparison = comparison,
                       outfile){
  
  seurat_object@active.ident = seurat_object@meta.data[, cluster]
  Idents(seurat_object) = cluster
  
  cluster_sample_tbl = table(seurat_object$seqID2, seurat_object@meta.data[, cluster])
  cluster_sample_sum(cluster_sample_tbl, paste0(outfile, ".cluster_sample_sum.xlsx"))
  
  markers = FindAllMarkers(seurat_object, group.by=cluster, only.pos = TRUE, min.pct = 0.25)
  markers$pct.diff = markers$pct.1 - markers$pct.2
  markers = markers %>% group_by(cluster) %>%
    dplyr::arrange(desc(pct.diff), .by_group=TRUE)
  write.table(markers,file=paste0(outfile, ".allmarker.0.25.long.txt"), sep="\t")
  
  marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
    dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
    dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
    select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                          v.names="gene", sep=" ", direction = "wide")
  colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
  write.table(marker_tbl,file=paste0(outfile, ".allmarker.0.25.wide.txt"), sep="\t", row.names = F)
  
  panglao_anno(marker_tbl, paste0(outfile, ".allmarker.0.25.wide.xlsx"))
  
  
}
