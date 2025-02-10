dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(SeuratDisk)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)
library(DoubletFinder)
library(parallel)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder")


## process each sample
sample_list = c( paste0("MLK", c(1:6)),  paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92",
                 paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)))

sample_list = c("RMCA7SN", "RMCA8SN")

para_list = data.frame(sample_ID = character(),
                       n_cell = numeric(),
                       n_doublet = numeric(),
                       doublet_rate = numeric(),
                       PC = numeric(),
                       pK = numeric())

for(sample_id in sample_list){
  ## use cellranger results
  # h5_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/", sample_id, "/outs/filtered_feature_bc_matrix.h5")
  
  ## use cellbender results
  h5_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/", sample_id, "/", sample_id, "_cellbender_output_filtered.h5")
  
  # fragment_file = paste0("/scratch/g/mliang/snRNA_vs_multi/analysis/CellRanger/", sample_id, "/outs/atac_fragments.tsv.gz")
  # outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder/", sample_id, "_cb.pc40.DoubletFinder.rds")
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder/", sample_id, "_doubletfinder.rds")
  
  
  if(!file.exists(outfile)){
    counts = Read10X_h5(h5_file)
    # if(grepl('LK', sample_id)){
    #   seurat_object = CreateSeuratObject(counts = counts$`Gene Expression`,
    #                                      project = sample_id)
    # }else{
      seurat_object = CreateSeuratObject(counts = counts, # counts$`Gene Expression`, 
                                         project = sample_id,
                                         min.cells = 3,
                                         min.features = 200)
    # }
    
    seurat_object[["percent.mt"]] = PercentageFeatureSet(seurat_object, pattern = "^[Mm]t-")
    
    seurat_object = NormalizeData(seurat_object)
    seurat_object = FindVariableFeatures(seurat_object)
    seurat_object = ScaleData(seurat_object)
    seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
    seurat_object = RunUMAP(seurat_object, dims=1:40)
    
    stdv = seurat_object[["pca"]]@stdev
    sum_stdv = sum(seurat_object[["pca"]]@stdev)
    percent_stdv = (stdv / sum_stdv) * 100
    cumulative = cumsum(percent_stdv)
    co1 = which(cumulative > 90 & percent_stdv < 5)[1]
    co2 = sort(which((percent_stdv[1:length(percent_stdv) - 1] -
                        percent_stdv[2:length(percent_stdv)]) > 0.1),
               decreasing = T)[1] + 1
    pc = min(co1, co2)
    pc_list = c(pc, 20, 30, 40)
    # pc_list = c(40)
    
    ## try different pcs
    for(pci in 1:length(pc_list)){
      
      pc = pc_list[pci]
      ## pK identification
      sweep.list = paramSweep_v3(seurat_object, PCs = 1:pc, num.cores = detectCores() - 1)
      sweep.stats = summarizeSweep(sweep.list)
      bcmvn = find.pK(sweep.stats)
      
      pK=as.numeric(as.character(bcmvn$pK))
      BCmetric=bcmvn$BCmetric
      pK_choose = pK[which(BCmetric %in% max(BCmetric))]
      
      bcmvn.max = bcmvn[which.max(bcmvn$BCmetric),]
      optimal.pk = bcmvn.max$pK
      optimal.pk = as.numeric(levels(optimal.pk))[optimal.pk]
      
      doublet_rate = ncol(seurat_object)*8*1e-6
      if((doublet_rate > 0.3) & (doublet_rate > optimal.pk)){
        doublet_rate=optimal.pk
      }
      
      ## Homotypic doublet proportion estimate
      annotations = seurat_object@meta.data$seurat_clusters
      homotypic.prop = modelHomotypic(annotations) 
      nExp.poi = round(doublet_rate * nrow(seurat_object@meta.data)) ## Assuming 7.5% doublet formation rate - tailor for your dataset
      nExp.poi.adj = round(nExp.poi * (1 - homotypic.prop))
      
      # run DoubletFinder
      seurat_object = doubletFinder_v3(seu = seurat_object, 
                                       PCs = 1:pc, 
                                       pK = optimal.pk,
                                       nExp = nExp.poi.adj)
      
      metadata = seurat_object@meta.data
      if(pci==1){
        colnames(metadata)[ncol(metadata)] = "doublet_pc.optm"
      }else{
        colnames(metadata)[ncol(metadata)] = paste0("doublet_pc.", pci)
      }
      
      seurat_object@meta.data = metadata
      
      seurat_object@meta.data = seurat_object@meta.data[, -which(grepl("pANN", colnames(seurat_object@meta.data)))]
      
      
      para_list[nrow(para_list)+1, ] = c(sample_id, nrow(metadata),
                                         table(metadata[,ncol(metadata)])['Doublet'],
                                         doublet_rate,
                                         pc, optimal.pk)
    }
    
    saveRDS(seurat_object, outfile)
    
  }
  
}


write.table(para_list, "doubletFinder.parameter_list.v2.out", sep='\t', quote=F, col.names=F, row.names=F,
            append = T)
