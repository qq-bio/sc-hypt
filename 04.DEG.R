library(Seurat)
library(dplyr)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")


################################################################################
### treatment vs control
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

cluster = "subclass_level2"

for(i in input_file){

  deg_merged = c()

  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                   gsub("anno.rds", "DEG_all.out", basename(i)))

  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics|EC]+.anno.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))

  # seurat_object@active.ident = seurat_object@meta.data[, cluster]
  Idents(seurat_object) = cluster

  meta_table = seurat_object@meta.data

  project_list = unique(meta_table$project)
  for(pi in project_list){

    strain_list = unique(meta_table[meta_table$project==pi, ]$strain)
    for(si in strain_list){

      cell_list = unique(meta_table[meta_table$project==pi &
                                      meta_table$strain==si, cluster])

      treatment = meta_table[meta_table$project==pi &
                               meta_table$strain==si, ]$treatment
      treatment = intersect(levels(treatment), unique(treatment))
      control = treatment[1]

      for(cell in cell_list){

        cell.1 = rownames(meta_table[meta_table$project==pi &
                                       meta_table[,cluster]==cell &
                                       meta_table$treatment==control &
                                       meta_table$strain==si, ])

        if(length(cell.1)>3){

          for(ti in treatment[-1]){

            cell.2 = rownames(meta_table[meta_table$project==pi &
                                           meta_table[,cluster]==cell &
                                           meta_table$treatment==ti &
                                           meta_table$strain==si, ])

            if(length(cell.2)>3){

              deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)

              # deg$gene_orig = rownames(deg)
              # deg$gene = convert_to_h_gene(rownames(deg))

              deg$gene_name = rownames(deg)
              deg$cell_type = cell
              deg$pct.diff = deg$pct.1 - deg$pct.2
              deg$project = pi
              deg$strain = si
              deg$tissue = tissue
              deg$control = control
              deg$treatment = ti
              deg$control_size = length(cell.1)
              deg$treatment_size = length(cell.2)
              deg$test_gene = nrow(deg)

              # deg = deg[deg$p_val_adj<0.05, ]

              deg_merged = rbind(deg_merged, deg)


            }

          }

        }

      }

    }

  }

  write.table(deg_merged, outfile)

  print(i)

}




################################################################################
### strain_wise: ss vs sd & shr vs wky
input_file = c(
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
  
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

cluster = "subclass_level1"

for(i in input_file){
  
  deg_merged = c()
  
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                   gsub("anno.rds", "strain_wise.DEG_all.out", basename(i)))
  
  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics]+.anno.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  # seurat_object@active.ident = seurat_object@meta.data[, cluster]
  Idents(seurat_object) = cluster
  
  meta_table = seurat_object@meta.data
  
  project_list = unique(meta_table$project)
  
  for(pi in project_list){
    
    treatment = meta_table[meta_table$project==pi, ]$treatment
    treatment = intersect(levels(treatment), unique(treatment))
    control = treatment[1]
    
    cell_list = unique(meta_table[meta_table$project==pi &
                                    meta_table$treatment==control, cluster])
    
    for(cell in cell_list){
      
      strain = meta_table[meta_table$project==pi &
                                         meta_table$treatment==control &
                                         meta_table[,cluster]==cell, ]$strain
      strain_list = unique(strain)
      
      if(length(strain_list)==2){
        
        responsive_strain = intersect(strain_list, c("SS", "SHR"))
        control_strain = intersect(strain_list, c("SD", "WKY"))
        
        cell.1 = rownames(meta_table[meta_table$project==pi &
                                       meta_table[,cluster]==cell &
                                       meta_table$treatment==control &
                                       meta_table$strain==control_strain, ])
        
        cell.2 = rownames(meta_table[meta_table$project==pi &
                                       meta_table[,cluster]==cell &
                                       meta_table$treatment==control &
                                       meta_table$strain==responsive_strain, ])
        
        if(length(cell.1)>3 & length(cell.2)>3){
          
          deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)
          
          # deg$gene_orig = rownames(deg)
          # deg$gene = convert_to_h_gene(rownames(deg))
          
          deg$gene_name = rownames(deg)
          deg$cell_type = cell
          deg$pct.diff = deg$pct.1 - deg$pct.2
          deg$project = pi
          deg$strain = pi
          deg$tissue = tissue
          deg$control = paste0(control_strain, "-", control)
          deg$treatment = paste0(responsive_strain, "-", control)
          deg$control_size = length(cell.1)
          deg$treatment_size = length(cell.2)
          deg$test_gene = nrow(deg)
          
          # deg = deg[deg$p_val_adj<0.05, ]
          
          deg_merged = rbind(deg_merged, deg)
          
          
        }
        
      }
      
    }
    
  }
  
  write.table(deg_merged, outfile)
  
  print(i)
  
}




################################################################################
### merge DEG results
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*DEG_all*")

deg_merged = c()

for(i in input_file){

  dat_tmp = try(read.table(i, header = T), silent = T)

  if(class(dat_tmp) != "try-error"){

    deg_merged = rbind(deg_merged, dat_tmp)

  }

}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', col.names = T, row.names = F)

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=T)
deg_merged_reshape = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5,] %>%
  `rownames<-`( NULL ) %>%
  # group_by(strain, treatment, tissue, cell_type) %>%
  dplyr::mutate(group = paste(strain, treatment, tissue, cell_type, sep=".")) %>%
  # as.data.frame() %>%
  group_by(group) %>%
  dplyr::arrange(p_val_adj, .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, group, gene_name) %>%
  reshape(., idvar = "id", timevar = "group", v.names="gene_name", direction = "wide")
colnames(deg_merged_reshape) = sub("gene_name.", "", colnames(deg_merged_reshape))

write.table(deg_merged_reshape, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.wide.out", sep='\t', col.names = T, row.names = F)




### merge DEG results for strain-wise results
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*strain_wise.DEG_all*")

deg_merged = c()

for(i in input_file){
  
  dat_tmp = try(read.table(i, header = T), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep='\t', col.names = T, row.names = F)

deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.all.out", sep='\t', header=T)
deg_merged_reshape = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5,] %>%
  `rownames<-`( NULL ) %>%
  # group_by(strain, treatment, tissue, cell_type) %>%
  dplyr::mutate(group = paste(strain, treatment, tissue, cell_type, sep=".")) %>%
  # as.data.frame() %>%
  group_by(group) %>%
  dplyr::arrange(p_val_adj, .by_group=TRUE) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  dplyr::select(id, group, gene_name) %>%
  reshape(., idvar = "id", timevar = "group", v.names="gene_name", direction = "wide")
colnames(deg_merged_reshape) = sub("gene_name.", "", colnames(deg_merged_reshape))

write.table(deg_merged_reshape, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.DEG.merged.wide.out", sep='\t', col.names = T, row.names = F)








################################################################################
# 
# input_file = c(
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.immune_cell.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.immune_cell.anno.rds",
#   "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.immune_cell.anno.rds"
# )
# 
# cluster = "subclass_level2"
# 
# for(i in input_file){
# 
#   deg_merged = c()
# 
#   outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
#                    gsub("anno.rds", "DEG_all.out", basename(i)))
# 
#   seurat_object = readRDS(i)
#   # dataset = gsub(".[RNA|multiomics].+", "", basename(i), perl = T)
# 
#   seurat_object@active.ident = as.factor(seurat_object@meta.data[, cluster])
#   Idents(seurat_object) = cluster
# 
#   # markers = FindAllMarkers(seurat_object, group.by=subclass_level2, only.pos = TRUE, min.pct = 0.25)
#   # markers$pct.diff = markers$pct.1 - markers$pct.2
#   # markers = markers %>% group_by(cluster) %>%
#   #   dplyr::arrange(desc(pct.diff), .by_group=TRUE)
#   # write.table(markers,file=paste0(outfile, ".allmarker.0.25.long.txt"), sep="\t")
# 
# 
#   meta_table = seurat_object@meta.data
# 
#   project_list = unique(meta_table$project)
#   for(pi in project_list){
# 
#     species_list = unique(meta_table[meta_table$project==pi, ]$species)
#     for(si in species_list){
# 
#       tissue_list = unique(meta_table[meta_table$project==pi & meta_table$species==si, ]$tissue)
#       for(tissue in tissue_list){
# 
#         cell_list = unique(meta_table[meta_table$project==pi &
#                                         meta_table$species==si &
#                                         meta_table$tissue==tissue, cluster])
# 
#         treatment = meta_table[meta_table$project==pi &
#                                  meta_table$species==si &
#                                  meta_table$tissue==tissue, ]$treatment
#         treatment = intersect(levels(treatment), unique(treatment))
#         control = treatment[1]
# 
# 
#         for(cell in cell_list){
# 
#           cell.1 = rownames(meta_table[meta_table$project==pi &
#                                          meta_table[,cluster]==cell &
#                                          meta_table$treatment==control &
#                                          meta_table$species==si &
#                                          meta_table$tissue==tissue, ])
# 
#           if(length(cell.1)>3){
# 
#             for(ti in treatment[-1]){
# 
#               cell.2 = rownames(meta_table[meta_table$project==pi &
#                                              meta_table[,cluster]==cell &
#                                              meta_table$treatment==ti &
#                                              meta_table$species==si &
#                                              meta_table$tissue==tissue, ])
# 
#               if(length(cell.2)>3){
# 
#                 deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)
# 
#                 # deg$gene_orig = rownames(deg)
#                 # deg$gene = convert_to_h_gene(rownames(deg))
# 
#                 deg$gene_name = rownames(deg)
#                 deg$cell_type = cell
#                 deg$pct.diff = deg$pct.1 - deg$pct.2
#                 deg$project = pi
#                 deg$species = si
#                 deg$tissue = tissue
#                 deg$control = control
#                 deg$treatment = ti
#                 deg$control_size = length(cell.1)
#                 deg$treatment_size = length(cell.2)
#                 deg$test_gene = nrow(deg)
# 
#                 # deg = deg[deg$p_val_adj<0.05, ]
# 
#                 deg_merged = rbind(deg_merged, deg)
# 
# 
#               }
# 
#             }
# 
#           }
# 
#         }
# 
#       }
# 
#     }
# 
#   }
# 
#   write.table(deg_merged, outfile)
# 
#   print(i)
# 
# }
# 







################################################################################
# input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*RNA.DEG.out")
# 
# deg_merged = c()
# 
# for(i in input_file){
#   
#   dat_tmp = try(read.table(i, header = T), silent = T)
#   
#   if(class(dat_tmp) != "try-error"){
#     
#     if(!(grepl("immune_cell", i))){
#       dat_tmp = dat_tmp[!(dat_tmp$cell_type %in% c("Microglia", "Activated microglia", "IMM")), ]
#     }
#     
#     deg_merged = rbind(deg_merged, dat_tmp)
#     
#   }
#   
# }
# 
# write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", sep='\t', col.names = T, row.names = F)
# 
# 




# ################################################################################
# ### pseudo
input_file = c(
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

cluster = "subclass_level1"

for(i in input_file){

  deg_merged = c()

  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                   gsub("anno.rds", "pseudo.DEG_all.out", basename(i)))

  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics]+.anno.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))

  seurat_object@active.ident = seurat_object@meta.data[, cluster]
  Idents(seurat_object) = cluster

  seurat_object$seqID2 = gsub(" ", "", seurat_object$seqID2)
  
  if(dataset=="rat.ss.MSA"){
    seurat_object$seqID2_pseudo = paste0(seurat_object$seqID2, sample(1:2, size=ncol(seurat_object), replace = T))
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2_pseudo", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2_pseudo", "strain", "treatment", cluster)], sep="_"))
  }else if(dataset=="rat.ss.LV"){
    seurat_object$seqID2_pseudo = seurat_object$seqID2
    ss_ls_n = nrow(seurat_object@meta.data[seurat_object$orig.ident=="RLV5SN", ])
    seurat_object@meta.data[seurat_object$orig.ident=="RLV5SN", ]$seqID2_pseudo = paste0("RLV5SN", sample(1:2, size=ss_ls_n, replace = T))
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2_pseudo", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2_pseudo", "strain", "treatment", cluster)], sep="_"))
  }else{
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2", "strain", "treatment", cluster)], sep="_"))
  }
  pseudo_object$stc <- gsub("^\\w*?_", "", colnames(pseudo_object), perl = T)
  Idents(pseudo_object) <- "stc"

  pseudo_stc <- AverageExpression(seurat_object, assays = "RNA", return.seurat = F, slot = "data",
                                       group.by = c("strain", "treatment", cluster))

  project_list = unique(meta_table$project)
  for(pi in project_list){

    species_list = unique(meta_table[meta_table$project==pi, ]$strain)
    for(si in species_list){

      tissue_list = unique(meta_table[meta_table$project==pi & meta_table$strain==si, ]$tissue)
      for(tissue in tissue_list){

        cell_list = unique(meta_table[meta_table$project==pi &
                                        meta_table$strain==si &
                                        meta_table$tissue==tissue, cluster])

        treatment = meta_table[meta_table$project==pi &
                                 meta_table$strain==si &
                                 meta_table$tissue==tissue, ]$treatment
        treatment = intersect(levels(treatment), unique(treatment))
        control = treatment[1]

        for(cell in cell_list){

          cell.1 = rownames(meta_table[meta_table$project==pi &
                                         meta_table[,cluster]==cell &
                                         meta_table$treatment==control &
                                         meta_table$strain==si &
                                         meta_table$tissue==tissue, ])

          ident.1 = paste(si, control, cell, sep = "_")
          sample.1 = colnames(pseudo_object)[pseudo_object$stc %in% c(ident.1)]

          if(length(cell.1)>3 & length(sample.1)>1){

            for(ti in treatment[-1]){

              cell.2 = rownames(meta_table[meta_table$project==pi &
                                             meta_table[,cluster]==cell &
                                             meta_table$treatment==ti &
                                             meta_table$strain==si &
                                             meta_table$tissue==tissue, ])

              ident.2 = paste(si, ti, cell, sep = "_")
              sample.2 = colnames(pseudo_object)[pseudo_object$stc %in% c(ident.2)]

              if(length(cell.2)>3 & length(sample.2)>1){

                deg = FindMarkers(pseudo_object, ident.1=sample.1, ident.2=sample.2,
                                  logfc.threshold = 0, min.cells.group = 1, test.use = "DESeq2")

                deg$pct.1 = c(rowSums(seurat_object@assays$RNA@counts[,cell.1]>0)/length(cell.1))[rownames(deg)]
                deg$pct.2 = c(rowSums(seurat_object@assays$RNA@counts[,cell.2]>0)/length(cell.2))[rownames(deg)]
                deg$avg_expr.1 = pseudo_stc$RNA[rownames(deg), ident.1]
                deg$avg_expr.2 = pseudo_stc$RNA[rownames(deg), ident.2]

                deg$gene_name = rownames(deg)
                deg$cell_type = cell
                deg$pct.diff = deg$pct.1 - deg$pct.2
                deg$project = pi
                deg$strain = si
                deg$tissue = tissue
                deg$control = control
                deg$treatment = ti
                deg$control_size = length(cell.1)
                deg$treatment_size = length(cell.2)
                deg$test_gene = nrow(deg)

                # deg = deg[deg$p_val_adj<0.05, ]

                deg_merged = rbind(deg_merged, deg)


              }

            }

          }

        }

      }

    }

  }

  write.table(deg_merged, outfile)

  print(i)

}




input_file = c(
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

cluster = "subclass_level1"

for(i in input_file){
  
  deg_merged = c()
  
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                   gsub("anno.rds", "strain_wise.pseudo.DEG_all.out", basename(i)))
  
  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics]+.anno.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  seurat_object@active.ident = seurat_object@meta.data[, cluster]
  Idents(seurat_object) = cluster
  
  seurat_object$seqID2 = gsub(" ", "", seurat_object$seqID2)
  
  if(dataset=="rat.ss.MSA"){
    seurat_object$seqID2_pseudo = paste0(seurat_object$seqID2, sample(1:2, size=ncol(seurat_object), replace = T))
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2_pseudo", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2_pseudo", "strain", "treatment", cluster)], sep="_"))
  }else if(dataset=="rat.ss.LV"){
    seurat_object$seqID2_pseudo = seurat_object$seqID2
    ss_ls_n = nrow(seurat_object@meta.data[seurat_object$orig.ident=="RLV5SN", ])
    seurat_object@meta.data[seurat_object$orig.ident=="RLV5SN", ]$seqID2_pseudo = paste0("RLV5SN", sample(1:2, size=ss_ls_n, replace = T))
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2_pseudo", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2_pseudo", "strain", "treatment", cluster)], sep="_"))
  }else{
    meta_table = seurat_object@meta.data
    pseudo_object <- AggregateExpression(seurat_object, assays = "RNA", return.seurat = T, slot = "count",
                                         group.by = c("seqID2", "strain", "treatment", cluster))
    meta_table$sstc <- do.call(paste, c(meta_table[c("seqID2", "strain", "treatment", cluster)], sep="_"))
  }
  pseudo_object$stc <- gsub("^\\w*?_", "", colnames(pseudo_object), perl = T)
  Idents(pseudo_object) <- "stc"
  
  pseudo_stc <- AverageExpression(seurat_object, assays = "RNA", return.seurat = F, slot = "data",
                                  group.by = c("strain", "treatment", cluster))
  
  project_list = unique(meta_table$project)
  for(pi in project_list){
    
    treatment = meta_table[meta_table$project==pi, ]$treatment
    treatment = intersect(levels(treatment), unique(treatment))
    control = treatment[1]
    
    cell_list = unique(meta_table[meta_table$project==pi &
                                    meta_table$treatment==control, cluster])
    
    for(cell in cell_list){
      
      strain = meta_table[meta_table$project==pi &
                            meta_table$treatment==control &
                            meta_table[,cluster]==cell, ]$strain
      strain_list = unique(strain)
      # for(si in species_list){
    
      if(length(strain_list)==2){
        
      responsive_strain = intersect(strain_list, c("SS", "SHR"))
      control_strain = intersect(strain_list, c("SD", "WKY"))
    
      cell.1 = rownames(meta_table[meta_table$project==pi &
                                     meta_table[,cluster]==cell &
                                     meta_table$treatment==control &
                                     meta_table$strain==control_strain, ])
          
          ident.1 = paste(control_strain, control, cell, sep = "_")
          sample.1 = colnames(pseudo_object)[pseudo_object$stc %in% c(ident.1)]
          
          if(length(cell.1)>3 & length(sample.1)>1){
            
            # for(ti in treatment[-1]){
              
              cell.2 = rownames(meta_table[meta_table$project==pi &
                                             meta_table[,cluster]==cell &
                                             meta_table$treatment==control &
                                             meta_table$strain==responsive_strain, ])
              
              ident.2 = paste(responsive_strain, control, cell, sep = "_")
              sample.2 = colnames(pseudo_object)[pseudo_object$stc %in% c(ident.2)]
              
              if(length(cell.2)>3 & length(sample.2)>1){
                
                deg = FindMarkers(pseudo_object, ident.1=sample.1, ident.2=sample.2,
                                  logfc.threshold = 0, min.cells.group = 1, test.use = "DESeq2")
                
                deg$pct.1 = c(rowSums(seurat_object@assays$RNA@counts[,cell.1]>0)/length(cell.1))[rownames(deg)]
                deg$pct.2 = c(rowSums(seurat_object@assays$RNA@counts[,cell.2]>0)/length(cell.2))[rownames(deg)]
                deg$avg_expr.1 = pseudo_stc$RNA[rownames(deg), ident.1]
                deg$avg_expr.2 = pseudo_stc$RNA[rownames(deg), ident.2]
                
                deg$gene_name = rownames(deg)
                deg$cell_type = cell
                deg$pct.diff = deg$pct.1 - deg$pct.2
                deg$project = pi
                deg$strain = pi
                deg$tissue = tissue
                deg$control = paste0(control_strain, "-", control)
                deg$treatment = paste0(responsive_strain, "-", control)
                deg$control_size = length(cell.1)
                deg$treatment_size = length(cell.2)
                deg$test_gene = nrow(deg)
                
                # deg = deg[deg$p_val_adj<0.05, ]
                
                deg_merged = rbind(deg_merged, deg)
                
                
              }
              
            # }
            
          }
          
        }
        
      # }
      
    # }
    
  }
  
  write.table(deg_merged, outfile)
  
  print(i)
  
}
}




setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*[RNA|multiomics].pseudo.DEG_all*")

deg_merged = c()

for(i in input_file){
  
  dat_tmp = try(read.table(i, header = T), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/pseudo.DEG.all.out", sep='\t', col.names = T, row.names = F)




setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "*strain_wise.pseudo.DEG_all*")

deg_merged = c()

for(i in input_file){
  
  dat_tmp = try(read.table(i, header = T), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/strain_wise.pseudo.DEG.all.out", sep='\t', col.names = T, row.names = F)



# 
# deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.pseudo.merged_all.out", sep='\t', header=T)
# 
# 
# deg_merged_reshape = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.25,] %>% 
#   `rownames<-`( NULL ) %>% 
#   # group_by(species, treatment, tissue, cell_type) %>%
#   dplyr::mutate(group = paste(species, treatment, tissue, cell_type, sep=".")) %>% 
#   # as.data.frame() %>%
#   group_by(group) %>%
#   dplyr::arrange(p_val_adj, .by_group=TRUE) %>%
#   dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
#   dplyr::select(id, group, gene_name) %>% 
#   reshape(., idvar = "id", timevar = "group", v.names="gene_name", direction = "wide")
# colnames(deg_merged_reshape) = sub("gene_name.", "", colnames(deg_merged_reshape))
# 
# write.table(deg_merged_reshape, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.wide.out", sep='\t', col.names = T, row.names = F)
# 
# 
# 
# 
