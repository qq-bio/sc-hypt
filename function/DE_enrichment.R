library(Seurat)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(msigdbr)
library(biomaRt)
# library(org.Rn.eg.db)
# library(org.Mm.eg.db)


r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out", sep='\t', header=T)
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out", sep='\t', header=T)

r2h_mod = data.frame(id=c(r2h$Gene.stable.ID, r2h$Gene.name), name=c(r2h$Human.gene.name, r2h$Human.gene.name))
r2h_mod = r2h_mod %>% dplyr::filter(name != "") %>% distinct() %>% group_by(id) %>% top_n(1) %>% as.data.frame()
rownames(r2h_mod) = r2h_mod$id

m2h_mod = data.frame(id=c(m2h$Gene.stable.ID, m2h$Gene.name), name=c(m2h$Human.gene.name, m2h$Human.gene.name))
m2h_mod = m2h_mod %>% dplyr::filter(name != "") %>% distinct() %>% group_by(id) %>% top_n(1) %>% as.data.frame()
rownames(m2h_mod) = m2h_mod$id


convert_to_h_gene <- function(gene_list){

  if(dim(table(grepl("ENSRNOG", gene_list)))>1){
    dat_mod = r2h_mod
  }else{
    dat_mod = m2h_mod
  }
  
  gene_list_new = lapply(gene_list, function(x) ifelse(x %in% dat_mod$id, dat_mod[x, ]$name, x)) %>% unlist()
  
  return(gene_list_new)

}

# enrichment_analysis = function(seurat_object,
#                        cluster = "new.cluster.ids_umap",
#                        cell = NULL,
#                        outfile){
#   
#   enrich_res_merged = data.frame()
#   deg_merged = data.frame()
#   
#   # m_t2g = msigdbr(species = "Mus musculus") %>% 
#   #   filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
#   #   dplyr::select(gs_name, gene_symbol) %>% unique()
#   # 
#   # r_t2g = msigdbr(species = "Rattus norvegicus") %>% 
#   #   filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
#   #   dplyr::select(gs_name, gene_symbol) %>% unique()
#   
#   h_t2g = msigdbr(species = "Homo sapiens") %>% 
#     filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
#     dplyr::select(gs_name, gene_symbol) %>% unique()
#   
#   meta_table = seurat_object@meta.data
#   
#   for(pi in unique(meta_table$project)){
#     
#     # if(pi=="AngII"){
#     #   t2g=m_t2g
#     #   name_id = m_name_id
#     # }else{
#     #   t2g=r_t2g
#     #   name_id = r_name_id
#     # }
#     
#     for(si in unique(meta_table[meta_table$project==pi, ]$species)){
#       
#       treatment = meta_table[meta_table$project==pi & 
#                                meta_table$species==si, ]$treatment
#       treatment = intersect(levels(treatment), unique(treatment))
#       
#       for(ti in treatment[-1]){
#       
#         control_sample = treatment[1]; treatment_sample = ti
#         cell.1 = rownames(meta_table[meta_table[,cluster]==cell &
#                                        meta_table$treatment==control_sample &
#                                        meta_table$species==si, ])
#         cell.2 = rownames(meta_table[meta_table[,cluster]==cell &
#                                        meta_table$treatment==treatment_sample &
#                                        meta_table$species==si, ])
#         deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0.25, min.pct=0.1)
#         deg$gene_orig = rownames(deg)
#         deg$gene = convert_to_h_gene(rownames(deg))
#         
#         deg_use = deg[deg$p_val_adj<0.05, ]
#         
#         deg_order = deg$avg_log2FC * -log(deg$p_val, 10)
#         names(deg_order) = rownames(deg)
#         deg_order= deg_order[order(deg_order, decreasing = T)]
#         
#         enrich_res = enricher(deg_use$gene, universe = convert_to_h_gene(rownames(seurat_object)),
#                               TERM2GENE=h_t2g)@result
#         enrich_res$cell = cell; deg_use$cell = cell
#         enrich_res$control = treatment[1]; deg_use$control = treatment[1]
#         enrich_res$treatment = ti; deg_use$treatment = ti
#         enrich_res$species = si; deg_use$species = si
#         enrich_res$project = pi; deg_use$project = pi
#         
#         
#         enrich_res_merged = rbind(enrich_res_merged, enrich_res)
#         deg_merged = rbind(deg_merged, deg_use)
#         
#       }
# 
#     }
#     
#   }
#   
#   write.table(enrich_res_merged, paste0(outfile, ".DE_enrichment.out"))
#   write.table(deg_merged, paste0(outfile, ".DEG.out"))
#   
# }
# 





enrichment_analysis = function(deg_df,
                               cluster = "cell_type",
                               cell = NULL,
                               outfile){
  
  enrich_res_merged = data.frame()
  
  # m_t2g = msigdbr(species = "Mus musculus") %>% 
  #   filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
  #   dplyr::select(gs_name, gene_symbol) %>% unique()
  # 
  # r_t2g = msigdbr(species = "Rattus norvegicus") %>% 
  #   filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
  #   dplyr::select(gs_name, gene_symbol) %>% unique()
  
  h_t2g = msigdbr(species = "Homo sapiens") %>% 
    filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
    dplyr::select(gs_name, gene_symbol) %>% unique()
  
  
  project_list = unique(deg_df$project)
  
  for(pi in project_list){
    
    # if(pi=="AngII"){
    #   t2g=m_t2g
    #   name_id = m_name_id
    # }else{
    #   t2g=r_t2g
    #   name_id = r_name_id
    # }
    
    species_list = unique(deg_df[deg_df$project==pi,]$strain)
    
    for(si in species_list){
      
      tissue_list = unique(deg_df[deg_df$project==pi & deg_df$strain==si, ]$tissue)
      
      for(tissue in tissue_list){
        
        cell_list = unique(deg_df[deg_df$project==pi & 
                                  deg_df$strain==si &
                                  deg_df$tissue==tissue, cluster])
        
        for(ci in cell_list){
          
          treatment_list = deg_df[deg_df$project==pi &
                                  deg_df$strain==si &
                                  deg_df$tissue==tissue &
                                  deg_df[,cluster]==ci, ]$treatment
          
          for(ti in treatment_list){
            
            
            
            
            deg_use = deg_df[deg_df$project==pi &
                               deg_df$strain==si &
                               deg_df$tissue==tissue &
                               deg_df[,cluster]==ci &
                               deg_df$treatment==ti, ]
            
            
            # deg_order = deg$avg_log2FC * -log(deg$p_val, 10)d
            # names(deg_order) = rownames(deg)
            # deg_order= deg_order[order(deg_order, decreasing = T)]
            # 
            # enrich_res = enricher(deg_use$gene, universe = convert_to_h_gene(rownames(seurat_object)),
            #                       TERM2GENE=h_t2g)@result
            
            enrich_res = enricher(convert_to_h_gene(deg_use$gene_name), 
                                  TERM2GENE=h_t2g)@result
            
            
            enrich_res$project = pi
            enrich_res$strain = si
            enrich_res$tissue = tissue
            enrich_res$cell_type = ci
            enrich_res$treatment = ti
            
            
            enrich_res_merged = rbind(enrich_res_merged, enrich_res)
            
            
            }
          
          }
        
        }
      
      }
      
    }
    
  # write.table(enrich_res_merged, paste0(outfile, ".DE_enrichment.out"))
  
  }
  








gsea_analysis = function(seurat_object,
                               cluster = "new.cluster.ids_umap",
                               cell = NULL,
                               treatment_group = NULL,
                               outfile){
  
  enrich_res_merged = data.frame()
  deg_merged = data.frame()
  
  m_t2g = msigdbr(species = "Mus musculus") %>% 
    filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
    dplyr::select(gs_name, gene_symbol) %>% unique()
  
  r_t2g = msigdbr(species = "Rattus norvegicus") %>% 
    filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
    dplyr::select(gs_name, gene_symbol) %>% unique()

  meta_table = seurat_object@meta.data
  
  for(pi in unique(meta_table$project)){
    
    if(pi=="C57BL/6"){t2g=m_tg2}else{t2g=r_t2g}
    
    for(si in unique(meta_table[meta_table$project==pi, ]$species)){
      
      treatment = meta_table[meta_table$project==pi & 
                               meta_table$species==si, ]$treatment
      treatment = intersect(levels(treatment), unique(treatment))
      
      for(ti in treatment[-1]){
        
        control_sample = treatment[1]; treatment_sample = ti
        cell.1 = rownames(meta_table[meta_table[,cluster]==cell &
                                       meta_table$treatment==control_sample &
                                       meta_table$species==si, ])
        cell.2 = rownames(meta_table[meta_table[,cluster]==cell &
                                       meta_table$treatment==treatment_sample &
                                       meta_table$species==si, ])
        deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0, min.pct=0.1)
        
        deg_use = rownames(deg[deg$p_val_adj<0.05, ])
        
        deg_order = deg$avg_log2FC * -log(deg$p_val, 10)
        names(deg_order) = rownames(deg)
        deg_order= deg_order[order(deg_order, decreasing = T)]
        
        enrich_res = GSEA(deg_order, TERM2GENE=t2g)@result
        enrich_res$cell = cell; deg$cell = cell
        enrich_res$control = treatment[1]; deg$control = treatment[1]
        enrich_res$treatment = ti; deg$treatment = ti
        enrich_res$species = si; deg$species = si
        enrich_res$project = pi; deg$project = pi
        
        enrich_res_merged = rbind(enrich_res_merged, enrich_res)
        deg_merged = rbind(deg_merged, deg)
      }
      
    }
    
  }
  
  write.table(enrich_res_merged, paste0(outfile, "DE_enrichment.out"))
  write.table(deg_merged, paste0(outfile, "DEG.out"))
  
}
