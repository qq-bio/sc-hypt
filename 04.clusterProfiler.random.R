#library(Seurat)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(msigdbr)
library(biomaRt)
# library(org.Rn.eg.db)
# library(org.Mm.eg.db)


r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep='\t', header=T)
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep='\t', header=T)
r2h = r2h[r2h$Human.orthology.confidence..0.low..1.high.==1, ]
m2h = m2h[m2h$Human.orthology.confidence..0.low..1.high.==1, ]

# r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out", sep='\t', header=T)
# m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out", sep='\t', header=T)

r2h_mod = data.frame(id=c(r2h$Gene.stable.ID, r2h$Gene.name), name=c(r2h$Human.gene.name, r2h$Human.gene.name))
r2h_mod = r2h_mod %>% dplyr::filter(name != "") %>% distinct() %>% group_by(id) %>% top_n(1) %>% as.data.frame()
rownames(r2h_mod) = r2h_mod$id

m2h_mod = data.frame(id=c(m2h$Gene.stable.ID, m2h$Gene.name), name=c(m2h$Human.gene.name, m2h$Human.gene.name))
m2h_mod = m2h_mod %>% dplyr::filter(name != "") %>% distinct() %>% group_by(id) %>% top_n(1) %>% as.data.frame()
rownames(m2h_mod) = m2h_mod$id


convert_to_h_gene <- function(gene_list, species){
  
  if(species=="rat"){
    dat_mod = r2h_mod
  }else{
    dat_mod = m2h_mod
  }
  
  gene_list_new = lapply(gene_list, function(x) ifelse(x %in% dat_mod$id, dat_mod[x, ]$name, x)) %>% unlist()
  
  gene_list_new = unique(gene_list_new)
  
  return(gene_list_new)
  
}


enrichment_analysis = function(deg_df, all_genes,
                               cluster = "cell_type",
                               cell = NULL){
  
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
    species = ifelse(pi=="AngII", "mouse", "rat")
    
    species_list = unique(deg_df[deg_df$project==pi,]$species)
    
    for(si in species_list){
      
      tissue_list = unique(deg_df[deg_df$project==pi & deg_df$species==si, ]$tissue)
      
      for(tissue in tissue_list){
        
        cell_list = unique(deg_df[deg_df$project==pi & 
                                    deg_df$species==si &
                                    deg_df$tissue==tissue, cluster])
        
        for(ci in cell_list){
          
          treatment_list = unique(deg_df[deg_df$project==pi &
                                    deg_df$species==si &
                                    deg_df$tissue==tissue &
                                    deg_df[,cluster]==ci, ]$treatment)
          
          for(ti in treatment_list){
            
            ngene = nrow(deg_df[deg_df$project==pi &
                               deg_df$species==si &
                               deg_df$tissue==tissue &
                               deg_df[,cluster]==ci &
                               deg_df$treatment==ti, ])
            
            deg_use = all_genes[all_genes$project==pi &
                                  all_genes$species==si &
                                  all_genes$tissue==tissue &
                                  all_genes[,cluster]==ci &
                                  all_genes$treatment==ti, ]
            # deg_order = deg$avg_log2FC * -log(deg$p_val, 10)d
            # names(deg_order) = rownames(deg)
            # deg_order= deg_order[order(deg_order, decreasing = T)]
            # 
            # enrich_res = enricher(deg_use$gene, universe = convert_to_h_gene(rownames(seurat_object)),
            #                       TERM2GENE=h_t2g)@result
            
            if(ngene>=5){
              
              enrich_res = enricher(gene = convert_to_h_gene(sample(deg_use$gene_name, ngene), species), 
                                    TERM2GENE=h_t2g,
                                    universe = convert_to_h_gene(deg_use$gene_name, species))
              if(!is.null(enrich_res)){
                enrich_res = enrich_res@result
                
                enrich_res$project = pi
                enrich_res$species = si
                enrich_res$tissue = tissue
                enrich_res$cell_type = ci
                enrich_res$treatment = ti
                
                enrich_res_merged = rbind(enrich_res_merged, enrich_res)
              }
                
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  return(enrich_res_merged)
  
}


deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header=T)
deg_df = deg_df[deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>0.25, ]

all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header=T, sep='\t')

enrich_res_all = c()
for(i in 1:10){
  
  enrich_res_merged = enrichment_analysis(deg_df, all_genes)
  enrich_res_all = rbind(enrich_res_all, enrich_res_merged)
  
}

write.table(enrich_res_all, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.ortho.random.out")

