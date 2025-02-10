#library(Seurat)
library(dplyr)
library(stringr)
library(clusterProfiler)
library(msigdbr)
library(biomaRt)
library(reshape2)
# library(org.Rn.eg.db)
# library(org.Mm.eg.db)


cell_order = c(
  c("Inhibitory neurons", "Excitatory neurons", "Avp+ neurons", "Neurons",
    "Astrocytes", # "Microglia", "Activated microglia", 
    "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycytes", "Ependymal cells", "Pars tuberalis cells"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "DCT/CT", "CT", "CD", "IC"),
  c("EC", "Venous EC", "VSMC", "VLMC", "Pericytes", "FIB", "ABC", "Adipocytes"),
  c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
    "NK cells", "NKT", "T cells", "B cells")
)


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


enrichment_analysis = function(deg_df,
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
    
    species_list = unique(deg_df[deg_df$project==pi,]$species)
    
    for(si in species_list){
      
      tissue_list = unique(deg_df[deg_df$project==pi & deg_df$species==si, ]$tissue)
      
      for(tissue in tissue_list){
        
        cell_list = unique(deg_df[deg_df$project==pi & 
                                    deg_df$species==si &
                                    deg_df$tissue==tissue, cluster])
        
        for(ci in cell_list){
          
          treatment_list = deg_df[deg_df$project==pi &
                                    deg_df$species==si &
                                    deg_df$tissue==tissue &
                                    deg_df[,cluster]==ci, ]$treatment
          
          for(ti in treatment_list){
            
            deg_use = deg_df[deg_df$project==pi &
                               deg_df$species==si &
                               deg_df$tissue==tissue &
                               deg_df[,cluster]==ci &
                               deg_df$treatment==ti, ]
            
            
            # deg_order = deg$avg_log2FC * -log(deg$p_val, 10)d
            # names(deg_order) = rownames(deg)
            # deg_order= deg_order[order(deg_order, decreasing = T)]
            # 
            # enrich_res = enricher(deg_use$gene, universe = convert_to_h_gene(rownames(seurat_object)),
            #                       TERM2GENE=h_t2g)@result
            
            if(nrow(deg_use)>=5){
              
              enrich_res = enricher(convert_to_h_gene(deg_use$gene_name), 
                                    TERM2GENE=h_t2g)@result
              
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
  
  return(enrich_res_merged)
  
}


categorize_within_cell_type <- function(data, model) {
  if(model=="ss"){
    ss_specific <- setdiff(
      unique(data$deg_name[data$sxt == "SS-HS 3d" | data$sxt == "SS-HS 21d"]),
      unique(data$deg_name[data$sxt == "SD-HS 3d"])
    )
    sd_specific <- setdiff(
      unique(data$deg_name[data$sxt == "SD-HS 3d"]),
      unique(data$deg_name[data$sxt == "SS-HS 3d" | data$sxt == "SS-HS 21d"])
    )
    overlapped <- intersect(
      unique(data$deg_name[data$sxt == "SS-HS 3d" | data$sxt == "SS-HS 21d"]),
      unique(data$deg_name[data$sxt == "SD-HS 3d"])
    )
    
    data$category <- case_when(
      data$deg_name %in% ss_specific ~ "SS-Specific",
      data$deg_name %in% sd_specific ~ "SD-Specific",
      data$deg_name %in% overlapped ~ "Overlapped",
      TRUE ~ NA_character_  # For genes not in either "SS-HS 3d" or "SD-HS 3d"
    )
  }

  if(model=="sp"){
    shr_specific <- setdiff(
      unique(data$deg_name[data$sxt == "SHR-26w"]),
      unique(data$deg_name[data$sxt == "WKY-26w"])
    )
    wky_specific <- setdiff(
      unique(data$deg_name[data$sxt == "WKY-26w"]),
      unique(data$deg_name[data$sxt == "SHR-26w"])
    )
    overlapped <- intersect(
      unique(data$deg_name[data$sxt == "SHR-26w"]),
      unique(data$deg_name[data$sxt == "WKY-26w"])
    )
    
    data$category <- case_when(
      data$deg_name %in% shr_specific ~ "SHR-Specific",
      data$deg_name %in% wky_specific ~ "WKY-Specific",
      data$deg_name %in% overlapped ~ "Overlapped",
      TRUE ~ NA_character_
    )
  }
  
  return(data)
  
}



deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header=T)
deg_df$deg_name = ifelse(deg_df$avg_log2FC > 0, paste(deg_df$gene_name, "↓", sep=""), paste(deg_df$gene_name, "↑", sep=""))
deg_df$sxt = paste0(deg_df$species, "-", deg_df$treatment)

deg_df_ss = deg_df[deg_df$project=="Salt-sensitive" & deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>0.25, ]
deg_df_ss = deg_df_ss %>%
  group_by(cell_type, tissue) %>%
  do(categorize_within_cell_type(., "ss")) %>% as.data.frame()
deg_df_count <- deg_df_ss %>%
  select(cell_type, tissue, category, deg_name) %>%
  unique %>%
  filter(!is.na(category)) %>%
  group_by(cell_type, tissue, category) %>%
  summarise(number_of_genes = n_distinct(deg_name)) %>% 
  as.data.frame()

species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(4, 5)])

deg_df_count$cell_type = factor(deg_df_count$cell_type, levels=cell_order)
deg_df_count$tissue = factor(deg_df_count$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
ggplot(mapping=aes(x = cell_type, y = number_of_genes)) +
  geom_col(data=deg_df_count[deg_df_count$category!="Overlapped",], aes(fill = category), position = "dodge") +
  geom_col(data=deg_df_count[deg_df_count$category=="Overlapped",], aes(fill = category), position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(species_col, "black"),
                    breaks = c("SS-Specific", "SD-Specific", "Overlapped"),
                    labels = c("SS-Specific", "SD-Specific", "Overlap")) +
  theme(axis.text.y = element_text(size=12, colour = 'black'),
        axis.text.x = element_text(size=10, colour = 'black'),
        legend.text = element_text(size=12, colour = 'black')) +
  labs(y="Number of DEGs", x="", fill="") +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")


deg_df_sp = deg_df[deg_df$project=="Spontaneous" & deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>1, ]
deg_df_sp = deg_df_sp %>%
  group_by(cell_type, tissue) %>%
  do(categorize_within_cell_type(., "sp")) %>% as.data.frame()
deg_df_count <- deg_df_sp %>%
  select(cell_type, tissue, category, deg_name) %>%
  unique %>%
  filter(!is.na(category)) %>%
  group_by(cell_type, tissue, category) %>%
  summarise(number_of_genes = n_distinct(deg_name)) %>% 
  as.data.frame()

species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(2, 3)])

deg_df_count$cell_type = factor(deg_df_count$cell_type, levels=cell_order)
deg_df_count$tissue = factor(deg_df_count$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
ggplot(mapping=aes(x = cell_type, y = number_of_genes)) +
  geom_col(data=deg_df_count[deg_df_count$category!="Overlapped",], aes(fill = category), position = "dodge") +
  geom_col(data=deg_df_count[deg_df_count$category=="Overlapped",], aes(fill = category), position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = c(species_col, "black"),
                    breaks = c("SHR-Specific", "WKY-Specific", "Overlapped"),
                    labels = c("SHR-Specific", "WKY-Specific", "Overlap")) +
  theme(axis.text.y = element_text(size=12, colour = 'black'),
        axis.text.x = element_text(size=10, colour = 'black'),
        legend.text = element_text(size=12, colour = 'black')) +
  labs(y="Number of DEGs", x="", fill="") +
  scale_x_discrete(limits=rev) +
  coord_flip() +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")




enrich_res_merged = enrichment_analysis(deg_df_ss[deg_df_ss$category=="SS-Specific", ])
enrich_res_merged = enrichment_analysis(deg_df_ss[deg_df_ss$category=="SHR-Specific", ])

# write.table(enrich_res_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.out")






