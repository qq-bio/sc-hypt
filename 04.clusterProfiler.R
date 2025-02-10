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


convert_to_h_gene <- function(gene_list, strain){
  
  if(strain=="rat"){
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
  
  m_t2g = msigdbr(species = "Mus musculus") %>%
    filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>%
    dplyr::select(gs_name, gene_symbol) %>% unique()

  r_t2g = msigdbr(species = "Rattus norvegicus") %>%
    filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>%
    dplyr::select(gs_name, gene_symbol) %>% unique()
  
  # h_t2g = msigdbr(species = "Homo sapiens") %>% 
  #   filter(gs_subcat %in% c("CP:KEGG", "CP:REACTOME", "GO:BP")) %>% 
  #   dplyr::select(gs_name, gene_symbol) %>% unique()
  
  project_list = unique(deg_df$project)
  
  for(pi in project_list){
    
    if(pi=="AngII"){
      t2g=m_t2g
      OrgDb="org.Mm.eg.db"
    }else{
      t2g=r_t2g
      OrgDb="org.Rn.eg.db"
    }
    
    strain = ifelse(pi=="AngII", "mouse", "rat")
    
    strain_list = unique(deg_df[deg_df$project==pi,]$strain)
    
    for(si in strain_list){
      
      tissue_list = unique(deg_df[deg_df$project==pi & deg_df$strain==si, ]$tissue)
      
      for(tissue in tissue_list){
        
        cell_list = unique(deg_df[deg_df$project==pi & 
                                    deg_df$strain==si &
                                    deg_df$tissue==tissue, cluster])
        
        for(ci in cell_list){
          
          treatment_list = unique(deg_df[deg_df$project==pi &
                                    deg_df$strain==si &
                                    deg_df$tissue==tissue &
                                    deg_df[,cluster]==ci, ]$treatment)
          
          for(ti in treatment_list){
            
            deg_use = deg_df[deg_df$project==pi &
                               deg_df$strain==si &
                               deg_df$tissue==tissue &
                               deg_df[,cluster]==ci &
                               deg_df$treatment==ti, ]
            
            all_gene_use = all_genes[all_genes$project==pi &
                                       all_genes$strain==si &
                                       all_genes$tissue==tissue &
                                       all_genes[,cluster]==ci &
                                       all_genes$treatment==ti, ]
            # deg_order = deg$avg_log2FC * -log(deg$p_val, 10)
            # names(deg_order) = rownames(deg)
            # deg_order= deg_order[order(deg_order, decreasing = T)]
            # 
            # enrich_res = enricher(deg_use$gene, universe = convert_to_h_gene(rownames(seurat_object)),
            #                       TERM2GENE=h_t2g)@result
            
            if(nrow(deg_use)>=5){
              
              # enrich_res = enricher(gene = deg_use$gene_name,
              #   #gene = convert_to_h_gene(deg_use$gene_name, strain),
              #                       TERM2GENE=t2g
              #                       # universe = all_gene_use$gene_name
              #                       # universe = convert_to_h_gene(all_gene_use$gene_name, strain)
              #   )
              
              enrich_res = enrichGO(deg_use$gene_name, OrgDb=OrgDb,
                                     keyType="SYMBOL", ont="BP")
              
              if(!is.null(enrich_res)){
                enrich_res = enrich_res@result
                
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
      
    }
    
  }
  
  return(enrich_res_merged)
  
}


deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header=T)
deg_df_use = deg_df[deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>0.5, ]

# all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header=T, sep='\t')

enrich_res_merged = enrichment_analysis(deg_df_use, deg_df)
write.table(enrich_res_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG.0.5.withot_background.out")













################################################################################
combine_pvalues <- function(row) {
  p_values <- na.omit(as.numeric(row))
  # chi_sq_statistic <- -2 * sum(log(p_values))
  chi_sq_statistic <- -2 * sum(p_values/log10(2))
  combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p_value)
}

merge_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG.0.5.withot_background.out", header = T, check.names = F)
merge_all$cell_type = factor(merge_all$cell_type, levels = cell_order)
merge_all$tissue = factor(merge_all$tissue, levels = tissue_order)
merge_all$path_cell = paste(merge_all$tissue, merge_all$cell_type,  merge_all$Description, sep="_")
merge_all$sxt = paste(merge_all$strain, merge_all$treatment, sep="-")

merge_background = merge_all[merge_all$strain %in% c("SD", "WKY"), ]
path_rm_list = merge_background[merge_background$pvalue<0.05, ]$path_cell
merge_filter = merge_all[! merge_all$path_cell %in% path_rm_list, ]


# merge_summary = merge_filter[grepl("Member", merge_filter$GroupID), ]
merge_summary = merge_filter; outfile = "clusterprofiler.merge_all.reshape.out"
top10 <- merge_summary[merge_summary$p.adjust<0.05,]
top20 <- merge_summary[merge_summary$p.adjust<0.05,] %>% group_by(strain, treatment, tissue, cell_type) %>% top_n(n = 20) %>% as.data.frame()
top_list <- top10

merge_filter2 = merge_filter[merge_filter$path_cell %in% top_list$path_cell, ]
merge_reshape = reshape(merge_filter2[,c("sxt", "path_cell", "p.adjust")], idvar=c("path_cell"), timevar = "sxt", v.names = "p.adjust", direction = "wide")
merge_reshape$`Log(comb.qval)` = log10(apply(merge_reshape[,2:6], 1, combine_pvalues))

merge_reshape[,grepl("Log", colnames(merge_reshape))] = -1 * merge_reshape[,grepl("Log", colnames(merge_reshape))]
merge_reshape[, c("tissue", "cell_type", "Description")] = do.call(rbind, strsplit(merge_reshape$path_cell, "_"))

merge_reshape = merge_reshape <- merge_reshape %>%
  arrange(desc(`Log(comb.qval)`)) %>%
  group_by(tissue, cell_type) %>%
  dplyr::mutate(id = row_number(),
                highlight="no") %>%
  ungroup() %>%
  arrange(tissue, cell_type) %>%  # Order by tissue and cell_type
  as.data.frame()

merge_reshape[merge_reshape=="NA"] <- NA
merge_reshape = na.replace(merge_reshape, 1)
merge_reshape$AngII = apply(merge_reshape[,2:3], 1, function(x) ifelse(any(as.numeric(x) < 0.05), 1, 0))
merge_reshape$SS = apply(merge_reshape[,6:7], 1, function(x) ifelse(any(as.numeric(x) < 0.05), 1, 0))
merge_reshape$SHR = ifelse(as.numeric(merge_reshape[,4]) < 0.05, 1, 0)

merge_reshape <- merge_reshape %>%
  mutate(category = case_when(
    AngII == 1 & SS == 0 & SHR == 0 ~ "Only AngII",
    AngII == 0 & SS == 1 & SHR == 0 ~ "Only SS",
    AngII == 0 & SS == 0 & SHR == 1 ~ "Only SHR",
    AngII == 1 & SS == 1 & SHR == 0 ~ "AngII & SS",
    AngII == 1 & SS == 0 & SHR == 1 ~ "AngII & SHR",
    AngII == 0 & SS == 1 & SHR == 1 ~ "SS & SHR",
    AngII == 1 & SS == 1 & SHR == 1 ~ "AngII & SS & SHR",
    TRUE ~ "Other"
  ))

merge_reshape$model_count = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
merge_reshape$cell_count = sapply(merge_reshape$Description, function(x) sum(merge_reshape$Description==x))
merge_reshape$pathway_sum = sapply(merge_reshape$Description, function(x) sum(merge_reshape$model_count[merge_reshape$Description==x]))

write.table(merge_reshape, outfile, sep='\t', quote=F, col.names=T, row.names=F)



################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_summary.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")

path_count = merge_reshape %>% filter(category != "Other") %>% 
  group_by(tissue, cell_type, category) %>%
  dplyr::summarise(count = n()) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)
path_count$category = factor(path_count$category, levels = c("Only AngII", "Only SS", "Only SHR", "AngII & SS",
                                                             "AngII & SHR", "SS & SHR", "AngII & SS & SHR"))
path_count$cell_type_mod = factor(paste0(path_count$cell_type, " •"), levels = paste0(cell_order, " ●"))
cell_col_mod = cell_col
names(cell_col_mod) = paste0(names(cell_col), " ●")
ggplot(path_count, aes(y = cell_type, x = category, label = count)) +
  geom_point(aes(size=count/3), alpha = .5, color="red") + geom_text() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  scale_y_discrete(limits=rev) +
  scale_size_continuous(breaks = c(1, 3, 6),
                        labels = c("3", "9", "18")) +
  labs(x="", y="", size="Number of\npathways") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free")


