library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggupset)
# library(VennDiagram)
# library(ggvenn)
# library(ggtern)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))


# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DEG_enrichment.R")
mouse2rat = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2rat.out.txt", header = T, sep = "\t")
mouse2rat = unique(mouse2rat[mouse2rat[,2]==1, c("Rat.gene.name", "Gene.name")])

calculateJaccardIndex <- function(set1, set2) {
  intersectionSize <- length(intersect(set1, set2))
  unionSize <- length(union(set1, set2))
  jaccardIndex <- intersectionSize / unionSize
  return(jaccardIndex)
}

calculateDiceCoefficient <- function(set1, set2) {
  # Convert input to sets if they are not already (to ensure uniqueness)
  set1 <- unique(set1)
  set2 <- unique(set2)
  
  # Calculate intersection and sizes
  intersectionSize <- length(intersect(set1, set2))
  setSizeSum <- length(set1) + length(set2)
  
  # Calculate Dice coefficient
  diceCoefficient <- (2 * intersectionSize) / setSizeSum
  
  return(diceCoefficient)
}

findUniqueDEGs <- function(targetList, allLists) {
  # Remove the target list from the comparison pool
  comparisonLists <- allLists[names(allLists) != targetList]
  # Combine all DEGs in the comparison pool
  comparisonDegs <- unlist(comparisonLists, use.names = FALSE)
  # Identify DEGs unique to the target list
  uniqueDegs <- setdiff(allLists[[targetList]], comparisonDegs)
  return(uniqueDegs)
}

transformGeneNames <- function(df, gene_name) {
  transformedNames <- ifelse(df$avg_log2FC > 0, paste(df[,gene_name], "↓", sep=""), paste(df[,gene_name], "↑", sep=""))
  transformedNames = transformedNames[!(is.na(transformedNames))]
  return(transformedNames)
}


performPermutationTest <- function(deg_c1, deg_c2, all_gene_c1, all_gene_c2, n_perm = 10000) {
  observed_similarity = calculateJaccardIndex(deg_c1, deg_c2)
  perm_similarities = numeric(n_perm)
  
  for (i in 1:n_perm) {
    perm_c1 = sample(all_gene_c1, length(deg_c1))
    perm_c2 = sample(all_gene_c2, length(deg_c2))
    perm_similarities[i] = calculateJaccardIndex(perm_c1, perm_c2)
  }
  
  p_value = mean(perm_similarities >= observed_similarity)
  return(p_value)
}


process_df_with_orthologs <- function(df, orthologs_df=mouse2rat) {
  # Ensure orthologs_df has the correct column names
  orthologs_df <- orthologs_df %>%
    dplyr::rename(gene_name = Rat.gene.name, gene_name_ortho2m = Gene.name)
  
  # orthologs_df = unique(orthologs_df[c("gene_name", "gene_name_ortho2m")])
  
  # Add 'gene_name_ortho2m' column to df for mouse records directly and as NA for rat records
  df <- df %>%
    mutate(gene_name_ortho2m = case_when(
      project == "AngII" ~ gene_name, # Copy gene_name for mouse records
      TRUE ~ NA_character_ # Set as NA for rat records, to be filled with orthologous names
    ))
  
  # Merge df with orthologs_df to fill in orthologous mouse gene names for rat records
  df <- merge(df, orthologs_df, by = "gene_name", all.x = TRUE) %>%
    mutate(gene_name_ortho2m = if_else(is.na(gene_name_ortho2m.x), gene_name_ortho2m.y, gene_name_ortho2m.x)) %>%
    dplyr::select(-gene_name_ortho2m.x, -gene_name_ortho2m.y) # Remove temporary columns
  
  # Filter to keep only rows where 'gene_name_ortho2m' exists in the orthologous table
  # df_filtered <- df %>%
  #   filter(gene_name_ortho2m %in% orthologs_df$gene_name_ortho2m)
  
  return(df)
}




all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T, sep='\t')
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

all_genes$ji_comp = paste(all_genes$species, all_genes$treatment,  sep="-")

mouse2rat=mouse2rat[mouse2rat$Gene.name %in% all_genes$gene_name,]
all_genes = process_df_with_orthologs(all_genes)

# deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header = T, row.names = 1)
# deg_df = process_df_with_orthologs(deg_df)

deg_df = all_genes[all_genes$p_val_adj < 0.05 & abs(all_genes$avg_log2FC) > 0.5, ]
# deg_df$ji_comp = paste(deg_df$strain, deg_df$treatment,  sep="-")

# comp_list1 = list(c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d"),
#                   c("SS-HS 3d", "SS-HS 21d"),
#                  c("SS-HS 3d", "SD-HS 3d"),
#                  c("SHR-26w", "WKY-26w"))

comp_list2 = list(c("AngII", "Salt-sensitive"),
                  c("Salt-sensitive", "Spontaneous"),
                  c("AngII", "Spontaneous"))

tissue_list = unique(deg_df$tissue)

ji_df = c()

for(ti in tissue_list){
  
  cell_list = unique(deg_df[deg_df$tissue==ti, ]$cell_type)
  
  for(ci in cell_list){
    
    deg_use = deg_df[deg_df$tissue==ti & deg_df$cell_type==ci, ]
    all_gene_use = all_genes[all_genes$tissue==ti & all_genes$cell_type==ci, ]
    
    # for(j in comp_list1){
    #   
    #   comp1 = j[1]
    #   comp2 = j[2]
    #   similarity = ""
    #   p_value = ""
    #   
    #   # deg_c1 = unique(deg_use[deg_use$ji_comp==comp1, ]$gene_name)
    #   # deg_c2 = unique(deg_use[deg_use$ji_comp==comp2, ]$gene_name)
    #   
    #   if(grepl("Ang", comp1) & !(grepl("Ang", comp2))){
    #     gene_name = "gene_name_ortho2m"
    #   }else{
    #     gene_name = "gene_name"
    #   }
    #   
    #   deg_c1 = unique(transformGeneNames(deg_use[deg_use$ji_comp==comp1, ], gene_name))
    #   deg_c2 = unique(transformGeneNames(deg_use[deg_use$ji_comp==comp2, ], gene_name))
    #   
    #   
    #   if(length(deg_c1)>0 & length(deg_c2)>0){
    #     similarity = calculateJaccardIndex(deg_c1, deg_c2)
    #     # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
    #     p_value = performPermutationTest(deg_c1, deg_c2, 
    #                                      transformGeneNames(all_gene_use[all_gene_use$ji_comp==comp1, ], gene_name),
    #                                      transformGeneNames(all_gene_use[all_gene_use$ji_comp==comp2, ], gene_name))
    #   }
    #   
    #   ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
    #   
    # }
    
    
    deg_use = deg_df[deg_df$tissue==ti & deg_df$cell_type==ci & 
                       deg_df$strain %in% c("C57BL/6", "SS", "SHR"), ]
    
    for(j in comp_list2){
      
      comp1 = j[1]
      comp2 = j[2]
      similarity = ""
      p_value = ""
      
      if(grepl("Ang", comp1) & !(grepl("Ang", comp2))){
        gene_name = "gene_name_ortho2m"
      }else{
        gene_name = "gene_name"
      }
      
      deg_c1 = unique(transformGeneNames(deg_use[deg_use$project==comp1, ], gene_name))
      deg_c2 = unique(transformGeneNames(deg_use[deg_use$project==comp2, ], gene_name))
      
      if(length(deg_c1)>0 & length(deg_c2)>0){
        similarity = calculateJaccardIndex(deg_c1, deg_c2)
        # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
        p_value = performPermutationTest(deg_c1, deg_c2, 
                                         transformGeneNames(all_gene_use[all_gene_use$project==comp1, ], gene_name),
                                         transformGeneNames(all_gene_use[all_gene_use$project==comp2, ], gene_name))
      }
      
      ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
      
    }
    
  }
  
}

write.table(ji_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.ji.out", col.names = T, row.names = F, sep = '\t', quote = F)




ji_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.ji.out", header=T, sep = '\t')
ji_df = as.data.frame(ji_df)
colnames(ji_df) = c("tissue", "cell_type", "comp1", "comp2", "ji_score", "pval")
ji_df = ji_df[ji_df$tissue!="MCA", ]
ji_df$cell_type = factor(ji_df$cell_type, levels=unique(cell_order))
ji_df$tissue = factor(ji_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
ji_df$ji_score = as.numeric(ji_df$ji_score)
ji_df$pval = as.numeric(ji_df$pval)
ji_df$padj = p.adjust(ji_df$pval, method = "BH")
ji_df$annotation <- ifelse(ji_df$padj < 0.05, "*", "")
ji_df[is.na(ji_df$annotation), ]$annotation=""
ji_df$comp = paste(ji_df[,3], ji_df[,4], sep="-")
# ji_df$comp = factor(ji_df$comp, c("C57BL/6-AngII 3d-C57BL/6-AngII 28d", "SS-HS 3d-SS-HS 21d", 
#                                   "SS-HS 3d-SD-HS 3d", "SHR-26w-WKY-26w", "AngII-Salt-sensitive",
#                                   "AngII-Spontaneous", "Salt-sensitive-Spontaneous"))

ji_df_use = ji_df[ji_df$comp %in% c("AngII-Salt-sensitive", "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),]

ggplot(ji_df_use) +
  geom_tile(mapping = aes(x = comp, y = cell_type, fill=ji_score)) +
  geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low="white", high="blue") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Jaccard index") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.ji.png", width=361/96, height=718/96, dpi=300)

ggplot(ji_df_use, aes(x = ji_score, y = cell_type)) +
  geom_vline(xintercept = c(0, 0.1, 0.2, 0.3), color="grey", alpha=0.3)+
  geom_boxplot(color="darkgrey", fill="white") +
  geom_point(color="blue", alpha=0.8) +
  scale_y_discrete(limits=rev) +
    # scale_color_gradient(low="white", high="blue") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="Dice coefficient", y="") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")










################################################################################
### correlation plot

all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T, sep='\t')
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

all_genes$comp = paste(all_genes$strain, all_genes$treatment,  sep="-")

mouse2rat=mouse2rat[mouse2rat$Gene.name %in% all_genes$gene_name,]
all_genes = process_df_with_orthologs(all_genes)

comp_list = list(c("C57BL/6-AngII 3d", "SS-HS 3d"),
                 c("C57BL/6-AngII 28d", "SS-HS 21d"),
                 c("C57BL/6-AngII 28d", "SHR-26w"),
                 c("SS-HS 21d", "SHR-26w"))

cor_df = c()
for( i in comp_list){
  
  comp1 = i[1]
  comp2 = i[2]
  
  tissue_list = unique(all_genes[all_genes$comp %in% c(comp1, comp2), ]$tissue)
  
  for( ti in tissue_list ){

    cell_list = unique(all_genes[all_genes$comp %in% c(comp1, comp2) & all_genes$tissue==ti, ]$cell_type)
    
    for( ci in cell_list ){
      
      dat1 = all_genes[all_genes$comp==comp1 & all_genes$tissue==ti & all_genes$cell_type==ci, ]
      dat2 = all_genes[all_genes$comp==comp2 & all_genes$tissue==ti & all_genes$cell_type==ci, ]
      
      if(nrow(dat1)>0 & nrow(dat2)>0){
        
        dat_tmp = merge(dat1, dat2, by = "gene_name_ortho2m")
        cor_res = cor.test(dat_tmp$avg_log2FC.x, dat_tmp$avg_log2FC.y)
        
        cor_df = rbind(cor_df, c(ti, ci, comp1, comp2, cor_res$estimate, cor_res$p.value))
        
      }
      
    }
    
  }
  
  
}

write.table(cor_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.cor.out", col.names = T, row.names = F, sep = '\t', quote = F)




cor_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.cor.out", header=T, sep = '\t')
cor_df = as.data.frame(cor_df)
colnames(cor_df) = c("tissue", "cell_type", "comp1", "comp2", "ji_score", "pval")
cor_df$cell_type = factor(cor_df$cell_type, levels=unique(cell_order))
cor_df$tissue = factor(cor_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
cor_df$ji_score = as.numeric(cor_df$ji_score)
cor_df$pval = as.numeric(cor_df$pval)
cor_df$padj = p.adjust(cor_df$pval, method = "BH")
cor_df$annotation <- ifelse(cor_df$padj < 0.05, "*", "")
cor_df$comp = paste(cor_df[,3], cor_df[,4], sep="-")
# cor_df$comp = factor(cor_df$comp, c("C57BL/6-AngII 3d-C57BL/6-AngII 28d", "SS-HS 3d-SS-HS 21d", 
#                                   "SS-HS 3d-SD-HS 3d", "SHR-26w-WKY-26w", "AngII-Salt-sensitive",
#                                   "AngII-Spontaneous", "Salt-sensitive-Spontaneous"))

# cor_df_use = cor_df[cor_df$comp %in% c("AngII-Salt-sensitive", "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),]

ggplot(cor_df) +
  geom_tile(mapping = aes(x = comp, y = cell_type, fill=ji_score)) +
  geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", mid = "white", high="red", midpoint = 0) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Correlation\ncoefficient") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")










################################################################################
### process yong's list
yong_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.txt", header = T, sep='\t')
yong_list = str_to_sentence(yong_list$Gene.product.name)
yong_list = gsub("_human", "", yong_list)
yong_list = setdiff(yong_list, "")

m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
m2h = unique(m2h[, c("Gene.name", "Human.gene.name")])

r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")
r2h = unique(r2h[, c("Gene.name", "Human.gene.name")])

yong_mod1 = rbind(m2h[m2h$Gene.name %in% yong_list, ],
                 r2h[r2h$Gene.name %in% yong_list, ]) %>% unique()
human_list_p = str_to_upper(yong_list[!(yong_list %in% yong_mod1$Gene.name)])
yong_mod2 = rbind(m2h[m2h$Human.gene.name %in% human_list_p, ],
      r2h[r2h$Human.gene.name %in% human_list_p, ]) %>% unique()
yong_final = rbind(yong_mod1, yong_mod2) %>% unique()

write.table(yong_final, "/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.mod.txt", col.names = T, sep='\t')

################################################################################
### enrichment analysis
library(stringr)

m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
m2h = unique(m2h[m2h[,4]==1, c("Gene.name", "Human.gene.name")])

r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")
r2h = unique(r2h[r2h[,5]==1, c("Gene.name", "Human.gene.name")])

all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T, sep='\t')
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]
all_genes = all_genes[!(all_genes$tissue=="MCA" & all_genes$project %in% c("AngII", "Salt-sensitive")),]

all_genes$comp = paste(all_genes$strain, all_genes$treatment,  sep="-")

# mouse2rat=mouse2rat[mouse2rat$Gene.name %in% all_genes$gene_name,]
# all_genes = process_df_with_orthologs(all_genes)

yong_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.mod.txt", header = T, sep='\t')
yong_list = unique(yong_list$Gene.name)
helen_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/NG_2024_Helen_bp_pred_gene_list.txt", header = T, sep='\t')
helen_list = unique(helen_list$Gene)

bp_rat = unique(r2h[r2h$Human.gene.name %in% helen_list, ]$Gene.name)
bp_mouse = unique(m2h[m2h$Human.gene.name %in% helen_list, ]$Gene.name)
helen_list = unique(c(bp_rat, bp_mouse))

bp_list = helen_list
enrich_df = c()
comp_list = unique(all_genes$comp)
for( i in comp_list){
  
  tissue_list = unique(all_genes[all_genes$comp==i, ]$tissue)
  
  for( ti in tissue_list ){
    
    cell_list = unique(all_genes[all_genes$comp==i & all_genes$tissue==ti, ]$cell_type)
    
    for( ci in cell_list ){
      
      dat_tmp = all_genes[all_genes$comp==i & all_genes$tissue==ti & all_genes$cell_type==ci, ]
      deg_list = unique(dat_tmp[dat_tmp$p_val_adj<0.05 & abs(dat_tmp$avg_log2FC) > 0.5, ]$gene_name)
      bg_list = unique(dat_tmp$gene_name)
      
      # if(unique(dat_tmp$strain)=="C57BL/6"){
      #   bp_list = bp_mouse
      # }else{
      #   bp_list = bp_rat
      # }
      
      x <- sum(deg_list %in% bp_list)  # Number of DEGs in the reference list
      m <- sum(bp_list %in% bg_list)  # Total number of genes in the reference list
      n <- length(bg_list) - m  # Total number of genes not in the reference list
      k <- length(deg_list)  # Total number of DEGs
      
      print(c(x/k, m/(n+m)))
      
      p_value <- phyper(x - 1, m, n, k, lower.tail = FALSE)
      enrich_df = rbind(enrich_df, c(ti, ci, i, x/k, m/(n+m), x/length(bp_list), p_value))
      
    }
    
  }
  
}

write.table(enrich_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.enrich.yong.out", col.names = T, row.names = F, sep = '\t', quote = F)
# write.table(enrich_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.enrich.helen.out", col.names = T, row.names = F, sep = '\t', quote = F)




yong_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.mod.txt", header = T, sep='\t')
yong_list = unique(yong_list$Gene.name)
helen_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/NG_2024_Helen_bp_pred_gene_list.txt", header = T, sep='\t')
helen_list = unique(helen_list$Gene)

bp_rat = unique(r2h[r2h$Human.gene.name %in% helen_list, ]$Gene.name)
bp_mouse = unique(m2h[m2h$Human.gene.name %in% helen_list, ]$Gene.name)

bp_list = yong_list

enrich_df = c()
comp_list = unique(all_genes$strain)
for( i in comp_list){
  
  tissue_list = unique(all_genes[all_genes$strain==i, ]$tissue)
  
  for( ti in tissue_list ){
    
    cell_list = unique(all_genes[all_genes$strain==i & all_genes$tissue==ti, ]$cell_type)
    
    for( ci in cell_list ){
      
      dat_tmp = all_genes[all_genes$strain==i & all_genes$tissue==ti & all_genes$cell_type==ci, ]
      deg_list = unique(dat_tmp[dat_tmp$p_val_adj<0.05 & abs(dat_tmp$avg_log2FC) > 0.5, ]$gene_name)
      bg_list = unique(dat_tmp$gene_name)
      
      # if(unique(dat_tmp$strain)=="C57BL/6"){
      #   bp_list = bp_mouse
      # }else{
      #   bp_list = bp_rat
      # }
      
      x <- sum(deg_list %in% bp_list)  # Number of DEGs in the reference list
      m <- sum(bp_list %in% bg_list)  # Total number of genes in the reference list
      n <- length(bg_list) - m  # Total number of genes not in the reference list
      k <- length(deg_list)  # Total number of DEGs
      
      print(c(x/k, m/(n+m)))
      
      p_value <- phyper(x - 1, m, n, k, lower.tail = FALSE)
      enrich_df = rbind(enrich_df, c(ti, ci, i, x/k, m/(n+m), x/length(bp_list), p_value))
      
    }
    
  }
  
}

write.table(enrich_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.strain_wise.enrich.yong.out", col.names = T, row.names = F, sep = '\t', quote = F)
# write.table(enrich_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.strain_wise.enrich.helen.out", col.names = T, row.names = F, sep = '\t', quote = F)



# enrich_df_helen = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.enrich.helen.out", header=T, sep = '\t')
# enrich_df_yong = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.enrich.yong.out", header=T, sep = '\t')
# comp_level = c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d", "SD-HS 3d", "SS-HS 3d", "SS-HS 21d", "WKY-26w", "SHR-26w")
enrich_df_helen = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.strain_wise.enrich.helen.out", header=T, sep = '\t')
enrich_df_yong = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.strain_wise.enrich.yong.out", header=T, sep = '\t')
comp_level = c("C57BL/6", "SS", "SD", "SHR", "WKY")
enrich_df_helen$gene_list = "Keaton JM, et al.\n(n=1,873)"
enrich_df_yong$gene_list = "Manoj KM, et al.\n(n=251)"
enrich_df = as.data.frame(rbind(enrich_df_helen, enrich_df_yong))
colnames(enrich_df) = c("tissue", "cell_type", "comp", "prop.deg", "prop.all", "deg.prop", "pval", "gene_list")
enrich_df = enrich_df[!is.na(enrich_df$prop.deg),]
enrich_df$cell_type = factor(enrich_df$cell_type, levels=unique(cell_order))
enrich_df$tissue = factor(enrich_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
enrich_df$ji_score = as.numeric(enrich_df$prop.deg)
enrich_df$pval = as.numeric(enrich_df$pval)
enrich_df$padj = p.adjust(enrich_df$pval, method = "BH")
enrich_df$annotation <- ifelse(enrich_df$padj < 0.05, "*", "")
enrich_df$comp = factor(enrich_df$comp, levels = comp_level)

# enrich_df_use = enrich_df[enrich_df$comp %in% c("AngII-Salt-sensitive", "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),]

# enrich_df[enrich_df$ji_score>0.1, ]$ji_score = 0.1
ggplot(enrich_df[enrich_df$comp %in% c("C57BL/6", "SS", "SHR"),]) +
  geom_tile(mapping = aes(x = comp, y = cell_type, fill=ji_score)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low="white", high="red",
                      breaks = c(0, 0.1, 0.2),
                      labels = c("0%", "10%", "20%")) +
  # scale_fill_gradient(low="white", high="red") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="DEG\npercentage") +
  facet_grid(tissue ~ gene_list, 
             scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.genset_enrich.png", width=492/96, height=818/96, dpi=300)



library(patchwork)
dat_split <- 
  enrich_df %>% mutate(color_high = "red") %>%
  group_split(gene_list)
purrr::map(dat_split, ~{
  ggplot(.x) + geom_tile(mapping = aes(x = comp, y = cell_type, fill=ji_score)) +
             # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
             scale_y_discrete(limits=rev) +
             scale_fill_gradient(low="white", high="red") +
             theme(
               panel.grid.major.y = element_blank(),   # No horizontal grid lines
               # legend.position = c(1, 0.55),           # Put legend inside plot area
               legend.justification = c(1, 0.5),
               axis.text.y = element_text(colour = 'black'),
               axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
               # text = element_text(colour = 'black'),
               strip.text = element_text(colour = 'black')
             ) +
             labs(x="", y="", fill="DEG\npercentage") +
             facet_grid(tissue ~ ., 
                        scales = "free", space = "free_y")
}) %>% wrap_plots()



helen_list = union(bp_rat, bp_mouse)
helen_prior_list = c(
  "Gstm1", "Casq2", "Mef2d", "Btn2a1", "Myl12a", "Ccdc97", "Ckb", 
  "Foxn3", "Actn4", "Amz1", "Pcnx", "Fubp1", "Adra1a", "Grb10", 
  "Notch4", "Arid3b", "Ube3c", "Fgfr2", "Lnpep", "Tmem51", "Gpc6", 
  "Scaf11", "Trim33", "Col9a2", "Klhl23", "Slc39a10", "Il20rb", 
  "Clip2", "Abcc8", "Gtf2ird1", "Slc15a2", "Dnajc13", "Ankh", 
  "Bin1", "Prrx2", "Naglu"
)

broad_deg_list = unique(all_genes[all_genes$p_val_adj<0.05 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)
deg_list = unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)
deg_count = table(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)
deg_tissue_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "tissue")])$gene_name)
deg_strain_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "strain")])$gene_name)

overlapped = intersect(yong_list, helen_list)
length(overlapped)

length(intersect(deg_list, helen_list))
length(intersect(broad_deg_list, helen_list))
length(intersect(deg_list, yong_list))
length(intersect(broad_deg_list, yong_list))
length(intersect(deg_list, helen_prior_list))
length(intersect(broad_deg_list, helen_prior_list))

deg_count[deg_count>10 & names(deg_count) %in% yong_list]
deg_count[deg_count>10 & names(deg_count) %in% helen_list]
deg_count[deg_count>10 & names(deg_count) %in% helen_prior_list]

deg_tissue_count[deg_tissue_count>2 & names(deg_tissue_count) %in% yong_list]
# deg_tissue_count[deg_tissue_count>2 & names(deg_tissue_count) %in% helen_list]
deg_tissue_count[deg_tissue_count>2 & names(deg_tissue_count) %in% helen_prior_list]

deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% yong_list]
deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% helen_prior_list]

unique(all_genes[all_genes$p_val_adj<0.05 & all_genes$gene_name %in% overlapped, ]$gene_name)
unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$gene_name %in% overlapped & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)


library(ggVennDiagram)
# library(ggvenn)
# 
x = list("Keaton JM, et al.\n(prioritized)"=helen_prior_list, 
         #"Keaton JM, et al."=helen_list, 
         "Manoj KM, et al."=yong_list,
         "Tissue shared DEG\n(n>=2)"=names(deg_tissue_count[deg_tissue_count>2]),
         "Strain shared DEG\n(n=3)"=names(deg_strain_count[deg_strain_count>2]))
# ggvenn(
#   x, show_percentage = F,
#   fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
#   stroke_size = 0.5, set_name_size = 4, text_size = 6
# )

ggVennDiagram(x, label="count", label_alpha = 0) +
  scale_fill_gradient(low="white",high = "red")



deg_count = table(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)
deg_tissue_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "tissue")])$gene_name)
deg_strain_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), 
                                          c("gene_name", "strain")])$gene_name)

data <- as.data.frame(table(deg_count))
names(data) <- c("Gene", "Count")
data <- data[order(-data$Count),]
shared_genes_count <- sum(deg_tissue_count > 1)
p <- ggplot(data, aes(x = factor(Gene, levels = Gene), y = Count)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  geom_text(aes(label = Count), vjust = -0.1, color = "black", size = 4) +
  labs(x = "Number of shared cell types", y = "Number of DEGs") +
  theme(axis.text = element_text(color="black"))
mid_x <- length(unique(data$Gene)) / 2 + 1
xmax <- length(unique(data$Gene))
p1 = p + annotate("text", x = mid_x, y = max(data$Count)*0.9, 
             label = paste0("Cell type shared\n(n = ",format(shared_genes_count, big.mark = ","), ")"), size = 4) +
  annotate("segment", x = 2, xend = xmax, y = max(data$Count)*0.9, yend = max(data$Count)*0.9, colour = "black", size = 0.5)

data <- as.data.frame(table(deg_tissue_count))
names(data) <- c("Gene", "Count")
data <- data[order(-data$Count),]
shared_genes_count <- sum(deg_tissue_count > 1)
p <- ggplot(data, aes(x = factor(Gene, levels = Gene), y = Count)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  geom_text(aes(label = Count), vjust = -0.1, color = "black", size = 4) +
  labs(x = "Number of shared tissues", y = "Number of DEGs") +
  theme(axis.text = element_text(color="black"))
mid_x <- length(unique(data$Gene)) / 2 + 1
xmax <- length(unique(data$Gene))
p2 = p + annotate("text", x = mid_x, y = max(data$Count)*0.9, 
             label = paste0("Tissue shared\n(n = ",format(shared_genes_count, big.mark = ","), ")"), size = 4) +
  annotate("segment", x = 2, xend = xmax, y = max(data$Count)*0.9, yend = max(data$Count)*0.9, colour = "black", size = 0.5)

data <- as.data.frame(table(deg_strain_count))
names(data) <- c("Gene", "Count")
data <- data[order(-data$Count),]
shared_genes_count <- sum(deg_strain_count > 1)
p <- ggplot(data, aes(x = factor(Gene, levels = Gene), y = Count)) +
  geom_bar(stat = "identity", fill = "#F8766D") +
  geom_text(aes(label = Count), vjust = -0.1, color = "black", size = 4) +
  labs(x = "Number of shared strains", y = "Number of DEGs") +
  theme(axis.text = element_text(color="black"))
mid_x <- length(unique(data$Gene)) / 2 + 1
xmax <- length(unique(data$Gene))
p3 = p + annotate("text", x = mid_x, y = max(data$Count)*0.9, 
             label = paste0("Strain shared\n(n = ",format(shared_genes_count, big.mark = ","), ")"), size = 4) +
  annotate("segment", x = 2, xend = xmax, y = max(data$Count)*0.9, yend = max(data$Count)*0.9, colour = "black", size = 0.5)

p2+p3

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.png", width=600/96, height=330/96, dpi=300)



length(deg_tissue_count[deg_tissue_count==1 & names(deg_tissue_count) %in% yong_list])/length(deg_tissue_count[deg_tissue_count==1])
length(deg_tissue_count[deg_tissue_count>1 & names(deg_tissue_count) %in% yong_list])/length(deg_tissue_count[deg_tissue_count>1])
length(deg_tissue_count[deg_tissue_count>2 & names(deg_tissue_count) %in% yong_list])/length(deg_tissue_count[deg_tissue_count>2])
# length(deg_tissue_count[deg_tissue_count>3 & names(deg_tissue_count) %in% yong_list])/length(deg_tissue_count[deg_tissue_count>3])

length(deg_tissue_count[deg_tissue_count==1 & names(deg_tissue_count) %in% helen_list])/length(deg_tissue_count[deg_tissue_count==1])
length(deg_tissue_count[deg_tissue_count>1 & names(deg_tissue_count) %in% helen_list])/length(deg_tissue_count[deg_tissue_count>1])
length(deg_tissue_count[deg_tissue_count>2 & names(deg_tissue_count) %in% helen_list])/length(deg_tissue_count[deg_tissue_count>2])
# length(deg_tissue_count[deg_tissue_count>3 & names(deg_tissue_count) %in% helen_list])/length(deg_tissue_count[deg_tissue_count>3])

deg_tissue_count[deg_tissue_count>1 & names(deg_tissue_count) %in% yong_list]
deg_tissue_count[deg_tissue_count>1 & names(deg_tissue_count) %in% helen_prior_list]

length(deg_strain_count[deg_strain_count==1 & names(deg_strain_count) %in% yong_list])/length(deg_strain_count[deg_strain_count==1])
length(deg_strain_count[deg_strain_count>1 & names(deg_strain_count) %in% yong_list])/length(deg_strain_count[deg_strain_count>1])
# length(deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% yong_list])/length(deg_strain_count[deg_strain_count>2])

length(deg_strain_count[deg_strain_count==1 & names(deg_strain_count) %in% helen_list])/length(deg_strain_count[deg_strain_count==1])
length(deg_strain_count[deg_strain_count>1 & names(deg_strain_count) %in% helen_list])/length(deg_strain_count[deg_strain_count>1])
# length(deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% helen_list])/length(deg_strain_count[deg_strain_count>2])

deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% yong_list]
deg_strain_count[deg_strain_count>2 & names(deg_strain_count) %in% helen_prior_list]

deg_strain_count[deg_strain_count>1 & names(deg_strain_count) %in% yong_list]
deg_strain_count[deg_strain_count>1 & names(deg_strain_count) %in% helen_prior_list]

intersect(deg_list, intersect(yong_list, helen_list))

deg_tissue_count[deg_tissue_count>5]


################################################################################
### barplot


# Calculate overlap proportions
calculate_overlap <- function(deg_count, ref_list, count_label) {
  if (count_label == "Tissue") {
    df <- data.frame(
      Count_Label = factor(rep(c("=1", ">1"), each = 2), levels = c("=1", ">1")),
      Reference = factor(rep(c("Mishra MK, et al.\n(n = 152)", "Keaton JM, et al.\n(n = 1,128)"), 2)),
      Proportion = c(
        length(deg_count[deg_count == 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count == 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count > 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count > 1])
        # length(deg_count[deg_count > 2 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count > 2]),
        # length(deg_count[deg_count > 2 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count > 2])
      )* 100,
      Count = factor(rep(count_label, 4))
    )
  } else {
    df <- data.frame(
      Count_Label = factor(rep(c("=1", ">1"), each = 2), levels = c("=1", ">1")),
      Reference = factor(rep(c("Mishra MK, et al.\n(n = 152)", "Keaton JM, et al.\n(n = 1,128)"), 2)),
      Proportion = c(
        length(deg_count[deg_count == 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count == 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count > 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count > 1])
      )* 100,
      Count = factor(rep(count_label, 4))
    )
  }
  return(df)
}

# Create data frames for tissue and strain counts
deg_tissue_df <- calculate_overlap(deg_tissue_count, yong_list, "Tissue")
deg_strain_df <- calculate_overlap(deg_strain_count, yong_list, "Strain")

# Combine data frames
df <- bind_rows(deg_tissue_df, deg_strain_df)


ggplot(df, aes(x = Count_Label, y = Proportion, fill = Count_Label)) +
  geom_bar(stat = "identity", position = "stack") +
  # facet_wrap(~Count, scales = "free_x") +
  labs(x = "Number of shared tissues/strains", y = "Overlapped DEG proportion (%)") +
  geom_text(aes(label = paste0(round(Proportion, 2), "%")), position = position_stack(vjust = 0.5)) +
  theme(axis.text = element_text(color="black"),
    axis.title.x = element_text(hjust = 0.5),
    legend.position = "None") + 
  coord_flip() + 
  facet_grid(Count ~ Reference, scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.prop.png", width=323/96, height=335/96, dpi=300)


################################################################################
### pie chart

library(plotly)
library(htmlwidgets)
op <- options()
options(viewer = NULL)

col_tmp = c("#CD534CFF", "lightgrey")

deg_list = unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), ]$gene_name)

ref_df <- data.frame(Gene = yong_list[yong_list %in% all_genes$gene_name])
result_df <- ref_df %>%
  mutate(Status = ifelse(Gene %in% deg_list, "DEG", "Non-DEG")) %>%
  count(Status) %>%
  mutate(Proportion = n / sum(n))

fig1 <- plot_ly(result_df, labels = ~Status, values = ~Proportion, marker = list(colors=col_tmp), type = 'pie', pull = c(0, 0.2), rotation = 180)
fig1 <- fig1 %>% layout(#title = list(text='Manoj KM, et al.\n(n=251)', font=list(family="Arial", size = 15)),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      legend = list(orientation = 'h'))

fig1 %>% htmlwidgets::onRender(
  "function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'png', width: 210, height: 210, filename: 'deg.pie-chart-yong', scale: 5});
  }"
)

ref_df <- data.frame(Gene = helen_list[helen_list %in% all_genes$gene_name])
result_df <- ref_df %>%
  mutate(Status = ifelse(Gene %in% deg_list, "DEG", "Non-DEG")) %>%
  count(Status) %>%
  mutate(Proportion = n / sum(n))

fig2 <- plot_ly(result_df, labels = ~Status, values = ~Proportion, marker = list(colors=col_tmp), type = 'pie', pull = c(0, 0.2), rotation = 90)
fig2 <- fig2 %>% layout(#title = list(text='Keaton JM, et al.\n(n=1,873)', font=list(family="Arial", size = 15)),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      legend = list(orientation = 'h'))

fig2 %>% htmlwidgets::onRender(
  "function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'png', width: 210, height: 210, filename: 'deg.pie-chart-helen', scale: 5});
  }"
)


ref_df <- data.frame(Gene = helen_prior_list)
result_df <- ref_df %>%
  mutate(Status = ifelse(Gene %in% deg_list, "DEG", "Non-DEG")) %>%
  count(Status) %>%
  mutate(Proportion = n / sum(n))

fig3 <- plot_ly(result_df, labels = ~Status, values = ~Proportion, marker = list(colors=col_tmp), type = 'pie', rotation = 90)
fig3 <- fig3 %>% layout(#title = list(text='Keaton JM, et al.\n(prioritized, n=36)', font=list(family="Arial", size = 15)),
                      xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                      legend = list(orientation = 'h'))

fig3 %>% htmlwidgets::onRender(
  "function(el, x) {
  var gd = document.getElementById(el.id); 
  Plotly.downloadImage(gd, {format: 'png', width: 180, height: 180, filename: 'deg.pie-chart-3', scale: 5});
  }"
)

options(viewer = op$viewer)



# > length(helen_list[helen_list %in% all_genes$gene_name])
# [1] 1128
# > length(yong_list[yong_list %in% all_genes$gene_name])
# [1] 152
# > length(deg_list)
# [1] 3328
# > length(unique(all_genes$gene_name))
# [1] 15107
# length(intersect(helen_list, deg_list))
# [1] 270
# > length(intersect(yong_list, deg_list))
# [1] 50

phyper(50, 152, 15107-152, 3328, lower.tail = F)
phyper(270, 1128, 15107-1128, 3328, lower.tail = F)


################################################################################
### 3d scatter plot
library(scatterplot3d)

deg_lk_ec = all_genes[all_genes$tissue=="LK" & all_genes$cell_type=="EC" & all_genes$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_lk_ec$ss = -1 * log10(deg_lk_ec$p_val_adj) * deg_lk_ec$avg_log2FC
deg_lk_ec$ss = -1 * log10(deg_lk_ec$p_val_adj)

deg_lk_ec_mod = deg_lk_ec[, c("gene_name", "strain", "ss")]  %>% group_by(gene_name, strain) %>% 
  slice(which.max(ss)) %>% ungroup() %>% as.data.frame() %>%
  reshape(., idvar = "gene_name", timevar = "strain", direction = "wide")
deg_lk_ec_mod$overlap = ifelse(deg_lk_ec_mod$gene_name %in% c(yong_list, helen_list), "Yes", "No")
deg_lk_ec_mod = deg_lk_ec_mod[rowSums(is.na(deg_lk_ec_mod))==0,]
colors <- c("No"="black", "Yes"="red")[deg_lk_ec_mod$overlap]
s3d<-scatterplot3d(deg_lk_ec_mod[,2:4], pch = 16, color = colors, angle = 55, 
                   type="h", grid=TRUE, box=FALSE)
text(s3d$xyz.convert(deg_lk_ec_mod[,2:4]), labels = deg_lk_ec_mod$gene_name,
     cex= 1, col = "steelblue")


deg_mg = all_genes[all_genes$tissue=="HYP" & all_genes$cell_type=="Microglia" & all_genes$strain %in% c("C57BL/6", "SHR", "SS") & all_genes$treatment %in% c("HS 3d", "AngII 3d", "26w"), ]
deg_lk_ec$ss = -1 * log10(deg_lk_ec$p_val_adj) * deg_lk_ec$avg_log2FC
deg_lk_ec$ss = -1 * log10(deg_lk_ec$p_val_adj)

deg_lk_ec_mod = reshape(deg_lk_ec[, c("gene_name", "strain", "ss")], idvar = "gene_name", timevar = "strain", direction = "wide")
deg_lk_ec_mod$overlap = ifelse(deg_lk_ec_mod$gene_name %in% c(yong_list, helen_list), "Yes", "No")
deg_lk_ec_mod = deg_lk_ec_mod[rowSums(is.na(deg_lk_ec_mod))==0,]
colors <- c("No"="black", "Yes"="red")[deg_lk_ec_mod$overlap]
s3d<-scatterplot3d(deg_lk_ec_mod[,2:4], pch = 16, color = colors, angle = 55, 
                   type="h", grid=TRUE, box=FALSE)
text(s3d$xyz.convert(deg_lk_ec_mod[,2:4]), labels = deg_lk_ec_mod$gene_name,
     cex= 1, col = "steelblue")


################################################################################
# interaction matrix
library(reshape2)

data <- data.frame(
  "AngII" = c(1, 1, 0),
  "Salt-sensitive" = c(1, 0, 1),
  "Spontaneous" = c(0, 1, 1),
  Combination = c("AngII-Salt-sensitive",
                  "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),
  stringsAsFactors = FALSE
)

data$Combination = factor(data$Combination, data$Combination)
melted_data <- melt(data, id.vars = "Combination")
melted_data <- subset(melted_data, value == 1)
melted_data$variable = gsub("\\.", "-", melted_data$variable)

species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(1, 4, 2)])

ggplot(melted_data, aes(x = Combination, y = variable)) +
  geom_point(aes(color = variable), size = 3) +
  scale_colour_manual(values = species_col) +
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size=16, colour = 'black'),
        legend.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "None") +
  scale_y_discrete(limits = rev(c("AngII", "Salt-sensitive", "Spontaneous")), position = "right")






################################################################################
# 
# # ternary plot
# all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = T, sep='\t')
# mouse2rat=mouse2rat[mouse2rat$Gene.name %in% all_genes$gene_name,]
# 
# all_genes = process_df_with_orthologs(all_genes)
# all_genes = all_genes[!(is.na(all_genes$gene_name_ortho2m)), ]
# all_genes$cell_gene = paste(all_genes$tissue, all_genes$cell_type,  all_genes$gene_name_ortho2m, sep="-")
# 
# deg_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header = T, row.names = 1)
# deg_df = process_df_with_orthologs(deg_df)
# deg_df = deg_df[!(is.na(deg_df$gene_name_ortho2m)), ]
# deg_df$cell_gene = paste(deg_df$tissue, deg_df$cell_type,  deg_df$gene_name_ortho2m, sep="-")
# cell_gene = unique(deg_df[deg_df$p_val_adj<0.05 & abs(deg_df$avg_log2FC)>0.5, ]$cell_gene)
# all_genes_use = all_genes[all_genes$cell_gene %in% cell_gene & all_genes$species %in% c("C57BL/6", "SS", "SHR"), ]
# 
# all_genes_use = all_genes_use %>%
#   group_by(species, cell_gene) %>%
#   slice(which.max(abs(avg_log2FC))) %>%
#   ungroup() %>% as.data.frame()
# 
# all_genes_reshape = reshape(all_genes_use[,c("species", "cell_gene", "avg_log2FC")], idvar="cell_gene", timevar = "species", direction = "wide") %>%
#   filter(complete.cases(.))
# 
# 
# log2FC_min <- min(all_genes_reshape$avg_log2FC.C57BL_6, all_genes_reshape$avg_log2FC.SS, all_genes_reshape$avg_log2FC.SHR, na.rm = TRUE)
# log2FC_max <- max(all_genes_reshape$avg_log2FC.C57BL_6, all_genes_reshape$avg_log2FC.SS, all_genes_reshape$avg_log2FC.SHR, na.rm = TRUE)
# 
# breaks <- seq(log2FC_min, log2FC_max, length.out = 5)
# labels <- sprintf("%.2f", breaks) # Customize labels as needed
# quantiles <- quantile(unlist(all_genes_reshape[c("avg_log2FC.C57BL/6", "avg_log2FC.SS", "avg_log2FC.SHR")]), probs = seq(0, 1, 0.2))
# 
# p = ggtern(data = all_genes_reshape, mapping = aes(x = "avg_log2FC.C57BL/6", y = avg_log2FC.SS, z = avg_log2FC.SHR)) +
#   geom_mask()+
#   # stat_density_tern(geom = 'polygon', n = 400, aes(fill  = ..level.., alpha = ..level..)) +
#   geom_point() +
#   # scale_fill_gradient(low = "blue", high = "red", name = "", breaks = breaks, labels = labels) +
#   scale_L_continuous(breaks = quantiles, labels = scales::percent(quantiles)) +
#   scale_R_continuous(breaks = quantiles, labels = scales::percent(quantiles)) +
#   scale_T_continuous(breaks = quantiles, labels = scales::percent(quantiles))
#   
#   # guides(fill = guide_colorbar(order = 1), alpha = guide_none())
#   # theme_rgbw() +
#   # theme_noarrows()
# 



# venn plot for selected cell types
ti = "LV"; ci = "EC"
title = paste0(ti, "-", ci)
deg_use = deg_df[deg_df$tissue==ti & deg_df$cell_type==ci, ]

set1 <- unique(transformGeneNames(deg_use[deg_use$project=="AngII", ], "gene_name"))
set2 <- unique(transformGeneNames(deg_use[deg_use$project=="Salt-sensitive", ], "gene_name"))
set3 <- unique(transformGeneNames(deg_use[deg_use$project=="Spontaneous", ], "gene_name"))

# Create a list of the gene sets
geneSets <- list("AngII" = set1, "Salt-sensitive" = set2, "Spontaneous" = set3)

species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(1,4,2)])

overlappedGenes <- Reduce(intersect, geneSets)
# overlappedGenesText <- paste("Overlapped Genes:\n", paste(overlappedGenes, collapse = ", "))

ggvenn(geneSets, c("AngII", "Salt-sensitive", "Spontaneous"), show_percentage = F,
       fill_color = species_col, fill_alpha=0.3, set_name_size = 0) + labs(title = title)

cat(paste(overlappedGenes, collapse = ", "))


# fill_colors <- alpha(species_col, 0.3)
# # Draw the Venn diagram
# venn.plot <- venn.diagram(
#   x = geneSets,
#   category.names = c("AngII", "Salt-sensitive", "Spontaneous"),
#   col=species_col,
#   fill = fill_colors,
#   # cex = 0.5,
#   fontfamily = "Arial",
#   # main.fontfamily="serif",
#   
#   filename = NULL,
#   output = NULL # This will plot directly to the R plotting window
#   
# )
# 
# grid.draw(venn.plot)









# Highlight genes overlapped in all three sets by printing them
























# shannon index for deg results
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.merged.out", header = T, row.names = 1)

DEG_df = deg_merged %>% 
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>% 
  group_by(project, species, treatment, tissue, cell_type) %>%
  dplyr::select(project, species, treatment, tissue, cell_type, control_size, treatment_size, test_gene) %>%
  mutate(
    DEG_num = n(),
    mean_cell_size = (mean(control_size, na.rm = TRUE) + mean(treatment_size, na.rm = TRUE)) / 2,
    min_cell_size = pmin(control_size, treatment_size),
    max_cell_size = pmax(control_size, treatment_size),
    cell_size_diff = abs(control_size - treatment_size),
    combined = factor(interaction(species, treatment, sep = "-"), 
                      levels = c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d",
                                 "SS-HS 3d", "SS-HS 21d", "SD-HS 3d",
                                 "SHR-26w", "WKY-26w"))
  ) %>%
  unique() %>%
  as.data.frame()

DEG_df$cell_type = factor(DEG_df$cell_type, levels=unique(cell_order))
DEG_df$tissue = factor(DEG_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))


total_DEGs_per_tissue <- DEG_df %>%
  group_by(tissue, combined) %>%
  summarize(Total_DEGs = sum(DEG_num))

# Join the total DEGs back to the original data
deg_data <- merge(DEG_df, total_DEGs_per_tissue, by = "tissue")

# Calculate the proportion of DEGs for each cell type within each tissue
deg_data$Proportion <- with(deg_data, DEG_num / Total_DEGs)

# Calculate Shannon entropy for each tissue
entropy_per_tissue <- deg_data %>%
  group_by(tissue, combined.y) %>%
  summarize(Entropy = -sum(Proportion * log2(Proportion), na.rm = TRUE))

# View the entropy values to assess tissue responsiveness
print(entropy_per_tissue)




























# data <- data.frame(
#   "C57BL/6-AngII 3d" = c(1, 0, 0, 0, 1, 1, 0),
#   "C57BL/6-AngII 28d" = c(1, 0, 0, 0, 1, 1, 0),
#   "SS-HS 3d" = c(0, 1, 1, 0, 1, 0, 1),
#   "SS-HS 21d" = c(0, 1, 0, 0, 1, 0, 1),
#   "SD-HS 3d" = c(0, 0, 1, 0, 0, 0, 0),
#   "SHR-26w" = c(0, 0, 0, 1, 0, 1, 1),
#   "SWKY-26w" = c(0, 0, 0, 1, 0, 0, 0),
#   Combination = c("C57BL/6-AngII 3d-C57BL/6-AngII 28d", "SS-HS 3d-SS-HS 21d", 
#                   "SS-HS 3d-SD-HS 3d", "SHR-26w-WKY-26w", "AngII-Salt-sensitive",
#                   "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),
#   stringsAsFactors = FALSE
# )
# 
# data$Combination = factor(data$Combination, data$Combination)
# melted_data <- melt(data, id.vars = "Combination")
# melted_data <- subset(melted_data, value == 1)
# 
# species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
# alpha <- setNames(c(0.5, 1, 0.5, 1, 1), unique(deg_df$treatment))
# 
# combined = unique(ji_df$comp)
# combined_colors <- setNames(species_col[gsub("\\-.*", "", combined)], combined)
# combined_alphas <- setNames(alpha[gsub(".*\\-", "", combined)], combined)
# 
# 






















































































































################################################################################
### pathway 

merge_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", sep='\t', header = T, quote = "")
merge_all$ji_comp = paste(merge_all$species, merge_all$treatment,  sep="-")
merge_all$project = "AngII"
merge_all[merge_all$species %in% c("SD", "SS"), ]$project = "Salt-sensitive"
merge_all[merge_all$species %in% c("SHR", "WKY"), ]$project = "Spontaneous"
merge_all = merge_all[! grepl("Summary", merge_all$GroupID), ]

merge_sig = merge_all[merge_all$Log.q.value. < log10(0.05), ]

comp_list1 = list(c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d"),
                  c("SS-HS 3d", "SS-HS 21d"),
                  c("SS-HS 3d", "SD-HS 3d"),
                  c("SHR-26w", "WKY-26w"))

comp_list2 = list(c("AngII", "Salt-sensitive"),
                  c("Salt-sensitive", "Spontaneous"),
                  c("AngII", "Spontaneous"))

tissue_list = unique(merge_all$tissue)

ji_df = c()

for(ti in tissue_list){
  
  cell_list = unique(merge_sig[merge_sig$tissue==ti, ]$cell_type)
  
  for(ci in cell_list){
    
    deg_use = merge_sig[merge_sig$tissue==ti & merge_sig$cell_type==ci, ]
    all_gene_use = merge_all[merge_all$tissue==ti & merge_all$cell_type==ci, ]
    
    for(j in comp_list1){
      
      comp1 = j[1]
      comp2 = j[2]
      similarity = ""
      p_value = ""
      
      deg_c1 = unique(deg_use[deg_use$ji_comp==comp1, ]$pathway)
      deg_c2 = unique(deg_use[deg_use$ji_comp==comp2, ]$pathway)
      
      
      if(length(deg_c1)>0 & length(deg_c2)>0){
        similarity = calculateJaccardIndex(deg_c1, deg_c2)
        # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
        p_value = performPermutationTest(deg_c1, deg_c2, 
                                         all_gene_use$pathway,
                                         all_gene_use$pathway)
      }
      
      ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
      
    }
    
    
    deg_use = merge_sig[merge_sig$tissue==ti & merge_sig$cell_type==ci & 
                       merge_sig$species %in% c("C57BL/6", "SS", "SHR"), ]
    
    for(j in comp_list2){
      
      comp1 = j[1]
      comp2 = j[2]
      similarity = ""
      p_value = ""
      
      if(grepl("Ang", comp1) & !(grepl("Ang", comp2))){
        gene_name = "gene_name_ortho2m"
      }else{
        gene_name = "gene_name"
      }
      
      deg_c1 = unique(deg_use[deg_use$project==comp1, ]$pathway)
      deg_c2 = unique(deg_use[deg_use$project==comp2, ]$pathway)
      
      if(length(deg_c1)>0 & length(deg_c2)>0){
        similarity = calculateJaccardIndex(deg_c1, deg_c2)
        # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
        p_value = performPermutationTest(deg_c1, deg_c2, 
                                         all_gene_use$pathway,
                                         all_gene_use$pathway)
      }
      
      ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
      
    }
    
  }
  
}

ji_df = as.data.frame(ji_df)
colnames(ji_df) = c("tissue", "cell_type", "comp1", "comp2", "ji_score", "pval")
ji_df$cell_type = factor(ji_df$cell_type, levels=unique(cell_order))
ji_df$tissue = factor(ji_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
ji_df$ji_score = as.numeric(ji_df$ji_score)
ji_df$pval = as.numeric(ji_df$pval)
ji_df$padj = p.adjust(ji_df$pval, method = "BH")
ji_df$annotation <- ifelse(ji_df$padj < 0.05, "*", "")
ji_df[is.na(ji_df$annotation), ]$annotation=""
ji_df$comp = paste(ji_df[,3], ji_df[,4], sep="-")
ji_df$comp = factor(ji_df$comp, c("C57BL/6-AngII 3d-C57BL/6-AngII 28d", "SS-HS 3d-SS-HS 21d", 
                                  "SS-HS 3d-SD-HS 3d", "SHR-26w-WKY-26w", "AngII-Salt-sensitive",
                                  "AngII-Spontaneous", "Salt-sensitive-Spontaneous"))

ji_df_use = ji_df[ji_df$comp %in% c("AngII-Salt-sensitive", "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),]

ggplot(ji_df_use, aes(x = comp, y = cell_type, fill=ji_score)) +
  geom_tile() +
  geom_text(aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low="white", high="blue") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Jaccard index") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")

















































