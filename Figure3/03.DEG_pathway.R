library(Seurat)
library(dplyr)
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/")



################################################################################
### DEG
################################################################################
i = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.rds"
outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/",
                 gsub("rds", "DEG_all.out", basename(i)))

seurat_object <- readRDS(i)
cluster = "seurat_clusters"
Idents(seurat_object) = cluster
meta_table = seurat_object@meta.data

project_list = unique(meta_table$project)

deg_merged = c()
for(pi in project_list){
  
  strain_list = unique(meta_table[meta_table$project==pi, ]$strain)
  for(si in strain_list){
    
    cell_list = unique(meta_table[meta_table$project==pi &
                                    meta_table$strain==si, cluster])
    
    treatment = meta_table[meta_table$project==pi &
                             meta_table$strain==si, ]$treatment
    treatment = intersect(treatment_order, unique(treatment))
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








################################################################################
### Prepare input files for metascape
################################################################################
# functions
map_genes_to_ensembl <- function(genes_str,
                                 species = c("mouse", "rat"),
                                 return_mapping = FALSE) {
  species <- match.arg(species)
  
  if (species == "mouse") {
    OrgDb <- org.Mm.eg.db::org.Mm.eg.db
    ens_pat <- "^ENSMUSG[0-9]+$"
  } else {
    OrgDb <- org.Rn.eg.db::org.Rn.eg.db
    ens_pat <- "^ENSRNOG[0-9]+$"
  }
  
  # parse input
  genes <- unique(trimws(unlist(strsplit(genes_str, ","))))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return("")
  
  # separate Ensembl vs SYMBOL
  is_ens <- grepl(ens_pat, toupper(genes))
  ens_ids <- genes[is_ens]
  symbols <- genes[!is_ens]
  
  # map symbol to ensembl
  mapped <- if (length(symbols)) {
    AnnotationDbi::select(OrgDb,
                          keys = symbols,
                          keytype = "SYMBOL",
                          columns = "ENSEMBL")
  } else data.frame(SYMBOL = character(), ENSEMBL = character())
  
  # combine with original Ensembl IDs
  all_ensembl <- unique(c(ens_ids, mapped$ENSEMBL))
  all_ensembl <- all_ensembl[!is.na(all_ensembl)]
  
  result <- paste(all_ensembl, collapse = ",")
  
  if (return_mapping) {
    return(list(ensembl = result, mapping = mapped))
  } else {
    return(result)
  }
}


input_file = c("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.DEG_all.out")
logFC_threshold = 0.5
deg_for_metascape <- c()
for( i in input_file ){
  
  deg <- read.table(i, header = T, sep = " ")
  species <- "mouse"
  deg_sum <- deg %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) >= logFC_threshold) %>%
    dplyr::mutate(
      cell_type = paste0("C", cell_type),
      name = paste(strain, cell_type, sep = "_"),
      name = gsub(" ", "_", name),
      name = gsub("/", "", name)
    ) %>%
    group_by(name) %>%
    dplyr::mutate(n_sig = n()) %>%
    filter(n_sig > 20) %>%
    slice_max(order_by = abs(avg_log2FC), n = 3000, with_ties = FALSE) %>%
    distinct(gene_name, .keep_all = TRUE) %>%
    dplyr::summarise(
      n_sig = dplyr::first(n_sig),
      gene  = paste(gene_name, collapse = ","),
      gene_id = map_genes_to_ensembl(gene, species),
      .groups = "drop"
    ) %>% 
    dplyr::select(name, gene_id)
  
  deg_for_metascape <- rbind(deg_for_metascape, deg_sum)
  
}

colnames(deg_for_metascape) <- c("#Name", "Genes")

write.table(deg_for_metascape, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.pre_filter.scvi.strain_covar.gene_nb.hvg_2k.refined.merged.DEG.for_metascape.ensembl.tsv", col.names = T, row.names = F, sep = "\t", quote = F)











################################################################################
### Pathway enrichment results (Fig. S12c)
################################################################################
library(openxlsx)
library(forcats)

### summary across ec subtypes
cell_type_list = c(1:15)
pathway_res_merged = c()
for(i in cell_type_list){
  
  path_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/metascape/output/C", i, "/metascape_result.xlsx")
  if(file.exists(path_file)){
    pathway_res <- read.xlsx(path_file, sheet = 2)
    pathway_res <- pathway_res[pathway_res$`Log(q-value)`< log10(0.05), ]
    if(nrow(pathway_res)>0){
      pathway_res$category <- sapply(strsplit(as.character(pathway_res$Category), " "), `[`, 1)
      pathway_res$pathway <- paste0(pathway_res$Description, " (", pathway_res$category, ")")
      pathway_res$cluster_idx <- i
      pathway_res$cluster_name <- cluster_names[as.character(i)]
      pathway_res_merged <- rbind(pathway_res_merged, pathway_res)
    }
  }
}
pathway_res_merged$cluster_name = factor(pathway_res_merged$cluster_name, levels = cluster_names)

summary_list <- pathway_res_merged %>% group_by(cluster_name) %>%
  filter(grepl("Summary", GroupID), `Log(q-value)`<log(0.05)) %>% ungroup() %>% dplyr::select(pathway)

pathway_res_use <- pathway_res_merged %>%
  filter(
    pathway %in% summary_list$pathway,
    grepl("Member", GroupID),
    `Log(q-value)` < -log10(0.05)
  ) %>%
  mutate(
    `Log(q-value)_raw` = `Log(q-value)`,
    `Log(q-value)` = pmax(`Log(q-value)`, -30)
  ) %>%
  arrange(cluster_name, desc(`Log(q-value)`)) %>% 
  mutate(pathway = forcats::fct_inorder(pathway))

ggplot(pathway_res_use, aes(x = pathway, y = cluster_name, fill = -1 * `Log(q-value)`)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "", y = "", fill = "-log10(q-value)") +
  theme_classic() +
  theme(
    plot.margin = unit(c(1,1,1,2), "cm"),
    text = element_text(family = "Arial"), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black', size = 12),
    axis.text.y = element_text(colour = 'black', size = 12)
  ) +
  coord_flip() +
  scale_x_discrete(position = "top")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/EC.metascape.summary.png", width = 1330 / 96, height = 1150 / 96, dpi = 300)













