library(Seurat)
library(fgsea)
library(msigdbr)
library(dplyr)

set.seed(42)

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.HYP.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LV.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.multiomics.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.MCA.RNA.DEG_all.out",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.HYP.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.LV.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.LK.multiomics.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.MSA.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.MCA.RNA.DEG_all.out",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.HYP.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.LV.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.LK.multiomics.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.MSA.RNA.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.MCA.RNA.DEG_all.out",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.immune_cell.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.ss.immune_cell.DEG_all.out",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/rat.sp.immune_cell.DEG_all.out"
)


# hkr_sets <- c("hallmark", "kegg", "reactome")
# hkr_sets <- c("hallmark")
# msigdbr("Mus musculus", category = "C2") %>%
#   filter(grepl(paste(hkr_sets, collapse = "|"), gs_name, ignore.case = T)) %>% 
#   split(x = .$gene_symbol, f = .$gs_name) %>% length()
pathways_mouse = msigdbr("Mus musculus", category = "C2") %>% filter(gs_subcat != "CGP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

pathways_rat = msigdbr("Rattus norvegicus", category = "C2") %>% filter(gs_subcat != "CGP") %>% 
  split(x = .$gene_symbol, f = .$gs_name)

for( i in input_file[10:17] ){
  
  DEG = try(read.table(i), silent = T)
  
  if(class(DEG) != "try-error"){
    
    
    outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/fgsea/",
                     gsub("DEG_all.out", "GSEA.out", basename(i)))
    
    fgseaRes_merged = c()
    
    species_list = unique(DEG$species)
    
    for(si in species_list){
      
      if(si=="C57BL/6"){
        pathways = pathways_mouse
      }else{
        pathways = pathways_rat
      }
      
      tissue_list = unique(DEG[DEG$species==si, ]$tissue)
      
      
      for(tissue in tissue_list){
        
        cell_list = unique(DEG[DEG$species==si & DEG$tissue==tissue, ]$cell_type)
        
        for( ci in cell_list ){
          
          treatment = unique(DEG[DEG$species==si & DEG$tissue==tissue & DEG$cell_type==ci, ]$treatment)
          
          for( ti in treatment ){
            
            DEG_tmp = DEG[DEG$species==si & DEG$tissue==tissue & DEG$cell_type==ci & DEG$treatment==ti , ]
            if(any(DEG_tmp$p_val==0)){
              DEG_tmp[DEG_tmp$p_val==0, ]$p_val = min(DEG_tmp[DEG_tmp$p_val!=0, ]$p_val)
            }
            # prerank = -log10(DEG_tmp$p_val) * sign(DEG_tmp$avg_log2FC)
            # prerank = -log10(DEG_tmp$p_val) * DEG_tmp$avg_log2FC
            prerank = DEG_tmp$avg_log2FC
            
            names(prerank) = rownames(DEG_tmp)
            prerank = prerank[!is.na(prerank)]
            prerank = sort(prerank, decreasing = T)
            
            fgseaRes = fgsea(pathways = pathways, 
                             stats    = prerank,
                             minSize  = 10)
            
            fgseaRes$project = unique(DEG_tmp$project)
            fgseaRes$species = si
            fgseaRes$tissue = tissue
            fgseaRes$cell_type = ci
            fgseaRes$control = unique(DEG_tmp$control)
            fgseaRes$treatment = ti
            
            fgseaRes_mod = data.frame(lapply(fgseaRes, as.character), stringsAsFactors = F)
            fgseaRes_mod$leadingEdge = do.call(rbind, lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse = ",")))
            
            fgseaRes_merged = rbind(fgseaRes_merged, fgseaRes_mod)
            
          }
          
        }
        
        
      }
      
      
    }
    
    fgseaRes_merged = fgseaRes_merged[order(fgseaRes_merged$padj),]
    
    write.table(fgseaRes_merged, outfile, quote = F, sep = '\t')
    
    print(i)
  
    
  }
  
}






# tail(fgseaRes_merged[fgseaRes_merged$padj<0.05,])
# table(fgseaRes_merged[fgseaRes_merged$padj<0.05,]$cell_type)




