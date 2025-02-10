library(plyr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(openxlsx)
getPalette = colorRampPalette(brewer.pal(9, "Set1"))




gene_m_h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out", header=T, sep='\t')
gene_r_h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out", header=T, sep='\t')

gene_rm_h = unique(rbind(gene_m_h[gene_m_h$Human.gene.name!="", c(2, 4)],
                         gene_r_h[gene_r_h$Human.gene.name!="", c(2, 4)]))




panglao_anno = function(df, outfile, anno_file="panglao"){
  
  if(anno_file=="panglao"){
    anno = read.table("/xdisk/mliang1/qqiu/reference/panglaoDB/PanglaoDB_markers_27_Mar_2020.tsv", header=T, sep='\t', fill=T)
    anno_mm = anno[anno$species %in% c("Mm"), ]
    anno_hs = anno[anno$species %in% c("Mm Hs", "Hs"), ]
    
    anno_mm$official.gene.symbol = str_to_title(anno_mm$official.gene.symbol)
    anno_hs = anno_hs %>% left_join(gene_rm_h, by=c("official.gene.symbol" = "Human.gene.name"), 
                                    relationship = "many-to-many") %>% 
      mutate(official.gene.symbol.old = official.gene.symbol) %>% 
      mutate(official.gene.symbol = Gene.name)
    anno = rbind(anno_mm, anno_hs[,1:14])
  }else{
    anno = read.table(anno_file, header=T, sep='\t', fill=T)
  }
  
  ## df is wide cluster-marker table
  gene_list = as.vector(unlist(df[, -1])); gene_list = unique(gene_list[!is.na(gene_list)])
  cell_count = table(anno[anno$official.gene.symbol %in% gene_list, ]$cell.type)
  cell_count = cell_count[order(cell_count, decreasing = T)]
  if(anno_file=="panglao"){
      nCell = length(cell_count > 5)
  }else{
      nCell = cell_count
    }
  myColors <- getPalette(nCell)
  names(myColors) <- names(cell_count)[1:nCell]
  anno_use = anno[anno$cell.type %in% names(myColors), ]
  
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName="reference")
  writeData(wb, sheet="reference", x=anno_use)
  addWorksheet(wb, sheetName="Merged")
  writeData(wb, "Merged", data.frame(t(names(myColors))), startCol = 1, startRow = 1, rowNames = F, colNames = F)
  writeData(wb, "Merged", startCol = 1, startRow = 3, x=df, colNames = T)
  
  for(i in 1:length(myColors)){
    cell = substr(names(myColors)[i], 1, 31)
    addWorksheet(wb, sheetName=cell)
    writeData(wb, cell, data.frame(t(names(myColors))), startCol = 1, startRow = 1, rowNames = F, colNames = F)
    writeData(wb, sheet=cell, startCol = 1, startRow = 3, x=df, colNames = T)
    
    # define style
    style = createStyle(fgFill=myColors[i])
    geneList = unique(anno[anno$cell.type == names(myColors)[i], ]$official.gene.symbol)
    
    for(y in 2:ncol(df)){
      x <- which(df[,y] %in% geneList)
      if(length(x)>0){
        addStyle(wb, sheet=cell, style=style, rows=1, cols=i, gridExpand=TRUE) # +1 for header line
        addStyle(wb, sheet=cell, style=style, rows=x+3, cols=y, gridExpand=TRUE) # +1 for header line
        
        addStyle(wb, sheet="Merged", style=style, rows=1, cols=i, gridExpand=TRUE) # +1 for header line
        addStyle(wb, sheet="Merged", style=style, rows=x+3, cols=y, gridExpand=TRUE) # +1 for header line
        
        x = which((anno_use$cell.type == names(myColors)[i]) & (anno_use$official.gene.symbol %in% df[,y]))
        addStyle(wb, sheet="reference", style=style, rows=x+1, cols=1:ncol(anno_use), gridExpand=TRUE)
      }
    }
  }
  
  
  tbl=table(anno[anno$cell.type %in% names(myColors), ]$official.gene.symbol)
  multi_geneList = names(tbl[tbl>1])
  style = createStyle(textDecoration="bold", borderStyle="thick", fontColour="red")

  for(y in 2:ncol(df)){
    x <- which(df[,y] %in% multi_geneList)
    if(length(x)>0){
      addStyle(wb, sheet="Merged", style=style, rows=x+3, cols=y, gridExpand=TRUE) # +1 for header line
    }
  }

  saveWorkbook(wb, file = outfile, overwrite = T)
  
}

