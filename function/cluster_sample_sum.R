library(plyr)
library(dplyr)
library(stringr)
# library(RColorBrewer)
library(openxlsx)
# getPalette = colorRampPalette(brewer.pal(9, "Set1"))

cluster_sample_sum = function(df, outfile){
  
  ## df is wide cluster-sample table
  colnames(df) = paste("Cluster", colnames(df))
  
  wb <- createWorkbook()
  addWorksheet(wb, sheetName="count")
  writeData(wb, sheet="count", x=df)
  
  addWorksheet(wb, sheetName="prop in sample")
  writeData(wb, "prop in sample", df/rowSums(df))
  
  addWorksheet(wb, sheetName="prop in cluster")
  writeData(wb, "prop in cluster", t(t(df)/colSums(df)))

  saveWorkbook(wb, file = outfile, overwrite = T)
  
}

