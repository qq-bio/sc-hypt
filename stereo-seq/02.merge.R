dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(tidyverse)
py_bin <- reticulate::conda_list() %>% 
  filter(name == "stereopy") %>% 
  pull(python)

Sys.setenv(RETICULATE_PYTHON = py_bin)
library(reticulate)
library(SCP)
library(Seurat)
# library(SeuratDisk)
# library(harmony)
library(ggplot2)
library(patchwork)
library(dplyr)
library(purrr)
library(mistyR)
library(gridExtra)
library(ggthemes)
library(RColorBrewer)


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/")
sc <- import("scanpy")

# merge spatial data
sp_sample_list = c("952BR_B4-R", "952BR_E4-L", "952BR_E4-R", "954BR_F4-D", "954BR_F4-U")
sp_section_list = c()
for(i in sp_sample_list){
  
  infile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/",i,".bin50.QC.h5ad")
  adata <- sc$read_h5ad(infile)
  space_object <- adata_to_srt(adata)
  
  space_object$orig.ident = i
  space_object = NormalizeData(space_object)
  space_object = FindVariableFeatures(space_object) 
  
  sample_name = gsub("-", "_", i)
  assign(sample_name, space_object)
  sp_section_list = c(sp_section_list, sample_name)

}

### merge all samples together
space_merge = merge(get(sp_section_list[1]), y = map(sp_section_list[-1], get), 
                    add.cell.ids = sp_section_list, 
                    project = "hypt-stereo-seq")

### stereo-seq sample info
sample_info = data.frame(sp_id = c("952BR_B4-R", "952BR_E4-L", "952BR_E4-R", "954BR_F4-D", "954BR_F4-U"),
                         sample_id = c("SD#2", "SS#48", "SD#2", "SD#2", "SS#48"),
                         strain = c("SD"),
                         treatment = c("LS"),
                         tissue = rep("kidney", 5)
                         )








