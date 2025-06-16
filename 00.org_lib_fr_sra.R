
lib_list = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_lib_fr_sra.txt", header = T)


setwd("/xdisk/mliang1/qqiu/data/scHTN/")


R1_files = list.files(".", ".*R1.*gz", full.names = T, recursive = T)
length(R1_files)
R2_files = list.files(".", ".*R2.*gz", full.names = T, recursive = T)
length(R2_files)

