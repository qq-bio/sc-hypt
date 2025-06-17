
lib_list <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_lib_fr_sra.txt", header = T)



################################################################################

setwd("/xdisk/mliang1/qqiu/data/scHTN/")

R1_files <- list.files(".", ".*R1.*gz", full.names = T, recursive = T)
length(R1_files)
R2_files <- list.files(".", ".*R2.*gz", full.names = T, recursive = T)
length(R2_files)

all(gsub("_R1_", "_R2_", R1_files)==R2_files)
# [1] TRUE

file_info <- list()

for (i in seq_along(R1_files)) {
  R1_path_parts <- strsplit(R1_files[i], "/")[[1]]
  R2_path_parts <- strsplit(R2_files[i], "/")[[1]]
  
  R1_file <- R1_path_parts[length(R1_path_parts)]
  R2_file <- R2_path_parts[length(R2_path_parts)]
  
  sample_id <- strsplit(R1_file, "_")[[1]][1]
  
  file_info[[i]] <- data.frame(sample_id = sample_id,
                               R1_file = R1_file,
                               R2_file = R2_file,
                               stringsAsFactors = FALSE)
}

file_info_df <- do.call(rbind, file_info)

table(toupper(lib_list$library_ID) %in% toupper(file_info_df$sample_id))

file_info_df <- file_info_df[file_info_df$sample_id %in% lib_list$library_ID, ]
file_info_df <- file_info_df[match(lib_list$library_ID, file_info_df$sample_id), ]

write.table(file_info_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_lib_fr_sra.txt")





################################################################################

file_info_df <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_lib_fr_sra.txt", header = T)

transfer_file_list <- c(file_info_df$R1_file, file_info_df$R2_file)
for(i in transfer_file_list){
  
  file_folder <- list.files(".", i, full.names = T, recursive = T)
  print(file_folder)
  
  mv_command <- paste0("mv /xdisk/mliang1/qqiu/data/scHTN/", file_folder, " /xdisk/mliang1/qqiu/data/scHTN/to_upload/")
  system(mv_command)
  
}












