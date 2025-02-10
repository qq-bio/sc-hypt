library(Seurat)

seurat_files <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",

  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

cell_order <- c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", # "Microglia", "Activated microglia", 
    "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("EC", "VSMC", "E/P transition cell", "Pericyte", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
    "NK cells", "NKT", "T cells", "B cells")
)


qc_summary_list <- list()

required_columns <- c("origID", "seqID", "model", "strain", "treatment", "tissue", "num_cells", "num_features", "num_umi", "num_fragments", "num_peaks", "tss_enrichment", cell_order)

for (file_path in seurat_files) {
  
  seurat_obj <- readRDS(file_path)
  
  seqID_groups <- split(seurat_obj@meta.data, seurat_obj@meta.data$seqID2)
  
  group_summary <- lapply(names(seqID_groups), function(seqID) {
    
    group_data <- seurat_obj@meta.data[seurat_obj@meta.data$seqID2 == seqID, ]
    
    origID = unique(group_data$orig.ident)
    model = unique(group_data$project)
    strain = unique(group_data$strain)
    treatment = unique(group_data$treatment)
    tissue = unique(group_data$tissue)
    
    num_cells <- nrow(group_data) 
    num_features <- median(group_data$nFeature_RNA, na.rm = TRUE) 
    num_umi <- median(group_data$nCount_RNA, na.rm = TRUE) 
    
    cell_type_counts <- as.data.frame(table(group_data$subclass_level1)) 
    colnames(cell_type_counts) <- c("cell_type", "count")
    
    cell_type_counts_wide <- pivot_wider(cell_type_counts, names_from = cell_type, values_from = count, values_fill = 0)
    
    combined_summary <- cbind(
      data.frame(origID = origID, seqID = seqID, model = model, strain = strain, treatment = treatment, tissue = tissue, 
                 num_cells = num_cells, num_features = num_features, num_umi = num_umi),
      cell_type_counts_wide
    )
    
    if ("ATAC" %in% Assays(seurat_obj)) {
      num_fragments <- median(group_data$nCount_ATAC, na.rm = TRUE) 
      num_peaks <- median(group_data$nFeature_ATAC, na.rm = TRUE) 
      tss_enrichment <- median(group_data$TSSEnrichment, na.rm = TRUE) 
      
      combined_summary <- cbind(combined_summary, data.frame(
        num_fragments = num_fragments, num_peaks = num_peaks, tss_enrichment = tss_enrichment
      ))
    } else {
      combined_summary <- cbind(combined_summary, data.frame(
        num_fragments = NA, num_peaks = NA, tss_enrichment = NA
      ))
    }
    
    return(combined_summary)
  })
  
  qc_summary <- do.call(rbind, group_summary)
  
  missing_cols <- setdiff(required_columns, colnames(qc_summary))
  if (length(missing_cols) > 0) {
    qc_summary[missing_cols] <- NA 
  }
  
  qc_summary_list[[file_path]] <- qc_summary[required_columns]
  
  rm(seurat_obj)
  gc() 
}

qc_summary_df <- do.call(rbind, lapply(qc_summary_list, function(x) {
  missing_cols <- setdiff(required_columns, colnames(x))
  if (length(missing_cols) > 0) {
    x[missing_cols] <- NA
  }
  return(x)
}))

qc_summary_df[cell_order] <- lapply(qc_summary_df[cell_order], function(x) replace(x, is.na(x), 0))

write.table(qc_summary_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/QC.summary.out", sep = "\t", quote = F, col.names = T, row.names = F)





###
qc_summary_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/QC.summary.out", sep = "\t", header = T)
qc_summary_df$origID = gsub(" ", "", qc_summary_df$origID)
qc_summary_df$seqID = gsub(" ", "", qc_summary_df$seqID)

sample_seq = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/sample_seq.summary.out", header = T, sep = '\t', comment.char = "")
sample_seq$Seq.ID = gsub(" ", "", sample_seq$Seq.ID)

table(sample_seq[sample_seq$Assay=="snRNA-seq", ]$Seq.ID) %>% mean()
table(sample_seq[sample_seq$Assay=="snRNA-seq", ]$Seq.ID) %>% median()

table(sample_seq[sample_seq$Assay=="bulk mRNA-seq", ]$Seq.ID) %>% median()

length(unique(sample_seq$Sample.ID))


setdiff(qc_summary_df$origID, sample_seq$Seq.ID)
setdiff(sample_seq[sample_seq$Assay!="bulk mRNA-seq", ]$Seq.ID, qc_summary_df$origID)


setdiff(qc_summary_df$seqID, sample_seq$Seq.ID)
setdiff(sample_seq[sample_seq$Assay!="bulk mRNA-seq", ]$Seq.ID, qc_summary_df$seqID)








