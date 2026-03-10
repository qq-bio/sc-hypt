
library(Seurat)
library(harmony)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")



################################################################################
### Annotate Major Clusters
################################################################################
input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.PBMC.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.cluster.rds"
)

# Load manually curated cluster annotation info and sample info
anno_info <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/major_cluster.anno.txt", header=TRUE, sep="\t", fill=TRUE)
sample_info <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_sample_info.txt", header=TRUE, sep='\t')
rownames(sample_info) <- sample_info$seqID

project_list <- setNames(c("AngII", "Salt-sensitive", "Spontaneous"), c("mouse", "rat.ss", "rat.sp"))

for (i in input_file) {
  
  # Extract metadata
  dataset <- gsub("\\.(RNA|multiomics)+.cluster.rds", "", basename(i), perl=TRUE)
  project <- gsub("\\.[A-Z]+", "", dataset)
  outfile <- gsub("cluster.rds", "anno.rds", i)
  
  seurat_object <- readRDS(i)
  
  # Map sample info
  sample_data <- sample_info[as.character(seurat_object$orig.ident),]
  seurat_object$cell_id <- colnames(seurat_object)
  seurat_object$seqID2 <- sample_data$seqID2
  seurat_object$sample <- sample_data$sample
  seurat_object$strain <- factor(sample_data$strain, levels=strain_order)
  seurat_object$tissue <- factor(sample_data$tissue.abbr, levels=tissue_order)
  seurat_object$treatment <- factor(sample_data$treatment, levels=treatment_order)
  seurat_object$project <- factor(project_list[project], levels=c("AngII", "Salt-sensitive", "Spontaneous"))
  
  # Determine cell type levels for annotation
  tissue <- unlist(strsplit(dataset, "\\."))[length(unlist(strsplit(dataset, "\\.")))]
  merged_levels <- cell_order
  
  # Retrieve annotation info for the current dataset
  anno_info_use <- anno_info[anno_info$dataset == dataset,]
  reso <- unique(anno_info_use$resolution)
  cluster_use <- anno_info_use[!(grepl("remove", anno_info_use$note)),]$cluster
  subclass_level1 <- anno_info_use[!(grepl("remove", anno_info_use$note)),]$subclass_level1
  class <- anno_info_use[!(grepl("remove", anno_info_use$note)),]$class
  
  # Filter Seurat object by cluster use
  Idents(seurat_object) <- reso
  seurat_object$seurat_clusters <- seurat_object@meta.data[, reso]
  seurat_object <- subset(seurat_object, seurat_clusters %in% cluster_use)
  seurat_object$seurat_clusters <- droplevels(seurat_object$seurat_clusters)
  
  # Annotate clusters
  seurat_object@meta.data$subclass_level1 <- factor(subclass_level1[match(seurat_object$seurat_clusters, cluster_use)], 
                                                    levels=merged_levels)
  seurat_object <- subset(seurat_object, subclass_level1 %in% merged_levels)
  seurat_object@meta.data$class <- factor(class[match(seurat_object$seurat_clusters, cluster_use)], 
                                          levels=class_order)
  
  saveRDS(seurat_object, outfile)
  print(c(dataset, all(unique(subclass_level1) %in% merged_levels)))
}








################################################################################
### Immune Cell Annotation Using SingleR 
################################################################################

#### Merge immune cells
### Function to load data and set the default assay
load_and_set_default_assay <- function(file_path, assay = "RNA") {
  obj <- readRDS(file_path)
  DefaultAssay(obj) <- assay
  return(obj)
}

### Function to subset immune cells based on resolution
subset_immune_cells <- function(object, resolution, clusters) {
  object@active.ident <- factor(object@active.ident)
  names(object@active.ident) <- colnames(object)
  subset(object, get(resolution) %in% clusters)
}

### Function to prepare reference data with Harmony and PCA
prepare_reference_data <- function(reference_list, project_name) {
  ref_object <- merge(reference_list[[1]], y = reference_list[-1], project = project_name)
  ref_object <- NormalizeData(ref_object) %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = c("percent.mt")) %>%
    RunPCA(features = VariableFeatures(ref_object), npcs = 30) %>%
    RunHarmony(group.by.vars = "orig.ident") %>%
    RunUMAP(reduction = "harmony", dims = 1:30, return.model = TRUE) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30)
  return(ref_object)
}

### Function to create and normalize query object
prepare_query_data <- function(query_object) {
  query_object <- CreateSeuratObject(counts = query_object@assays$RNA@counts,
                                     min.cells = 3, min.features = 200, meta.data = query_object@meta.data)
  query_object[["percent.mt"]] <- PercentageFeatureSet(query_object, pattern = "^[Mm]t-")
  NormalizeData(query_object)
}

### Function for anchor finding and query mapping
map_query_to_reference <- function(ref_object, query_object) {
  transfer_anchors <- FindTransferAnchors(reference = ref_object, query = query_object, k.anchor = 30,
                                          k.filter = NA, reference.reduction = 'pca', dims = 1:30)
  MapQuery(anchorset = transfer_anchors, query = query_object, reference = ref_object,
           reference.reduction = "pca", reduction.model = "umap")
}

### Function to merge and harmonize mapped query with reference
finalize_merged_object <- function(ref_object, query_object) {
  merged_object <- merge(ref_object, query_object)
  merged_object[["pca"]] <- merge(ref_object[["pca"]], query_object[["ref.pca"]])
  merged_object <- RunHarmony(merged_object, group.by.vars = "orig.ident", project.dim = FALSE) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony", dims = 1:30) %>%
    FindClusters(resolution = seq(0.5, 3, 0.5))
  return(merged_object)
}

### Function to process species data
process_species <- function(LK_infile, MCA_infile, LV_infile, clusters_LK, clusters_LV, clusters_MCA, output_file) {
  # Load files
  LK_object <- load_and_set_default_assay(LK_infile)
  MCA_object <- load_and_set_default_assay(MCA_infile)
  LV_object <- load_and_set_default_assay(LV_infile)
  
  # Subset immune cells
  LK_immune <- subset_immune_cells(LK_object, "wsnn_res.0.4", clusters_LK)
  LV_immune <- subset_immune_cells(LV_object, "RNA_snn_res.0.1", clusters_LV)
  MCA_immune <- subset_immune_cells(MCA_object, "RNA_snn_res.0.4", clusters_MCA)
  
  # Prepare reference and query data
  ref_object <- prepare_reference_data(list(MCA_immune, LV_immune), project_name = gsub("\\..*", "", LK_infile))
  query_object <- prepare_query_data(LK_immune)
  
  # Map query to reference and finalize merged object
  query_object <- map_query_to_reference(ref_object, query_object)
  merged_object <- finalize_merged_object(ref_object, query_object)
  
  # Save result
  saveRDS(merged_object, output_file)
}

### Process for Mouse Data
process_species(
  LK_infile = "mouse.LK.multiomics.anno.rds",
  MCA_infile = "mouse.MCA.RNA.anno.rds",
  LV_infile = "mouse.LV.RNA.anno.rds",
  clusters_LK = c(12, 17, 18),
  clusters_LV = c(2),
  clusters_MCA = c(18),
  output_file = "mouse.immune_cell.cluster.rds"
)

### Process for Rat SS Data
process_species(
  LK_infile = "rat.ss.LK.multiomics.anno.rds",
  MCA_infile = "rat.ss.MCA.RNA.anno.rds",
  LV_infile = "rat.ss.LV.RNA.anno.rds",
  clusters_LK = c(4, 13, 15),
  clusters_LV = c(4, 5),
  clusters_MCA = c(2),
  output_file = "rat.ss.immune_cell.cluster.rds"
)

### Process for Rat SP Data
process_species(
  LK_infile = "rat.sp.LK.multiomics.anno.rds",
  MCA_infile = "rat.sp.MCA.RNA.anno.rds",
  LV_infile = "rat.sp.LV.RNA.anno.rds",
  clusters_LK = c(12),
  clusters_LV = c(4),
  clusters_MCA = c(5, 18),
  output_file = "rat.sp.immune_cell.cluster.rds"
)



#### Annotate immune cells using singleR
library(SingleR)
library(Seurat)
library(celldex)

input_file = c("mouse.immune_cell.cluster.rds",
               "rat.ss.immune_cell.cluster.rds",
               "rat.sp.immune_cell.cluster.rds")

mimd.sc <- celldex::ImmGenData()

i = input_file[1]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = c("RNA_snn_res.1"))
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.1
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
main.group$labels = c(main.group$labels)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
saveRDS(seurat_object, outfile)


i = input_file[2]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
saveRDS(seurat_object, outfile)


i = input_file[3]
seurat_object = readRDS(i)
outfile = gsub("cluster.rds", "anno.rds", i)
DimPlot(seurat_object, label = T, group.by = "RNA_snn_res.0.5")
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.0.5
main.group <- SingleR(method = "cluster", sc_data = seurat_object@assays$RNA@data, ref = mimd.sc@assays@data$logcounts, types = mimd.sc$label.main, clusters=seurat_object$seurat_clusters)
seurat_object$subclass_level1 = main.group$labels[seurat_object$seurat_clusters]
saveRDS(seurat_object, outfile)







################################################################################
### Adding level-2 annotation
################################################################################
### functions
clean_format <- function(x){
  if (is.character(x)) {
    return(gsub("\\s+$", "", x))
  }
  x
}

add_sxt <- function(df){
  required <- c("strain", "treatment")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  df$sxt <- paste0(df$strain, "-", df$treatment)
  return(df)
}

clean_cols <- function(df, L2_label = "subclass_level2", ATAC=FALSE){
  meta_map <- c(
    Sample_ID          = "seqID2",
    Sample_ID_original = "orig.ident",
    nCount_RNA         = "nCount_RNA",
    nFeature_RNA       = "nFeature_RNA",
    Project            = "project",
    Strain             = "strain",
    Treatment          = "treatment",
    SxT                = "sxt",
    Tissue             = "tissue",
    Cell_ID            = "cell_id",
    Cell_type_L1       = "subclass_level1",
    Cell_type_L2       = L2_label
  )
  
  if(ATAC){
    meta_map <- c(meta_map, 
                  nCount_ATAC = "nCount_ATAC",
                  nFeature_ATAC = "nFeature_ATAC",
                  TSSEnrichment = "TSSEnrichment")
  }
  
  missing <- setdiff(unname(meta_map), colnames(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  
  out <- df[, unname(meta_map), drop = FALSE]
  colnames(out) <- names(meta_map)
  return(out)
}

clean_metadata <- function(df, L2_label = "subclass_level2"){
  df <- add_sxt(df)
  df[] <- lapply(df, clean_format)
  df <- clean_cols(df, L2_label = L2_label)
  return(df)
}

correct_labels <- function(df){
  
  stopifnot("popv_prediction" %in% names(df), "subclass_level1" %in% names(df))
  df[df$popv_prediction=="", ]$popv_prediction <- df[df$popv_prediction=="", ]$subclass_level1
  
  return(df)
}

update_seurat_with_meta <- function(
    so,
    meta_merged,
    cell_col = "Cell_ID",
    verbose = TRUE
) {
  stopifnot(is.data.frame(meta_merged), cell_col %in% colnames(meta_merged))
  
  ids_obj  <- colnames(so)
  ids_meta <- meta_merged[[cell_col]]
  keep <- intersect(ids_obj, ids_meta)
  
  if (length(keep) == 0L) {
    stop("No overlapping cells between Seurat object and meta_merged.")
  }
  if (verbose) {
    message("Keeping ", length(keep), " / ", length(ids_obj), " cells (",
            length(ids_obj) - length(keep), " dropped).")
  }
  
  so <- subset(so, cells = keep)
  
  new_md <- meta_merged[match(keep, ids_meta), , drop = FALSE]
  rownames(new_md) <- keep
  
  so@meta.data <- new_md
  return(so)
}


### files
so_popv_list <- data.frame(
  
  seurat_object_file = c(
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
    
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
    
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds",
    
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.PBMC.RNA.anno.rds"
  ),
  
  popv_label_file = c(
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/mouse.HYP.RNA.anno.v2.hypomap.c3.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/mouse.LV.RNA.anno.hca.heart.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/mouse.LK.multiomics.anno.v2.kpmp.sn_sc.retrain.popv_labels.tsv",
    
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.ss.HYP.RNA.anno.v2.hypomap.c3.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.ss.LV.RNA.anno.hca.heart.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.ss.LK.multiomics.anno.v2.kpmp.sn_sc.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.ss.MSA.RNA.anno.v2.hca.vascular.retrain.popv_labels.tsv",
    
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.sp.HYP.RNA.anno.v2.hypomap.c3.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.sp.LV.RNA.anno.hca.heart.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.sp.LK.multiomics.anno.v2.kpmp.sn_sc.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.sp.MSA.RNA.anno.v2.hca.vascular.retrain.popv_labels.tsv",
    "/xdisk/mliang1/qqiu/project/multiomics-hypertension/popv/rat.sp.MCA.RNA.anno.v2.hca.vascular.retrain.popv_labels.tsv",
    
    NA
  )
  
)


### Get metadata from immune cells
mouse_immune <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.immune_cell.anno.rds")
ss_immune <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.immune_cell.anno.rds")
sp_immune <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.immune_cell.anno.rds")

mouse_immune_meta <- clean_metadata(mouse_immune@meta.data, L2_label = "subclass_level1")
ss_immune_meta <- clean_metadata(ss_immune@meta.data, L2_label = "subclass_level1")
sp_immune_meta <- clean_metadata(sp_immune@meta.data, L2_label = "subclass_level1")

mouse_immune_meta$Cell_type_L1 <- "Immune cell"
ss_immune_meta$Cell_type_L1 <- "Immune cell"
sp_immune_meta$Cell_type_L1 <- "Immune cell"

immune_meta <- rbind(mouse_immune_meta, ss_immune_meta, sp_immune_meta)


### Combine and orgnalize metadata
all_meta <- list()
for (i in seq_len(nrow(so_popv_list))) {
  seurat_object <- readRDS(so_popv_list$seurat_object_file[i])
  obj_cells <- colnames(seurat_object)
  md <- seurat_object@meta.data
  
  if(!is.na(so_popv_list$popv_label_file[i])){
    
    md$subclass_level2 <- as.character(md$subclass_level1)
    
    popv_label <- read.table(so_popv_list$popv_label_file[i],
                             header = TRUE, sep = "\t",
                             comment.char = "", quote = "")
    
    if (!all(c("cell","popv_prediction") %in% colnames(popv_label))) {
      stop("popv_label must contain columns: 'cell' and 'popv_prediction'")
    }
    
    popv_label <- correct_labels(popv_label)
    
    overlap <- intersect(popv_label$cell, obj_cells)
    if (length(overlap) == 0L) {
      message(sprintf("[i=%d] No overlapping cells between Seurat object and popv labels.", i))
    } else {
      idx_obj   <- match(overlap, rownames(md))
      idx_popv  <- match(overlap, popv_label$cell)
      md$subclass_level2[idx_obj] <- popv_label$popv_prediction[idx_popv]
      
      if (nrow(popv_label) != length(obj_cells)) {
        message(sprintf("[i=%d] Popv labels (%d) != Seurat cells (%d). Using overlap = %d.",
                        i, nrow(popv_label), length(obj_cells), length(overlap)))
      }
    }
    
  }
  
  tissue <- unique(md$tissue)
  if(tissue!="PBMC"){
    
    if(tissue=="HYP"){
      immune_types <- c("Monocytes", "Macrophages", "DC", "Neutrophils",
                        "NK cells", "NKT", "T cells", "B cells", "IMM")
    }else{
      immune_types <- c("Microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
                        "NK cells", "NKT", "T cells", "B cells", "IMM")
    }
    
    keep <- !(md$subclass_level1 %in% immune_types)
    md <- md[keep, , drop = FALSE]
    
  }
  
  md_clean <- clean_metadata(md, L2_label = "subclass_level2")
  
  all_meta[[i]] <- md_clean
}

all_meta <- do.call(rbind, all_meta)
meta_merged <- rbind(all_meta, immune_meta[!(immune_meta$Cell_ID %in% all_meta$Cell_ID), ])



### Update metadata in seurat objects
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.PBMC.RNA.anno.rds"
)



for (f in input_file) {
  outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/", gsub("anno.rds", "anno.L2.rds", basename(f)))
  so <- readRDS(f)
  so <- update_seurat_with_meta(
    so,
    meta_merged = meta_merged,
    cell_col = "Cell_ID"
  )
  saveRDS(so, outfile)
}







################################################################################
### Visualize Cell Types Using UMAP (Figure 1e; Fig. S8a)
################################################################################

# Input files for each dataset
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LK.multiomics.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.MSA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.PBMC.RNA.anno.L2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MSA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MCA.RNA.anno.L2.rds"
)

# Project list for labeling
project_list <- c("AngII", "Salt-sensitive", "Spontaneous")
names(project_list) <- c("mouse", "rat.ss", "rat.sp")

# Define PDF output for UMAP plots
pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.umap.L1.pdf", width=350/96, height=350/96)

for (i in input_file) {
  
  seurat_object <- readRDS(i)
  
  dataset <- gsub("\\.(RNA|multiomics)+.*rds", "", basename(i), perl = TRUE)
  project <- gsub("\\.[A-Z]+", "", dataset)
  tissue <- sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
  title <- paste0(project_list[project], "-", tissue, "(N=",
                  prettyNum(ncol(seurat_object), big.mark = ',', scientific = FALSE), ")")
  
  # Choose UMAP reduction based on tissue type
  reduction <- ifelse(tissue == "LK", "wnn.umap.harmony", "umap")
  
  p <- DimPlot(seurat_object, reduction = reduction, group.by = "subclass_level1", label = FALSE, repel = FALSE) +
    blank_theme +
    scale_color_manual(values = cell_col) +
    ggtitle(title)
  
  print(p)
  
}

dev.off()


### Create and save a legend for cell types
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_type.legend.png", width=3500, height=2000, res=300)
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =names(cell_col), pch=19, pt.cex=2, cex=1, bty='n',
       col = cell_col, ncol=6)
mtext("Major cell type", at=0.1, cex=1.2)
dev.off()





