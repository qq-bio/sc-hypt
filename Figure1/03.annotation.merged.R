
library(Seurat)
library(harmony)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")



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

# Define cell type levels based on tissue
get_levels_by_tissue <- function(tissue) {
  switch(
    tissue,
    "HYP" = c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", "Astrocyte", "Microglia", 
              "Activated microglia", "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
              "Tanycyte", "Ependymal cell", "Pars tuberalis cell", "EC", "E/P transition cell", 
              "Pericyte", "Fibroblast"),
    "LV" = c("CM", "EC", "Pericyte", "Fibroblast", "IMM"),
    "LK" = c("POD", "PT", "TL", "TAL", "DCT", "DCT/CT", "CT", "CD", "IC", "EC", "VSMC", 
             "Fibroblast", "IMM"),
    "MCA" = c("Microglia", "EC", "VSMC", "Pericyte", "Fibroblast", "IMM"),
    "MSA" = c("EC", "VSMC", "Fibroblast", "IMM", "Adipocyte")
  )
}

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
  merged_levels <- get_levels_by_tissue(tissue)
  
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
### Mapping EC Results Back with Annotations
################################################################################

input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.MCA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.HYP.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MCA.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.HYP.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LV.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MCA.EC.cluster.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.MSA.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.HYP.EC.cluster.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.cluster.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.cluster.rds"
)

anno_info <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/EC.anno.txt", header=T, sep='\t')
for(i in input_file){
  
  dataset <- gsub(".cluster.rds", "", basename(filepath))
  EC_object <- readRDS(filepath)
  outfile <- gsub("cluster.rds", "anno.rds", filepath)
  
  # Filter annotation info for this dataset
  anno_info_use <- anno_data[anno_data$dataset == dataset, ]
  reso <- unique(anno_info_use$resolution)
  cluster_use <- anno_info_use$cluster
  subclass_level2 <- anno_info_use$subclass_level2
  
  # Map annotations
  Idents(EC_object) <- reso
  EC_object$seurat_clusters <- EC_object@meta.data[[reso]]
  EC_object@meta.data$"subclass_level2" <- subclass_level2[match(EC_object$seurat_clusters, cluster_use)]
  
  saveRDS(EC_object, outfile)
  
}






################################################################################
### Mapping Immune Cell Results Back to Main Seurat Objects
################################################################################

# Load immune cell objects
ang_immune_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.immune_cell.anno.rds")
ss_immune_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.immune_cell.anno.rds")
sp_immune_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.immune_cell.anno.rds")

# Create a list of immune cell subclass annotations for replacement
subcluster_replace_list <- c(as.character(ang_immune_object$subclass_level1),
                             as.character(ss_immune_object$subclass_level1),
                             as.character(sp_immune_object$subclass_level1))
names(subcluster_replace_list) <- c(colnames(ang_immune_object), 
                                    colnames(ss_immune_object), 
                                    colnames(sp_immune_object))

input_files <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

project_list <- c("AngII", "Salt-sensitive", "Spontaneous")
names(project_list) <- c("mouse", "rat.ss", "rat.sp")

for (i in input_files) {
  
  seurat_object <- readRDS(i)
  
  subclass_level1 <- as.character(seurat_object$subclass_level1)
  names(subclass_level1) <- colnames(seurat_object)
  
  # Identify cells for which annotations will be replaced
  replace_list <- intersect(names(subclass_level1), names(subcluster_replace_list))
  subclass_level1[replace_list] <- subcluster_replace_list[replace_list]
  
  # Update and subset Seurat object based on new annotations
  seurat_object$subclass_level1 <- factor(subclass_level1, levels = cell_order)
  seurat_object <- subset(seurat_object, subclass_level1 %in% cell_order)
  
  saveRDS(seurat_object, i)
}







################################################################################
### Visualize Cell Types Using UMAP (Figure 1e)
################################################################################

# Define a blank theme for UMAP plots
blank_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)

# Input files for each dataset
input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.v2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.v2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.v2.rds"
)

# Project list for labeling
project_list <- c("AngII", "Salt-sensitive", "Spontaneous")
names(project_list) <- c("mouse", "rat.ss", "rat.sp")

# Define PDF output for UMAP plots
pdf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1e.umap.pdf", width = 350 / 96, height = 350 / 96)

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
png("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_type.legend.png", width = 2247, height = 1896, res = 300)
plot(NULL, xaxt = 'n', yaxt = 'n', bty = 'n', ylab = '', xlab = '', xlim = 0:1, ylim = 0:1)
legend("topleft", legend = names(cell_col), pch = 19, pt.cex = 2, cex = 1, bty = 'n',
       col = cell_col, ncol = 2)
mtext("Cell type", at = 0.1, cex = 1.2)
dev.off()






