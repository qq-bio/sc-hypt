
library(Seurat)
library(DoubletFinder)
library(parallel)


################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder")

# Example sample lists
sample_list <- c("RMCA7SN", "RMCA8SN")

# Initialize parameter list to store results
para_list <- data.frame(
  sample_ID = character(),
  n_cell = numeric(),
  n_doublet = numeric(),
  doublet_rate = numeric(),
  PC = numeric(),
  pK = numeric()
)

# Process each sample
for (sample_id in sample_list) {
  
  # Define file paths
  h5_file <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellbender/", sample_id, "/", sample_id, "_cellbender_output_filtered.h5")
  outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DoubletFinder/", sample_id, "_doubletfinder.rds")
  
  # Skip processing if the output file already exists
  if (!file.exists(outfile)) {
    
    # Load counts
    counts <- Read10X_h5(h5_file)
    seurat_object <- CreateSeuratObject(
      counts = counts,
      project = sample_id,
      min.cells = 3,
      min.features = 200
    )
    
    # Quality control and normalization
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^[Mm]t-")
    seurat_object <- NormalizeData(seurat_object)
    seurat_object <- FindVariableFeatures(seurat_object)
    seurat_object <- ScaleData(seurat_object)
    seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
    seurat_object <- RunUMAP(seurat_object, dims = 1:40)
    
    # Determine optimal PCs based on cumulative variance
    stdv <- seurat_object[["pca"]]@stdev
    percent_stdv <- (stdv / sum(stdv)) * 100
    cumulative <- cumsum(percent_stdv)
    co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
    co2 <- sort(which(diff(percent_stdv) > 0.1), decreasing = TRUE)[1] + 1
    optimal_pc <- min(co1, co2)
    pc_list <- c(optimal_pc, 20, 30, 40)
    
    # Run DoubletFinder with multiple PC choices
    for (pci in seq_along(pc_list)) {
      pc <- pc_list[pci]
      
      # Identify optimal pK
      sweep.list <- paramSweep_v3(seurat_object, PCs = 1:pc, num.cores = detectCores() - 1)
      sweep.stats <- summarizeSweep(sweep.list)
      bcmvn <- find.pK(sweep.stats)
      optimal.pk <- as.numeric(as.character(bcmvn[which.max(bcmvn$BCmetric), "pK"]))
      
      # Define doublet rate and adjust for homotypic doublets
      doublet_rate <- min(ncol(seurat_object) * 8 * 1e-6, 0.3)
      homotypic.prop <- modelHomotypic(seurat_object@meta.data$seurat_clusters)
      nExp.poi <- round(doublet_rate * nrow(seurat_object@meta.data))
      nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
      
      # Run DoubletFinder
      seurat_object <- doubletFinder_v3(
        seu = seurat_object,
        PCs = 1:pc,
        pK = optimal.pk,
        nExp = nExp.poi.adj
      )
      
      # Update metadata and remove intermediate columns
      col_name <- if (pci == 1) "doublet_pc.optm" else paste0("doublet_pc.", pci)
      colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)] <- col_name
      seurat_object@meta.data <- seurat_object@meta.data[, !grepl("pANN", colnames(seurat_object@meta.data))]
      
      # Save parameters
      para_list <- rbind(para_list, data.frame(
        sample_ID = sample_id,
        n_cell = nrow(seurat_object@meta.data),
        n_doublet = table(seurat_object@meta.data[[col_name]])["Doublet"],
        doublet_rate = doublet_rate,
        PC = pc,
        pK = optimal.pk
      ))
    }
    
    # Save processed Seurat object
    saveRDS(seurat_object, outfile)
  }
}

write.table(para_list, "doubletFinder.parameter_list.v2.out", sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE, append = TRUE)
