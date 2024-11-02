
library(Seurat)
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TFBSTools)
load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')


################################################################################
### https://github.com/GreenleafLab/ArchR/discussions/880
### https://github.com/GreenleafLab/ArchR/blob/c44b50c56323d5995fb3d2fb28a9072255fec644/R/AnnotationPeaks.R#L1053-L1146

### Functions for ArchR Enrichment Computation

computeEnrichment=function (matches = NULL, compare = NULL, background = NULL)
{
  matches <- .getAssay(matches, grep("matches", names(assays(matches)),
                                     value = TRUE, ignore.case = TRUE))
  matchCompare <- matches[compare, , drop = FALSE]
  matchBackground <- matches[background, , drop = FALSE]
  matchCompareTotal <- Matrix::colSums(matchCompare)
  matchBackgroundTotal <- Matrix::colSums(matchBackground)
  pOut <- data.frame(feature = colnames(matches), CompareFrequency = matchCompareTotal,
                     nCompare = nrow(matchCompare), CompareProportion = matchCompareTotal/nrow(matchCompare),
                     BackgroundFrequency = matchBackgroundTotal, nBackground = nrow(matchBackground),
                     BackgroundProporition = matchBackgroundTotal/nrow(matchBackground))
  pOut$Enrichment <- pOut$CompareProportion/pOut$BackgroundProporition
  pOut$mlog10p <- lapply(seq_len(nrow(pOut)), function(x) {
    p <- -phyper(pOut$CompareFrequency[x] - 1, pOut$BackgroundFrequency[x],
                 pOut$nBackground[x] - pOut$BackgroundFrequency[x],
                 pOut$nCompare[x], lower.tail = FALSE, log.p = TRUE)
    return(p/log(10))
  }) %>% unlist %>% round(4)
  pOut$mlog10Padj <- pmax(pOut$mlog10p - log10(ncol(pOut)),
                          0)
  pOut <- pOut[order(pOut$mlog10p, decreasing = TRUE), , drop = FALSE]
  pOut
}

.getAssay=function (se = NULL, assayName = NULL)
{
  .assayNames <- function(se) {
    names(SummarizedExperiment::assays(se))
  }
  if (is.null(assayName)) {
    o <- SummarizedExperiment::assay(se)
  }
  else if (assayName %in% .assayNames(se)) {
    o <- SummarizedExperiment::assays(se)[[assayName]]
  }
  else {
    stop(sprintf("assayName '%s' is not in assayNames of se : %s",
                 assayName, paste(.assayNames(se), collapse = ", ")))
  }
  return(o)
}




################################################################################
### Process snMultiome-seq Data using ArchR
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR")

# Sample and annotation setup
sample_lists <- list(
  angii_sample_list = paste0("MLK", c(1:6)),
  ss_sample_list = c(paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92"),
  shr_sample_list = c(paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)))
)

model_info <- list(
  angii = list(seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.v2.rds", outputDirectory = "angii"),
  ss = list(seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.v2.rds", outputDirectory = "ss"),
  sp = list(seurat_file = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.v2.rds", outputDirectory = "sp")
)

# Define reusable functions
initialize_project <- function(sample_list, geneAnnotation, genomeAnnotation, seurat_file, outputDirectory) {
  ArrowFiles <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", sample_list, ".arrow")
  proj <- ArchRProject(
    ArrowFiles = ArrowFiles,
    outputDirectory = outputDirectory,
    geneAnnotation = geneAnnotation,
    genomeAnnotation = genomeAnnotation,
    copyArrows = TRUE
  )
  seurat_object <- readRDS(seurat_file)
  proj <- subsetArchRProject(proj, cells = gsub("_", "#", colnames(seurat_object)), outputDirectory = outputDirectory, force = TRUE)
  proj$cellNames_mod <- gsub("#", "_", proj$cellNames)
  return(list(proj = proj, seurat_object = seurat_object))
}

add_metadata_and_umap <- function(proj, seurat_object) {
  # Add UMAP
  seurat_umap <- seurat_object@reductions$wnn.umap.harmony@cell.embeddings
  df <- DataFrame(row.names = proj$cellNames,
                  "seurat#wnn.umap.harmony1" = seurat_umap[proj$cellNames_mod, "wnnUMAPHarmony_1"],
                  "seurat#wnn.umap.harmony2" = seurat_umap[proj$cellNames_mod, "wnnUMAPHarmony_2"],
                  check.names = FALSE)
  proj@embeddings$wnn.umap.harmony <- SimpleList(df = df, params = list())
  
  # Add metadata
  seurat_meta <- seurat_object@meta.data[proj$cellNames_mod, ]
  for (j in colnames(seurat_meta)) {
    proj@cellColData[, j] <- seurat_meta[, j]
  }
  return(proj)
}

# Loop over each model in the model_list
for (model_name in names(sample_lists)) {
  sample_list <- sample_lists[[model_name]]
  model_settings <- model_info[[model_name]]
  
  if (grepl("MLK", sample_list[1])) {
    geneAnnotation <- geneAnnoMm10
    genomeAnnotation <- genomeAnnoMm10
  } else {
    load('/xdisk/mliang1/qqiu/reference/ArchR/rn7/rn7.ArchR_annotations.rda')
  }
  
  # Initialize project and add metadata/UMAP
  res <- initialize_project(sample_list, geneAnnotation, genomeAnnotation, model_settings$seurat_file, model_settings$outputDirectory)
  proj <- add_metadata_and_umap(res$proj, res$seurat_object)
  
  # Process project with ArchR functions
  proj <- addGroupCoverages(proj, groupBy = "subclass_level2", minCells = 15, threads = 1)
  proj <- addReproduciblePeakSet(proj, groupBy = "subclass_level2", pathToMacs2 = "/home/u1/qqiu/.conda/envs/macs2/bin/macs2", genomeSize = 1.86e9)
  proj <- addPeakMatrix(proj)
  proj <- addMotifAnnotations(proj, motifSet = "cisbp", name = "Motif", force = TRUE, species = "mus musculus")
  
  # Gene expression matrix setup
  counts <- res$seurat_object@assays$RNA@counts
  counts <- counts[rownames(counts) %in% geneAnnotation$genes$symbol, ]
  colnames(counts) <- gsub("[SN]*_", "#", colnames(counts))
  metadata <- res$seurat_object@meta.data
  rownames(metadata) <- gsub("[SN]*_", "#", rownames(metadata))
  index <- match(rownames(counts), geneAnnotation$genes$symbol)
  ordered_gr <- geneAnnotation$genes[index]; names(ordered_gr) <- NULL
  seRNA <- SummarizedExperiment(counts, metadata = metadata, rowRanges = ordered_gr)
  seRNA$Group <- paste0(metadata[gsub("#", "_", colnames(seRNA)), "orig.ident"])
  
  # Add Gene Expression Matrix
  proj <- addGeneExpressionMatrix(proj, seRNA = seRNA, force = TRUE)
  
  # Perform Iterative LSI on ATAC and RNA matrices
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
    saveIterations = FALSE,
    useMatrix = "TileMatrix",
    depthCol = "nFrags",
    name = "LSI_ATAC"
  )
  
  proj <- addIterativeLSI(
    ArchRProj = proj,
    clusterParams = list(resolution = 0.2, sampleCells = 10000, n.start = 10),
    saveIterations = FALSE,
    useMatrix = "GeneExpressionMatrix",
    depthCol = "Gex_nUMI",
    varFeatures = 2500,
    firstSelection = "variable",
    binarize = FALSE,
    name = "LSI_RNA"
  )
  
  # Combine LSI Dimensions and Perform UMAP and Clustering
  proj <- addCombinedDims(proj, reducedDims = c("LSI_ATAC", "LSI_RNA"), name = "LSI_Combined")
  proj <- addUMAP(proj, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)
  proj <- addClusters(proj, reducedDims = "LSI_Combined", name = "Clusters", resolution = 0.4, force = TRUE)
  
  # Integrate Gene Expression
  proj <- addGeneIntegrationMatrix(
    ArchRProj = proj,
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "LSI_Combined",
    seRNA = seRNA,
    addToArrow = TRUE,
    force = TRUE,
    groupRNA = "Group",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
  )
  
  # Add Peak to Gene Links
  proj <- addPeak2GeneLinks(
    ArchRProj = proj,
    reducedDims = "LSI_Combined"
  )
  
  saveArchRProject(ArchRProj = proj, outputDirectory = model_settings$outputDirectory, load = FALSE)
  saveRDS(proj, model_settings$outputDirectory)
}









################################################################################
### Major Cell Type Peak-to-Gene Linkage Analysis
################################################################################
archr_files <- list(
  mouse = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds",
  rat_ss = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds",
  rat_sp = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
)

for (infile in archr_files) {
  
  proj <- readRDS(infile)
  
  strain <- gsub(".LK.archr.rds", "", basename(infile))
  outputDirectory <- paste0("ArchRSubset_", strain)
  
  # Initialize storage for merged p2g results across cell types
  p2geneDF_merge <- data.frame()
  
  cell_type_list <- unique(proj@cellColData$subclass_level1)
  
  # Iterate over each major cell type
  for (i in cell_type_list) {
    cell_list <- rownames(proj@cellColData[proj@cellColData$subclass_level1 == i, ])
    proj_tmp <- subsetArchRProject(proj, cells = cell_list, outputDirectory = outputDirectory, force = TRUE)
    
    # Add peak-to-gene links
    proj_tmp <- addPeak2GeneLinks(
      ArchRProj = proj_tmp,
      k = 50, 
      reducedDims = "LSI_Combined"
    )
    
    # Retrieve peak-to-gene links with filtering
    p2geneDF <- metadata(proj_tmp@peakSet)$Peak2GeneLinks
    p2geneDF$gene_symbol <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
    p2geneDF$peakName <- paste0(seqnames(metadata(p2geneDF)$peakSet), ":", start(metadata(p2geneDF)$peakSet), "-", end(metadata(p2geneDF)$peakSet))[p2geneDF$idxATAC]
    p2geneDF <- as.data.frame(p2geneDF)
    p2geneDF$cell_type <- i
    
    # Save individual cell type results
    P2G_outfile <- gsub("rds", paste0(i, ".p2g.out"), infile)
    write.table(p2geneDF, P2G_outfile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
    
    # Accumulate results for merged output
    p2geneDF_merge <- rbind(p2geneDF_merge, p2geneDF)
  }
  
  # Save merged output for all cell types in the model
  merged_outfile <- gsub("archr.rds", "merged.p2g.out", infile)
  write.table(p2geneDF_merge, merged_outfile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
}




################################################################################
### EC subtype Peak-to-Gene Linkage Analysis
################################################################################
archr_files <- list(
  mouse = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds",
  ss = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.rds",
  sp = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.rds"
)

for (infile in archr_files) {
  
  proj <- readRDS(infile)
  
  # Extract strain information and define output directory
  strain <- gsub(".LK.archr.rds", "", basename(infile))
  outputDirectory <- paste0("ArchRSubset_", strain)
  
  # Identify unique EC subtypes
  ec_subtypes <- unique(proj@cellColData[proj@cellColData$subclass_level1 == "EC", ]$subclass_level2)
  
  # Loop through each EC subtype to perform peak-to-gene linkage analysis
  for (subtype in ec_subtypes) {
    
    # Subset project for the specific EC subtype
    cell_subset <- rownames(proj@cellColData[proj@cellColData$subclass_level2 == subtype, ])
    proj_sub <- subsetArchRProject(proj, cells = cell_subset, outputDirectory = outputDirectory, force = TRUE, dropCells = FALSE)
    
    p2g_success <- tryCatch({
      proj_sub <- addPeak2GeneLinks(
        ArchRProj = proj_sub,
        k = 50,  
        reducedDims = "LSI_Combined"
      )
      TRUE
    }, error = function(e) {
      message("addPeak2GeneLinks failed for subtype ", subtype, ": ", e$message)
      FALSE
    })
    
    if (p2g_success) {
      
      retrieve_success <- tryCatch({
        p2g <- getPeak2GeneLinks(
          ArchRProj = proj_sub,
          corCutOff = 0.45,
          resolution = 10000,
          returnLoops = FALSE
        )
        TRUE
      }, error = function(e) {
        message("getPeak2GeneLinks failed for subtype ", subtype, ": ", e$message)
        FALSE
      })
      
      if (retrieve_success) {
        p2g_df <- metadata(proj_sub@peakSet)$Peak2GeneLinks
        p2g_df$gene_symbol <- mcols(metadata(p2g_df)$geneSet)$name[p2g_df$idxRNA]
        p2g_df$peakName <- (metadata(p2g_df)$peakSet %>% {paste0(seqnames(.), ":", start(.), "-", end(.))})[p2g_df$idxATAC]
        p2g_df <- as.data.frame(p2g_df)
        p2g_df$cell_type <- subtype
        
        p2g_outfile <- gsub("rds", paste0(subtype, ".EC_sub.p2g.out"), infile)
        write.table(p2g_df, p2g_outfile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
      }
    }
  }
}







################################################################################
### Motif Enrichment for Peaks Linked to Marker Genes in EC and EC Subtypes
################################################################################
output_files <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.EC_sub.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
)

for (outfile in output_files) {
  
  proj <- readRDS(outfile)
  
  strain <- gsub(".LK.archr.EC_sub.rds", "", basename(outfile))
  
  ### Gather Peak-to-Gene (P2G) Linkage Data for EC and EC Subtypes
  p2g_list <- list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", pattern = paste0(strain, ".*EC.*p2g.out"))
  p2geneDF <- do.call(rbind, lapply(p2g_list, function(i) {
    read.table(i, header = TRUE, sep = '\t')
  }))
  
  ### Filter DEGs Based on Given Criteria
  DEG_file <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", strain, ".LK.EC.allmarker.0.25.long.txt")
  DEG <- read.table(DEG_file, header = TRUE, sep = "\t")
  DEG_use <- DEG %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
  
  ### Perform Motif Enrichment Analysis
  motif_merge <- data.frame()
  matches <- getMatches(proj, "Motif")
  cell_types <- unique(DEG_use$cluster)
  
  # Loop through each cell type to analyze motif enrichment
  for (cell_type in cell_types) {
    
    # Filter DEGs for the specific cell type
    DEG_filtered <- DEG_use %>% filter(cluster == cell_type)
    
    # Select unique peaks from P2G links that are significantly associated with DEGs
    selected_peaks <- p2geneDF %>%
      filter(!is.na(FDR) & FDR < 0.05 & gene_symbol %in% DEG_filtered$gene) %>%
      pull(idxATAC) %>%
      unique()
    
    # Compute motif enrichment if selected peaks are available
    if (length(selected_peaks) > 0) {
      enrichment_result <- computeEnrichment(matches, selected_peaks, seq_len(nrow(matches)))
      enrichment_result$cell_type <- cell_type  # Annotate results with cell type
      
      # Merge results across all cell types
      motif_merge <- rbind(motif_merge, enrichment_result)
    }
  }
  
  motif_outfile <- gsub("rds", "marker_gene_EC_p2g_motif.out", outfile)
  write.table(motif_merge, motif_outfile, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
  
}






################################################################################
### Compare Different Motif Enrichment Results (Figure S12b)
################################################################################
ss_motif_ec <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')
shr_motif_ec <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header=T, sep='\t')

for(i in c("ss_motif_ec", "shr_motif_ec")){
  
  p <- get(i) %>%
    filter(mlog10Padj > 2) %>%
    group_by(cell_type) %>%
    slice_max(order_by = mlog10Padj, n = 10) %>%  # Select top 10 features by mlog10Padj within each cell_type
    ungroup() %>%  # Ungroup after selection
    mutate(feature = factor(feature, levels = unique(feature))) %>%  # Convert feature to a factor
    ggplot(aes(x = cell_type, y = feature, fill = mlog10Padj)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(get(i)$mlog10Padj, na.rm = TRUE)) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(color = "black"),
      axis.text.x = element_text(color = "black", angle = 45, hjust = 1),
      panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # Add a black border
      panel.background = element_rect(fill = NA, color = NA),  # Ensure no background color
      panel.grid = element_blank()  # Remove grid lines for a cleaner look
    ) +
    labs(title = i,
         subtitle = "Top 10 motifs enriched in each cell type",
         x = "", y = "",
         fill = "-log10(p-adj)")
  print(p)
  
}





################################################################################
### Subtype-Specific Motif Analysis through Max Deviation (Figure S12c)
################################################################################
outfile <- list(
  "rat.ss.LK.archr.EC_sub.rds" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds",
  "rat.sp.LK.archr.EC_sub.rds" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
)

for (model in names(outfile)) {
  
  proj <- readRDS(outfile[[model]])
  strain <- gsub(".LK.archr.EC_sub.rds", "", basename(outfile[[model]]))
  
  # Extract Z-score for motif matrix for EC subtypes
  seGroupMotif <- getGroupSE(ArchRProj = proj, useMatrix = "MotifMatrix", groupBy = "subclass_level2")
  seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames == "z", ]
  
  # Calculate delta values for each cell type by subtracting mean values from other subtypes
  delta_list <- list()
  for (i in seq_len(ncol(seZ))) {
    subtype_name <- colnames(seZ)[i]
    delta_values <- assay(seZ)[, i] - rowMeans(assay(seZ)[, -i])
    delta_df <- data.frame(
      feature = rowData(seZ)$name,
      delta = delta_values,
      cell_type = subtype_name
    )
    delta_list[[subtype_name]] <- delta_df
  }
  delta_df_combined <- do.call(rbind, delta_list)
  delta_file <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", strain, ".LK.archr.EC_sub.motif_delta.out")
  write.table(delta_df_combined, delta_file, col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
}

# Load motif enrichment results and delta values for comparison
ss_motif_ec <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header = TRUE, sep = '\t')
shr_motif_ec <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.marker_gene_EC_p2g_motif.out", header = TRUE, sep = '\t')
delta_ss <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.motif_delta.out", header = TRUE, sep = '\t')
delta_shr <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.motif_delta.out", header = TRUE, sep = '\t')

# Identify overlapping motifs between salt-sensitive and spontaneous models
overlap_motif <- intersect(
  ss_motif_ec[ss_motif_ec$mlog10Padj > 1.3 & ss_motif_ec$cell_type == "EC_Mecom", ]$feature,
  shr_motif_ec[shr_motif_ec$mlog10Padj > 1.3 & shr_motif_ec$cell_type == "EC_Mecom", ]$feature
)

# Merge delta values and filter for overlapping motifs
merged_data <- merge(delta_ss, delta_shr, by = c("feature", "cell_type")) %>%
  filter(cell_type == "Ang EC") %>%
  filter(feature %in% overlap_motif) %>%
  mutate(TFRegulator = ifelse(delta.x > quantile(delta.x, 0.75) & delta.y > quantile(delta.y, 0.75), "Yes", "No")) %>%
  arrange(TFRegulator)

# Plot the results
ggplot(merged_data, aes(delta.x, delta.y, color = TFRegulator)) +
  geom_point() +
  ggrepel::geom_label_repel(
    data = merged_data[merged_data$TFRegulator == "Yes", ],
    aes(x = delta.x, y = delta.y, label = feature),
    size = 4,
    nudge_x = 0,
    color = "black",
    max.overlaps = Inf
  ) +
  geom_vline(xintercept = quantile(merged_data$delta.x, 0.75), lty = "dashed") +
  geom_hline(yintercept = quantile(merged_data$delta.y, 0.75), lty = "dashed") +
  theme_classic() +
  scale_color_manual(values = c("No" = "darkgrey", "Yes" = "firebrick3")) +
  scale_y_continuous(
    expand = c(0, 0),
    limits = c(min(merged_data$delta.y),
               max(merged_data$delta.y) * 2)) +
  scale_x_continuous(
    expand = c(0, 0),
    limits = c(min(merged_data$delta.x),
               max(merged_data$delta.x) * 2)) +
  xlab("Z-score difference between subtypes\n(in salt-sensitive model)") +
  ylab("Z-score difference between subtypes\n(in spontaneous model)") +
  labs(color = "Subtype specific\nTF regulator") +
  facet_wrap(~ cell_type)







################################################################################
### Visualize Selected Motif Activity (Figure 3h)
################################################################################
proj_files <- list(
  "rat.ss" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds",
  "rat.sp" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
)

seurat_files <- list(
  "rat.ss" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds",
  "rat.sp" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds"
)

# Function to update ArchR embeddings with Seurat UMAP coordinates
update_embeddings <- function(proj, seurat_object) {
  cell_list <- gsub("#", "_", rownames(proj@embeddings$wnn.umap.harmony$df))
  proj@embeddings$wnn.umap.harmony$df[, 1] <- seurat_object@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 1]
  proj@embeddings$wnn.umap.harmony$df[, 2] <- seurat_object@reductions$wnn.umap.harmony@cell.embeddings[cell_list, 2]
  return(proj)
}

# Load data and update embeddings
proj_ss <- update_embeddings(readRDS(proj_files[["rat.ss"]]), readRDS(seurat_files[["rat.ss"]]))
proj_shr <- update_embeddings(readRDS(proj_files[["rat.sp"]]), readRDS(seurat_files[["rat.sp"]]))

# Define the motif list
motifs <- c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", 
            "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", 
            "Prrx1_455", "Vax1_413", "Vax2_497")

# Function to visualize motif activity in UMAP
visualize_motifs <- function(proj, motifs, strain, output_dir) {
  for (motif in motifs) {
    try({
      # Select motif features
      markerMotifs <- getFeatures(proj, select = paste(motif, collapse = "|"), useMatrix = "MotifMatrix")
      markerMotifs <- grep("z:", markerMotifs, value = TRUE)
      
      # Plot motif activity in UMAP
      p <- plotEmbedding(
        ArchRProj = proj,
        colorBy = "MotifMatrix",
        name = sort(markerMotifs),
        embedding = "wnn.umap.harmony",
        imputeWeights = getImputeWeights(proj),
        size = 2,
        plotAs = "points"
      ) + labs(x = "UMAP 1", y = "UMAP 2", title = markerMotifs, color = "Motif activity")
      
      # Save the plot
      outfile <- file.path(output_dir, paste0(strain, ".", motif, ".umap.png"))
      ggsave(outfile, plot = p, width = 523 / 96, height = 377 / 96, dpi = 300)
      print(paste("Saved:", outfile))
      
    }, error = function(e) {
      message(paste("Error processing motif:", motif, "-", e$message))
    })
  }
}

output_dir <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/TFregulatorPlots"

visualize_motifs(proj_ss, motifs, "rat.ss", output_dir)
visualize_motifs(proj_shr, motifs, "rat.sp", output_dir)







##########################################################################################
### Identify Regulatory Targets of TFs (Figure 3i and S12d)
#### reference: https://github.com/GreenleafLab/scScalpChromatin/blob/main/Figure_3_TFregulators.R
##########################################################################################
# library(RcppAlgos)
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/archr_helpers.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/matrix_helpers.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/ArchR_helpers/misc_helpers.R")


proj_files <- list(
  "rat.ss" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.ss.LK.archr.EC_sub.rds",
  "rat.sp" = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/rat.sp.LK.archr.EC_sub.rds"
)


for (strain in names(proj_files)) {
  
  atac_proj <- readRDS(proj_files[[strain]])
  p2g_list <- list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/", pattern=paste0(strain, ".*EC.*p2g.out"))
  p2geneDF <- do.call(rbind, lapply(p2g_list, read.table, header=TRUE, sep='\t'))
  
  
  # ChromVAR deviations matrix: (rows motif names x cols cell names)
  motifMatrix <- getMatrixFromProject(atac_proj, useMatrix="MotifMatrix")
  deviationsMatrix <- assays(motifMatrix)$deviations
  
  # GeneIntegration Matrix: (rows gene names x cols cell names)
  GIMatrix <- getMatrixFromProject(atac_proj, useMatrix="GeneIntegrationMatrix")
  GImat <- assays(GIMatrix)$GeneIntegrationMatrix
  rownames(GImat) <- rowData(GIMatrix)$name
  GImat <- as(GImat[Matrix::rowSums(GImat) > 0,], "sparseMatrix") # Remove unexpressed genes
  
  regulators <- c("Hlx_518", "Hoxa1_479", "Hoxa4_396", "Hoxa5_510", "Hoxb4_511", "Hoxc4_581", "Hoxc6_406", "Lhx3_459", "Lmx1b_515", "Prop1_539", "Prrx1_455", "Vax1_413", "Vax2_497" )
  deviationsMatrix <- deviationsMatrix[regulators,]
  
  # Identify pseudobulks for performing matrix correlations
  knn_groups <- getLowOverlapAggregates(atac_proj, target.agg=500, k=100, 
                                        overlapCutoff=0.8, dimReduc="LSI_Combined")
  
  kgrps <- unique(knn_groups$group)
  
  # GeneIntegrationMatrix
  GIMatPsB <- lapply(kgrps, function(x){
    use_cells <- knn_groups[knn_groups$group==x,]$cell_name
    Matrix::rowMeans(GImat[,use_cells])
  }) %>% do.call(cbind,.)
  colnames(GIMatPsB) <- kgrps
  
  # In rare instances, we can get pseudo-bulked genes that have zero averages
  GIMatPsB <- GIMatPsB[Matrix::rowSums(GIMatPsB) > 0,]
  
  # DeviationsMatrix
  DevMatPsB <- lapply(kgrps, function(x){
    use_cells <- knn_groups[knn_groups$group==x,]$cell_name
    Matrix::rowMeans(deviationsMatrix[,use_cells])
  }) %>% do.call(cbind,.)
  colnames(DevMatPsB) <- kgrps
  
  # Perform chromVAR deviations to Integrated RNA correlation analysis:
  start <- Sys.time()
  geneCorMat <- cor2Matrices(DevMatPsB, GIMatPsB)
  colnames(geneCorMat) <- c("motifName", "symbol", "Correlation", "FDR")
  end <- Sys.time()
  message(sprintf("Finished correlations in %s minutes.", round((end  - start)/60.0, 2)))
  
  allGenes <- rownames(GIMatPsB) %>% sort() # Already filtered to only expressed genes
  
  
  # Get locations of motifs of interest:
  motifPositions <- getPositions(atac_proj, name="Motif")
  motifGR <- stack(motifPositions, index.var="motifName")
  
  # Get peak to gene GR
  corrCutoff <- 0.45
  metadata(atac_proj@peakSet)$Peak2GeneLinks = p2geneDF
  p2gDF <- metadata(proj@peakSet)$Peak2GeneLinks
  p2gDF$symbol <- mcols(metadata(p2gDF)$geneSet)$name[p2gDF$idxRNA] %>% as.character()
  p2gDF$peakName <- (metadata(p2gDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2gDF$idxATAC]
  p2gDF <- p2gDF[!is.na(p2gDF$Correlation),]
  p2gDF <- p2gDF[(p2gDF$Correlation > corrCutoff),]
  p2gDF <- p2gDF[which(p2gDF$VarQATAC > 0.25 & p2gDF$VarQRNA > 0.25),]
  p2gGR <- metadata(p2gDF)$peakSet[p2gDF$idxATAC]
  mcols(p2gGR) <- p2gDF
  p2gGR
  
  calculateLinkageScore <- function(motifLocs, p2gGR){
    # Calculate Linkage Score (LS) for each gene in p2gGR with regards to a motif location GR
    ###################################
    # For a given gene, the LS = sum(corr peak R2 * motifScore)
    ol <- findOverlaps(motifLocs, p2gGR, maxgap=0, type=c("any"), ignore.strand=TRUE)
    olGenes <- p2gGR[to(ol)]
    olGenes$motifScore <- motifLocs[from(ol)]$score
    olGenes$R2 <- olGenes$Correlation**2 # All p2g links here are already filtered to only be positively correlated
    LSdf <- mcols(olGenes) %>% as.data.frame() %>% group_by(symbol) %>% summarise(LS=sum(R2*motifScore)) %>% as.data.frame()
    LSdf <- LSdf[order(LSdf$LS, decreasing=TRUE),]
    LSdf$rank <- 1:nrow(LSdf)
    return(LSdf)
  }
  
  calculateMotifEnrichment <- function(motifLocs, p2gGR){
    # Calculate Motif enrichment per gene
    ###################################
    # For a given gene, calculate the hypergeometric enrichment of motifs in 
    # linked peaks (generally will be underpowered)
    motifP2G <- p2gGR[overlapsAny(p2gGR, motifLocs, maxgap=0, type=c("any"), ignore.strand=TRUE)]
    m <- length(motifP2G) # Number of possible successes in background
    n <- length(p2gGR) - m # Number of non-successes in background
    
    motifLinks <- motifP2G$symbol %>% getFreqs()
    allLinks <- p2gGR$symbol %>% getFreqs()
    df <- data.frame(allLinks, motifLinks=motifLinks[names(allLinks)])
    df$motifLinks[is.na(df$motifLinks)] <- 0
    df$mLog10pval <- apply(df, 1, function(x) -phyper(x[2]-1, m, n, x[1], lower.tail=FALSE, log.p=TRUE)/log(10))
    df <- df[order(df$mLog10pval, decreasing=TRUE),]
    df$symbol <- rownames(df)
    return(df)
  }
  
  regPlotDir <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/TFregulatorPlots"
  dir.create(regPlotDir, showWarnings = FALSE, recursive = TRUE)
  
  # plot all TF regulators
  markerGenes  <- c("Atp1b1", "Chrm3", "Col4a3", "Lrp2", "Nedd4l", "Nox4", "Slc16a12", "Slc4a4", "Slit2",
                    "Arhgap24", "Efna5", "Lgf1r", "Pkhd1", "Ptprd", "Vegfa", "Mecom") %>% unique()
  
  res_list <- list()
  ##########################################################################################
  regulators <- c("Hoxb4_511")
  
  for(motif in regulators){
    motif_short <- strsplit(motif,"_")[[1]][1]
    # First get motif positions
    motifLocs <- motifGR[motifGR$motifName == motif]
    # Calculate Linkage Score for motif
    LS <- calculateLinkageScore(motifLocs, p2gGR)
    # Get just genes correlated to motif
    motifGeneCorDF <- unique(geneCorMat[geneCorMat$motifName == motif,])
    plot_df <- merge(LS, motifGeneCorDF, by="symbol", all.x=TRUE)
    # Calculate motif enrichment per gene
    ME <- calculateMotifEnrichment(motifLocs, p2gGR)
    plot_df <- merge(plot_df, ME, by="symbol", all.x=TRUE)
    plot_df = plot_df[rowSums(is.na(plot_df))==0,]
    plot_df <- plot_df[,c("symbol", "LS", "Correlation", "FDR", "mLog10pval")]
    plot_df$toLabel <- "NO"
    topN <- 5
    plot_df <- plot_df[order(plot_df$LS, decreasing=TRUE),]
    plot_df$rank_LS <- 1:nrow(plot_df)
    plot_df$toLabel[1:topN] <- "YES"
    plot_df <- plot_df[order(plot_df$Correlation, decreasing=TRUE),]
    plot_df$rank_Corr <- 1:nrow(plot_df)
    plot_df$toLabel[1:topN] <- "YES"
    plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=TRUE),]
    plot_df$rank_Pval <- 1:nrow(plot_df)
    plot_df$toLabel[1:10] <- "YES"
    # plot_df$meanRank <- apply(plot_df[,c("rank_LS", "rank_Corr", "rank_Pval")], 1, mean)
    # plot_df <- plot_df[order(plot_df$meanRank, decreasing=FALSE),]
    # plot_df$toLabel[1:topN] <- "YES"
    # Label any marker genes in window of interest
    LS_window <- quantile(plot_df$LS, 0.75)
    corr_window <- 0.25
    pos_top_genes <- plot_df[plot_df$LS > LS_window & plot_df$Correlation > corr_window,]$symbol
    neg_top_genes <- plot_df[plot_df$LS > LS_window & -plot_df$Correlation > corr_window,]$symbol
    if(nrow(plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]) > 0){
      plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes) & plot_df$symbol %in% markerGenes,]$toLabel <- "YES"
    }
    plot_df[plot_df$LS < LS_window | abs(plot_df$Correlation) < corr_window,]$toLabel <- "NO"
    res_list[[motif_short]] <- pos_top_genes # Save regulatory targets
    # Save dataframe of results
    save_df <- plot_df[plot_df$symbol %in% c(pos_top_genes, neg_top_genes),c(1:5)]
    save_df <- save_df[order(save_df$Correlation, decreasing=TRUE),]
    saveRDS(save_df, paste0(regPlotDir, sprintf("/%s.regulatory_targets_%s.rds", strain, motif_short)))
    plot_df <- plot_df[order(plot_df$mLog10pval, decreasing=FALSE),]
    # Label motif as well
    plot_df$toLabel[which(plot_df$symbol == motif_short)] <- "YES"
    plot_df$symbol[which(plot_df$toLabel == "NO")] <- ""
    # Threshold pvalue for plotting
    maxPval <- 5
    plot_df$mLog10pval <- ifelse(plot_df$mLog10pval > maxPval, maxPval, plot_df$mLog10pval)
    #Plot results
    p <- (
      ggplot(plot_df, aes(x=Correlation, y=LS, color=mLog10pval)) 
      #+ geom_point(size = 2)
      + ggrastr::geom_point_rast(size=2)
      + ggrepel::geom_text_repel(
        data=plot_df[plot_df$toLabel=="YES",], aes(x=Correlation, y=LS, label=symbol), 
        #data = plot_df, aes(x=Correlation, y=LS, label=symbol), #(don't do this, or the file will still be huge...)
        size=4,
        point.padding=0, # additional pading around each point
        box.padding=0.5,
        min.segment.length=0, # draw all line segments
        max.overlaps=Inf, # draw all labels
        #nudge_x = 2,
        color="black"
      ) 
      + geom_vline(xintercept=0, lty="dashed") 
      + geom_vline(xintercept=corr_window, lty="dashed", color="red")
      + geom_vline(xintercept=-corr_window, lty="dashed", color="red")
      + geom_hline(yintercept=LS_window, lty="dashed", color="red")
      + theme_classic()
      + theme(panel.grid.major=element_blank(), 
              panel.grid.minor= element_blank(), 
              plot.margin = unit(c(0.25,1,0.25,1), "cm"), 
              aspect.ratio=1.0,
              #legend.position = "none", # Remove legend
              axis.text.x = element_text(angle=90, hjust=1))
      + ylab("Linkage Score") 
      + xlab("Motif Correlation to Gene") 
      + scale_color_gradientn(colors=wesanderson::wes_palette("Zissou1", type = "continuous"), limits=c(0, maxPval), name = "-log10(p-val)")
      + scale_y_continuous(expand = expansion(mult=c(0,0.05)))
      + scale_x_continuous(limits = c(-0.85, 0.955)) # Force plot limits
      + ggtitle(sprintf("%s putative targets", motif_short))
    )
    
    pdf(paste0(regPlotDir, sprintf("/%s.%s_LS.pdf", strain, motif_short)), width=8, height=6)
    print(p)
    dev.off()
    
    outfile = paste0(regPlotDir, sprintf("/%s.%s_LS.png", strain, motif_short))
    print(p)
    ggsave(outfile, width=678/96, height=420/96, dpi=300)
    
    print(c(motif, length(unique(plot_df[plot_df$toLabel=="YES",]$symbol)), 
            length(unique(plot_df[plot_df$toLabel=="YES" & plot_df$symbol %in% markerGenes,]$symbol))))
  }
  
}

