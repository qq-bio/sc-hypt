
library(Seurat)
library(miloR)
library(ggh4x)
library(ggbeeswarm)
library(ggtext)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")




################################################################################
input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

# Define parameters
cluster <- "subclass_level1" # Define main clustering identifier

# Loop through each input file
for(i in input_file) {
  
  # Extract tissue information from filename
  tissue <- sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
  
  # Load Seurat object and set default assay
  seurat_object <- readRDS(i)
  DefaultAssay(seurat_object) <- "RNA"
  
  seurat_object@active.ident <- seurat_object@meta.data[, cluster]
  Idents(seurat_object) <- cluster
  
  # Clean strain names and extract metadata
  seurat_object$strain <- gsub(" ", "", seurat_object$strain)
  meta_table <- seurat_object@meta.data
  
  # Define list of unique species (strains) present in the dataset
  species_list <- unique(meta_table$strain)
  
  # Loop through each species
  for (si in species_list) {
    
    # Subset Seurat object by strain
    seurat_object_use <- subset(seurat_object, strain == si)
    
    # Adjust strain label for mouse
    if (si == "C57BL/6") { si <- "mouse" }
    
    # Set output file path
    outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/",
                      si, ".", tissue, ".miloR.rds")
    
    # Set dimensionality reduction method based on tissue type
    reduced.dim <- ifelse(tissue == "LK", "HARMONY.RNA", "HARMONY")
    
    # Construct KNN graph and identify representative neighborhoods
    sce <- as.SingleCellExperiment(seurat_object_use, assay = "RNA")
    milo_object <- Milo(sce)
    milo_object <- buildGraph(milo_object, k = 30, d = 30, reduced.dim = reduced.dim)
    milo_object <- makeNhoods(milo_object, prop = 0.2, k = 30, d = 30, refined = TRUE, 
                              reduced_dims = reduced.dim, refinement_scheme = "graph")
    
    # Visualize neighborhood size histogram
    plotNhoodSizeHist(milo_object) # Ideal average size: 5 x N_samples / 50-100
    
    # Count cells within each neighborhood and create design matrix
    milo_object <- countCells(milo_object, meta.data = seurat_object_use@meta.data, sample = "orig.ident")
    
    # Create and format design matrix for Milo analysis
    milo_design <- distinct(seurat_object_use@meta.data[, c("orig.ident", "treatment")])
    rownames(milo_design) <- milo_design$orig.ident
    milo_design$orig.ident <- as.factor(milo_design$orig.ident)
    milo_design$treatment <- as.factor(milo_design$treatment)
    
    # Calculate neighborhood connectivity
    milo_object <- calcNhoodDistance(milo_object, d = 30, reduced.dim = reduced.dim)
    
    # Save the Milo object for downstream analysis
    saveRDS(milo_object, outfile)
  }
}









################################################################################
# Load MiloR files and perform differential abundance testing
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR")


input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MSA.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MCA.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.MSA.miloR.rds"
)


da_merge <- data.frame()

# Loop through each MiloR file for differential abundance testing
for (i in input_file) {
  
  # Extract tissue type from the filename and set reduction method
  tissue <- sub(".*\\.([^.]*)\\.(miloR).*", "\\1", basename(i))
  reduction <- ifelse(tissue == "LK", "HARMONY.RNA", "HARMONY")
  
  milo_object <- readRDS(i)
  
  milo_object <- countCells(milo_object, meta.data = as.data.frame(colData(milo_object)), sample = "orig.ident")
  milo_design <- distinct(data.frame(colData(milo_object))[, c("treatment", "strain", "orig.ident")])
  rownames(milo_design) <- milo_design$orig.ident
  milo_design$treatment.new <- as.factor(gsub(" ", "", milo_design$treatment))
  milo_design$strain <- as.factor(milo_design$strain)
  
  # Define species and treatment lists
  species_list <- unique(milo_design$strain)
  treatment_list <- intersect(levels(colData(milo_object)$treatment), unique(milo_design$treatment))
  treatment.new_list <- gsub(" ", "", treatment_list)
  
  # Loop through each species and perform differential abundance testing for each treatment
  for (si in species_list) {
    for (j in seq_along(treatment.new_list)[-1]) {
      
      # Define control and treatment groups for contrast
      control.new <- treatment.new_list[1]
      treatment.new <- treatment.new_list[j + 1]
      control <- treatment_list[1]
      treatment <- treatment_list[j + 1]
      
      # Define model contrast
      model.contrasts <- paste0("treatment.new", control.new, " - ", "treatment.new", treatment.new)
      
      # Run differential abundance test with error handling
      da_tmp <- try(testNhoods(milo_object, design = ~ 0 + treatment.new, design.df = milo_design,
                               fdr.weighting = "graph-overlap", reduced.dim = reduction,
                               model.contrasts = model.contrasts), silent = TRUE)
      
      if (class(da_tmp) != "try-error") {
        da_tmp <- annotateNhoods(milo_object, da_tmp, coldata_col = "subclass_level1")
        da_tmp$subcluster <- factor(da_tmp$subclass_level1, levels = cell_order)
        da_tmp$species <- si
        da_tmp$tissue <- tissue
        da_tmp$control <- control
        da_tmp$treatment <- treatment
        da_merge <- rbind(da_merge, da_tmp)
        
      } else {
        print(c(i, treatment))
      }
    }
  }
}

write.table(da_merge, "milo.da_result.out", sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)










################################################################################
# Load and Process Differential Abundance Data for Visualization (Figure 1f)
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/")

da_merge <- read.table("milo.da_result.out", sep='\t', header=TRUE)

da_merge$subcluster <- factor(da_merge$subcluster, levels = cell_order)
da_merge$species <- factor(da_merge$species, levels = strain_order)
da_merge$tissue <- factor(da_merge$tissue, levels = tissue_order)
da_merge$treatment <- factor(da_merge$treatment, levels = treatment_order)


# Create a plotting order based on significance level
da_merge_plot <- da_merge %>%
  mutate(order = case_when(
    SpatialFDR > 0.05 ~ 1,
    SpatialFDR <= 0.05 ~ 2
  )) %>%
  arrange(order)


# Define grouping and significance level
group.by <- "subcluster"
alpha <- 0.1

# Create plot data, marking significant points and coloring by log Fold Change
da.res <- da_merge_plot %>%
  mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0),
         logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
  arrange(get(group.by)) %>%
  mutate(Nhood = factor(Nhood, levels = unique(Nhood)))


ggplot(da.res, aes_string(x = group.by, y = "-logFC", color = "logFC_color")) +
  guides(color = "none") +
  xlab(group.by) + ylab("Log Fold Change") +
  geom_quasirandom(alpha = 1) +
  coord_flip() +
  scale_x_discrete(limits = rev) +
  theme_bw(base_size = 22) +
  theme(
    strip.text.y = element_text(angle = 0),
    plot.title = element_markdown(lineheight = 1.1),
    legend.text = element_markdown(size = 11),
    text = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, colour = 'black'),
    axis.text.x = element_text(size = 9, colour = 'black'),
    axis.title = element_text(size = 11, colour = 'black'),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(
    x = "", y = "Log Fold Change",
    title = "Strain - treatment (vs. control)<br>
    <span style='font-size:11pt'>
    <span style='color:#832424;'>Enriched in treatment<br></span>
    <span style='color:#3A3A98;'>Depleted in treatment</span>
    </span>"
  ) +
  facet_nested(tissue ~ species + treatment, scales = "free", space = "free") +
  geom_hline(yintercept = c(-5, 5), linetype = "solid", color = "grey") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun = median, geom = "point", color = "black", size = 2) +
  scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
  scale_color_gradient2(na.value = "lightgrey")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1f.milo.png", width = 731/96, height = 957/96, dpi = 300)

