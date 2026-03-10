
library(Seurat)
library(miloR)
library(ggh4x)
library(ggbeeswarm)
library(ggtext)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt//utils/00.initial_setting.R")




################################################################################
### Run MiloR
################################################################################
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LV.RNA.anno.L2.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LV.RNA.anno.L2.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.MSA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MCA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MSA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LK.multiomics.anno.L2.rds", 
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.PBMC.RNA.anno.L2.rds"
)

cluster = "Cell_type_L1"

# Loop through each input file
for(i in input_file) {
  
  tissue <- sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
  
  seurat_object <- readRDS(i)
  DefaultAssay(seurat_object) <- "RNA"
  
  seurat_object@active.ident <- seurat_object@meta.data[, cluster]
  Idents(seurat_object) <- cluster
  
  # Clean strain names and extract metadata
  seurat_object$strain <- gsub(" ", "", seurat_object$strain)
  meta_table <- seurat_object@meta.data
  
  species_list <- unique(meta_table$strain)
  
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
### Load MiloR files and perform differential abundance testing
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR")


input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LK.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MSA.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MCA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.PBMC.miloR.rds"
)


da_merge <- data.frame()

# Loop through each MiloR file for differential abundance testing
for (i in input_file) {
  
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
      
      # positive LFCs = control > treatment
      model.contrasts <- paste0("treatment.new", control.new, " - ", "treatment.new", treatment.new)
      
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
### Load and Process Differential Abundance Data for Visualization (Figure 1f)
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







################################################################################
### Identify over-represented L2 cell types in each L1 neighborhoods
################################################################################

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/mouse.LK.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SD.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MSA.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SHR.MCA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.HYP.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LV.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.LK.miloR.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/WKY.MSA.miloR.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/SS.PBMC.miloR.rds"
)

da_merge = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/milo.da_result.out", sep='\t', header=T)

da_l2_enrich_list <- list()
da_merge_l2_enrich <- c()
for(i in input_file){
  
  QR_mat <- c()
  P_mat <- c()
  Q_global <- c()
  
  milo_obj <- readRDS(i)
  meta <- as.data.frame(colData(milo_obj))
  strain <- unique(meta$Strain)
  tissue <- unique(meta$Tissue)
  
  da_res <- da_merge[da_merge$strain == strain & da_merge$tissue == tissue, ]
  da_res$nhood <- as.character(da_res$nhood_id)
  
  M <- miloR::nhoods(milo_obj)  # cells x nhoods
  stopifnot(identical(rownames(M), rownames(meta)))
  
  L1_ind  <- model.matrix(~ 0 + meta$Cell_type_L1)
  tab_L1  <- as.matrix(Matrix::t(M) %*% L1_ind)
  colnames(tab_L1) <- sub("^meta\\$Cell_type_L1", "", colnames(tab_L1))
  L2_ind  <- model.matrix(~ 0 + meta$Cell_type_L2)
  tab_L2 <- as.matrix(Matrix::t(M) %*% L2_ind)
  colnames(tab_L2) <- sub("^meta\\$Cell_type_L2", "", colnames(tab_L2))
  nh_ids <- rownames(tab_L2)
  L2_levels <- colnames(tab_L2)
  
  nh_size <- Matrix::colSums(M)
  names(nh_size) <- nh_ids
  
  row_has_L1 <- rowSums(tab_L1) > 0
  row_has_L2 <- rowSums(tab_L2) > 0
  
  nh_major_L1 <- rep(NA_character_, nrow(tab_L1))
  nh_major_L1[row_has_L1] <- colnames(tab_L1)[max.col(tab_L1[row_has_L1, , drop = FALSE], ties.method = "first")]
  
  nh_major_L2 <- rep(NA_character_, nrow(tab_L2))
  nh_major_L2[row_has_L2] <- colnames(tab_L2)[max.col(tab_L2[row_has_L2, , drop = FALSE], ties.method = "first")]
  
  a_mat <- tab_L2
  b_mat <- nh_size[nh_ids] - a_mat
  tot_L2 <- colSums(L2_ind)
  Ncells <- nrow(M)
  
  c_mat <- matrix(rep(tot_L2, each = nrow(a_mat)), nrow = nrow(a_mat)) - a_mat
  d_mat <- (Ncells - nh_size[nh_ids]) - c_mat
  
  OR_mat <- (a_mat * d_mat + 0.5) / (b_mat * c_mat + 0.5)
  dimnames(OR_mat) <- list(nh_ids, L2_levels)
  
  P_mat <- matrix(NA_real_, nrow = nrow(a_mat), ncol = ncol(a_mat),
                  dimnames = list(nh_ids, L2_levels))
  
  for (j in seq_len(ncol(a_mat))) {
    K <- tot_L2[j]
    n <- nh_size[nh_ids]
    a <- a_mat[, j]
    # invalid tests: empty nhood, or K==0 (no such L2 anywhere) -> set p=1
    invalid <- (n == 0L) | (K == 0L)
    pj <- rep(1, length(a))
    ok <- !invalid
    if (any(ok)) {
      pj[ok] <- phyper(q = pmax(a[ok] - 1, 0),
                       m = K,
                       n = Ncells - K,
                       k = n[ok],
                       lower.tail = FALSE)
    }
    P_mat[, j] <- pj
  }
  
  Q_rowwise <- t(apply(P_mat, 1, function(p) p.adjust(replace(p, !is.finite(p), 1), method = "BH")))
  dimnames(Q_rowwise) <- dimnames(P_mat)
  
  Q_global <- matrix(p.adjust(replace(as.vector(P_mat), !is.finite(as.vector(P_mat)), 1), method = "BH"),
                     nrow = nrow(P_mat), ncol = ncol(P_mat),
                     dimnames = dimnames(P_mat))
  
  da_l2_enrich_list[i]['QR_mat'] <- QR_mat
  da_l2_enrich_list[i]['P_mat'] <- P_mat
  da_l2_enrich_list[i]['Q_global'] <- Q_global
  
  top_idx <- max.col(-Q_rowwise, ties.method = "first")
  L2_top  <- L2_levels[top_idx]
  L2_top[!row_has_L2] <- NA_character_
  
  top_tbl <- tibble(
    nhood    = nh_ids,
    L1_major = nh_major_L1, 
    L2_top   = L2_top,  
    L2_abund = nh_major_L2,               # most abundant (majority)
    OR       = OR_mat[cbind(seq_len(nrow(OR_mat)), top_idx)],
    p        = P_mat[cbind(seq_len(nrow(P_mat)),  top_idx)],
    q        = Q_rowwise[cbind(seq_len(nrow(Q_rowwise)), top_idx)]
  ) %>%
    arrange(q)
  
  da_res_l2_enrich <- da_res %>%
    left_join(top_tbl, by = "nhood")
  
  da_merge_l2_enrich <- rbind(da_merge_l2_enrich, da_res_l2_enrich)
  
}

write.table(da_merge_l2_enrich, "milo.da_result.l2_enrich.out", sep='\t', quote=F, col.names = T, row.names = F)
saveRDS(da_l2_enrich_list, "milo.da_result.l2_enrich.rds")











################################################################################
### Visualize over-represented L2 cell types in each L1 neighborhoods (Fig. S9b-e)
################################################################################
# function
l2_prop_vis <- function(l2_enrich_df, l1_target = NULL, l2_use = "L2_top", alpha = 0.1){
  l2_enrich_df %>%
    dplyr::mutate(L2_use = .data[[l2_use]],
                  L1_L2 = paste(Cell_type_L1, "-", L2_use)) %>%
    filter(L1_L2 %in% l1_l2_cor) %>%
    filter(Cell_type_L1 == l1_target) %>%
    dplyr::mutate(is_signif = ifelse(SpatialFDR < alpha, 1, 0)) %>%
    dplyr::mutate(logFC_color = ifelse(is_signif == 1, logFC, NA)) %>%
    # arrange(L2_use) %>% mutate(Nhood = factor(Nhood, levels = unique(Nhood))) %>%
    ggplot(aes(L2_use, -logFC, color = logFC_color)) + guides(color = "none") +
    xlab("") + ylab("Log Fold Change") + geom_quasirandom(alpha = 1) +
    coord_flip() + theme_bw(base_size = 22) + theme(strip.text.y = element_text(angle = 0)) + scale_x_discrete(limits=rev) +
    labs(
      x = "", y = "Log Fold Change",
      title = "Strain - treatment (vs. control)<br>
    <span style='font-size:11pt'>
    <span style='color:#832424;'>Enriched in treatment<br></span>
    <span style='color:#3A3A98;'>Depleted in treatment</span>
    </span>"
    ) +
    theme_bw() +
    theme(plot.title = element_markdown(lineheight = 1.1),
          legend.text = element_markdown(size = 11),
          text=element_text(size=11, color="black", family = "Arial"),
          axis.text.y = element_text(size=11, colour = 'black'),
          axis.text.x = element_text(size=9, colour = 'black'),
          axis.title = element_text(size=11, colour = 'black'),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(colour = "black", fill = NA)) +
    facet_nested(tissue ~ strain + treatment , scales = "free", space = "free") +
    geom_hline(yintercept = c(-5,  5), linetype = "solid", color = "grey") +
    geom_hline(yintercept = c(0), linetype = "dashed", color = "black") +
    stat_summary(fun = median, geom = "point", color = "black", size = 2) +
    scale_y_continuous(breaks = c(-10, -5, 0, 5, 10), oob = scales::oob_squish) +
    scale_color_gradient2( na.value = "lightgrey" )
}


meta_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/merged_metadata.tsv", sep = "\t", header = T)
l1_l2_cor <- unique(paste(meta_merged$Cell_type_L1, "-", meta_merged$Cell_type_L2))

da_merge_l2_enrich <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/miloR/milo.da_result.l2_enrich.out", sep = "\t", header = T)
da_merge_l2_enrich = da_merge_l2_enrich[!(da_merge_l2_enrich$tissue=="MCA" & da_merge_l2_enrich$strain %in% c("C57BL/6", "SS")), ]
da_merge_l2_enrich$strain = factor(da_merge_l2_enrich$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
da_merge_l2_enrich$tissue = factor(da_merge_l2_enrich$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
da_merge_l2_enrich$treatment = factor(da_merge_l2_enrich$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))

l2_prop_vis(da_merge_l2_enrich, "Excitatory neuron")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/Excitatory_neuron.subtype.milo.png", width=975/96, height=1105/96, dpi=300)

l2_prop_vis(da_merge_l2_enrich, "Inhibitory neuron")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/Inhibitory_neuron.subtype.milo.png", width=975/96, height=1105/96, dpi=300)

l2_prop_vis(da_merge_l2_enrich, "Astrocyte")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/Astrocyte.subtype.milo.png", width=975/96, height=261/96, dpi=300)

l2_prop_vis(da_merge_l2_enrich, "Myelinating OL")
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/Myelinating_OL.subtype.milo.png", width=975/96, height=261/96, dpi=300)






















