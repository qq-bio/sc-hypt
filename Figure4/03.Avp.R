library(Seurat)
library(dplyr)
library(ggplot2)


################################################################################
### Avp+ neuron annotation (Fig. S16a-b)
################################################################################
# load files with subclustering results of neurons
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.cluster.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.cluster.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.cluster.rds')

neuron_list = c("Inhibitory neurons", "Excitatory neurons", "Avp+ neurons")
hyp_m %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_m$project))
hyp_ss %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_ss$project))
hyp_shr %>% subset(new.cluster.ids_umap %in% neuron_list) %>% DotPlot(., features = "Avp", group.by = "subcluster") + labs(title = unique(hyp_shr$project))

avp_m <- rownames(hyp_m@meta.data[hyp_m$new.cluster.ids_umap=="Avp+ neurons", ])
avp_ss_gaba <- rownames(hyp_ss@meta.data[hyp_ss$subcluster=="GABA-10", ])
avp_ss_glu <- rownames(hyp_ss@meta.data[hyp_ss$subcluster=="Glu-15", ])
avp_shr_gaba <- rownames(hyp_shr@meta.data[hyp_shr$subcluster=="GABA-5", ])
avp_shr_glu <- rownames(hyp_shr@meta.data[hyp_shr$subcluster=="Glu-15", ])

avp_neuron_df <- data.frame(
  cell_id = c(avp_m, avp_ss_gaba, avp_ss_glu, avp_shr_gaba, avp_shr_glu),
  Avp_anno = c(rep("Avp+ neurons", length(avp_m)),
               rep("Avp+ neurons, GABAergic", length(avp_ss_gaba)),
               rep("Avp+ neurons, glutamatergic", length(avp_ss_glu)),
               rep("Avp+ neurons, GABAergic", length(avp_shr_gaba)),
               rep("Avp+ neurons, glutamatergic", length(avp_shr_glu))
  )
)

write.table(avp_neuron_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/avp_neuron.anno.out", col.names = T, row.names = F, sep = "\t")


# map Avp neuron flag to seurat object metatable
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds')

avp_neuron_df <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/avp_neuron.anno.out", header = T)
rownames(avp_neuron_df) <- avp_neuron_df$cell_id

avp_cells <- intersect(colnames(hyp_m), avp_neuron_df$cell_id)
hyp_m@meta.data$flag_Avp_neuron <- FALSE
hyp_m@meta.data[avp_cells, ]$flag_Avp_neuron <- avp_neuron_df[avp_cells, ]$Avp_anno

avp_cells <- intersect(colnames(hyp_ss), avp_neuron_df$cell_id)
hyp_ss@meta.data$flag_Avp_neuron <- FALSE
hyp_ss@meta.data[as.character(avp_cells), ]$flag_Avp_neuron <- avp_neuron_df[avp_cells, ]$Avp_anno

avp_cells <- intersect(colnames(hyp_shr), avp_neuron_df$cell_id)
hyp_shr@meta.data$flag_Avp_neuron <- FALSE
hyp_shr@meta.data[as.character(avp_cells), ]$flag_Avp_neuron <- avp_neuron_df[avp_cells, ]$Avp_anno


hyp_m %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_m$Project))
hyp_ss %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_ss$Project))
hyp_shr %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_shr$Project))

hyp_m$SxT = factor(hyp_m$SxT, levels = sxt_order)
hyp_ss$SxT = factor(hyp_ss$SxT, levels = sxt_order)
hyp_shr$SxT = factor(hyp_shr$SxT, levels = sxt_order)

saveRDS(hyp_m, '/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds')
saveRDS(hyp_ss, '/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds')
saveRDS(hyp_shr, '/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds')



# visualize Avp expression 
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds')
hyp_shr = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds')


p1 = hyp_m %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_m$Project), y = "Avp neuron label (yes/no)", x = "")
p2 = hyp_ss %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_ss$Project), y = "Avp neuron label (yes/no)", x = "")
p3 = hyp_shr %>% DotPlot(., features = "Avp", group.by = "flag_Avp_neuron") + labs(title = unique(hyp_shr$Project), y = "Avp neuron label (yes/no)", x = "")
p1 + p2 + p3 & theme(text = element_text(family="Arial"))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.avp.dotplot.png", width=1393/96, height=335/96, dpi=300)


subset_data = hyp_m %>% subset(flag_Avp_neuron == "Avp+ neurons")
p1 = VlnPlot(subset_data, features = "Avp", group.by = "SxT", cols = sxt_col[unique(subset_data$SxT)]) + theme(legend.position = 'None') + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons\n(AngII)")
print(p1)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.angii.avp.expr.png", width=356/96, height=307/96, dpi=300)

subset_data = hyp_ss %>% subset(flag_Avp_neuron == "Avp+ neurons, glutamatergic")
p2 = VlnPlot(subset_data, features = "Avp", group.by = "SxT", cols = sxt_col[unique(subset_data$SxT)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p2)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.ss.avp.expr-1.png", width=420/96, height=270/96, dpi=300)
subset_data = hyp_ss %>% subset(flag_Avp_neuron == "Avp+ neurons, GABAergic")
p3 = VlnPlot(subset_data, features = "Avp", group.by = "SxT", cols = sxt_col[unique(subset_data$SxT)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") + theme(legend.position = 'None')
print(p3)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.ss.avp.expr-2.png", width=420/96, height=270/96, dpi=300)

subset_data = hyp_shr %>% subset(flag_Avp_neuron == "Avp+ neurons, glutamatergic")
p4 = VlnPlot(subset_data, features = "Avp", group.by = "SxT", cols = sxt_col[unique(subset_data$SxT)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, glutamatergic\n(Spontaneous)") + theme(legend.position = 'None')
print(p4)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.shr.avp.expr-1.png", width=420/96, height=270/96, dpi=300)

subset_data = hyp_shr %>% subset(flag_Avp_neuron == "Avp+ neurons, GABAergic")
p5 = VlnPlot(subset_data, features = "Avp", group.by = "SxT", cols = sxt_col[unique(subset_data$SxT)]) + 
  labs(x = "", y = "Avp\nexpression", title = "Avp+ neurons, GABAergic\n(Spontaneous)") + theme(legend.position = 'None')
print(p5)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.shr.avp.expr-2.png", width=420/96, height=270/96, dpi=300)




################################################################################
### Visualize SCENIC results for Ep300 activity in Avp+ neurons (Fig. S16c)
################################################################################
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds')
hyp_ss = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds')

m_scenic <- fread("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/mhyp_all_regulon_AUC.csv")
ss_scenic <- fread("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_ssall_regulon_AUC.csv")
sd_scenic <- fread("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scenic_mhyp/rhyp_sdall_regulon_AUC.csv")

colnames(ss_scenic) <- gsub("\\.","-", colnames(ss_scenic))
colnames(sd_scenic) <- gsub("\\.","-", colnames(sd_scenic))

length(intersect(colnames(ss_scenic), hyp_ss$Cell_ID))

m_avp_cell <- hyp_m@meta.data %>% filter(flag_Avp_neuron!=FALSE)
ss_avp_cell <- hyp_ss@meta.data %>% filter(flag_Avp_neuron!=FALSE)

m_avp_scenic <- m_scenic %>%
  tidyr::pivot_longer(
    cols = -Regulon,
    names_to = "Cell_ID",
    values_to = "AUC"
  ) %>%
  filter(grepl("Ep300", Regulon)) %>%
  inner_join(m_avp_cell, by = "Cell_ID")

ss_avp_scenic <- ss_scenic %>%
  tidyr::pivot_longer(
    cols = -Regulon,
    names_to = "Cell_ID",
    values_to = "AUC"
  ) %>%
  filter(grepl("Ep300", Regulon)) %>%
  inner_join(ss_avp_cell, by = "Cell_ID")

sd_avp_scenic <- sd_scenic %>%
  tidyr::pivot_longer(
    cols = -Regulon,
    names_to = "Cell_ID",
    values_to = "AUC"
  ) %>%
  filter(grepl("Ep300", Regulon)) %>%
  inner_join(ss_avp_cell, by = "Cell_ID")


p1 <- m_avp_scenic %>%
  mutate(Treatment = factor(Treatment, levels = treatment_order)) %>%
  ggplot(aes(x = Treatment, y = AUC)) +
  geom_violin(aes(fill = Strain)) +
  geom_jitter(color = "grey", size = 0.5) +
  scale_fill_manual(values = strain_col) +
  labs(x = "", y = "Ep300(+)\nactivity", title = "Avp+ neurons\n(AngII)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "black", fill = NA)
  )

p1
ggsave(paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/angii.avp_ep300.png"), 
       dpi = 300,        
       width = 505/96,  
       height = 250/96)


p2 <- rbind(sd_avp_scenic, ss_avp_scenic) %>%
  filter(flag_Avp_neuron=="Avp+ neurons, glutamatergic") %>%
  ggplot(aes(x = SxT, y = AUC)) +
  geom_violin(aes(fill = Strain)) +
  geom_jitter(color = "grey", size = 0.5) +
  scale_fill_manual(values = strain_col) +
  labs(x = "", y = "Ep300(+)\nactivity", title = "Avp+ neurons, glutamatergic\n(Salt-sensitive)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "black", fill = NA)
  )

p3 <- rbind(sd_avp_scenic, ss_avp_scenic) %>%
  filter(flag_Avp_neuron=="Avp+ neurons, GABAergic") %>%
  ggplot(aes(x = SxT, y = AUC)) +
  geom_violin(aes(fill = Strain)) +
  geom_jitter(color = "grey", size = 0.5) +
  scale_fill_manual(values = strain_col) +
  labs(x = "", y = "Ep300(+)\nactivity", title = "Avp+ neurons, GABAergic\n(Salt-sensitive)") +
  theme_classic() +
  theme(
    text = element_text(family = "Arial"),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    strip.text = element_text(colour = "black"),
    strip.background = element_rect(colour = "black", fill = NA)
  )

plot_grid(p2, p3, ncol = 1, align = "v")
ggsave(paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/ss.avp_ep300.png"), 
       dpi = 300,        
       width = 505/96,  
       height = 403/96)







################################################################################
### ATAC-seq deviation matrix
################################################################################
subset_archr_by_seurat_celltype_add_meta <- function(
    seurat_obj,
    archr_obj,
    cell_type_col,
    cell_types,
    seurat_cell_col = NULL,       # optional: meta column in Seurat storing ArchR cell IDs
    seurat_to_archr = TRUE,       # convert "_" -> "#"
    outputDirectory = NULL,      
    verbose = TRUE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(inherits(archr_obj, "ArchRProject"))
  
  md <- seurat_obj@meta.data
  
  if (!cell_type_col %in% colnames(md)) {
    stop("cell_type_col not found in seurat_obj@meta.data: ", cell_type_col)
  }
  if (!is.null(seurat_cell_col) && !seurat_cell_col %in% colnames(md)) {
    stop("seurat_cell_col not found in seurat_obj@meta.data: ", seurat_cell_col)
  }
  
  if (is.null(outputDirectory) || !nzchar(outputDirectory)) {
    stop("save_subset=TRUE but outputDirectory is NULL/empty.")
  }
  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
  
  # 1) pick Seurat cells of interest
  keep_seu_cells <- rownames(md)[md[[cell_type_col]] %in% cell_types]
  if (length(keep_seu_cells) == 0) stop("No Seurat cells found for requested cell_types.")
  if (verbose) message("Seurat cells selected: ", length(keep_seu_cells))
  
  # 2) decide ArchR cell names corresponding to these Seurat cells
  if (!is.null(seurat_cell_col)) {
    archr_ids_from_seu <- as.character(md[keep_seu_cells, seurat_cell_col])
    archr_ids_from_seu <- archr_ids_from_seu[!is.na(archr_ids_from_seu)]
  } else {
    archr_ids_from_seu <- keep_seu_cells
  }
  archr_ids_from_seu <- unique(archr_ids_from_seu)
  
  if (seurat_to_archr) {
    archr_ids_from_seu <- gsub("_", "#", archr_ids_from_seu)
  }
  
  # 3) intersect with ArchR project
  archr_cellnames <- archr_obj$cellNames
  common_cells <- intersect(archr_ids_from_seu, archr_cellnames)
  
  if (length(common_cells) == 0) {
    stop(
      "No overlapping cells between Seurat and ArchR after name normalization.\n",
      "Tip: verify cell IDs or provide seurat_cell_col."
    )
  }
  if (verbose) {
    message(
      "Matched ArchR cells: ", length(common_cells),
      " (", round(100 * length(common_cells) / length(archr_ids_from_seu), 1), "%)"
    )
  }
  
  # 4) subset ArchR project
  archr_sub <- ArchR::subsetArchRProject(
    ArchRProj = archr_obj,
    cells = common_cells,
    dropCells = TRUE,
    force = TRUE,
    outputDirectory = outputDirectory
  )
  
  return(archr_sub)
}





################################################################################
### Subset archr object to reduce computing time
################################################################################
# functions
subset_archr_by_seurat_celltype_add_meta <- function(
    seurat_obj,
    archr_obj,
    cell_type_col,
    cell_types,
    seurat_cell_col = NULL,       # optional: meta column in Seurat storing ArchR cell IDs
    seurat_to_archr = TRUE,       # convert "_" -> "#"
    outputDirectory = NULL,      
    verbose = TRUE
) {
  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(inherits(archr_obj, "ArchRProject"))
  
  md <- seurat_obj@meta.data
  
  if (!cell_type_col %in% colnames(md)) {
    stop("cell_type_col not found in seurat_obj@meta.data: ", cell_type_col)
  }
  if (!is.null(seurat_cell_col) && !seurat_cell_col %in% colnames(md)) {
    stop("seurat_cell_col not found in seurat_obj@meta.data: ", seurat_cell_col)
  }
  
  if (is.null(outputDirectory) || !nzchar(outputDirectory)) {
    stop("save_subset=TRUE but outputDirectory is NULL/empty.")
  }
  dir.create(outputDirectory, recursive = TRUE, showWarnings = FALSE)
  
  # 1) pick Seurat cells of interest
  keep_seu_cells <- rownames(md)[md[[cell_type_col]] %in% cell_types]
  if (length(keep_seu_cells) == 0) stop("No Seurat cells found for requested cell_types.")
  if (verbose) message("Seurat cells selected: ", length(keep_seu_cells))
  
  # 2) decide ArchR cell names corresponding to these Seurat cells
  if (!is.null(seurat_cell_col)) {
    archr_ids_from_seu <- as.character(md[keep_seu_cells, seurat_cell_col])
    archr_ids_from_seu <- archr_ids_from_seu[!is.na(archr_ids_from_seu)]
  } else {
    archr_ids_from_seu <- keep_seu_cells
  }
  archr_ids_from_seu <- unique(archr_ids_from_seu)
  
  if (seurat_to_archr) {
    archr_ids_from_seu <- gsub("_", "#", archr_ids_from_seu)
  }
  
  # 3) intersect with ArchR project
  archr_cellnames <- archr_obj$cellNames
  common_cells <- intersect(archr_ids_from_seu, archr_cellnames)
  
  if (length(common_cells) == 0) {
    stop(
      "No overlapping cells between Seurat and ArchR after name normalization.\n",
      "Tip: verify cell IDs or provide seurat_cell_col."
    )
  }
  if (verbose) {
    message(
      "Matched ArchR cells: ", length(common_cells),
      " (", round(100 * length(common_cells) / length(archr_ids_from_seu), 1), "%)"
    )
  }
  
  # 4) subset ArchR project
  archr_sub <- ArchR::subsetArchRProject(
    ArchRProj = archr_obj,
    cells = common_cells,
    dropCells = TRUE,
    force = TRUE,
    outputDirectory = outputDirectory
  )
  
  return(archr_sub)
}


LK_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.rds')
LK_m_archr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse.LK.archr.rds")

LK_m_avpr_archr <- subset_archr_by_seurat_celltype_add_meta(
  seurat_obj = LK_m,
  archr_obj = LK_m_archr,
  cell_type_col = "Cell_type_L1",
  cell_types = c("IC", "Fibroblast", "CT", "CD"),
  outputDirectory = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse_avpr_subset/"
)

saveRDS(LK_m_avpr_archr, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse_avpr_subset/mouse.LK.avpr.archr.rds")





################################################################################
### Wilcoxon test on ATAC-seq deviation matrix in Avp receptor gene-expressing cells
################################################################################
input_archr <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse_avpr_subset/mouse.LK.avpr.archr.rds"
proj = readRDS(input_archr)

proj$cell_grp = paste0(proj$Cell_type_L1, "_", proj$Strain, "_", proj$Treatment)
cell_grp_count = table(proj$cell_grp)
cell_list = unique(proj$Cell_type_L1)
strain_list = unique(proj$Strain)
treatment_list = intersect(treatment_order, unique(proj$Treatment))

marker_merge = c()
for(ci in cell_list){
  
  for(si in strain_list){
    
    for(ti in treatment_list[-1]){
      
      bgdGroups = paste0(ci, "_", si, "_", treatment_list[1])
      useGroups = paste0(ci, "_", si, "_", ti)
      
      bgd_count <- ifelse(bgdGroups %in% names(cell_grp_count), cell_grp_count[bgdGroups], 0)
      use_count <- ifelse(useGroups %in% names(cell_grp_count), cell_grp_count[useGroups], 0)
      
      if(bgd_count > 10 & use_count > 10){
        
        markerTest <- getMarkerFeatures(
          ArchRProj = proj,
          useMatrix = "MotifMatrix",
          groupBy = "cell_grp",
          testMethod = "wilcoxon",
          useGroups = useGroups,
          bgdGroups = bgdGroups,
          # bias = c("TSSEnrichment", "log10(nFrags)"),
          threads = 1
        )
        
        markerList <- getMarkers(markerTest, cutOff = "FDR <= 1")
        markerList = markerList[[1]]
        
        if(nrow(markerList)>0){
          
          markerList$cell = ci
          markerList$species = si
          markerList$control = treatment_list[1]
          markerList$treatment = ti
          
          marker_merge = rbind(marker_merge, markerList)
          
        }
      }
      
    }
    
  }
  
}

DEM_outfile = gsub("rds", "DEM.out", input_archr)
write.table(marker_merge, DEM_outfile, col.names = T, row.names = T, sep='\t', quote = F)





################################################################################
### Summarize differential deviation score results (Fig. S16e)
################################################################################
input_archr <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/ArchR/mouse_avpr_subset/mouse.LK.avpr.archr.rds"
DEM_outfile = gsub("rds", "DEM.out", input_archr)
DEM = read.table(DEM_outfile, header=T, sep='\t')
DEM$treatment = factor(DEM$treatment, 
                       levels = treatment_order)

camp_cre <- c("Creb1", "Atf1", "Crem")
ap1_ie <- c("Jund", "Jun", "Junb", "Junc", "Fos", "Fosl1", "Fosl2", "Fosb")
nfat_calcineurin <- c("Nfatc1", "Nfatc2", "Nfatc3", "Nfatc4")
cebp <- c("Cebpa", "Cebpb", "Cebpd", "Cebpe", "Cebpg")
cd_identity <- c("Elf2", "Elf3", "Elf4", "Elf5", "Hnf1b", "Gata3", "Tfap2a", "Tfap2b")
hormonal_modulators <- c("Nr3c1", "Pparg", "Esr1")

ets_mapk <- c("Elk1", "Ets1", "Ets2", "Etv1", "Etv4", "Etv5")
nuclear_modulators <- c("Nr3c1", "Pparg", "Esr1")


tf_all_avpr <- list(
  `cAMP-CREB` = camp_cre,
  `AP-1 immediate early` = ap1_ie,
  `NFAT-calcineurin` = nfat_calcineurin,
  `C/EBP` = cebp,
  `Renal epithelial lineage` = cd_identity,
  `ETS/MAPK` = ets_mapk,
  nuclear_modulators = nuclear_modulators
)
tf_avpr_vec <- setNames(
  rep(names(tf_all_avpr), lengths(tf_all_avpr)),
  unlist(tf_all_avpr)
)

df_plot <- DEM %>%
  as.data.frame() %>%
  dplyr::filter(cell %in% c("CD", "Fibroblast", "IC", "CT")) %>%
  dplyr::group_by(treatment, cell) %>%
  dplyr::arrange(MeanDiff, .by_group = TRUE) %>%
  dplyr::mutate(
    rank = dplyr::row_number(),
    sig  = ifelse(FDR < 0.05, "Yes", "No"),
    cat  = dplyr::coalesce(tf_avpr_vec[gene_name], "Others"),
    receptor = ifelse(cell %in% c("CD", "CT"), "Avpr2-expressing", "Avpr1a-expressing")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    cat = factor(cat, levels = c("Others", setdiff(sort(unique(cat)), "Others")))
  ) %>%
  dplyr::arrange(cat != "Others")  # Others first => bottom layer

df_top1 <- df_plot %>%
  dplyr::filter(cat != "Others", sig == "Yes") %>%   
  dplyr::group_by(treatment, cell, receptor, cat) %>%
  dplyr::slice_max(order_by = MeanDiff, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

cat_levels <- levels(df_plot$cat)
cat_non_others <- setdiff(cat_levels, "Others")
legend_breaks <- c(setdiff(cat_levels, "Others"), "Others")

bmj_pal_fun <- scale_color_bmj()$palette
bmj_cols <- bmj_pal_fun(length(cat_non_others))

color_map <- c(
  "Others" = "grey70",
  setNames(bmj_cols, cat_non_others)
)

ggplot(df_plot, aes(rank, MeanDiff)) +
  geom_point(aes(size = sig, color = cat)) +
  
  ggrepel::geom_text_repel(
    data = df_top1,
    aes(label = gene_name),
    segment.color = "black",
    segment.size  = 0.3,
    box.padding   = 0.6,
    point.padding = 0.6,
    force         = 2,
    min.segment.length = 0,
    size = 3.5,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  
  theme_classic() +
  theme(text = element_text(family = "Arial"), size = 12) +
  scale_color_manual(values = color_map, breaks = legend_breaks) +
  scale_size_manual(values = c("No" = 1, "Yes" = 3)) +
  labs(color = "TF category", size = "Significant", 
       x = "Rank", y = "Motif deviation difference") +
  facet_nested(
    receptor + cell ~ treatment,
    scales = "free",
    independent = "x"
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/avp_tf.motif_dev.w_CT.mouse.dotplot.png", width=518/96, height=558/96, dpi = 300)






################################################################################
### UMAP visualization of Creb1/Crem/Alf1 motif score (Fig. S16f)
################################################################################
motifs <- c("Creb1", "Atf1", "Crem")
markerMotifs <- getFeatures(proj, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")

markerMotifs <- grep("z:", markerMotifs, value = TRUE)

p <- plotEmbedding(
  ArchRProj = proj, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "wnn.umap.harmony",
  imputeWeights = getImputeWeights(proj)
)

p2 <- lapply(p, function(x){
  x + # guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 12) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))






