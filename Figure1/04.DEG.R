library(Seurat)
library(dplyr)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src_pub/utils/00.initial_setting.R")




################################################################################
# Differential Gene Expression (DGE) Analysis - Treatment vs Control
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file <- c(
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
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

cluster <- "subclass_level1"

# Loop through each input file to perform DGE analysis
for (i in input_file) {
  
  # Initialize container for merged differential expression results
  deg_merged <- c()
  
  # Define output file path
  outfile <- paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                    gsub("anno.rds", "DEG_all.out", basename(i)))
  
  # Load the Seurat object
  seurat_object <- readRDS(i)
  
  # Extract dataset and tissue information
  dataset <- gsub("\\.[RNA|multiomics|EC]+.anno.rds", "", basename(i), perl = TRUE)
  tissue <- unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  Idents(seurat_object) <- cluster
  
  meta_table <- seurat_object@meta.data
  
  # Loop through each project
  project_list <- unique(meta_table$project)
  for (pi in project_list) {
    
    # Loop through each strain within the project
    strain_list <- unique(meta_table[meta_table$project == pi, ]$strain)
    for (si in strain_list) {
      
      # Loop through each cell type in the cluster
      cell_list <- unique(meta_table[meta_table$project == pi & meta_table$strain == si, cluster])
      
      # Determine control and treatment groups
      treatment <- meta_table[meta_table$project == pi & meta_table$strain == si, ]$treatment
      treatment <- intersect(levels(treatment), unique(treatment))
      control <- treatment[1]
      
      for (cell in cell_list) {
        
        # Select control cells
        cell.1 <- rownames(meta_table[meta_table$project == pi &
                                        meta_table[, cluster] == cell &
                                        meta_table$treatment == control &
                                        meta_table$strain == si, ])
        
        # Only proceed if enough control cells are available
        if (length(cell.1) > 3) {
          
          for (ti in treatment[-1]) {
            
            cell.2 <- rownames(meta_table[meta_table$project == pi &
                                            meta_table[, cluster] == cell &
                                            meta_table$treatment == ti &
                                            meta_table$strain == si, ])
            
            # Only proceed if enough treatment cells are available
            if (length(cell.2) > 3) {
              
              deg <- FindMarkers(seurat_object, ident.1 = cell.1, ident.2 = cell.2, logfc.threshold = 0)
              
              # Add metadata to DEG results
              deg$gene_name <- rownames(deg)
              deg$cell_type <- cell
              deg$pct.diff <- deg$pct.1 - deg$pct.2
              deg$project <- pi
              deg$strain <- si
              deg$tissue <- tissue
              deg$control <- control
              deg$treatment <- ti
              deg$control_size <- length(cell.1)
              deg$treatment_size <- length(cell.2)
              deg$test_gene <- nrow(deg)
              
              # Merge DEG results
              deg_merged <- rbind(deg_merged, deg)
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  write.table(deg_merged, outfile)
  
}








################################################################################
# Merge DEG Results
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file <- list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", pattern = "*DEG_all*")

deg_merged <- c()

# Loop through each input file to read and merge the data
for (i in input_file) {
  
  dat_tmp <- try(read.table(i, header = TRUE), silent = TRUE)
  
  if (class(dat_tmp) != "try-error") {
    deg_merged <- rbind(deg_merged, dat_tmp)
  }
  
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", 
            sep = '\t', col.names = TRUE, row.names = FALSE)







################################################################################
### Visualize DEG Results (Figure 1g)
################################################################################
library(lemon)
library(ggplot2)
library(dplyr)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep='\t', header=TRUE)

DEG_df <- deg_merged %>% 
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5) %>% 
  group_by(project, strain, treatment, tissue, cell_type) %>%
  select(project, strain, treatment, tissue, cell_type, control_size, treatment_size, test_gene) %>%
  mutate(
    DEG_num = n(),
    mean_cell_size = (mean(control_size, na.rm = TRUE) + mean(treatment_size, na.rm = TRUE)) / 2,
    min_cell_size = pmin(control_size, treatment_size),
    max_cell_size = pmax(control_size, treatment_size),
    cell_size_diff = abs(control_size - treatment_size),
    combined = factor(interaction(strain, treatment, sep = "-"), 
                      levels = c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d",
                                 "SS-HS 3d", "SS-HS 21d", "SD-HS 3d",
                                 "SHR-26w", "WKY-26w"))
  ) %>%
  unique() %>%
  as.data.frame()

DEG_df$cell_type <- factor(DEG_df$cell_type, levels = cell_order)
DEG_df$tissue <- factor(DEG_df$tissue, levels = tissue_order)

# Define color and alpha mappings for strains and treatments
treatment_alpha <- setNames(c(0.5, 1, 1, 0.5, 1), unique(DEG_df$treatment))
combined <- unique(DEG_df$combined)
combined_colors <- setNames(species_col[gsub("\\-.*", "", combined)], combined)
combined_alphas <- setNames(treatment_alpha[gsub(".*\\-", "", combined)], combined)

ggplot(DEG_df, aes(x = DEG_num, y = cell_type)) +
  geom_segment(aes(yend = cell_type), xend = 0, colour = "grey50") +
  geom_point(aes(colour = combined, alpha = combined, size = log10(max_cell_size) * 3)) +
  scale_y_discrete(limits = rev) +
  scale_colour_manual(values = combined_colors) +
  scale_alpha_manual(values = combined_alphas) +
  scale_size_continuous(breaks = c(6, 9, 12), labels = c("100", "1,000", "10,000")) +
  scale_x_continuous(breaks = seq(0, 500, 100)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Number of DEGs", y = "", color = "Strain x Treatment", alpha = "Strain x Treatment", size = "Number of Cells") +
  facet_grid2(cols = vars(project), rows = vars(tissue), scales = "free", space = "free_y") +
  coord_capped_cart(xlim = c(0, 500))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1g.deg_count.png", 
       width = 577/96, height = 913/96, dpi = 300)








################################################################################
### Visualize Overlap of DEG with Known BP Gene Lists (Figure 1h)
################################################################################
library(plotly)
library(htmlwidgets)
library(dplyr)

options(viewer = NULL)


# Load data and filter for unique gene mappings to human
m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = TRUE, sep = "\t") %>%
  filter(V4 == 1) %>%
  select(Gene.name, Human.gene.name) %>%
  unique()

r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", header = TRUE, sep = "\t") %>%
  filter(V5 == 1) %>%
  select(Gene.name, Human.gene.name) %>%
  unique()


# Load DEG data and known BP gene lists
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = TRUE, sep = '\t')
yong_list <- read.table("/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.mod.txt", header = TRUE, sep = '\t')$Gene.name %>% unique()
helen_list <- read.table("/xdisk/mliang1/qqiu/data/gene_list/NG_2024_Helen_bp_pred_gene_list.txt", header = TRUE, sep = '\t')$Gene %>% unique()

# Convert Helen's list to mouse and rat genes
bp_rat <- r2h %>% filter(Human.gene.name %in% helen_list) %>% pull(Gene.name) %>% unique()
bp_mouse <- m2h %>% filter(Human.gene.name %in% helen_list) %>% pull(Gene.name) %>% unique()
helen_list <- unique(c(bp_rat, bp_mouse))

deg_list <- deg_merged %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5, strain %in% c("C57BL/6", "SS", "SHR")) %>%
  pull(gene_name) %>%
  unique()

generate_pie_chart <- function(gene_list, deg_list, color, filename, rotation = 180) {
  ref_df <- data.frame(Gene = gene_list[gene_list %in% deg_merged$gene_name]) %>%
    mutate(Status = ifelse(Gene %in% deg_list, "DEG", "Non-DEG")) %>%
    count(Status) %>%
    mutate(Proportion = n / sum(n))
  
  fig <- plot_ly(ref_df, labels = ~Status, values = ~Proportion, 
                 marker = list(colors = color), type = 'pie', pull = c(0, 0.2), rotation = rotation) %>%
    layout(xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
           legend = list(orientation = 'h')) %>%
    htmlwidgets::onRender(sprintf("function(el, x) {
          var gd = document.getElementById(el.id); 
          Plotly.downloadImage(gd, {format: 'png', width: 210, height: 210, filename: '%s', scale: 5});
      }", filename))
  fig
}

color_scheme <- c("#CD534CFF", "lightgrey")
fig1 <- generate_pie_chart(yong_list, deg_list, color_scheme, "deg.pie-chart-yong")
fig2 <- generate_pie_chart(helen_list, deg_list, color_scheme, "deg.pie-chart-helen", rotation = 90)

fig1
fig2








################################################################################
### Visualize Shared DEGs (Figure 1i)
################################################################################
library(ggplot2)
library(dplyr)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", header = TRUE, sep = '\t')

deg_filter <- deg_merged %>%
  filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.5, strain %in% c("C57BL/6", "SS", "SHR"))

deg_tissue_count <- table(unique(deg_filter[, c("gene_name", "tissue")])$gene_name)
deg_strain_count <- table(unique(deg_filter[, c("gene_name", "strain")])$gene_name)

plot_shared_DEGs <- function(deg_count, x_label, shared_label) {
  data <- as.data.frame(table(deg_count))
  colnames(data) <- c("Gene", "Count")
  data <- data[order(-data$Count), ]
  shared_genes_count <- sum(deg_count > 1)
  
  mid_x <- length(unique(data$Gene)) / 2 + 1
  xmax <- length(unique(data$Gene))
  
  ggplot(data, aes(x = factor(Gene, levels = Gene), y = Count)) +
    geom_bar(stat = "identity", fill = "#F8766D") +
    geom_text(aes(label = Count), vjust = -0.1, color = "black", size = 4) +
    labs(x = x_label, y = "Number of DEGs") +
    theme(axis.text = element_text(color = "black")) +
    annotate("text", x = mid_x, y = max(data$Count) * 0.9,
             label = paste0(shared_label, "\n(n = ", format(shared_genes_count, big.mark = ","), ")"), size = 4) +
    annotate("segment", x = 2, xend = xmax, y = max(data$Count) * 0.9,
             yend = max(data$Count) * 0.9, colour = "black", size = 0.5)
}

p1 <- plot_shared_DEGs(deg_tissue_count, "Number of Shared Tissues", "Tissue shared")
p2 <- plot_shared_DEGs(deg_strain_count, "Number of Shared Strains", "Strain shared")

p1 + p2
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig1j.deg.shared.png", width = 600 / 96, height = 330 / 96, dpi = 300)











