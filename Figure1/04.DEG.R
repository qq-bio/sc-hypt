library(Seurat)
library(dplyr)

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")




################################################################################
### DEG Analysis - Treatment vs Control
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")

input_file <- c(
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

treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")

for(i in input_file){
  
  deg_merged = c()
  
  outfile <- file.path("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                       gsub("anno.L2.rds", "DEG.L1.out", basename(i)))
  
  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics|EC]+.anno.L2.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  Idents(seurat_object) = cluster
  
  meta_table = seurat_object@meta.data
  
  project_list = unique(meta_table$Project)
  for(pi in project_list){
    
    strain_list = unique(meta_table[meta_table$Project==pi, ]$Strain)
    for(si in strain_list){
      
      cell_list = unique(meta_table[meta_table$Project==pi &
                                      meta_table$Strain==si, cluster])
      
      treatment = meta_table[meta_table$Project==pi &
                               meta_table$Strain==si, ]$Treatment
      treatment = intersect(treatment_order, unique(treatment))
      control = treatment[1]
      
      for(cell in cell_list){
        
        cell.1 = rownames(meta_table[meta_table$Project==pi &
                                       meta_table[,cluster]==cell &
                                       meta_table$Treatment==control &
                                       meta_table$Strain==si, ])
        
        if(length(cell.1)>3){
          
          for(ti in treatment[-1]){
            
            cell.2 = rownames(meta_table[meta_table$Project==pi &
                                           meta_table[,cluster]==cell &
                                           meta_table$Treatment==ti &
                                           meta_table$Strain==si, ])
            
            if(length(cell.2)>3){
              
              deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)
              
              deg$gene_name = rownames(deg)
              deg$cell_type = cell
              deg$pct.diff = deg$pct.1 - deg$pct.2
              deg$project = pi
              deg$strain = si
              deg$tissue = tissue
              deg$control = control
              deg$treatment = ti
              deg$control_size = length(cell.1)
              deg$treatment_size = length(cell.2)
              deg$test_gene = nrow(deg)
              
              deg_merged = rbind(deg_merged, deg)
              
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  write.table(deg_merged, outfile)
  
  print(i)
  
}




cluster = "Cell_type_L2"

treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")

for(i in input_file){
  
  deg_merged = c()
  
  outfile <- file.path("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/",
                       gsub("anno.L2.rds", "DEG.L2.out", basename(i)))
  
  seurat_object = readRDS(i)
  dataset = gsub("\\.[RNA|multiomics|EC]+.anno.L2.rds", "", basename(i), perl = T)
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  Idents(seurat_object) = cluster
  
  meta_table = seurat_object@meta.data
  
  project_list = unique(meta_table$Project)
  for(pi in project_list){
    
    strain_list = unique(meta_table[meta_table$Project==pi, ]$Strain)
    for(si in strain_list){
      
      cell_list = unique(meta_table[meta_table$Project==pi &
                                      meta_table$Strain==si, cluster])
      
      treatment = meta_table[meta_table$Project==pi &
                               meta_table$Strain==si, ]$Treatment
      treatment = intersect(treatment_order, unique(treatment))
      control = treatment[1]
      
      for(cell in cell_list){
        
        cell.1 = rownames(meta_table[meta_table$Project==pi &
                                       meta_table[,cluster]==cell &
                                       meta_table$Treatment==control &
                                       meta_table$Strain==si, ])
        
        if(length(cell.1)>3){
          
          for(ti in treatment[-1]){
            
            cell.2 = rownames(meta_table[meta_table$Project==pi &
                                           meta_table[,cluster]==cell &
                                           meta_table$Treatment==ti &
                                           meta_table$Strain==si, ])
            
            if(length(cell.2)>3){
              
              deg = FindMarkers(seurat_object, ident.1=cell.1, ident.2=cell.2, logfc.threshold = 0)
              
              deg$gene_name = rownames(deg)
              deg$cell_type = cell
              deg$pct.diff = deg$pct.1 - deg$pct.2
              deg$project = pi
              deg$strain = si
              deg$tissue = tissue
              deg$control = control
              deg$treatment = ti
              deg$control_size = length(cell.1)
              deg$treatment_size = length(cell.2)
              deg$test_gene = nrow(deg)
              
              deg_merged = rbind(deg_merged, deg)
              
              
            }
            
          }
          
        }
        
      }
      
    }
    
  }
  
  write.table(deg_merged, outfile)
  
  print(i)
  
}







################################################################################
### Merge and summarize DEG results
################################################################################
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/")
logFC_threshold <- 0.5

l1_deg_list <- list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "DEG.L1.out")

deg_count_merged <- c()
deg_merged <- c()
for( i in input_file ){

  deg <- read.table(i, header = T, sep = " ")
  deg_sum <- deg %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) >= logFC_threshold) %>%
    group_by(project, strain, tissue, control, treatment, cell_type, control_size, treatment_size) %>%
    dplyr::summarize(n_sig = n())

  deg_count_merged <- rbind(deg_count_merged, deg_sum)
  deg_merged <- rbind(deg_merged, deg)
}

write.table(deg_count_merged, "DEG_summary.L1.tsv", col.names = T, row.names = F, sep = "\t")
write.table(deg_merged, "DEG.L1.all.out", col.names = T, row.names = F, sep = "\t")




l2_deg_list <- list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/", "DEG.L2.out")

deg_count_merged <- c()
deg_merged <- c()
for( i in input_file ){
  
  deg <- read.table(i, header = T, sep = " ")
  deg_sum <- deg %>%
    filter(p_val_adj < 0.05 & abs(avg_log2FC) >= logFC_threshold) %>%
    group_by(project, strain, tissue, control, treatment, cell_type, control_size, treatment_size) %>%
    dplyr::summarize(n_sig = n())
  
  deg_count_merged <- rbind(deg_count_merged, deg_sum)
  deg_merged <- rbind(deg_merged, deg)
}

write.table(deg_count_merged, "DEG_summary.L2.tsv", col.names = T, row.names = F, sep = "\t")
write.table(deg_merged, "DEG.L2.all.out", col.names = T, row.names = F, sep = "\t")







################################################################################
### Visualize DEG Results (Figure 1h; Fig. S11a)
################################################################################
library(lemon)
library(ggplot2)
library(dplyr)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", sep='\t', header=TRUE)

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
                      levels = sxt_order)
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

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/DEG.L1.png", 
       width = 577/96, height = 913/96, dpi = 300)






# L2-DEG
deg_count_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG_summary.L2.tsv", header = T, sep = "\t")

deg_count_merged$Cell_type_L1 = factor(deg_count_merged$Cell_type_L1, levels = cell_order)
deg_count_merged$tissue = factor(deg_count_merged$tissue, levels = tissue_order)
deg_count_merged$strain = factor(deg_count_merged$strain, levels = strain_order)

### overall plot
deg_count_merged_top3 <- deg_count_merged %>%
  group_by(project, strain, tissue) %>%
  slice_max(n_sig, n = 3)

deg_count_merged_top3 <- deg_count_merged %>%
  group_by(project, strain, tissue, Cell_type_L2) %>%
  slice_max(order_by = n_sig, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  group_by(project, strain, tissue) %>%
  slice_max(order_by = n_sig, n = 3, with_ties = FALSE) %>%
  ungroup()

ggplot(deg_count_merged, aes(x = control_size + treatment_size, y = n_sig)) +
  geom_point(aes(color = Cell_type_L1)) +
  scale_color_manual(values=cell_col) +
  scale_x_continuous(trans='log10') +
  labs(x = "Cell type size", y = "Number of DEGs") +
  theme_classic(base_family = "Arial") +
  theme(
    axis.text.x = element_text(colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black', size=10),
    panel.border = element_rect(color = "black", fill = NA),
    legend.position = "none"
  ) +
  geom_text_repel(
    data = deg_count_merged_top3,
    aes(x = control_size + treatment_size, y = n_sig, label = Cell_type_L2),
    size = 3,
    nudge_y = 0.5 * max(deg_count_merged$n_sig),
    force = 10,
    max.overlaps = Inf,
    box.padding = 0.01,
    segment.color = "grey80"
  ) +
  facet_grid(strain ~ tissue, scales = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/DEG.L2.count_size.dotplot.png", width=1160/96, height=850/96, dpi=300)










################################################################################
### Visualize Overlap of DEG with Known BP Gene Lists (Fig. S11b, d)
################################################################################
library(plotly)
library(htmlwidgets)
library(dplyr)

options(viewer = NULL)


# Load data and filter for unique gene mappings to human
m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2human.out.txt", header = TRUE, sep = "\t") %>%
  filter(V4 == 1) %>%
  select(Gene.name, Human.gene.name) %>%
  unique()

r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.rat2human.out.txt", header = TRUE, sep = "\t") %>%
  filter(V5 == 1) %>%
  select(Gene.name, Human.gene.name) %>%
  unique()


# Load DEG data and known BP gene lists
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = TRUE, sep = '\t')
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




### Fig. S11d
calculate_overlap <- function(deg_count, count_label) {
  if (count_label == "Tissue") {
    df <- data.frame(
      Count_Label = factor(rep(c("=1", ">1"), each = 2), levels = c("=1", ">1")),
      Reference = factor(rep(c("Mishra MK, et al.\n(n = 210)", "Keaton JM, et al.\n(n = 1,371)"), 2)),
      Proportion = c(
        length(deg_count[deg_count == 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count == 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count > 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count > 1])
      )* 100,
      Count = factor(rep(count_label, 4))
    )
  } else {
    df <- data.frame(
      Count_Label = factor(rep(c("=1", ">1"), each = 2), levels = c("=1", ">1")),
      Reference = factor(rep(c("Mishra MK, et al.\n(n = 210)", "Keaton JM, et al.\n(n = 1,371)"), 2)),
      Proportion = c(
        length(deg_count[deg_count == 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count == 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count == 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% yong_list]) / length(deg_count[deg_count > 1]),
        length(deg_count[deg_count > 1 & names(deg_count) %in% helen_list]) / length(deg_count[deg_count > 1])
      )* 100,
      Count = factor(rep(count_label, 4))
    )
  }
  return(df)
}

# Create data frames for tissue and strain counts
deg_tissue_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "tissue")])$gene_name)
deg_strain_count = table(unique(all_genes[all_genes$p_val_adj<0.05 & abs(all_genes$avg_log2FC) > 0.5 & all_genes$strain %in% c("C57BL/6", "SS", "SHR"), c("gene_name", "strain")])$gene_name)

deg_tissue_df <- calculate_overlap(deg_tissue_count, "Tissue")
deg_strain_df <- calculate_overlap(deg_strain_count, "Strain")

# Combine data frames
df <- bind_rows(deg_tissue_df, deg_strain_df)

ggplot(df, aes(x = Count_Label, y = Proportion, fill = Count_Label)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Number of shared tissues/strains", y = "Overlapped DEG proportion (%)") +
  geom_text(aes(label = paste0(round(Proportion, 2), "%")), position = position_stack(vjust = 0.5)) +
  theme_classic() +
  theme(text = element_text(family = "Arial"),
        axis.text = element_text(color="black"),
        axis.title.x = element_text(hjust = 0.5),
        legend.position = "None") + 
  coord_flip() + 
  facet_grid(Count ~ Reference, scales = "free", space = "free_y")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.multi_st.overlap.prop.png", width=323/96, height=335/96, dpi=300)







################################################################################
### Permutation test for DEG (Figure 1i; Fig. S11c)
################################################################################
# function for permutation
analyze_deg_overlap_table <- function(data, group_var = "tissue", 
                                      strain_filter = c("C57BL/6", "SS", "SHR"), 
                                      pval_thresh = 0.05, logfc_thresh = 0.5, 
                                      n_sim = 1000, seed = 42) {
  set.seed(seed)
  group_var <- rlang::sym(group_var)
  
  # Filter DEGs
  deg_filtered <- data %>%
    filter(p_val_adj < pval_thresh,
           abs(avg_log2FC) > logfc_thresh,
           strain %in% strain_filter)
  
  # Observed: how many groups each gene is DEG in
  deg_gene_group <- deg_filtered %>%
    distinct(gene_name, !!group_var)
  
  gene_group_table_obs <- table(deg_gene_group$gene_name)
  freq_obs <- as.vector(table(gene_group_table_obs))  # how many genes appear in 1, 2, 3... groups
  names(freq_obs) <- names(table(gene_group_table_obs))
  
  # Prepare for simulation
  all_genes <- unique(data$gene_name)
  
  # DEG counts per group
  deg_count_per_group <- deg_gene_group %>%
    dplyr::count(!!group_var, name = "n_deg")
  
  # Simulate
  sim_freq_list <- replicate(n_sim, {
    sim_deg_gene_group <- deg_count_per_group %>%
      mutate(gene_name = map(n_deg, ~ sample(all_genes, .x))) %>%
      unnest(gene_name)
    
    gene_group_table_sim <- table(sim_deg_gene_group$gene_name)
    freq_sim <- table(gene_group_table_sim)
    
    # Ensure all group sizes are present
    out <- rep(0, max(length(freq_obs), as.numeric(names(freq_sim)) %>% max()))
    names(out) <- as.character(1:length(out))
    out[names(freq_sim)] <- as.numeric(freq_sim)
    out
  }, simplify = "matrix")
  
  sim_freq_mean <- rowMeans(sim_freq_list)
  sim_freq_sd <- apply(sim_freq_list, 1, sd)
  z_scores <- (freq_obs - sim_freq_mean) / sim_freq_sd
  
  result <- list(
    observed_counts = freq_obs,
    sim_mean = sim_freq_mean,
    sim_sd = sim_freq_sd,
    z_score = z_scores,
    sim_matrix = sim_freq_list
  )
  
  return(result)
}


all_genes = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = T)
all_genes$comp = paste(all_genes$strain, all_genes$treatment,  sep="-")

# permutation
result_tissue <- analyze_deg_overlap_table(all_genes, group_var = "tissue")
result_strain <- analyze_deg_overlap_table(all_genes, group_var = "strain")

# summarize permutation results
result = result_tissue
observed_kplus <- sum(result$observed_counts[names(result$observed_counts) >= 2])
sim_kplus <- colSums(result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 2, ])
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                        "\n(1,000 perms)")
result_tissue=result


result = result_strain
observed_kplus <- sum(result$observed_counts[names(result$observed_counts) >= 2])
sim_kplus <- colSums(result$sim_matrix[as.numeric(rownames(result$sim_matrix)) >= 2, ])
result$p_empirical <- mean(sim_kplus >= observed_kplus)
result$p_text <- paste0("p-value = ", signif(result$p_empirical, 2), 
                        "\n(1,000 perms)")
result_strain=result


# visualize results using bar plot
result = result_tissue; group_type="tissues"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_counts)),
  Observed = as.numeric(result$observed_counts),
  Expected = result$sim_mean[1:length(result$observed_counts)],
  SD = result$sim_sd[1:length(result$observed_counts)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
# write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.tissue.obs_exp.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p1<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  scale_x_continuous(
    breaks = seq(min(plot_df$GroupCount), max(plot_df$GroupCount), by = 1)
  ) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4.2) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Number of genes"
  )

result = result_strain; group_type="strains"
plot_df <- tibble(
  GroupCount = as.numeric(names(result$observed_counts)),
  Observed = as.numeric(result$observed_counts),
  Expected = result$sim_mean[1:length(result$observed_counts)],
  SD = result$sim_sd[1:length(result$observed_counts)]
) %>%
  mutate(
    Lower = Expected - SD,
    Upper = Expected + SD
  )
# write.table(plot_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/deg.count.strain.obs_exp.tsv", sep = "\t", col.names = T, row.names = F)
plot_long <- plot_df %>%
  dplyr::select(GroupCount, Observed, Expected) %>%
  pivot_longer(cols = c("Observed", "Expected"), names_to = "Type", values_to = "Value")
x_pos <- max(plot_df$GroupCount) + 0.5
y_pos <- max(plot_df$Expected) * 0.95
p_val <- result$p_text
p2<-ggplot(plot_df, aes(x = GroupCount)) +
  geom_col(aes(y = Observed), fill = "red", alpha = 0.3, width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, color = "gray40") +
  geom_line(data = plot_long, aes(y = Value, color = Type), size = 1.2) +
  scale_color_manual(values = c("Observed" = "red", "Expected" = "gray40")) +
  annotate("text", x = x_pos, y = y_pos, label = p_val,
           hjust = 1, vjust = 1, size = 4.2) +
  labs(
    x = paste("Number of shared", group_type),
    y = "Number of genes"
  )

p1+p2&theme(text = element_text(family = "Arial"))

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/deg.count.obs_exp.png", width=793/96, height=220/96, dpi=300)









