library(GenomicRanges)
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)
library(ggraph)
library(leiden)
library(cluster)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")


################################################################################
### load reference files
### h2m
r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) <- c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
r2h[r2h$gene_name_rat=="",]$gene_name_rat <- r2h[r2h$gene_name_rat=="",]$gene_id_rat

m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) <- c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")
m2h[m2h$gene_name_mouse=="",]$gene_name_mouse <- m2h[m2h$gene_name_mouse=="",]$gene_id_mouse

### grch38
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

gene_gr <- GRanges(
  seqnames = Rle(ensembl$Chromosome.scaffold.name),
  ranges = IRanges(start = ensembl$Gene.start..bp., end = ensembl$Gene.end..bp.),
  strand = Rle(ensembl$Strand),
  gene_id = ensembl$Gene.stable.ID,
  gene_name = ensembl$Gene.name,
  gene_type = ensembl$Gene.type
)





################################################################################
### Identify Cell Clusters Based on DE SNP-Gene
################################################################################
### load enrichment results
fisher_test_df <- read.table("trait.fisher.snp_gene.enrichment_results.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df %>%
  mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep = "_"))

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_filtered <- deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

snp_gene_df <- read.table("snp_gene.evi_org.out", sep = "\t", header = TRUE)
deg_use <- deg_filtered %>% filter(gene_name %in% snp_gene_df$combined_gene_name)

# Generate Binary Table
set.seed(42)
deg_use <- deg_use %>% mutate(condition = paste(strain, treatment, tissue, cell_type, sep = "_"))
binary_matrix <- deg_use %>%
  select(gene_name, condition) %>%
  mutate(value = 1) %>%
  spread(key = condition, value = value, fill = 0)
row.names(binary_matrix) <- binary_matrix$gene_name
binary_matrix <- binary_matrix %>% select(-gene_name) %>% as.matrix()

# Compute Jaccard Distance and Similarity
jaccard_dist_columns <- vegdist(t(binary_matrix), method = "jaccard")
similarity_matrix <- 1 - as.matrix(jaccard_dist_columns)
similarity_graph <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "undirected", weighted = TRUE)

# Consensus Clustering
n_runs <- 1000
all_clusterings <- matrix(0, nrow = length(V(similarity_graph)), ncol = n_runs)

set.seed(42)
for (i in 1:n_runs) {
  all_clusterings[, i] <- leiden(similarity_graph, resolution_parameter = 1)
}

# Create Consensus Matrix
consensus_matrix <- matrix(0, nrow = length(V(similarity_graph)), ncol = length(V(similarity_graph)))
for (i in 1:n_runs) {
  clustering <- all_clusterings[, i]
  for (j in 1:length(clustering)) {
    for (k in j:length(clustering)) {
      if (clustering[j] == clustering[k]) {
        consensus_matrix[j, k] <- consensus_matrix[j, k] + 1
        consensus_matrix[k, j] <- consensus_matrix[k, j] + 1
      }
    }
  }
}
consensus_matrix <- consensus_matrix / n_runs
hclust_result <- hclust(as.dist(1 - consensus_matrix), method = "complete")
consensus_clusters <- cutree(hclust_result, k = 9)
names(consensus_clusters) <- colnames(binary_matrix)

# Create Annotation Dataframe
annotation_df <- data.frame(
  cell_group = colnames(binary_matrix),
  strain = fisher_test_df$strain[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  tissue = fisher_test_df$tissue[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  cell_type = fisher_test_df$cell_type[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  cluster = consensus_clusters
)

# Define New Cluster Assignments
new_clusters <- factor(ifelse(annotation_df$cluster %in% c(6), "Cluster 4",
                              ifelse(annotation_df$cluster %in% c(8), "Cluster 5",
                                     ifelse(annotation_df$cluster %in% c(2), "Cluster 3",
                                            ifelse(annotation_df$cluster %in% c(3, 7), "Cluster 2",
                                                   ifelse(annotation_df$cluster %in% c(1, 4, 5), "Cluster 1",
                                                          ifelse(annotation_df$cluster == 9, "Cluster 6", as.character(annotation_df$cluster))))))))
annotation_df$new_cluster <- new_clusters

# Generate Heatmap Annotations
strain_colors <- species_col
tissue_colors <- tissue_col
cell_colors <- cell_col

top_annotation <- HeatmapAnnotation(
  Strain = annotation_df$strain,
  Tissue = annotation_df$tissue,
  Cell = annotation_df$cell_type,
  col = list(
    Strain = strain_colors,
    Tissue = tissue_colors,
    Cell = cell_colors
  )
)

# Generate and Draw Heatmap
heatmap <- Heatmap(as.matrix(consensus_matrix), 
                   name = "Similarity",
                   col = colorRamp2(c(0, 1), c("white", "red")),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_columns = "euclidean",
                   clustering_method_rows = "complete",
                   clustering_method_columns = "complete",
                   top_annotation = top_annotation,
                   row_split = annotation_df$new_cluster,
                   column_split = annotation_df$new_cluster, 
                   column_title_rot = 45,
                   row_title_rot = 360,
                   heatmap_legend_param = list(title = "Consensus Similarity", legend_direction = "vertical"))

draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")

fisher_test_df <- fisher_test_df %>%
  filter(cell_group %in% names(new_clusters)) %>%
  mutate(cluster = consensus_clusters[cell_group],
         consensus_cluster = new_clusters[cell_group])

write.table(fisher_test_df, "trait.fisher.cluster.out", col.names = T, row.names = F, sep = "\t", quote=F)





################################################################################
### Representative Cell Types in Each Cluster (Figure 4g)
################################################################################
fisher_test_df = read.table("trait.fisher.cluster.out", sep="\t", header=T)

top5_data <- fisher_test_df %>%
  filter(p.adj < 0.05) %>%
  group_by(consensus_cluster) %>%
  arrange(desc(NES)) %>%
  slice_head(n = 10) %>%
  ungroup()

# Generate the plot
p <- ggplot(top5_data, aes(x = cell_group, y = NES)) +
  geom_point(aes(fill = NES, size = -log10(p.value)), shape = 21) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(), 
    legend.justification = c(1, 0.5),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.spacing.y = unit(-0.2, 'cm'),
    
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(colour = 'black'),
    strip.text.x = element_text(colour = 'black', angle = 360),
    axis.title.x.top = element_text(), 
    axis.text.x.top = element_text()
  ) +
  labs(title = "Representative cell type in each cluster",
       x = "", y = "NES", fill = "NES", size = "-log10(p-value)") +
  facet_grid(. ~ consensus_cluster, scales = "free", space = "free", switch = "x") +
  scale_x_discrete(position = "top") +
  guides(
    fill = guide_legend(ncol = 2), 
    size = guide_legend(ncol = 2) 
  )

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4g.rep_cell_types_in_clusters.png", width=626/96, height=383/96, dpi=300)





################################################################################
### Enrichment of traits in Each Cluster (Figure 4h)
################################################################################
# Summarize the NES scores for each cluster-trait combination
nes_summary <- fisher_test_df %>%
  group_by(consensus_cluster, trait) %>%
  summarise(average_NES = mean(NES, na.rm = TRUE), 
            cluster_size = n(),
            NES_values = list(NES), .groups = 'drop') %>%
  ungroup()

# Calculate the p-value for NES values in each cluster compared to all other clusters for the same trait
nes_summary <- nes_summary %>%
  rowwise() %>%
  mutate(
    p_value = {
      cluster_nes <- unlist(NES_values)
      current_trait <- trait
      current_cluster <- consensus_cluster
      other_nes <- fisher_test_df %>%
        filter(trait == current_trait, consensus_cluster != current_cluster) %>%
        pull(NES)
      
      if (length(cluster_nes) > 1 && length(other_nes) > 1) {
        wilcox.test(cluster_nes, other_nes, alternative = "greater")$p.value
      } else if (length(cluster_nes) > 1 && length(other_nes) == 0) {
        NA  # If there are no other NES values to compare, return NA
      } else {
        1  # If not enough data for a meaningful test, default to non-significant
      }
    }
  ) %>%
  ungroup()

# Mark significant results
nes_summary <- nes_summary %>%
  mutate(significant = ifelse(p_value < 0.05, "*", ""),
         cluster = consensus_cluster) %>%
  filter(consensus_cluster!="Cluster 6")

# Plotting
p1 <- ggplot(nes_summary, aes(x = trait, y = factor(cluster), size = cluster_size, fill = average_NES)) +
  geom_point(shape = 21, stroke = 0) +
  geom_point(data = nes_summary[nes_summary$p_value < 0.05, ],
             aes(x = trait, y = factor(cluster), size = cluster_size, fill = average_NES),
             shape = 21, stroke = 1, color = "black") +
  scale_size_continuous(range = c(2, 6), breaks = c(10, 30, 50)) +
  scale_fill_gradient2(low = "lightblue", mid = "white", high = "red", midpoint = 2.5, name = "Average NES") +
  labs(title = " ",
       x = "Trait",
       y = "",
       size = "Cluster size",
       fill = "Average NES") +
  theme_minimal() +
  theme(axis.text.y = element_text(colour = 'black'),
        # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        axis.text.x = element_blank()) +
  coord_flip()
print(p1)



### Visualize pathway enrichment results in each cluster
merge_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/metascape.clusters.out", sep='\t', header=T)
top_list = merge_all %>%
  filter(Log.q.value. < -2) %>%
  group_by(Description) %>%
  mutate(pathway_count = n_distinct(cluster)) %>%
  filter(pathway_count<3) %>% ungroup() %>%
  filter(grepl("Summary", GroupID)) %>%
  group_by(cluster) %>% arrange(Log.q.value.) %>%
  slice_head(n=5)
merge_use = merge_all[merge_all$Description %in% top_list$Description & 
                        merge_all$Log.q.value. < -2 &
                        grepl("Member", merge_all$GroupID), c("pathway", "cluster", "Log.q.value.")]
sorted_df <- merge_use[order(merge_use$cluster, merge_use$pathway), ]
sorted_df[sorted_df$pathway=="Diseases of signal transduction by growth factor receptors and second messengers (Reactome)", ]$pathway = "Diseases of signal transduction by growth factor\nreceptors and second messengers (Reactome)"
sorted_df$pathway <- factor(sorted_df$pathway, levels = unique(sorted_df$pathway))

p2 = ggplot(sorted_df, aes(y = pathway, x = cluster, fill = -1*Log.q.value.)) +
  geom_tile(color="black") + 
  scale_fill_gradient(low="white", high="purple") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black', size=10)) +
  labs(x="", y="Pathway", fill="-log(q-value)")

design = "A
          B"
combined_plot <- (p1 + p2) + plot_layout(design = design, guides = 'collect', heights = c(0.12, 0.3)) &
  theme(legend.position = 'right', 
        legend.box = 'vertical',
        plot.margin = margin(t = -10, r = 0, b = -10, l = 10, unit = "pt"))
print(combined_plot)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig4h.cell_cluster.asso.png", width=540/96, height=597/96, dpi=300)










################################################################################
### Cluster-Based Multi-Functional SNPs (Figure S15b-d)
################################################################################
fisher_test_df <- read.table("trait.fisher.cluster.out", header = T, sep = "\t", comment.char = "")
snp_gene_df <- read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

all_traits <- c("diastolic_bp", "eGFR", "systolic_bp", "CAD", "albuminuria", "pulse_pressure", "stroke", "essential_hypertension", "BUN")
gwas_merge_use <- gwas_merge_use %>%
  group_by(SNPS) %>%
  complete(trait = all_traits, fill = list(P.VALUE = NA)) %>%
  ungroup()

# Separate expressional and e2g tissues/groups
snp_gene_df_expanded <- snp_gene_df %>%
  separate_rows(expressional, sep = ", ") %>%
  separate_rows(e2g, sep = ", ") %>%
  mutate(
    expressional = ifelse(expressional == "no", NA, expressional),
    e2g = ifelse(e2g == "no", NA, e2g)
  ) %>%
  pivot_longer(cols = c(expressional, e2g), names_to = "evidence_type", values_to = "evidence_value") %>%
  filter(!is.na(evidence_value))

# Fill missing tissues and groups
all_tissues <- unique(eqtl_merge$tissue)
all_groups <- unique(e2g_merge$group)
snp_gene_df_expanded <- snp_gene_df_expanded %>%
  complete(evidence_value = unique(c(all_tissues, all_groups)), fill = list(evidence_type = NA))

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
deg_filtered$cell_group = paste(deg_filtered$strain, deg_filtered$treatment, deg_filtered$tissue, deg_filtered$cell_type, sep = "_")
consensus_clusters = fisher_test_df$consensus_cluster
names(consensus_clusters) = fisher_test_df$cell_group
deg_filtered$consensus_cluster = consensus_clusters[deg_filtered$cell_group]



plot_snp_analysis <- function(snp_id, cluster_name, gwas_merge_use, snp_gene_df, fisher_test_df, deg_filtered) {
  
  # Filter the GWAS data for the selected SNP
  gwas_data <- gwas_merge_use %>% filter(SNPS %in% snp_id)
  
  # Plot GWAS results
  p1 <- ggplot(gwas_data, aes(x = SNPS, y = trait)) +
    geom_tile(aes(fill = -log10(P.VALUE)), color = "white") +
    scale_fill_gradient(low = "white", high = "red", na.value = "grey90", name = "-log10(P.VALUE)") +
    # theme_minimal() +
    theme_bw() +
    theme(
      legend.position = "bottom",
      # legend.direction = "vertical",
      legend.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      axis.text.y = element_text(colour = 'black')
    ) +
    labs(x = "SNP", y = "Trait")
  
  # Filter the SNP-gene data for the selected SNP
  snp_gene_data <- snp_gene_df_expanded %>% 
    filter(SNP %in% snp_id, !(is.na(combined_gene_name)))
  
  # Combine all rows, including those with only "proximal" evidence
  snp_gene_data <- bind_rows(
    snp_gene_data,
    snp_gene_df %>%
      filter(SNP %in% snp_id, evidence_summary == "proximal", !(is.na(combined_gene_name))) %>%
      mutate(evidence_type = "proximal", evidence_value = "proximal") %>%
      select(SNP, gene, Gene.name, combined_gene_name, proximal, evidence_summary, evidence_type, evidence_value)
  )
  
  # Plot evidence summary
  p2 <- ggplot(snp_gene_data, aes(x = combined_gene_name, y = evidence_value)) +
    geom_tile(aes(fill = evidence_type), color = "white") +
    scale_fill_manual(values = c("expressional" = "brown", "e2g" = "darkgreen", "proximal" = "orange"), na.value = "grey90", name = "Evidence Type") +
    theme_bw() +
    theme(
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.title = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      axis.text.y = element_text(colour = 'black')
    ) +
    labs(x = "Gene", y = "Tissue/Group")
  
  # Filter the DEGs for the selected cluster
  result <- fisher_test_df %>% filter(consensus_cluster %in% cluster_name)
  deg_filtered_matching <- deg_filtered %>%
    filter(cell_group %in% result$cell_group,
           gene_name %in% snp_gene_data$combined_gene_name)
  
  # Plot DEG results
  p3 <- NULL
  if(nrow(deg_filtered_matching) > 0) {
    p3 <- ggplot(deg_filtered_matching, aes(x = treatment, y = paste(tissue, cell_type, sep = "-"))) +
      geom_point(aes(fill = -avg_log2FC, size = -log10(p_val)), shape = 21) +
      scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab", name = "log2(FC)") +
      scale_size_continuous(name = "-log10(p-value)") +
      theme_bw() +
      theme(
        panel.grid.major.y = element_blank(),
        legend.justification = c(1, 0.5),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        strip.text = element_text(colour = 'black'),
        strip.text.y = element_text(colour = 'black', angle = 360),
        strip.background = element_rect(colour = "black", fill = NA)
      ) +
      labs(x = "", y = "", fill = "log2(FC)", size = "-log10(p-value)") +
      facet_nested(consensus_cluster + gene_name ~ project + strain, scales = "free", space = "free")
  }
  
  # Combine plots
  if (is.null(p3)) {
    combined_plot <- p1 + p2 + plot_layout(ncol = 2)
  } else {
    combined_plot <- (p1 + p2 + p3) + plot_layout(ncol = 3, widths = c(0.2, 0.2, 0.6)) &
      theme(legend.position = 'bottom', 
            legend.box = 'vertical')
  }
  
  print(combined_plot)
}



# Usage example
plot_snp_analysis(snp_id = "rs179993", cluster_name = "Cluster 4", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs11838776", cluster_name = "Cluster 4", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)


plot_snp_analysis(snp_id = c("rs77924615"), cluster_name = c("Cluster 2"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

