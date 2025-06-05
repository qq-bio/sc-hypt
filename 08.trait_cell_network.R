# library(biomaRt)
library(rtracklayer)
library(GenomicRanges)
# library(openxlsx)
# library(httr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(ggh4x)
library(patchwork)
library(stringr)
library(igraph)
library(vegan)
library(ggraph)
Sys.setenv(RETICULATE_PYTHON = "/home/u1/qqiu/.conda/envs/stereopy/bin/python")
library(reticulate)
# use_python("/home/u1/qqiu//Python/pyenv_v3.8/bin/python")
library(leiden)
library(lsa)
library(cluster)

# options(timeout = 300)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data")


################################################################################
### load reference files
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

### h2m
r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", sep = '\t', header=T)
colnames(r2h) = c("gene_id_rat", "gene_name_rat", "gene_id", "gene_name", "r2h_orthology_conf")
r2h[r2h$gene_name_rat=="",]$gene_name_rat = r2h[r2h$gene_name_rat=="",]$gene_id_rat
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", sep = '\t', header=T)
colnames(m2h) = c("gene_id_mouse", "gene_id", "gene_name", "m2h_orthology_conf", "gene_name_mouse")
m2h[m2h$gene_name_rat=="",]$gene_name_rat = m2h[m2h$gene_name_rat=="",]$gene_id_rat


gene_gr <- GRanges(
  seqnames = Rle(ensembl$Chromosome.scaffold.name),
  ranges = IRanges(start = ensembl$Gene.start..bp., end = ensembl$Gene.end..bp.),
  strand = Rle(ensembl$Strand),
  gene_id = ensembl$Gene.stable.ID,
  gene_name = ensembl$Gene.name,
  gene_type = ensembl$Gene.type
)












################################################################################
### load enrichment results
fisher_test_df = read.table("trait.fisher.snp_gene.0712.fc_0.25.permut_norm.out", header = T, sep = "\t", comment.char = "")
fisher_test_df <- fisher_test_df[!(fisher_test_df$tissue == "MCA" & fisher_test_df$strain %in% c("C57BL/6", "SS")), ]
fisher_test_df <- fisher_test_df %>%
  mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep = "_"))

#### process snp-gene table
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

proximal_merge <- read.table("gwas_snp_gene_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
e2g_merge <- read.table("e2g_merged.txt", header = TRUE, sep = '\t')

proximal_merge = proximal_merge[grepl("rs", proximal_merge$SNP),]
eqtl_merge = eqtl_merge[grepl("rs", eqtl_merge$SNP),]
e2g_merge = e2g_merge[grepl("rs", e2g_merge$SNP),]

snp_gene_prox <- data.frame(SNP = proximal_merge$SNP, gene = proximal_merge$Gene_ID, proximal = "yes")
snp_gene_eqtl <- eqtl_merge %>% select(SNP, gene_id_mod, tissue) %>%
  group_by(SNP, gene_id_mod) %>%
  summarise(expressional = paste(unique(tissue), collapse = ", ")) %>%
  rename(gene = gene_id_mod)
snp_gene_e2g <- e2g_merge %>% select(SNP, gene, group) %>%
  group_by(SNP, gene) %>%
  summarise(e2g = paste(unique(group), collapse = ", "))
snp_gene_df <- full_join(snp_gene_prox, snp_gene_eqtl, by = c("SNP", "gene")) %>%
  full_join(., snp_gene_e2g, by = c("SNP", "gene"))
snp_gene_df <- snp_gene_df %>%
  rowwise() %>%
  mutate(proximal = ifelse(is.na(proximal), "no", "yes"),
         expressional = ifelse(is.na(expressional), "no", expressional),
         e2g = ifelse(is.na(e2g), "no", e2g),
         evidence_summary = paste(na.omit(c(ifelse(proximal == "yes", "proximal", NA), 
                                            ifelse(expressional != "no", "expressional", NA), 
                                            ifelse(e2g != "no", "e2g", NA))), collapse = ", "))

r2h_agg <- r2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_rat = paste(unique(gene_name_rat), collapse = ";"))

m2h_agg <- m2h %>%
  group_by(gene_id) %>%
  summarise(gene_name_mouse = paste(unique(gene_name_mouse), collapse = ";"))

snp_gene_df <- snp_gene_df %>%
  mutate(pair = paste(SNP, gene, sep = "-")) %>%
  left_join(ensembl, by = c("gene" = "Gene.stable.ID")) %>%
  left_join(r2h_agg, by = c("gene" = "gene_id")) %>%
  left_join(m2h_agg, by = c("gene" = "gene_id")) %>%
  rowwise() %>%
  mutate(
    combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(SNP, gene, Gene.name, combined_gene_name, proximal, expressional, e2g, evidence_summary) %>%
  unique()

write.table(snp_gene_df, "snp_gene.evi_org.out", col.names = T, row.names = F, sep = "\t", quote=F)

snp_gene_df_mod <- read.table("snp_gene.evi_org.out", sep = "\t", header = T)
### process expr data
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
deg_use = deg_filtered[deg_filtered$gene_name %in% snp_gene_df_mod$combined_gene_name, ]


################################################################################
#### identify cell community based on DE SNP-gene
### generate binary table
set.seed(42)

deg_use <- deg_use %>%
  mutate(condition = paste(strain, treatment, tissue, cell_type, sep = "_"))
binary_matrix <- deg_use %>%
  dplyr::select(gene_name, condition) %>%
  mutate(value = 1) %>%
  spread(key = condition, value = value, fill = 0)
dim(binary_matrix)
# [1] 1320  148
row.names(binary_matrix) <- binary_matrix$gene_name
binary_matrix <- binary_matrix %>% dplyr::select(-gene_name) %>% as.matrix()

jaccard_dist_rows <- vegdist(binary_matrix, method = "jaccard")
jaccard_dist_columns <- vegdist(t(binary_matrix), method = "jaccard")
similarity_matrix <- 1 - as.matrix(jaccard_dist_columns)

# similarity_matrix <- cosine(binary_matrix)
similarity_graph <- graph_from_adjacency_matrix(as.matrix(similarity_matrix), mode = "undirected", weighted = TRUE)

### consensus clustering
n_runs <- 1000
all_clusterings <- matrix(0, nrow = length(V(similarity_graph)), ncol = n_runs)

set.seed(42)
for (i in 1:n_runs) {
  all_clusterings[, i] <- leiden(similarity_graph, resolution_parameter = 1)
}

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
consensus_clusters <- cutree(hclust_result, k = 9) # Adjust the number of clusters as needed
names(consensus_clusters) = colnames(binary_matrix)


annotation_df <- data.frame(
  cell_group = colnames(binary_matrix),
  strain = fisher_test_df$strain[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  tissue = fisher_test_df$tissue[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  cell_type = fisher_test_df$cell_type[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  community = consensus_clusters,
  new_community = new_clusters
)

strain_colors <- species_col
tissue_colors <- tissue_col
cell_colors <- cell_col
community_colors <- rainbow(length(unique(annotation_df$community)))

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

# Create heatmap and split by community
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
                   row_split = annotation_df$community,
                   column_split = annotation_df$community, 
                   column_title_rot = 45,
                   row_title_rot = 360,
                   heatmap_legend_param = list(title = "Consensus Similarity", legend_direction = "vertical"))

# Draw heatmap
draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")


new_clusters <- as.factor(ifelse(annotation_df$community %in% c(6), "Cluster 4",
                                 ifelse(annotation_df$community %in% c(8), "Cluster 5",
                                        ifelse(annotation_df$community %in% c(2), "Cluster 3",
                                               ifelse(annotation_df$community %in% c(3, 7), "Cluster 2",
                                                      ifelse(annotation_df$community %in% c(1, 4, 5), "Cluster 1",
                                                             ifelse(annotation_df$community == 9, "Cluster 6", as.character(annotation_df$community))))))))
names(new_clusters) = colnames(binary_matrix)

consensus_dist_matrix <- as.dist(1 - consensus_matrix)
sil <- silhouette(as.integer(new_clusters), consensus_dist_matrix)
plot(sil, main = "Silhouette plot for concensus clusters")
### test resolution
# resolution_range <- seq(0.1, 2.0, by = 0.1)
# num_communities <- numeric(length(resolution_range))
# for (i in seq_along(resolution_range)) {
#   resolution <- resolution_range[i]
#   communities <- leiden(similarity_graph, resolution_parameter = resolution)
#   num_communities[i] <- length(unique(communities))
# }
# plot(resolution_range, num_communities, type = "b", xlab = "Resolution Parameter", ylab = "Number of Communities",
#      main = "Effect of Resolution Parameter on Community Detection")


# # Map genes to clusters
# deg_use <- deg_use %>% mutate(cluster = new_clusters[deg_use$condition])
# gene_clusters <- deg_use %>% select(gene_name, cluster)
# head(gene_clusters)
# 
# # To see the genes in each cluster:
# genes_in_clusters <- split(gene_clusters$gene_name, gene_clusters$cluster)
# genes_in_clusters


### cell-wise similarity
annotation_df <- data.frame(
  cell_group = colnames(binary_matrix),
  strain = fisher_test_df$strain[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  tissue = fisher_test_df$tissue[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  cell_type = fisher_test_df$cell_type[match(colnames(binary_matrix), fisher_test_df$cell_group)],
  community = consensus_clusters,
  new_community = new_clusters
)

strain_colors <- species_col
tissue_colors <- tissue_col
cell_colors <- cell_col
community_colors <- rainbow(length(unique(annotation_df$community)))

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

# Create heatmap and split by community
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
                   row_split = annotation_df$new_community,
                   column_split = annotation_df$new_community, 
                   column_title_rot = 45,
                   row_title_rot = 360,
                   heatmap_legend_param = list(title = "Consensus Similarity", legend_direction = "vertical"))

# Draw heatmap
draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")



################################################################################
#### community - trait association
fisher_test_df <- fisher_test_df %>%
  # mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep = "_")) %>%
  filter(cell_group %in% names(new_clusters)) %>%
  mutate(cluster = consensus_clusters[cell_group],
         consensus_cluster = new_clusters[cell_group])

write.table(fisher_test_df, "trait.fisher.cluster.out", col.names = T, row.names = F, sep = "\t", quote=F)


fisher_test_df = read.table("trait.fisher.cluster.out", sep="\t", header=T)

gene_clusters <- fisher_test_df %>%
  separate_rows(gene_list, sep = ", ") %>%
  select(consensus_cluster, gene_list) %>%
  filter(gene_list!="") %>%
  distinct()

# Create a list of genes for each cluster
genes_in_clusters <- gene_clusters %>%
  group_by(consensus_cluster) %>%
  summarize(genes = list(gene_list))

# Convert the list into a dataframe with each column representing a cluster
max_genes <- max(sapply(genes_in_clusters$genes, length))
genes_df <- do.call(cbind, lapply(genes_in_clusters$genes, `length<-`, max_genes))
colnames(genes_df) <- genes_in_clusters$consensus_cluster

# Print the dataframe
write.table(genes_df, "trait.fisher.cluster_gene.out", col.names = T, row.names = F, sep = "\t", quote=F)

# > length(intersect(gene_clusters[gene_clusters$consensus_cluster=="Cluster 3",]$gene_list, gene_clusters[gene_clusters$consensus_cluster=="Cluster 4",]$gene_list))
# [1] 326
# > length(gene_clusters[gene_clusters$consensus_cluster=="Cluster 3",]$gene_list)
# [1] 1078
# > length(gene_clusters[gene_clusters$consensus_cluster=="Cluster 4",]$gene_list)
# [1] 548


################################################################################
### representative cell types in each community
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
    panel.grid.major.x = element_blank(),   # No vertical grid lines
    legend.justification = c(1, 0.5),
    # legend.direction = "horizontal",
    # legend.box = "vertical",
    # legend.box.just = "top",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10),
    legend.spacing.y = unit(-0.2, 'cm'),  # Reduce vertical spacing between legend items
    
    axis.text.x = element_text(angle = 45, hjust = 0, vjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(colour = 'black'),
    strip.text.x = element_text(colour = 'black', angle = 360),
    axis.title.x.top = element_text(),  # Move x-axis title to top
    axis.text.x.top = element_text()  # Move x-axis text to top
  ) +
  labs(title = "Representative cell type in each cluster",
       x = "", y = "NES", fill = "NES", size = "-log10(p-value)") +
  facet_grid(. ~ consensus_cluster, scales = "free", space = "free", switch = "x") +
  scale_x_discrete(position = "top") +
  guides(
    fill = guide_legend(ncol = 2),   # Set legend to have 2 columns
    size = guide_legend(ncol = 2)    # Set legend size to have 2 columns
  )

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/rep_cell_types_in_clusters.png", width=626/96, height=383/96, dpi=300)



# Summarize the NES scores for each community-trait combination
nes_summary <- fisher_test_df %>%
  group_by(consensus_cluster, trait) %>%
  summarise(average_NES = mean(NES, na.rm = TRUE), 
            community_size = n(),
            NES_values = list(NES), .groups = 'drop') %>%
  ungroup()

# Calculate the p-value for NES values in each community compared to all other communities for the same trait
nes_summary <- nes_summary %>%
  rowwise() %>%
  mutate(
    p_value = {
      community_nes <- unlist(NES_values)
      current_trait <- trait
      current_community <- consensus_cluster
      other_nes <- fisher_test_df %>%
        filter(trait == current_trait, consensus_cluster != current_community) %>%
        pull(NES)
      
      if (length(community_nes) > 1 && length(other_nes) > 1) {
        wilcox.test(community_nes, other_nes, alternative = "greater")$p.value
      } else if (length(community_nes) > 1 && length(other_nes) == 0) {
        NA  # If there are no other NES values to compare, return NA
      } else {
        1  # If not enough data for a meaningful test, default to non-significant
      }
    }
  ) %>%
  ungroup()

# Check the resulting p_values
print(nes_summary)

# Mark significant results
nes_summary <- nes_summary %>%
  mutate(significant = ifelse(p_value < 0.05, "*", ""),
         community = consensus_cluster) %>%
  filter(consensus_cluster!="Cluster 6")

# Plotting
p1 <- ggplot(nes_summary, aes(x = trait, y = factor(community), size = community_size, fill = average_NES)) +
  geom_point(shape = 21, stroke = 0) +
  # geom_point(data = nes_summary[nes_summary$average_NES > 5, ],
  #            aes(x = trait, y = factor(community), size = community_size, fill = average_NES),
  #            shape = 21, stroke = 1, color = "black") +
  geom_point(data = nes_summary[nes_summary$p_value < 0.05, ],
             aes(x = trait, y = factor(community), size = community_size, fill = average_NES),
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



### pathway enrichment in each community
# setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS")
# 
# merge_all = c()
# input_files = list.files(path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/", pattern = "zip")
# 
# for(i in input_files){
#   
#   unzip(zipfile=i, files = "metascape_result.xlsx", exdir="./tmp", overwrite = T)
#   data = read.xlsx("tmp/metascape_result.xlsx", 2)
#   
#   i_mod = gsub("(trait.fisher.|.zip)", "", i)
#   data$cluster = gsub("cluster", "Cluster ", i_mod)
#   
#   data$category_mod = lapply(strsplit(data$Category, " "), function(x) x[1]) %>% as.character()
#   data$pathway = paste0(data$Description, " (", data$category_mod, ")")
#   merge_all = rbind(merge_all, data)
#   
# }
# 
# write.table(merge_all, "metascape.clusters.out", sep='\t', col.names = T, row.names = F, quote = F)

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
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cell_cluster.asso.png", width=540/96, height=597/96, dpi=300)









################################################################################
#### find central (overlapped) genes within certain trait-community association
extract_intersecting_genes_and_cell_types <- function(fisher_test_df, traits_of_interest, community_of_interest) {
  
  # Extract the intersecting genes
  gene_lists <- fisher_test_df %>%
    filter(consensus_cluster == community_of_interest) %>%
    filter(trait %in% traits_of_interest) %>%
    select(trait, gene_list) %>%
    separate_rows(gene_list, sep = ", ") %>%
    group_by(trait) %>%
    summarise(genes = list(unique(gene_list)))
  
  intersecting_genes <- Reduce(intersect, gene_lists$genes)
  intersecting_genes <- intersecting_genes[intersecting_genes != ""]
  
  # Find the cell types associated with the intersecting genes
  detailed_info <- fisher_test_df %>%
    filter(consensus_cluster == community_of_interest) %>%
    filter(trait %in% traits_of_interest) %>%
    separate_rows(gene_list, sep = ", ") %>%
    filter(gene_list %in% intersecting_genes) %>%
    select(strain, treatment, tissue, cell_type, gene_list) %>%
    distinct() %>%
    arrange(trait, cell_type, gene_list)
  
  return(list(intersecting_genes = intersecting_genes, detailed_info = detailed_info))
}


traits_of_interest <- c("BUN", "eGFR", "systolic_bp")
community_of_interest <- "Cluster 2"
result <- extract_intersecting_genes_and_cell_types(fisher_test_df, traits_of_interest, community_of_interest)
print(result$intersecting_genes)
print(result$detailed_info)

traits_of_interest <- c("CAD", "diastolic_bp", "pulse_pressure", "systolic_bp")
community_of_interest <- 5
result <- extract_intersecting_genes_and_cell_types(fisher_test_df, traits_of_interest, community_of_interest)
print(result$intersecting_genes)
print(result$detailed_info)


traits_of_interest <- c("BUN")
community_of_interest <- 4
intersecting_genes <- extract_intersecting_genes(fisher_test_df, traits_of_interest, community_of_interest)
print(intersecting_genes)







deg_filtered_matching <- deg_filtered %>%
  filter(strain %in% result$detailed_info$strain,
         treatment %in% result$detailed_info$treatment,
         tissue %in% result$detailed_info$tissue,
         cell_type %in% result$detailed_info$cell_type,
         gene_name %in% result$intersecting_genes)
if(nrow(deg_filtered_matching)>0){
  p1=ggplot(deg_filtered_matching, aes(x = treatment, y = paste(tissue, cell_type, sep = "-"))) +
    geom_point(aes(fill = -avg_log2FC, size=-log10(p_val)), shape=21) +
    scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white",
                         high = "red", space = "Lab" ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),   # No horizontal grid lines
      legend.justification = c(1, 0.5),
      legend.position = "bottom",
      panel.spacing = unit(0.2, "lines"),
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      strip.text = element_text(colour = 'black'),
      strip.text.y = element_text(colour = 'black', angle = 360),
      strip.background = element_rect(colour = "black", fill = NA)
    ) +
    labs(x="", y="", fill="log2(FC)", size="-log10(p-value)") +
    facet_nested(gene_name ~ project + strain, scales = "free", space = "free")
  
  snp_trait_df = snp_gene_df_mod %>%
    merge(., gwas_merge_use, by.x="SNP", by.y="SNPS") %>%
    filter(combined_gene_name %in% result$intersecting_genes,
           trait %in% traits_of_interest) %>%
    group_by(Gene.name) %>% mutate(sum_p = sum(-log10(P.VALUE))) %>% 
    arrange(desc(sum_p))
  
  # top10_genes = head(unique(snp_trait_df$combined_gene_name), n=10)
  top10_genes = unique(snp_trait_df$combined_gene_name)
  
  top5_snps_to_genes = snp_trait_df %>%
    filter(combined_gene_name %in% top10_genes) %>% 
    group_by(SNP, combined_gene_name) %>% 
    mutate(snp_sum_p = sum(-log10(P.VALUE)), .groups = 'drop') %>% 
    ungroup() %>% 
    mutate(combined_gene_name = factor(combined_gene_name, levels = top10_genes)) %>% 
    group_by(combined_gene_name) %>% 
    arrange(combined_gene_name, desc(sum_p)) %>% 
    slice_head(n = 5)
  
  p2 = ggplot(top5_snps_to_genes, aes(x = trait, y = SNP)) +
    geom_tile(aes(fill = -log10(P.VALUE))) +
    scale_y_discrete(limits=rev, position = "right") +
    scale_fill_gradient(low = "white", high = "purple") +
    theme_bw() +
    labs(x="", y="", fill="-log10(p-value)", title = "Top SNPs") +
    facet_nested(Gene.name ~ ., scales = "free", space = "free", switch="y") +
    theme(
      # panel.grid.major.y = element_blank(),   # No horizontal grid lines
      # legend.justification = c(1, 0.5),
      # legend.title = element_text(size = 10),
      # legend.title.align = 0.5,
      legend.position = "bottom",
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      # strip.text = element_text(colour = 'black'),
      strip.background = element_rect(colour = "black", fill = NA),
      strip.text.y.left = element_text(colour = 'black', angle = 360)
    )

  combined_plot <- (p1 + p2) + plot_layout(guides = 'collect', widths = c(0.5, 0.2)) &
    theme(legend.position = 'bottom', 
          legend.box = 'vertical')
  print(combined_plot)
}else{
  p1=NULL
}







################################################################################
### community based multi-functional SNPs
fisher_test_df = read.table("trait.fisher.cluster.out", header = T, sep = "\t", comment.char = "")
snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
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
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)
# deg_filtered = deg_merged %>% filter(p_val_adj < 0.05)
deg_filtered$cell_group = paste(deg_filtered$strain, deg_filtered$treatment, deg_filtered$tissue, deg_filtered$cell_type, sep = "_")
consensus_clusters = fisher_test_df$consensus_cluster
names(consensus_clusters) = fisher_test_df$cell_group
deg_filtered$consensus_cluster = consensus_clusters[deg_filtered$cell_group]


# Usage example
plot_snp_analysis(snp_id = "rs11838776", cluster_name = "Cluster 4", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs179993", cluster_name = "Cluster 4", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs1032763", cluster_name = c("Cluster 1", "Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs1530440", cluster_name = "Cluster 2", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs7766596", cluster_name = "Cluster 2", 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = "rs3785880", cluster_name = c("Cluster 1", "Cluster 2"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs2953516"), cluster_name = c("Cluster 1", "Cluster 2"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs77924615"), cluster_name = c("Cluster 2"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)


plot_snp_analysis(snp_id = c("rs7920682"), cluster_name = c("Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs8102847"), cluster_name = c("Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)


plot_snp_analysis(snp_id = c("rs740406"), cluster_name = c("Cluster 2", "Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs1032763"), cluster_name = c("Cluster 2", "Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs4468717"), cluster_name = c("Cluster 2", "Cluster 3"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)



### experimental candidates
plot_snp_analysis(snp_id = c("rs28451064"), cluster_name = c("Cluster 2"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs2782980"), cluster_name = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs1882961"), cluster_name = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs1173771"), cluster_name = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)

plot_snp_analysis(snp_id = c("rs217727"), cluster_name = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6"), 
                  gwas_merge_use = gwas_merge_use, snp_gene_df = snp_gene_df, 
                  fisher_test_df = fisher_test_df, deg_filtered = deg_filtered)




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


################################################################################
### load LD information
LD_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/LD_infor.txt", header = T, sep = "\t")
snps_LD = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/SNPS_in_LD.txt", header = T, sep = "\t")

ld_1 = names(table(snps_LD$rsID1)[table(snps_LD$rsID1)==1])

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(p_val_adj < 0.05)

ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

snp_gene_df_mod = snp_gene_df %>%
  separate_rows(combined_gene_name , sep=";") %>%
  merge(., SNPS, by="SNP") %>%
  merge(., ensembl, by.x="gene", by.y="Gene.stable.ID") %>%
  mutate(Chromosome = sapply(strsplit(as.character(POS), "_"), `[`, 1),
         Position = as.numeric(sapply(strsplit(as.character(POS), "_"), `[`, 2)),
         TSS = ifelse(Strand == 1, Gene.start..bp., Gene.end..bp.),
         Distance_to_TSS = abs(TSS-Position)) %>%
  filter(! is.na(Distance_to_TSS)) %>%
  unique()

snp_gene_df_ld = snp_gene_df_mod %>%
  filter(SNP %in% ld_1, combined_gene_name %in% deg_filtered$gene_name)
length(unique(snp_gene_df_ld$SNP))
length(unique(snp_gene_df_ld$gene))


snp_gene_df_ld = snp_gene_df_mod %>%
  filter(SNP %in% ld_1, combined_gene_name %in% deg_filtered$gene_name, Distance_to_TSS>10000)
length(unique(snp_gene_df_ld$SNP))
length(unique(snp_gene_df_ld$gene))


# gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_filter = gwas_merge[gwas_merge$SNPS %in% snp_gene_df_ld$SNP, ] %>% 
  group_by(trait, SNPS) %>%
  arrange(P.VALUE) %>% slice_head(n=1)
gwas_merge_filter <- gwas_merge_filter[order(gwas_merge_filter$P.VALUE), ]
gwas_merge_filter$Rank <- rank(gwas_merge_filter$P.VALUE, ties.method = "min")
gwas_merge_filter[gwas_merge_filter$SNPS=="rs1882961", ]$Rank


gwas_merge_filter = gwas_merge %>% 
  filter(SNPS %in% snp_gene_df_ld$SNP, trait=="systolic_bp") %>% 
  group_by(trait, SNPS) %>%
  arrange(P.VALUE) %>% slice_head(n=1)
gwas_merge_filter <- gwas_merge_filter[order(gwas_merge_filter$P.VALUE), ]
gwas_merge_filter$Rank <- rank(gwas_merge_filter$P.VALUE, ties.method = "min")
gwas_merge_filter[gwas_merge_filter$SNPS=="rs1882961", ]$Rank






################################################################################
################################################################################
### load LD information
ensembl <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.GRCH38.out", sep = '\t', header=T)
rownames(ensembl) = ensembl$Gene.stable.ID

gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

LD_info = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/LD_infor.txt", header = T, sep = "\t")
snps_LD = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/SNPS_in_LD.txt", header = T, sep = "\t")
ld_1 = names(table(snps_LD$rsID1)[table(snps_LD$rsID1)==1])

snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
SNPS = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/gwas_catalog_bp_relevant.snp.txt", header=T, sep='\t')

### negative gene list
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("SD", "WKY"), ]
neg_list = deg_merged %>% filter(cell_type %in% c("EC", "VSMC") & p_val < 0.05)

deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]

snp_gene_df_mod = snp_gene_df %>%
  separate_rows(combined_gene_name , sep=";") %>%
  merge(., SNPS, by="SNP") %>%
  merge(., ensembl, by.x="gene", by.y="Gene.stable.ID") %>%
  mutate(Chromosome = sapply(strsplit(as.character(POS), "_"), `[`, 1),
         Position = as.numeric(sapply(strsplit(as.character(POS), "_"), `[`, 2)),
         TSS = ifelse(Strand == 1, Gene.start..bp., Gene.end..bp.),
         Distance_to_TSS = abs(TSS-Position)) %>%
  filter(! is.na(Distance_to_TSS))


deg_filtered = deg_merged %>% filter(cell_type %in% c("EC", "VSMC") & p_val_adj < 0.05) # 32 SNP; 25 gene
deg_filtered = deg_merged %>% filter(cell_type %in% c("EC", "VSMC") & p_val_adj < 0.05 & abs(avg_log2FC)>0.25) # 30 SNP; 23 gene

snp_gene_df_ld = snp_gene_df_mod %>%
  filter(SNP %in% ld_1, 
         SNP %in% gwas_merge_use[gwas_merge_use$trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"), ]$SNPS,
         combined_gene_name %in% deg_filtered$gene_name, 
         Distance_to_TSS>10000)
length(unique(snp_gene_df_ld$SNP))
length(unique(snp_gene_df_ld$gene))

gwas_merge_use[gwas_merge_use$SNPS %in% snp_gene_df_ld$SNP, ]



snp_gene_df_ld = snp_gene_df_mod %>%
  merge(., gwas_merge_use, by.x = "SNP", by.y="SNPS") %>%
  filter(SNP %in% ld_1, 
         trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"),
         Distance_to_TSS>10000) %>%
  merge(., enhanced_index, by.x = "combined_gene_name", by.y="gene_name")
  
snp_gene_df_ld = snp_gene_df_mod %>%
  merge(., gwas_merge_use, by.x = "SNP", by.y="SNPS") %>%
  filter(SNP %in% ld_1, 
         trait %in% c("diastolic_bp", "systolic_bp", "pulse_pressure"),
         Distance_to_TSS>10000) %>%
  merge(., combined_index, by.x = "combined_gene_name", by.y="gene_name") %>%
  group_by(SNP) %>% mutate(n_gene_per_SNP = n_distinct(gene)) %>% ungroup() %>%
  group_by(gene) %>% mutate(n_SNP_per_gene = n_distinct(SNP)) %>% ungroup()

length(unique(snp_gene_df_ld$SNP))
length(unique(snp_gene_df_ld$gene))

snp_gene_df_ld  = rbind(unique(snp_gene_df_ld[snp_gene_df_ld$combined_gene_name!="Nrip1",]), unique(snp_gene_df_ld[snp_gene_df_ld$combined_gene_name=="Nrip1",]))
snp_gene_df_ld %>% filter(n_SNP_per_gene==1, n_gene_per_SNP==1) %>%
  ggplot(., aes(x=-log10(P.VALUE), y=combined_score, color= combined_gene_name=="Nrip1")) + geom_point()

snp_gene_df_ld %>% filter(n_SNP_per_gene==1, n_gene_per_SNP==1) %>%
  ggplot(., aes(x=breadth_score, y=index_value, color= combined_gene_name=="Nrip1")) + geom_point()


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
                        grepl("Member", merge_all$GroupID), ] %>%
  filter(cluster %in% c("Cluster 3", "Cluster 5")) %>%
  separate_rows(Symbols , sep=",")

deg_filtered = deg_merged %>% filter(cell_type %in% c("EC", "VSMC") & p_val_adj < 0.05)
deg_filtered = deg_merged %>% filter(cell_type %in% c("EC", "VSMC") & p_val_adj < 0.05 & abs(avg_log2FC)>0.25)

snp_gene_df_ld = snp_gene_df_mod %>%
  filter(SNP %in% ld_1, 
         combined_gene_name %in% deg_filtered$gene_name, 
         !(combined_gene_name %in% neg_list$gene_name), 
         Distance_to_TSS>10000, 
         Gene.name.x %in% merge_use$Symbols)
length(unique(snp_gene_df_ld$SNP))
length(unique(snp_gene_df_ld$gene))











################################################################################


calculate_enhanced_index <- function(deg_data) {
  
  # Calculate the proportion of conditions with significant changes (consistency metric)
  consistency <- deg_data %>%
    group_by(gene_name) %>%
    summarise(consistency_metric = mean(p_val_adj < 0.05, na.rm = TRUE))

  # Normalize the fold changes across tissues/models using z-score transformation
  deg_data <- deg_data %>%
    group_by(tissue, strain) %>%
    mutate(z_score_fc = scale(-avg_log2FC)) %>%
    ungroup()
  
  # Calculate the standard deviation of z-scored fold changes across tissues/models
  deg_data <- deg_data %>%
    group_by(gene_name) %>%
    mutate(variance_fc = ifelse(n() == 1, 0, sd(z_score_fc, na.rm = TRUE))) %>%
    ungroup()
  
  # Join the consistency metric back to the main data
  deg_data <- deg_data %>%
    left_join(consistency, by = "gene_name")
  
  # Handling p_val_adj = 1 by replacing it with a small number slightly less than 1
  deg_data <- deg_data %>%
    mutate(adjusted_p_val = ifelse(p_val_adj == 1, 0.9999, p_val_adj))
  
  # Combine the normalized fold change, variance, and consistency into an overall index
  index <- deg_data %>%
    group_by(gene_name) %>%
    summarise(
      weighted_z_score_fc = sum(z_score_fc * -log10(adjusted_p_val), na.rm = TRUE) / sum(-log10(adjusted_p_val), na.rm = TRUE),
      robustness_penalty = 1 / (1 + variance_fc),
      consistency_metric = first(consistency_metric),
      # index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * robustness_penalty * consistency_metric)
      index_value = ifelse(is.na(weighted_z_score_fc) | is.na(robustness_penalty), NA, weighted_z_score_fc * robustness_penalty)
    ) %>%
    mutate(index_value = ifelse(is.nan(index_value), 0, index_value)) %>%
    unique()
  
  return(index)
}

snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
snp_gene_df_mod = snp_gene_df %>%
  separate_rows(combined_gene_name , sep=";")

# deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)
# deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
enhanced_index <- calculate_enhanced_index(deg_merged[deg_merged$cell_type %in% c("EC"), ])
enhanced_index <- calculate_enhanced_index(deg_merged[deg_merged$cell_type %in% c("VSMC"), ])
enhanced_index <- calculate_enhanced_index(deg_merged)

total_tissue_cell_combinations = nrow(unique(deg_merged[, c("tissue", "cell_type")]))
breadth_score <- deg_merged %>%
  group_by(gene_name) %>%
  summarise(breadth_score = n_distinct(tissue, cell_type) / total_tissue_cell_combinations)

combined_index <- enhanced_index %>%
  inner_join(breadth_score, by = "gene_name") %>%
  mutate(combined_score = index_value * breadth_score)

combined_index[combined_index$gene_name=="Nrip1",]


# View the results
head(enhanced_index)


# Create two subsets of the data
yong_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/HYPT_2020_Yong_bp_physiology_gene_list.mod.txt", header = T, sep='\t')
helen_list = read.table("/xdisk/mliang1/qqiu/data/gene_list/NG_2024_Helen_bp_pred_gene_list.txt", header = T, sep='\t')

overlap_genes <- enhanced_index[enhanced_index$gene_name %in% snp_gene_df_mod$combined_gene_name, ]$index_value
non_overlap_genes <- enhanced_index[!(enhanced_index$gene_name %in% snp_gene_df_mod$combined_gene_name), ]$index_value

overlap_genes <- enhanced_index[enhanced_index$gene_name %in% yong_list$Gene.name, ]$index_value
non_overlap_genes <- enhanced_index[!(enhanced_index$gene_name %in% yong_list$Gene.name), ]$index_value

# overlap_genes <- enhanced_index[enhanced_index$gene_name %in% helen_list$combined_gene_name, ]$index_value
# non_overlap_genes <- enhanced_index[!(enhanced_index$gene_name %in% helen_list$combined_gene_name), ]$index_value
# Remove NA values
overlap_genes <- na.omit(overlap_genes)
non_overlap_genes <- na.omit(non_overlap_genes)


plot_data <- data.frame(
  index_value = c(overlap_genes, non_overlap_genes),
  group = factor(c(rep("Overlap", length(overlap_genes)), rep("Non-Overlap", length(non_overlap_genes))))
)

# Plot cumulative density
ggplot(plot_data, aes(x = index_value, color = group)) +
  stat_ecdf(size = 1.5) +
  labs(title = "Cumulative Density Plot of Index Values",
       x = "Index Value",
       y = "Cumulative Density") +
  theme_minimal() +
  theme(legend.title = element_blank())

ks_test <- ks.test(overlap_genes, non_overlap_genes)
print(ks_test)







# Manually define bins to capture both extreme and near-zero values
enhanced_index <- enhanced_index %>%
  mutate(index_group = cut(index_value, 
                           breaks = c(-Inf, -0.5, -0.1, 0.1, 0.5, Inf), 
                           labels = c("Strongly Negative", "Moderately Negative", "Neutral", "Moderately Positive", "Strongly Positive")))

# Verify the distribution across bins
table(enhanced_index$index_group)
# Calculate Proportions of Overlap and Non-Overlap in Each Group
proportion_data <- enhanced_index %>%
  mutate(overlap = ifelse(gene_name %in% snp_gene_df_mod$combined_gene_name, "Overlap", "Non-Overlap")) %>%
  group_by(index_group, overlap) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()

# Visualize the Proportions
ggplot(proportion_data, aes(x = index_group, y = proportion, fill = overlap)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Proportion of Overlap and Non-Overlap Genes Across Index Value Groups",
       x = "Index Value Group",
       y = "Proportion") +
  theme_minimal()








################################################################################
### eqtl - allele specific info
fisher_test_df = read.table("trait.fisher.cluster.out", header = T, sep = "\t", comment.char = "")
snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
eqtl_merge <- read.table("eqtl_merged.txt", header = TRUE, sep = '\t')
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.all.out", sep = '\t', header = TRUE)

gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  mutate(risk_allele = str_extract(STRONGEST.SNP.RISK.ALLELE, "[A-Za-z]+$")) %>%
  filter(risk_allele %in% c("A", "T", "C" ,"G")) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE, STRONGEST.SNP.RISK.ALLELE, risk_allele, OR.or.BETA, X95..CI..TEXT.)

##
# Normalized effect size (NES), previously known as the effect size on the portal, 
# is defined as the slope of the linear regression, and is computed as the effect of the 
# alternative allele (ALT) relative to the reference allele (REF) in the human genome reference 
# (i.e., the eQTL effect allele is the ALT allele). NES are computed in a normalized space where magnitude has no direct biological interpretation.
eqtl_merge <- eqtl_merge %>%
  filter(grepl("rs", SNP)) %>%
  separate(variant_id, into = c("chr", "pos", "ref_allele", "alt_allele", "build"), sep = "_") %>%
  mutate(higher_expression_allele = ifelse(slope > 0, alt_allele, ref_allele),
         combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  select(gene_id_mod, POS, maf, pval_nominal, slope, gene_name.x, tissue, SNP, 
         gene_id_rat, gene_name_rat, r2h_orthology_conf, gene_id_mouse, gene_name_mouse, 
         m2h_orthology_conf, higher_expression_allele, ref_allele, alt_allele, combined_gene_name)

# deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(strain %in% c("C57BL/6", "SS", "SHR") & p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

merged_data <- gwas_merge_use %>%
  left_join(eqtl_merge, by = c("SNPS" = "SNP"), relationship = "many-to-many") %>%
  left_join(deg_filtered, by = c("combined_gene_name" = "gene_name"), relationship = "many-to-many") %>%
  filter(risk_allele %in% c(ref_allele, alt_allele))

merged_data <- merged_data %>% ungroup() %>%
  filter(complete.cases(.)) %>%
  mutate(effect_direction = ifelse(grepl("decrease", X95..CI..TEXT.), "Protective effect to trait", "Adverse effect to trait")) %>%
  mutate(consistency = ifelse((effect_direction == "Adverse effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC < 0) | 
                                                                   (risk_allele != higher_expression_allele & avg_log2FC > 0))) |
                                (effect_direction == "Protective effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC > 0) | 
                                                                     (risk_allele != higher_expression_allele & avg_log2FC < 0))), 
                              "Consistent", "Inconsistent"),
         adjusted_slope = ifelse(risk_allele == alt_allele, slope, -1 * slope),
         expr_effect = ifelse(risk_allele == higher_expression_allele, "Increase", "Decrease")) %>% unique()
merged_data = merged_data[rowSums(is.na(merged_data))==0,]
table(merged_data$consistency)
table(interaction(merged_data$tissue.x, merged_data$tissue.y), merged_data$consistency)

consistent_result = merged_data[rowSums(is.na(merged_data))==0 & merged_data$consistency=="Consistent",]
write.table(consistent_result, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/GWAS/eqtl.consistent_result.out", sep="\t", quote=F, col.names=T, row.names=T)




### filter by the most significant results > not good
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  mutate(risk_allele = str_extract(STRONGEST.SNP.RISK.ALLELE, "[A-Za-z]+$")) %>%
  filter(risk_allele %in% c("A", "T", "C" ,"G")) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE, STRONGEST.SNP.RISK.ALLELE, risk_allele, OR.or.BETA, X95..CI..TEXT.)

eqtl_merge <- eqtl_merge %>%
  filter(grepl("rs", SNP)) %>%
  mutate(variant_id_rep = variant_id) %>%
  separate(variant_id_rep, into = c("chr", "pos", "ref_allele", "alt_allele", "build"), sep = "_") %>%
  mutate(higher_expression_allele = ifelse(slope > 0, alt_allele, ref_allele),
         combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  group_by(variant_id, gene_id_mod) %>%
  arrange(pval_nominal) %>%
  slice_head(n = 1) %>%
  select(gene_id_mod, POS, maf, pval_nominal, slope, gene_name.x, tissue, SNP, variant_id, 
         gene_id_rat, gene_name_rat, r2h_orthology_conf, gene_id_mouse, gene_name_mouse, 
         m2h_orthology_conf, higher_expression_allele, ref_allele, alt_allele, combined_gene_name)

# deg_merged <- deg_merged[deg_merged$strain %in% c("C57BL/6", "SHR", "SS"), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
# deg_merged <- deg_merged[!(deg_merged$tissue == "MCA" & deg_merged$strain %in% c("C57BL/6", "SS")), ]
deg_filtered = deg_merged %>% filter(strain %in% c("C57BL/6", "SS", "SHR") & p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  group_by(gene_name) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 1)
  

merged_data <- gwas_merge_use %>%
  left_join(eqtl_merge, by = c("SNPS" = "SNP"), relationship = "many-to-many") %>%
  left_join(deg_filtered, by = c("combined_gene_name" = "gene_name"), relationship = "many-to-many") %>%
  filter(risk_allele %in% c(ref_allele, alt_allele))

merged_data <- merged_data %>% ungroup() %>%
  filter(complete.cases(.)) %>%
  mutate(effect_direction = ifelse(grepl("decrease", X95..CI..TEXT.), "Protective effect to trait", "Adverse effect to trait")) %>%
  mutate(consistency = ifelse((effect_direction == "Adverse effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC < 0) | 
                                                                                  (risk_allele != higher_expression_allele & avg_log2FC > 0))) |
                                (effect_direction == "Protective effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC > 0) | 
                                                                                       (risk_allele != higher_expression_allele & avg_log2FC < 0))), 
                              "Consistent", "Inconsistent"),
         adjusted_slope = ifelse(risk_allele == alt_allele, slope, -1 * slope),
         expr_effect = ifelse(risk_allele == higher_expression_allele, "Increase", "Decrease")) %>% unique()

table(merged_data$consistency)
table(merged_data$consistency, merged_data$trait)


### weighted methods
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  mutate(risk_allele = str_extract(STRONGEST.SNP.RISK.ALLELE, "[A-Za-z]+$")) %>%
  filter(risk_allele %in% c("A", "T", "C" ,"G")) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE, STRONGEST.SNP.RISK.ALLELE, risk_allele, OR.or.BETA, X95..CI..TEXT.)

eqtl_mod <- eqtl_merge %>%
  filter(grepl("rs", SNP)) %>%
  mutate(variant_id_rep = variant_id) %>%
  separate(variant_id_rep, into = c("chr", "pos", "ref_allele", "alt_allele", "build"), sep = "_") %>%
  mutate(higher_expression_allele = ifelse(slope > 0, alt_allele, ref_allele),
         combined_gene_name = coalesce(gene_name_rat, gene_name_mouse)) %>%
  # group_by(variant_id, gene_id_mod) %>%
  # mutate(weight = -log10(pval_nominal)) %>%
  # mutate(slope = sum(slope * weight) / sum(weight)) %>%
  select(gene_id_mod, POS, maf, pval_nominal, slope, gene_name.x, tissue, SNP, variant_id, 
         gene_id_rat, gene_name_rat, r2h_orthology_conf, gene_id_mouse, gene_name_mouse, 
         m2h_orthology_conf, higher_expression_allele, ref_allele, alt_allele, combined_gene_name)

deg_filtered = deg_merged %>% filter(strain %in% c("C57BL/6", "SS", "SHR") & p_val_adj < 0.05 & abs(avg_log2FC) > 0.25) %>%
  group_by(gene_name, tissue) %>%
  mutate(weight = -log10(p_val)) %>%
  summarise(avg_log2FC = sum(avg_log2FC * weight) / sum(weight))

merged_data <- gwas_merge_use %>%
  left_join(eqtl_mod, by = c("SNPS" = "SNP"), relationship = "many-to-many") %>%
  left_join(deg_filtered, by = c("combined_gene_name" = "gene_name"), relationship = "many-to-many") %>%
  filter(risk_allele %in% c(ref_allele, alt_allele))

merged_data <- merged_data %>% ungroup() %>%
  filter(complete.cases(.)) %>%
  mutate(effect_direction = ifelse(grepl("decrease", X95..CI..TEXT.), "Protective effect to trait", "Adverse effect to trait")) %>%
  mutate(consistency = ifelse((effect_direction == "Adverse effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC < 0) | 
                                                                                  (risk_allele != higher_expression_allele & avg_log2FC > 0))) |
                                (effect_direction == "Protective effect to trait" & ((risk_allele == higher_expression_allele & avg_log2FC > 0) | 
                                                                                       (risk_allele != higher_expression_allele & avg_log2FC < 0))), 
                              "Consistent", "Inconsistent"),
         adjusted_slope = ifelse(risk_allele == alt_allele, slope, -1 * slope),
         expr_effect = ifelse(risk_allele == higher_expression_allele, "Increase", "Decrease")) %>% unique()
merged_data = merged_data[rowSums(is.na(merged_data))==0,]

table(merged_data$consistency)
table(merged_data$consistency, merged_data$effect_direction)
table(merged_data$consistency, merged_data$trait)
tmp = table(interaction(merged_data$tissue.x, merged_data$tissue.y), merged_data$consistency)
tmp




consistent_result %>% filter(strain %in% c("SD", "WKY")) %>%
  select(SNPS, gene_id_mod) %>% unique() %>% summarise(snp_gene_pair = n())

consistent_result %>% filter(strain %in% c("SD", "WKY")) %>%
  select(SNPS, trait, gene_id_mod, effect_direction) %>% unique() %>% group_by(effect_direction) %>% summarise(snp_gene_pair = n())

consistent_result %>% filter(strain %in% c("SD", "WKY")) %>%
  group_by(trait) %>% summarize(eqtl_consistent_snp = n_distinct(SNPS),
                                eqtl_consistent_snp_gene_pair =  n_distinct(SNPS, gene_id_mod)) %>% 
  merge(., gwas_merge %>% group_by(trait) %>% summarize(all_snp = n_distinct(SNPS)), by="trait") %>% 
  mutate(eqtl_consistent_snp_prop = eqtl_consistent_snp/all_snp)


gene_effect_summary <- merged_data %>%
group_by(combined_gene_name) %>%
summarise(
  gene_effect = case_when(
    all(consistency == "Consistent" & effect_direction == "increase") ~ "Adverse",
    all(consistency == "Consistent" & effect_direction == "decrease") ~ "Protective",
    all(consistency == "Consistent" & effect_direction == "disease_risk") ~ "Adverse",
    TRUE ~ "Mixed"
  )
)

merged_data <- merged_data %>%
  left_join(gene_effect_summary, by = "combined_gene_name")




table(merged_data$consistency)
# 
# Consistent Inconsistent 
# 23794        18553 
table(merged_data[rowSums(is.na(merged_data))==0,]$gene_effect, merged_data[rowSums(is.na(merged_data))==0,]$effect_direction)
# 
# decrease disease_risk increase
# Adverse           0           31      433
# Mixed         18717         2159    19018
# Protective      253            0        0




eqtl_weighted_avg <- consistent_result %>%
  group_by(SNPS, gene_id_mod, risk_allele, trait) %>%
  unique() %>%
  summarise(weighted_avg_slope = sum(adjusted_slope * (1 / pval_nominal)) / sum(1 / pval_nominal)) %>%
  filter(! is.na(gene_id_mod))

merged_data <- consistent_result %>%
  left_join(eqtl_weighted_avg, by = c("SNPS", "gene_id_mod", "risk_allele", "trait"), relationship = "many-to-many")

expression_weighted_avg <- consistent_result %>%
  group_by(SNPS, gene_id_mod, risk_allele) %>%
  unique() %>%
  summarise(weighted_avg_log2FC = sum(avg_log2FC * (1 / p_val_adj)) / sum(1 / p_val_adj))

# expression_weighted_avg <- consistent_result %>%
#   mutate(cell_group = paste(strain, treatment, tissue, cell_type, sep="_")) %>%
#   merge(., fisher_test_df[c("cell_group", "consensus_cluster")], by = "cell_group") %>%
#   group_by(gene_name, consensus_cluster) %>%
#   summarise(weighted_avg_log2FC = sum(avg_log2FC * (1 / p_val_adj)) / sum(1 / p_val_adj)) %>%
#   rename(gene = gene_name)

# Integrate weighted average expression change with the main dataset
merged_data <- merged_data %>%
  left_join(expression_weighted_avg, by = c("SNPS", "gene_id_mod", "risk_allele"), relationship = "many-to-many")


min_nonzero_pval <- sort(unique(merged_data$P.VALUE[merged_data$P.VALUE > 0]))[1]
merged_data <- merged_data %>%
  mutate(P.VALUE = ifelse(P.VALUE == 0, min_nonzero_pval, P.VALUE))
cap_threshold <- 50  # Set a reasonable threshold for the cap
merged_data <- merged_data %>%
  mutate(capped_log10_pvalue = ifelse(-log10(P.VALUE) > cap_threshold, cap_threshold, -log10(P.VALUE)))
merged_data <- merged_data %>%
  filter(!is.na(P.VALUE) & !is.na(weighted_avg_log2FC))


# Create a new column to categorize the effect_direction for x-axis positioning
merged_data <- merged_data %>%
  mutate(effect_category = case_when(
    effect_direction == "increase" | effect_direction == "disease_risk" ~ "Adeverse SNP effect",
    effect_direction == "decrease" ~ "Favorable SNP effect",
    TRUE ~ "Other"
  ))
# 
# merged_data <- merged_data %>%
#   mutate(adjusted_x = ifelse(effect_category == "Right", log10(P.VALUE), -log10(P.VALUE)))

data_plot = merged_data[, c("SNPS", "risk_allele", "effect_category", "gene_id_mod", "P.VALUE", "weighted_avg_log2FC", "weighted_avg_slope")]
ggplot(data_plot, aes(x = weighted_avg_slope, y = weighted_avg_log2FC, fill = -log10(P.VALUE))) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  labs(x = "-log10(p-value)", y = "Weighted average log2FC", fill = "Weighted\naverage\nslope") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  ) +
  facet_grid( ~ effect_category, scales = "free")
  # scale_x_continuous(limits = c(-max(-log10(merged_data$P.VALUE)), max(-log10(merged_data$P.VALUE))),
  #                    breaks = seq(-max(-log10(merged_data$P.VALUE)), max(-log10(merged_data$P.VALUE)), by = 2))






table(merged_data[rowSums(is.na(merged_data))==0,]$gene_effect, merged_data[rowSums(is.na(merged_data))==0,]$effect_direction)




snp_eqtl_summary <- merged_data %>%
  group_by(trait, SNPS) %>%
  summarise(has_eqtl = any(!is.na(higher_expression_allele)),
            consistent = any(!is.na(consistency) & consistency == "Consistent")) %>%
  group_by(trait) %>%
  summarise(
    total_snps = n(),
    snps_with_eqtl = sum(has_eqtl),
    consistent_snps = sum(consistent)
  )

gene_effect_summary <- merged_data %>%
  filter(!is.na(combined_gene_name) | combined_gene_name != "") %>%
  group_by(trait, combined_gene_name) %>%
  summarise(gene_effect = unique(gene_effect)[!is.na(unique(gene_effect))]) %>%
  group_by(trait) %>%
  summarise(
    total_genes = n(),
    adverse_genes = sum(gene_effect == "Adverse", na.rm = TRUE),
    protective_genes = sum(gene_effect == "Protective", na.rm = TRUE)
  )

combined_summary <- left_join(snp_eqtl_summary, gene_effect_summary, by = "trait")
print(combined_summary)




genes_by_trait <- merged_data %>%
  filter(!is.na(gene_effect) & combined_gene_name != "") %>%
  group_by(trait, gene_effect) %>%
  summarise(genes = paste(unique(combined_gene_name), collapse = ", ")) %>%
  ungroup()

genes_by_trait_readable <- genes_by_trait %>%
  pivot_wider(names_from = gene_effect, values_from = genes, values_fill = list(genes = "None")) %>%
  arrange(trait)

print(genes_by_trait_readable)

genes_by_trait_readable[,-2]

################################################################################
### community based multi-functional SNPs
fisher_test_df = read.table("trait.fisher.cluster.out", header = T, sep = "\t", comment.char = "")
snp_gene_df = read.table("snp_gene.evi_org.out", header = T, sep = "\t", comment.char = "")
gwas_merge <- read.table("gwas_merged.txt", header = TRUE, sep = '\t', quote = "")
gwas_merge_use <- gwas_merge %>%
  filter(grepl("rs", SNPS)) %>%
  group_by(SNPS, trait) %>%
  arrange(P.VALUE) %>%
  slice_head(n = 1) %>%
  select(SNPS, trait, P.VALUE)

fisher_test_df_expanded <- fisher_test_df %>%
  separate_rows(gene_list, sep = ", ") %>%
  rename(gene = gene_list)

snp_gene_df_expanded <- snp_gene_df %>%
  separate_rows(combined_gene_name, sep = ";")

merged_df <- snp_gene_df_expanded %>%
  left_join(fisher_test_df_expanded %>% select(-trait) %>% unique(), by = c("combined_gene_name" = "gene"), keep = TRUE)

overall_proportion <- fisher_test_df_expanded %>%
  group_by(consensus_cluster) %>%
  summarise(cluster_genes = n()) %>%
  mutate(total_genes = sum(cluster_genes),
         proportion = cluster_genes / sum(cluster_genes)) %>%
  ungroup()

snp_gene_proportion <- merged_df %>%
  group_by(SNP) %>%
  mutate(linked_genes = n_distinct(combined_gene_name)) %>%
  ungroup() %>% filter(! is.na(strain)) %>%
  group_by(SNP, consensus_cluster) %>%
  mutate(degs = n_distinct(combined_gene_name),
         deg_prop = degs/linked_genes) %>%
  ungroup() %>% unique()

# Step 5: Perform statistical tests
snp_cluster_significance <- snp_gene_proportion %>%
  left_join(overall_proportion, by = "consensus_cluster") %>%
  mutate(expected_genes = linked_genes * proportion,
         p_value = mapply(function(obs, exp) chisq.test(c(obs, exp))$p.value, degs, expected_genes))

# Step 6: Compute weighted enrichment score
snp_cluster_significance <- snp_cluster_significance %>%
  mutate(weighted_enrichment_score = (linked_genes / total_genes) * (deg_genes / linked_genes))

# Step 7: Highlight significant SNP-cluster associations
significant_snp_cluster <- snp_cluster_significance %>%
  filter(p_value < 0.05)

# Display the results
print(significant_snp_cluster)


View(unique(snp_cluster_significance[snp_cluster_significance$deg_prop>0.3 & snp_cluster_significance$degs>2 & snp_cluster_significance$consensus_cluster=="Cluster 4", c("SNP", "consensus_cluster")]))

View(unique(snp_cluster_significance[snp_cluster_significance$deg_prop>0.5 & snp_cluster_significance$degs>2 & snp_cluster_significance$consensus_cluster=="Cluster 2", c("SNP", "consensus_cluster")]))

View(unique(snp_cluster_significance[snp_cluster_significance$deg_prop>0.5 & snp_cluster_significance$degs>2 & snp_cluster_significance$consensus_cluster=="Cluster 3", c("SNP", "consensus_cluster")]))
