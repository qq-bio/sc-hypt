
library(dplyr)
library(tidyr)
library(vegan)
library(ggplot2)
library(gtools)

combine_pvalues <- function(row, logged=TRUE) {
  p_values <- na.omit(as.numeric(row))

  if(logged){
    chi_sq_statistic <- -2 * sum(p_values/log10(2))
  }else{
    chi_sq_statistic <- -2 * sum(log(p_values))
  }
  
  combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p_value)
}


################################################################################
### merge clusterprofiler results

merge_all = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG.0.5.withot_background.out", header = T, check.names = F)

merge_all$cell_type = factor(merge_all$cell_type, levels = cell_order)
merge_all$tissue = factor(merge_all$tissue, levels = tissue_order)
merge_all$path_cell = paste(merge_all$tissue, merge_all$cell_type,  merge_all$Description, sep="_")
merge_all$sxt = paste(merge_all$strain, merge_all$treatment, sep="-")

merge_background = merge_all[merge_all$strain %in% c("SD", "WKY"), ]
path_rm_list = merge_background[merge_background$p.adjust<0.05, ]$path_cell
merge_filter = merge_all[! merge_all$path_cell %in% path_rm_list, ]

# merge_summary = merge_filter[grepl("Summary", merge_filter$GroupID), ]
# top10 <- merge_summary[merge_summary$`Log(q-value)`< log10(0.05),]
# top20 <- merge_summary[merge_summary$`Log(q-value)`< log10(0.05),] %>% group_by(species, treatment, tissue, cell_type) %>% top_n(n = 20) %>% as.data.frame()
# top_list <- top10
# merge_filter2 = merge_filter[merge_filter$path_cell %in% top_list$path_cell, ]

merge_filter2 = merge_filter
merge_filter2$`Log(q-value)` = -log10(merge_filter2$qvalue)

merge_reshape = reshape(merge_filter2[,c("sxt", "path_cell", "qvalue")], idvar=c("path_cell"), timevar = "sxt", v.names = "qvalue", direction = "wide")
merge_reshape$`Log(comb.qval)` = log10(apply(merge_reshape[,2:6], 1, function(x) combine_pvalues(x, logged=FALSE)))

merge_reshape[,grepl("Log", colnames(merge_reshape))] = -1 * merge_reshape[,grepl("Log", colnames(merge_reshape))]
merge_reshape[, c("tissue", "cell_type", "pathway")] = do.call(rbind, strsplit(merge_reshape$path_cell, "_"))

merge_reshape = merge_reshape <- merge_reshape %>%
  arrange(desc(`Log(comb.qval)`)) %>%
  group_by(tissue, cell_type) %>%
  dplyr::mutate(id = row_number(),
                highlight="no") %>%
  ungroup() %>%
  arrange(tissue, cell_type) %>%  # Order by tissue and cell_type
  as.data.frame()

merge_reshape[merge_reshape=="NA"] <- NA
merge_reshape = na.replace(merge_reshape, 1)
merge_reshape$AngII = apply(merge_reshape[,2:3], 1, function(x) ifelse(any(as.numeric(x) < 0.05), 1, 0))
merge_reshape$SS = apply(merge_reshape[,6:7], 1, function(x) ifelse(any(as.numeric(x) < 0.05), 1, 0))
merge_reshape$SHR = ifelse(as.numeric(merge_reshape[,5]) < 0.05, 1, 0)

merge_reshape <- merge_reshape %>%
  mutate(category = case_when(
    AngII == 1 & SS == 0 & SHR == 0 ~ "Only AngII",
    AngII == 0 & SS == 1 & SHR == 0 ~ "Only SS",
    AngII == 0 & SS == 0 & SHR == 1 ~ "Only SHR",
    AngII == 1 & SS == 1 & SHR == 0 ~ "AngII & SS",
    AngII == 1 & SS == 0 & SHR == 1 ~ "AngII & SHR",
    AngII == 0 & SS == 1 & SHR == 1 ~ "SS & SHR",
    AngII == 1 & SS == 1 & SHR == 1 ~ "AngII & SS & SHR",
    TRUE ~ "Other"
  ))

merge_reshape$module_count = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
merge_reshape$cell_count = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$pathway==x))
merge_reshape$pathway_sum = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$module_count[merge_reshape$pathway==x]))

write.table(merge_reshape, "clusterprofiler.merge.reshape.out", sep='\t', quote=F, col.names=T, row.names=F)



















# ################################################################################
# ### merge metascape results
# combine_pvalues <- function(row) {
#   p_values <- na.omit(as.numeric(row))
#   # chi_sq_statistic <- -2 * sum(log(p_values))
#   chi_sq_statistic <- -2 * sum(p_values/log10(2))
#   combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
#   return(combined_p_value)
# }
# 
# merge_all = read.table("metascape.merge_all.out", sep='\t', header = T, quote = "", check.names = F)
# merge_all$cell_type = factor(merge_all$cell_type, levels = cell_order)
# merge_all$tissue = factor(merge_all$tissue, levels = tissue_order)
# merge_all$path_cell = paste(merge_all$tissue, merge_all$cell_type,  merge_all$pathway, sep="_")
# merge_all$sxt = paste(merge_all$species, merge_all$treatment, sep="-")
# 
# merge_background = merge_all[merge_all$species %in% c("SD", "WKY"), ]
# path_rm_list = merge_background[merge_background$LogP<log10(0.05), ]$path_cell
# merge_filter = merge_all[! merge_all$path_cell %in% path_rm_list, ]
# 
# merge_summary = merge_filter[grepl("Summary", merge_filter$GroupID), ]
# merge_summary = merge_filter[grepl("Member", merge_filter$GroupID), ]
# top10 <- merge_summary[merge_summary$`Log(q-value)`< log10(0.05),]
# top20 <- merge_summary[merge_summary$`Log(q-value)`< log10(0.05),] %>% group_by(species, treatment, tissue, cell_type) %>% top_n(n = 20) %>% as.data.frame()
# top_list <- top10
# 
# merge_filter2 = merge_filter[merge_filter$path_cell %in% top_list$path_cell, ]
# merge_reshape = reshape(merge_filter2[,c("sxt", "path_cell", "Log(q-value)")], idvar=c("path_cell"), timevar = "sxt", v.names = "Log(q-value)", direction = "wide")
# merge_reshape$`Log(comb.qval)` = log10(apply(merge_reshape[,2:6], 1, combine_pvalues))
# 
# merge_reshape[,grepl("Log", colnames(merge_reshape))] = -1 * merge_reshape[,grepl("Log", colnames(merge_reshape))]
# merge_reshape[, c("tissue", "cell_type", "pathway")] = do.call(rbind, strsplit(merge_reshape$path_cell, "_"))
# 
# merge_reshape = merge_reshape <- merge_reshape %>%
#   arrange(desc(`Log(comb.qval)`)) %>%
#   group_by(tissue, cell_type) %>%
#   dplyr::mutate(id = row_number(),
#                 highlight="no") %>%
#   ungroup() %>%
#   arrange(tissue, cell_type) %>%  # Order by tissue and cell_type
#   as.data.frame()
# 
# merge_reshape[merge_reshape=="NA"] <- NA
# merge_reshape = na.replace(merge_reshape, 0)
# merge_reshape$AngII = apply(merge_reshape[,2:3], 1, function(x) ifelse(any(as.numeric(x) > -log10(0.05)), 1, 0))
# merge_reshape$SS = apply(merge_reshape[,5:6], 1, function(x) ifelse(any(as.numeric(x) > -log10(0.05)), 1, 0))
# merge_reshape$SHR = ifelse(as.numeric(merge_reshape[,4]) > -log10(0.05), 1, 0)
# 
# merge_reshape <- merge_reshape %>%
#   mutate(category = case_when(
#     AngII == 1 & SS == 0 & SHR == 0 ~ "Only AngII",
#     AngII == 0 & SS == 1 & SHR == 0 ~ "Only SS",
#     AngII == 0 & SS == 0 & SHR == 1 ~ "Only SHR",
#     AngII == 1 & SS == 1 & SHR == 0 ~ "AngII & SS",
#     AngII == 1 & SS == 0 & SHR == 1 ~ "AngII & SHR",
#     AngII == 0 & SS == 1 & SHR == 1 ~ "SS & SHR",
#     AngII == 1 & SS == 1 & SHR == 1 ~ "AngII & SS & SHR",
#     TRUE ~ "Other"
#   ))
# 
# merge_reshape$module_count = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
# merge_reshape$cell_count = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$pathway==x))
# merge_reshape$pathway_sum = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$module_count[merge_reshape$pathway==x]))
# 
# write.table(merge_reshape, "metascape.merge_summary.reshape.out", sep='\t', quote=F, col.names=T, row.names=F)
# 








################################################################################
### visualiza summarized results

merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/clusterprofiler.merge.reshape.out", header = T, sep = '\t', quote = "")

path_count = merge_reshape %>% filter(category != "Other") %>% 
  group_by(tissue, cell_type, category) %>%
  dplyr::summarise(count = n()) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)
path_count$category = factor(path_count$category, levels = c("Only AngII", "Only SS", "Only SHR", "AngII & SS",
                                                             "AngII & SHR", "SS & SHR", "AngII & SS & SHR"))
path_count$cell_type_mod = factor(paste0(path_count$cell_type, " •"), levels = paste0(cell_order, " ●"))
cell_col_mod = cell_col
names(cell_col_mod) = paste0(names(cell_col), " ●")
ggplot(path_count, aes(y = cell_type, x = category, label = count)) +
  geom_point(aes(size=count/3), alpha = .5, color="red") + geom_text() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  scale_y_discrete(limits=rev) +
  scale_size_continuous(breaks = c(1, 3, 6),
                        labels = c("3", "9", "18")) +
  labs(x="", y="", size="Number of\npathways") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free")

### dot for visualization
ggplot(path_count, aes(y = cell_type_mod, x = category, label = count)) +
  geom_text(aes(y=cell_type_mod, label = cell_type_mod, color = cell_type_mod), x=0.05, hjust = 1, size=16)+
  lims(x=c(-1,1)) +
  theme_void() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "None") +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=cell_col_mod) +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free")










































# enrich_res_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.without_background.out")
# enrich_res_merged_filter = enrich_res_merged[enrich_res_merged$p.adjust<0.05, ]
# binary_matrix <- enrich_res_merged_filter %>%
#   unite("combo", c(species, tissue, cell_type, treatment), remove = FALSE) %>%
#   mutate(presence = 1) %>%
#   dplyr::select(combo, Description, presence) %>%
#   distinct() %>%
#   spread(key = Description, value = presence, fill = 0)


enrich_res_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.out", sep='\t', header = T, quote = "", check.names = F)
enrich_res_merged = enrich_res_merged[!(enrich_res_merged$tissue == "MCA" & enrich_res_merged$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")), ]
enrich_res_merged = enrich_res_merged[!(enrich_res_merged$tissue == "MCA" & enrich_res_merged$species %in% c("C57BL/6", "SS")), ]
enrich_res_merged = enrich_res_merged[!(enrich_res_merged$species %in% c("Salt-sensitive", "Spontaneous")), ]
enrich_res_merged_filter = enrich_res_merged[enrich_res_merged$`Log(q-value)`<log10(0.05) & grepl("Member", enrich_res_merged$GroupID), ]
binary_matrix <- enrich_res_merged_filter %>%
  tidyr::unite("combo", c(species, tissue, cell_type, treatment), remove = FALSE) %>%
  mutate(presence = 1) %>%
  dplyr::select(combo, pathway, presence) %>%
  distinct() %>%
  tidyr::pivot_wider(names_from = pathway, values_from = presence, values_fill = 0)

jaccard_dist <- vegdist(binary_matrix[, -1], method = "jaccard")
mds_result <- metaMDS(jaccard_dist, k = 2)  # k = 2 for 2-dimensional MDS
plot(mds_result)

# outlier_indices <- which(abs(mds_result$points[, 2]) > 6000)
# outliers <- binary_matrix[outlier_indices, ]
# jaccard_dist <- vegdist(binary_matrix[-outlier_indices, -1], method = "jaccard")
# mds_result <- metaMDS(jaccard_dist, k = 2)  # k = 2 for 2-dimensional MDS
# plot(mds_result)
# 
# # Prepare data for plotting
# mds_points <- as.data.frame(mds_result$points)
# mds_points$combo <- sapply(binary_matrix[-outlier_indices, 1], as.character)
# mds_points[, c("species", "tissue", "cell_type", "treatment")] = do.call(rbind, strsplit(mds_points$combo, "_"))
# mds_points$cell_type = factor(mds_points$cell_type, levels = cell_order)
# cell_type_vis = names(table(mds_points$cell_type))[table(mds_points$cell_type)>1]
# # cell_type_vis = setdiff(mds_points$cell_type, c("IC", "Adipocytes", "POD", "NFO", "Ependymal cells", "DCT/CT", "DCT"))
# species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))

# Plot
ggplot(mds_points[mds_points$cell_type %in% cell_type_vis, ], aes(x = MDS1, y = MDS2, color=species, shape=tissue)) +
  geom_point(size=3) +
  # geom_text(vjust = 1.5, hjust = 0.5) +
  theme_bw() +
  scale_color_manual(values = species_col) + 
  labs(x = "MDS Dimension 1", y = "MDS Dimension 2", 
       color = "Strain", shape = "Tissue") + 
  facet_wrap(. ~ cell_type, ncol=6, scales = "fixed")


library(pheatmap)
pheatmap(as.matrix(1-jaccard_dist), clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")





anno_df = as.data.frame(sapply(binary_matrix[,1], as.character))
anno_df[, c("strain", "tissue", "cell_type", "treatment")] = do.call(rbind, strsplit(anno_df$combo, "_"))

strain_colors <- species_col
tissue_colors <- tissue_col
cell_colors <- cell_col

top_annotation <- HeatmapAnnotation(
  Strain = anno_df$strain,
  Tissue = anno_df$tissue,
  Cell = anno_df$cell_type,
  col = list(
    Strain = strain_colors,
    Tissue = tissue_colors,
    Cell = cell_colors
  )
)

# binary_matrix <- enrich_res_merged_filter %>%
#   tidyr::unite("combo", c(species, tissue, cell_type, treatment), remove = FALSE) %>%
#   dplyr::select(combo, pathway, `Log(q-value)`) %>%
#   distinct() %>%
#   tidyr::pivot_wider(names_from = pathway, values_from = `Log(q-value)`, values_fill = 0)
# binary_matrix = binary_matrix[, colSums(binary_matrix<log10(0.05))>1]


heatmap <- Heatmap(as.matrix(t(1-jaccard_dist)), 
                   name = "Similarity",
                   # col = colorRamp2(c(0, 1), c("white", "red")),
                   show_row_names = FALSE,
                   show_column_names = FALSE,
                   cluster_rows = TRUE,
                   cluster_columns = TRUE,
                   clustering_distance_rows = "euclidean",
                   clustering_distance_columns = "euclidean",
                   clustering_method_rows = "complete",
                   clustering_method_columns = "complete",
                   column_split = 3,
                   row_split = 3,
                   top_annotation = top_annotation,
                   column_title_rot = 45,
                   row_title_rot = 360,
                   heatmap_legend_param = list(title = "Jaccard\nsimilarity", legend_direction = "vertical"))

# Draw heatmap
draw(heatmap, heatmap_legend_side = "right", annotation_legend_side = "right")


column_clusters <- column_order(heatmap)
highlighted_cell_types = binary_matrix[,1][column_clusters[[2]],]

enriched_pathways_highlighted <- enrich_res_merged_filter %>%
  mutate(cell_group = paste(species, tissue, cell_type, treatment, sep = "_")) %>%
  filter(cell_group %in% highlighted_cell_types$combo)

pathway_contribution <- enriched_pathways_highlighted %>%
  group_by(pathway) %>%
  dplyr::summarise(frequency = n(), mean_enrichment_score = mean(`Log(q-value)`)) %>%
  arrange(desc(frequency))

ggplot(pathway_contribution[1:50,], aes(x = reorder(pathway, frequency), y = frequency, fill = -mean_enrichment_score)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient2(low="white", high="red") +
  coord_flip() +
  labs(title = "", x = "Pathway", y = "Frequency", fill="Average\nenrichment\nscore")

gene_list = unique(unlist(sapply(enriched_pathways_highlighted[enriched_pathways_highlighted$pathway %in% pathway_contribution$pathway[1:50],]$Symbols, function(x) strsplit(x, split=","))))


mouse2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
mouse2human = unique(mouse2human[, c("Human.gene.name", "Gene.name")])
rat2human = read.table("/xdisk/mliang1/qqiu/reference/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")
rat2human = unique(rat2human[, c("Human.gene.name", "Gene.name")])

pathway_gene_use = unique(c(mouse2human[mouse2human$Human.gene.name %in% gene_list,]$Gene.name,
                            rat2human[rat2human$Human.gene.name %in% gene_list,]$Gene.name))

Clic2, Cox6b1, Cox6b2









distances <- dist(mds_points[,c(1:2)], method = "euclidean")

# full_matrix <- as.matrix(jaccard_dist)
full_matrix <- as.matrix(distances)
rownames(full_matrix) <- binary_matrix$combo
colnames(full_matrix) <- binary_matrix$combo
long_format_df <- melt(full_matrix)
long_format_df = long_format_df[long_format_df$Var1 != long_format_df$Var2, ]
long_format_df[, c("species.x", "tissue.x", "cell_type.x", "treatment.x")] = do.call(rbind, strsplit(as.character(long_format_df$Var1), "_"))
long_format_df[, c("species.y", "tissue.y", "cell_type.y", "treatment.y")] = do.call(rbind, strsplit(as.character(long_format_df$Var2), "_"))
long_format_df = long_format_df[long_format_df$cell_type.x==long_format_df$cell_type.y,]
long_format_df = long_format_df[long_format_df$tissue.x==long_format_df$tissue.y,]
long_format_df$cell_type.x = factor(long_format_df$cell_type.x, levels=cell_order)
long_format_df = long_format_df[long_format_df$species.x %in% c("C57BL/6", "SS", "SHR") & long_format_df$species.y %in% c("C57BL/6", "SS", "SHR"), ]

ggplot(long_format_df, aes(x = value, y = reorder_within(cell_type.x, value, tissue.x, FUN=median))) +
  # geom_vline(xintercept = c(0, 0.1, 0.2, 0.3), color="grey", alpha=0.3)+
  geom_boxplot(color="darkgrey", fill="white") +
  geom_point(color="blue", alpha=0.8) +
  # scale_y_discrete(limits=rev) +
  scale_y_reordered() +
  # scale_color_gradient(low="white", high="blue") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  theme_classic() +
  labs(x="Dissimilarity distance", y="") +
  facet_grid(rows = vars(tissue.x), 
             scales = "free", space = "free_y")
















library(ComplexHeatmap)
library(circlize)

enrich_res_merged_filter = enrich_res_merged[enrich_res_merged$p.adjust<0.05, ]
binary_matrix <- enrich_res_merged_filter %>%
  unite("combo", c(species, tissue, cell_type, treatment), remove = FALSE) %>%
  mutate(presence = 1) %>%
  dplyr::select(combo, Description, presence) %>%
  distinct() %>%
  spread(key = Description, value = presence, fill = 0)

dist_matrix <- dist(t(binary_matrix[-1]))
hc <- hclust(dist_matrix)
cluster_order <- order(hc$order)

enrich_res_merged_filter$ID <- factor(enrich_res_merged_filter$ID, levels = colnames(binary_matrix)[cluster_order])



# Assuming 'df' is your data frame and 'data_matrix' is a matrix of the values you want to plot
# Convert your data frame to a matrix if necessary, making sure rows are genes/cell types and columns are conditions/samples

# Example conversion (adjust according to your actual data structure)
# data_matrix <- matrix(df$Expression, nrow = length(unique(df$Gene)), byrow = TRUE,
#                       dimnames = list(unique(df$Gene), unique(df$Condition)))

# Perform hierarchical clustering on rows and columns if you haven't already
row_hclust <- hclust(dist(t(binary_matrix[-1])))  # For rows
col_hclust <- hclust(dist(binary_matrix[-1]))  # For columns

# Make the heatmap with dendrograms
Heatmap(binary_matrix[-1],
        name = "Expression",
        row_title = "Genes", column_title = "Conditions",
        row_dend_reorder = TRUE, column_dend_reorder = TRUE,
        row_dend_side = "left", column_dend_side = "top",
        clustering_distance_rows = dist(t(binary_matrix[-1])),
        clustering_distance_columns = dist(binary_matrix[-1]),
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete",
        clustering_method_columns = "complete")

binary_matrix = binary_matrix[-1]
Heatmap(binary_matrix,
        name = "Presence", 
        show_row_names = TRUE, 
        show_column_names = TRUE,
        row_title = "Pathways", 
        column_title = "Samples",
        clustering_distance_rows = "euclidean", 
        clustering_distance_columns = "euclidean",
        clustering_method_rows = "complete", 
        clustering_method_columns = "complete",
        row_dend_side = "left",
        column_dend_side = "top",
        heatmap_legend_param = list(title = "Pathway Presence", at = c(0, 1), labels = c("Absent", "Present")),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (binary_matrix[i, j] == 1) {
            grid.rect(x, y, width, height, gp = gpar(fill = "black", col = NA))
          } else {
            grid.rect(x, y, width, height, gp = gpar(fill = "white", col = NA))
          }
        })
















ggplot(enrich_res_merged_filter, aes(x = ID, y = cell_type, fill=-log10(p.adjust)>0)) +
  geom_tile() +
  scale_y_discrete(limits=rev) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    # axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.x = element_blank(),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")




################################################################################

calculateJaccardIndex <- function(set1, set2) {
  intersectionSize <- length(intersect(set1, set2))
  unionSize <- length(union(set1, set2))
  jaccardIndex <- intersectionSize / unionSize
  return(jaccardIndex)
}

performPermutationTest <- function(deg_c1, deg_c2, all_gene_c1, all_gene_c2, n_perm = 1000) {
  observed_similarity = calculateJaccardIndex(deg_c1, deg_c2)
  perm_similarities = numeric(n_perm)
  
  for (i in 1:n_perm) {
    perm_c1 = sample(all_gene_c1, length(deg_c1))
    perm_c2 = sample(all_gene_c2, length(deg_c2))
    perm_similarities[i] = calculateJaccardIndex(perm_c1, perm_c2)
  }
  
  p_value = mean(perm_similarities >= observed_similarity)
  return(p_value)
}

enrich_res_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.ortho.out")
enrich_res_merged$ji_comp = paste(enrich_res_merged$species, enrich_res_merged$treatment,  sep="-")

enrich_res_merged_filter = enrich_res_merged[enrich_res_merged$p.adjust<0.05, ]

comp_list1 = list(c("C57BL/6-AngII 3d", "C57BL/6-AngII 28d"),
                  c("SS-HS 3d", "SS-HS 21d"),
                  c("SS-HS 3d", "SD-HS 3d"),
                  c("SHR-26w", "WKY-26w"))

comp_list2 = list(c("AngII", "Salt-sensitive"),
                  c("Salt-sensitive", "Spontaneous"),
                  c("AngII", "Spontaneous"))

tissue_list = unique(enrich_res_merged_filter$tissue)

ji_df = c()

for(ti in tissue_list){
  
  cell_list = unique(enrich_res_merged_filter[enrich_res_merged_filter$tissue==ti, ]$cell_type)
  
  for(ci in cell_list){
    
    deg_use = enrich_res_merged_filter[enrich_res_merged_filter$tissue==ti & enrich_res_merged_filter$cell_type==ci, ]
    all_gene_use = enrich_res_merged[enrich_res_merged$tissue==ti & enrich_res_merged$cell_type==ci, ]
    
    for(j in comp_list1){
      
      comp1 = j[1]
      comp2 = j[2]
      similarity = ""
      p_value = ""
      
      # deg_c1 = unique(deg_use[deg_use$ji_comp==comp1, ]$gene_name)
      # deg_c2 = unique(deg_use[deg_use$ji_comp==comp2, ]$gene_name)
      
      gene_name = "ID"
      
      deg_c1 = unique(deg_use[deg_use$ji_comp==comp1, ]$ID)
      deg_c2 = unique(deg_use[deg_use$ji_comp==comp2, ]$ID)
      
      
      if(length(deg_c1)>0 & length(deg_c2)>0){
        similarity = calculateJaccardIndex(deg_c1, deg_c2)
        # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
        p_value = performPermutationTest(deg_c1, deg_c2, 
                                         unique(all_gene_use[all_gene_use$ji_comp==comp1, ]$ID),
                                         unique(all_gene_use[all_gene_use$ji_comp==comp2, ]$ID))
      }
      
      ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
      
    }
    
    
    deg_use = enrich_res_merged_filter[enrich_res_merged_filter$tissue==ti & enrich_res_merged_filter$cell_type==ci & 
                                         enrich_res_merged_filter$species %in% c("C57BL/6", "SS", "SHR"), ]
    
    for(j in comp_list2){
      
      comp1 = j[1]
      comp2 = j[2]
      similarity = ""
      p_value = ""
      
      gene_name = "ID"
      
      deg_c1 = unique(deg_use[deg_use$project==comp1, ]$ID)
      deg_c2 = unique(deg_use[deg_use$project==comp2, ]$ID)
      
      if(length(deg_c1)>0 & length(deg_c2)>0){
        similarity = calculateJaccardIndex(deg_c1, deg_c2)
        # similarity = calculateDiceCoefficient(deg_c1, deg_c2)
        p_value = performPermutationTest(deg_c1, deg_c2, 
                                         unique(all_gene_use[all_gene_use$project==comp1, ]$ID),
                                         unique(all_gene_use[all_gene_use$project==comp2, ]$ID))
      }
      
      ji_df = rbind(ji_df, c(ti, ci, comp1, comp2, similarity, p_value))
      
    }
    
  }
  
}

ji_df = as.data.frame(ji_df)
colnames(ji_df) = c("tissue", "cell_type", "comp1", "comp2", "ji_score", "pval")
ji_df$cell_type = factor(ji_df$cell_type, levels=unique(cell_order))
ji_df$tissue = factor(ji_df$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
ji_df$ji_score = as.numeric(ji_df$ji_score)
ji_df$pval = as.numeric(ji_df$pval)
ji_df$padj = p.adjust(ji_df$pval, method = "BH")
ji_df$annotation <- ifelse(ji_df$padj < 0.05, "*", "")
ji_df[is.na(ji_df$annotation), ]$annotation=""
ji_df$comp = paste(ji_df[,3], ji_df[,4], sep="-")
ji_df$comp = factor(ji_df$comp, c("C57BL/6-AngII 3d-C57BL/6-AngII 28d", "SS-HS 3d-SS-HS 21d", 
                                  "SS-HS 3d-SD-HS 3d", "SHR-26w-WKY-26w", "AngII-Salt-sensitive",
                                  "AngII-Spontaneous", "Salt-sensitive-Spontaneous"))

ji_df_use = ji_df[ji_df$comp %in% c("AngII-Salt-sensitive", "AngII-Spontaneous", "Salt-sensitive-Spontaneous"),]

ggplot(ji_df_use, aes(x = comp, y = cell_type, fill=ji_score)) +
  geom_tile() +
  geom_text(aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low="white", high="blue") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Jaccard index") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free_y")









combine_pvalues <- function(row) {
  # Extract the relevant p-values from the row
  p_values <- na.omit(as.numeric(row[3:5]))  # Adjust indices as needed for your selected columns
  
  # Calculate the combined test statistic
  chi_sq_statistic <- -2 * sum(log(p_values))
  
  # Calculate combined p-value using the chi-squared distribution
  combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
  
  return(combined_p_value)
}

# enrich_res_background = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.random.out")
# enrich_res_background_filter = enrich_res_background[enrich_res_background$p.adjust<0.05 & enrich_res_background$species %in% c("C57BL/6", "SHR", "SS"), ]

enrich_res_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG.0.5.withot_background.out")
enrich_res_merged_filter = enrich_res_merged[enrich_res_merged$p.adjust<0.05 & enrich_res_merged$species %in% c("C57BL/6", "SHR", "SS"), ]
enrich_res_merged_filter$cell_gene = paste(enrich_res_merged_filter$tissue, enrich_res_merged_filter$cell_type,  enrich_res_merged_filter$ID, sep="-")
enrich_res_merged_filter = enrich_res_merged_filter %>%
  group_by(species, cell_gene) %>%
  slice(which.min(na.omit(p.adjust))) %>%
  ungroup() %>% as.data.frame()

enrich_res_reshape = reshape(enrich_res_merged_filter[,c("species", "cell_gene", "p.adjust", "Description")], idvar=c("cell_gene", "Description"), timevar = "species", direction = "wide")
enrich_res_reshape[, c("tissue", "cell_type", "ID")] = do.call(rbind, strsplit(enrich_res_reshape$cell_gene, "-"))
enrich_res_reshape$comb.pval = apply(enrich_res_reshape, 1, combine_pvalues)
enrich_res_reshape = enrich_res_reshape[order(enrich_res_reshape$comb.pval), ]

write.table(enrich_res_reshape, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/clusterProfiler/DEG_enrichment.without_background.wide.out")

# venn plot for selected cell types
species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(1,4,2)])

ti = "HYP"; ci = "Microglia"
title = paste0(ti, "-", ci)
deg_use = enrich_res_merged_filter[enrich_res_merged_filter$tissue==ti & enrich_res_merged_filter$cell_type==ci, ]

set1 <- unique(deg_use[deg_use$project=="AngII", ]$ID)
set2 <- unique(deg_use[deg_use$project=="Salt-sensitive", ]$ID)
set3 <- unique(deg_use[deg_use$project=="Spontaneous", ]$ID)

# Create a list of the gene sets
geneSets <- list("AngII" = set1, "Salt-sensitive" = set2, "Spontaneous" = set3)

overlappedGenes <- Reduce(intersect, geneSets)

ggvenn(geneSets, c("AngII", "Salt-sensitive", "Spontaneous"), show_percentage = F,
       fill_color = species_col, fill_alpha=0.3, set_name_size = 0) + labs(title = title)

enrich_res_barplot = enrich_res_reshape[enrich_res_reshape$tissue==ti & 
                                          enrich_res_reshape$cell_type==ci & 
                                          enrich_res_reshape$ID %in% overlappedGenes, ]

enrich_res_barplot = enrich_res_reshape[enrich_res_reshape$tissue==ti & 
                                          enrich_res_reshape$cell_type==ci, ]

enrich_res_barplot = enrich_res_barplot[rowSums(is.na(enrich_res_barplot[,2:4]))<2, ]
# enrich_res_barplot = enrich_res_barplot[grepl("REACTOME", enrich_res_barplot$ID), ]

ggplot(enrich_res_barplot[1:10, ]) +
  geom_col(aes(-log10(comb.pval), reorder(Description, -log10(comb.pval))), fill = "#076fa2", width = 0.6) +
  theme(
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="-log10(Combined p-value)", y="", title = title)

enrich_res_barplot[!(grepl("GO", enrich_res_barplot$ID)),]







































































combine_pvalues <- function(row) {
  # Extract the relevant p-values from the row
  p_values <- na.omit(as.numeric(row[2:4]))  # Adjust indices as needed for your selected columns
  
  # Calculate the combined test statistic
  chi_sq_statistic <- -2 * sum(log(p_values))
  
  # Calculate combined p-value using the chi-squared distribution
  combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
  
  return(combined_p_value)
}


setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/fgsea")
input_files = list.files()

deg_merged = c()
for(i in input_files){
  
  dat_tmp = try(read.table(i, header = T, sep='\t'), silent = T)
  
  if(class(dat_tmp) != "try-error"){
    
    # tissue = sub(".*\\.([^.]*)\\.(RNA|multiomics).*", "\\1", basename(i))
    # dat_tmp$tissue = tissue
    # 
    if(!(grepl("immune_cell", i))){
      dat_tmp = dat_tmp[!(dat_tmp$cell_type %in% c("Microglia", "Activated microglia", "IMM")), ]
    }
    
    print(table(dat_tmp$species))
    deg_merged = rbind(deg_merged, dat_tmp)
    
  }
}

write.table(deg_merged, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/fgsea/GSEA.merged.out")


gsea_res = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/fgsea/GSEA.merged.out")
gsea_res_filter = gsea_res[gsea_res$padj<0.1 & gsea_res$species %in% c("C57BL/6", "SHR", "SS"), ]
gsea_res_filter = gsea_res[gsea_res$species %in% c("C57BL/6", "SHR", "SS"), ]

gsea_res_filter$cell_gene = paste(gsea_res_filter$tissue, gsea_res_filter$cell_type,  gsea_res_filter$pathway, sep="-")
gsea_res_filter$sxt = paste(gsea_res_filter$species, gsea_res_filter$treatment, sep="-")
gsea_res_filter = gsea_res_filter %>%
  group_by(species, cell_gene) %>%
  slice(which.min(na.omit(padj))) %>%
  ungroup() %>% as.data.frame()

enrich_res_reshape = reshape(gsea_res_filter[,c("sxt", "cell_gene", "pval")], idvar="cell_gene", timevar = "sxt", direction = "wide")
enrich_res_reshape[, c("tissue", "cell_type", "ID")] = do.call(rbind, strsplit(enrich_res_reshape$cell_gene, "-"))
enrich_res_reshape$comb.pval = apply(enrich_res_reshape, 1, combine_pvalues)
enrich_res_reshape = enrich_res_reshape[order(enrich_res_reshape$comb.pval), ]



# venn plot for selected cell types
species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
species_col = as.character(species_col[c(1,4,2)])

ti = "LK"; ci = "PT"
title = paste0(ti, "-", ci)
deg_use = gsea_res_filter[gsea_res_filter$tissue==ti & gsea_res_filter$cell_type==ci, ]

set1 <- unique(deg_use[deg_use$species=="C57BL/6", ]$ID)
set2 <- unique(deg_use[deg_use$species=="SS", ]$ID)
set3 <- unique(deg_use[deg_use$species=="SHR", ]$ID)

# Create a list of the gene sets
geneSets <- list("AngII" = set1, "Salt-sensitive" = set2, "Spontaneous" = set3)

overlappedGenes <- Reduce(intersect, geneSets)

ggvenn(geneSets, c("AngII", "Salt-sensitive", "Spontaneous"), show_percentage = F,
       fill_color = species_col, fill_alpha=0.3, set_name_size = 0) + labs(title = title)

enrich_res_barplot = enrich_res_reshape[enrich_res_reshape$tissue==ti & 
                                          enrich_res_reshape$cell_type==ci & 
                                          enrich_res_reshape$ID %in% overlappedGenes, ]

enrich_res_barplot = enrich_res_reshape[enrich_res_reshape$tissue==ti & 
                                          enrich_res_reshape$cell_type==ci, ]
enrich_res_barplot
# enrich_res_barplot = enrich_res_barplot[rowSums(is.na(enrich_res_barplot[,2:6]))<2, ]
# enrich_res_barplot = enrich_res_barplot[grepl("REACTOME", enrich_res_barplot$ID), ]

ggplot(enrich_res_barplot[1:10, ]) +
  geom_col(aes(-log10(comb.pval), reorder(ID, -log10(comb.pval))), fill = "#076fa2", width = 0.6) +
  theme(
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="-log10(Combined p-value)", y="", title = title)

enrich_res_barplot[!(grepl("GO", enrich_res_barplot$ID)),]




