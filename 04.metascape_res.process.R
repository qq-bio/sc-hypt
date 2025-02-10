
library(openxlsx)
library(reshape2)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(stringr)
library(gtools)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape")


################################################################################
### merge enrichment result & enrichment level
merge_all = c()

input_files = list.files(path = "./", pattern = "zip")

for(i in input_files){

  unzip(zipfile=i, files = "metascape_result.xlsx", exdir="./tmp", overwrite = T)
  data = read.xlsx("tmp/metascape_result.xlsx", 2)

  i_mod = gsub("C57BL.6", "C57BL/6", i)
  i_mod = gsub("E.P transition cell", "E/P transition cell", i_mod)
  data[, c("species", "treatment", "tissue", "cell_type")] = strsplit(i_mod, "\\.")[[1]][1:4]
  print(strsplit(i_mod, "\\.")[[1]][1:4])
  data$species = strsplit(i_mod, "\\.")[[1]][1]
  data$treatment = strsplit(i_mod, "\\.")[[1]][2]
  data$tissue = strsplit(i_mod, "\\.")[[1]][3]
  data$cell_type = strsplit(i_mod, "\\.")[[1]][4]

  data$category_mod = lapply(strsplit(data$Category, " "), function(x) x[1]) %>% as.character()
  data$pathway = paste0(data$Description, " (", data$category_mod, ")")

  merge_all = rbind(merge_all, data)

}

merge_all = merge_all[!(merge_all$tissue=="MCA" & merge_all$cell_type %in% c("Neuron", "Myelinating OL", "OPC", "Astrocyte")), ]

write.table(merge_all, "metascape.merge_all.out", sep='\t', col.names = T, row.names = F, quote = F)





################################################################################

combine_pvalues <- function(row) {
  p_values <- na.omit(as.numeric(row))
  # chi_sq_statistic <- -2 * sum(log(p_values))
  chi_sq_statistic <- -2 * sum(p_values/log10(2))
  combined_p_value <- pchisq(chi_sq_statistic, df = 2 * length(p_values), lower.tail = FALSE)
  return(combined_p_value)
}

merge_all = read.table("metascape.merge_all.out", sep='\t', header = T, quote = "", check.names = F)

merge_all$cell_type = factor(merge_all$cell_type, levels = cell_order)
merge_all$tissue = factor(merge_all$tissue, levels = tissue_order)
merge_all$path_cell = paste(merge_all$tissue, merge_all$cell_type,  merge_all$pathway, sep="_")
merge_all$sxt = paste(merge_all$species, merge_all$treatment, sep="-")
merge_all[merge_all$species=="Salt-sensitive", ]$sxt = "SS-SD"
merge_all[merge_all$species=="Spontaneous", ]$sxt = "SHR-WKY"

merge_background = merge_all[merge_all$species %in% c("SD", "WKY"), ]
path_rm_list = merge_background[merge_background$LogP<log10(0.05), ]$path_cell
merge_filter = merge_all[! merge_all$path_cell %in% path_rm_list, ]

merge_summary = merge_filter[grepl("Summary", merge_filter$GroupID), ]$path_cell
# merge_summary = merge_filter[grepl("Member", merge_filter$GroupID), ]
# merge_summary = merge_filter; outfile = "metascape.merge_all.reshape.out"
all_sig <- merge_filter[merge_filter$`Log(q-value)`< log10(0.05),]
top20 <- merge_filter[merge_filter$`Log(q-value)`< log10(0.05),] %>% group_by(species, treatment, tissue, cell_type) %>% top_n(n = 20) %>% as.data.frame()
top_list <- all_sig

merge_filter2 = merge_filter[merge_filter$path_cell %in% top_list$path_cell, ]
merge_reshape = reshape(merge_filter2[,c("sxt", "path_cell", "Log(q-value)")], idvar=c("path_cell"), timevar = "sxt", v.names = "Log(q-value)", direction = "wide")
merge_reshape$`Log(comb.qval)` = log10(apply(merge_reshape[,2:6], 1, combine_pvalues))

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
merge_reshape = na.replace(merge_reshape, 0)
merge_reshape$AngII = apply(merge_reshape[,2:3], 1, function(x) ifelse(any(as.numeric(x) > -log10(0.05)), 1, 0))
merge_reshape$SS = apply(merge_reshape[,5:6], 1, function(x) ifelse(any(as.numeric(x) > -log10(0.05)), 1, 0))
merge_reshape$SHR = ifelse(as.numeric(merge_reshape[,4]) > -log10(0.05), 1, 0)

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

merge_reshape$model_count = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
merge_reshape$cell_count = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$pathway==x))
merge_reshape$pathway_sum = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$model_count[merge_reshape$pathway==x]))
merge_reshape$summary = ifelse(merge_reshape$path_cell %in% merge_summary, "Yes", "No")

write.table(merge_reshape, "metascape.merge_all.reshape.out", sep='\t', quote=F, col.names=T, row.names=F)



################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")

merge_reshape = merge_reshape[!(merge_reshape$tissue=="MCA" & merge_reshape$cell_type %in% c("Neuron", "Astrocyte", "OPC", "Myelinating OL")),]

path_count = merge_reshape %>% filter(category != "Other", summary == "Yes") %>% 
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



################################################################################
### model-wise difference

path_count = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`SS-LS vs. SD-LS` = sum(Log.q.value..SS.SD>1.3),
                   `SS-HS 3d vs. SS-LS` = sum(Log.q.value..SS.HS.3d>1.3),
                   `SS-HS 21d vs. SS-LS` = sum(Log.q.value..SS.HS.21d>1.3)
                   ) %>% 
  melt() %>% filter(value>0) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)

ggplot(path_count, aes(y = cell_type, x = variable, label = value)) +
  geom_point(aes(size=value/3), alpha = .5, color="red") + geom_text() +
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


path_count = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`SHR-10w vs. WKY-10w` = sum(Log.q.value..SHR.WKY>1.3),
                   `SHR-26w vs. SHR-10w` = sum(Log.q.value..SHR.26w>1.3)
  ) %>% 
  melt() %>% filter(value>0) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)

ggplot(path_count, aes(y = cell_type, x = variable, label = value)) +
  geom_point(aes(size=value/3), alpha = .5, color="red") + geom_text() +
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


path_count = merge_reshape %>% filter(summary == "Yes") %>% 
  group_by(tissue, cell_type) %>%
  dplyr::summarise(`C57BL/6-AngII 3d vs. Saline 3d` = sum(Log.q.value..C57BL.6.AngII.3d>1.3),
                   `C57BL/6-AngII 28d vs. Saline 3d` = sum(Log.q.value..C57BL.6.AngII.28d>1.3)
  ) %>% 
  melt() %>% filter(value>0) %>% as.data.frame()

path_count$cell_type = factor(path_count$cell_type, levels = cell_order)
path_count$tissue = factor(path_count$tissue, levels = tissue_order)

ggplot(path_count, aes(y = cell_type, x = variable, label = value)) +
  geom_point(aes(size=value/3), alpha = .5, color="red") + geom_text() +
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




length(intersect(merge_reshape[merge_reshape$Log.q.value..SS.SD>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                 merge_reshape[merge_reshape$Log.q.value..SS.HS.3d>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway))
length(intersect(merge_reshape[merge_reshape$Log.q.value..SS.SD>1.3 & merge_reshape$tissue=="HYP" & merge_reshape$cell_type=="Astrocyte" & merge_reshape$summary=="Yes", ]$pathway,
                 merge_reshape[merge_reshape$Log.q.value..SS.HS.3d>1.3 & merge_reshape$tissue=="HYP" & merge_reshape$cell_type=="Astrocyte" & merge_reshape$summary=="Yes", ]$pathway))

length(intersect(merge_reshape[merge_reshape$Log.q.value..SHR.WKY>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="PT" & merge_reshape$summary=="Yes", ]$pathway,
                 merge_reshape[merge_reshape$Log.q.value..SHR.26w>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="PT" & merge_reshape$summary=="Yes", ]$pathway))

length(intersect(merge_reshape[merge_reshape$Log.q.value..SS.SD>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                 intersect(merge_reshape[merge_reshape$Log.q.value..SS.HS.3d>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway,
                           merge_reshape[merge_reshape$Log.q.value..SHR.26w>1.3 & merge_reshape$tissue=="LK" & merge_reshape$cell_type=="EC" & merge_reshape$summary=="Yes", ]$pathway)))




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


### heatmap
library(tidyverse)

merge_reshape$cell_type = factor(merge_reshape$cell_type, levels = cell_order)
merge_reshape$tissue = factor(merge_reshape$tissue, levels = tissue_order)
merge_reshape$size = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
merge_reshape$pathway_count = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$pathway==x))
merge_reshape$pathway_sum = sapply(merge_reshape$pathway, function(x) sum(merge_reshape$size[merge_reshape$pathway==x]))

# merge_reshape_use = merge_reshape[grepl("AngII", merge_reshape$category),]
merge_reshape_use = merge_reshape[merge_reshape$size>1 & merge_reshape$summary=="Yes",]
# pathway_list = merge_reshape_use$pathway[merge_reshape_use$pathway_sum>2]
# pathway_list = merge_reshape_use$pathway[merge_reshape_use$size>1]
# merge_reshape_use = merge_reshape_use[merge_reshape_use$pathway %in% pathway_list, ]
ggplot(merge_reshape_use, 
       aes(y = cell_type, x = fct_inorder(pathway), fill = factor(as.integer(size)))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  # scale_fill_gradient(low = "peachpuff", high = "orangered4") +
  scale_fill_manual(values = c("peachpuff", "orangered4")) +
  scale_y_discrete(limits=rev) +
  labs(x="", y="", fill="Number of\nmodel") +
  facet_grid(cols = vars(tissue), 
             scales = "free", space = "free")+
  coord_flip()





merge_reshape$size = rowSums(merge_reshape[, c("AngII", "SS", "SHR")])
merge_reshape_use = merge_reshape[merge_reshape$size>1,]
# merge_reshape_use$pathway_count = sapply(merge_reshape_use$pathway, function(x) sum(merge_reshape_use$pathway==x))
merge_reshape_use$pathway_sum = sapply(merge_reshape_use$pathway, function(x) sum(merge_reshape_use$size[merge_reshape_use$pathway==x]))
pathway_list = merge_reshape_use$pathway[merge_reshape_use$pathway_sum>2]
merge_reshape_use = merge_reshape_use[merge_reshape_use$pathway %in% pathway_list, ]
ggplot(merge_reshape_use, 
       aes(y = cell_type, x = reorder(pathway, pathway_sum), fill = as.numeric(size))) +
  geom_tile() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = 'black'),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
        legend.text = element_text(colour = 'black')) +
  scale_fill_gradient(low = "peachpuff", high = "orangered4") +
  scale_y_discrete(limits=rev) +
  labs(x="", y="", fill="Number of\nmodel") +
  facet_grid(rows = vars(tissue), 
             scales = "free", space = "free")



Heatmap(merge_reshape[,2:6], 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE
)

ggplot(merge_reshape, aes(y = pathway, x = cell_type, fill=-log10(comb.qval)>2)) +
  geom_tile() +
  scale_x_discrete(limits=rev) +
  # scale_fill_gradient(low="white", high="red", na.value = "white") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="-log10(adjusted pvalue)") +
  facet_grid(cols = vars(tissue), 
             scales = "free", space = "free_x")











library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)

merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")

rownames(merge_reshape) = merge_reshape$path_cell
merge_reshape_ang = merge_reshape[, c(3, 2)]
merge_reshape_ss = merge_reshape[, c(7, 6, 5)]
merge_reshape_shr = merge_reshape[, c(8, 4)]

data = merge_reshape_ss
data = as.data.frame(ifelse(data > -log10(0.05), 1, 0))
data$condition = apply(data, 1, function(row) paste0(row, collapse = ""))
data$tissue_cell = sapply(rownames(data), function(x) {tmp=strsplit(x, "_")[[1]]; paste0(tmp[1], "-", tmp[2])})

filtered_df <- data_df %>% 
  filter(log10_p_value > -log10(0.05)) %>%
  group_by(tissue_cell, Condition) %>%
  summarise(Count = n(), .groups = "drop")

alluvial_data <- filtered_df %>%
  tidyr::pivot_longer(cols = Count, names_to = "Key", values_to = "Value")


ggplot(data = data,
       aes(axis1 = Log.q.value..SS.SD, axis2 = Log.q.value..SS.HS.3d, y = Weight)) +
  geom_alluvium(aes(fill = Pathway)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  ggtitle("Alluvial Diagram of Pathway Enrichment under Different Conditions")



















path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>15]=15

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[unique(top_list$pathway),]

write.table(path_all_pat_df_plot, "pathway.unique.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

library(circlize)
f1 = colorRamp2(seq(0, 15, length = 15), c(magma(n=15)))
select_idx = c(1,2,4,6,9,11,12,17,21,22,24,32,36,44,45,52,56,59,61,69,71,76)
ha = rowAnnotation(foo = anno_mark(at = select_idx, 
                                   labels = rownames(path_all_pat_df_plot)[select_idx]))
pdf("cluster.pathway.unique.pdf", height = 6, width = 8)
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           right_annotation = ha,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           # column_names_rot = 45,
           show_row_names = FALSE,
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")






























summary_list = unique(merge_all[grepl("Summary", merge_all$GroupID), ]$pathway)
metascape_use = merge_all[merge_all$pathway %in% summary_list, ]

ggplot(merge_reshape, aes(y = pathway, x = cell_type, fill=-log10(comb.qval)>2)) +
  geom_tile() +
  scale_x_discrete(limits=rev) +
  # scale_fill_gradient(low="white", high="red", na.value = "white") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="-log10(adjusted pvalue)") +
  facet_grid(cols = vars(tissue), 
             scales = "free", space = "free_x")




path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>10]=10

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[top10$pathway,]

write.table(path_all_pat_df_plot, "pathway.summary.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

library(circlize)
f1 = colorRamp2(seq(0, 10, length = 10), c(magma(n=10)))
pdf("cluster.pathway.summary.pdf")
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()






## code for pre-metastasis & primary
path_use = c(path_pat_df[path_pat_df$patient=="Cluster 1",]$pathway[4:7],
             path_pat_df[path_pat_df$patient=="Cluster 2",]$pathway[c(1,3,5)],
             path_pat_df[path_pat_df$patient=="Cluster 3",]$pathway[c(1,16,19)],
             path_pat_df[path_pat_df$patient=="Cluster 4",]$pathway[c(3,5,8)],
             path_pat_df[path_pat_df$patient=="Cluster 5",]$pathway[1],
             path_pat_df[path_pat_df$patient=="Cluster 6",]$pathway[c(1,2,7)],
             path_pat_df[path_pat_df$patient=="Cluster 7",]$pathway[c(1)],
             path_pat_df[path_pat_df$patient=="Cluster 8",]$pathway[c(5,20)])
path_use = unique(path_use)
path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% path_use, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)-1*as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>10]=10

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))

# library(circlize)
# f1 = colorRamp2(seq(min(path_all_pat_df_plot), max(path_all_pat_df_plot), length = 2), c("#EEEEEE", "red"))
pdf("P_pre_P.pathway_enrich.pdf")
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           row_names_max_width = max_text_width(
             rownames(path_all_pat_df_plot),
             gp = gpar(fontsize = 12)
           ),
           # show_column_names = FALSE
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()


### visualize top 10 heatmap & highlight selected pathway
path_pat_df_use = path_pat_df[path_pat_df$es>2,]
path_pat_df_use$count <- table(path_pat_df_use$pathway)[path_pat_df_use$pathway]
path_pat_df_use = path_pat_df_use[path_pat_df_use$count<3, -4]

top10 <- path_pat_df_use %>% group_by(patient) %>% top_n(n = 10, es) %>% as.data.frame()
path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% top10$pathway, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>15]=15

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[unique(top10$pathway),]

write.table(path_all_pat_df_plot, "pathway.summary.unique.filter.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

# library(circlize)
f1 = colorRamp2(seq(0, 15, length = 15), c(magma(n=15)))
pdf("cluster.pathway.summary.unique.filter.pdf")
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()

### visualize top 10 heatmap & highlight selected pathway
top10 <- path_all_pat_df %>% group_by(patient) %>% top_n(n = 10) %>% as.data.frame()
path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% top10$pathway, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>10]=10

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[top10$pathway,]

write.table(path_all_pat_df_plot, "pathway.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

# library(circlize)
f1 = colorRamp2(seq(0, 10, length = 10), c(magma(n=10)))
pdf("cluster.pathway.pdf")
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()


### visualize top 10 heatmap & highlight selected pathway
# path_all_pat_df_use = path_all_pat_df[path_all_pat_df$es>2,]
# path_all_pat_df_use$count <- table(path_all_pat_df_use$pathway)[path_all_pat_df_use$pathway]
# path_all_pat_df_use = path_all_pat_df_use[path_all_pat_df_use$count<3, ]

top10 <- path_all_pat_df[path_all_pat_df$es>2,] %>% group_by(patient) %>% top_n(n = 10) %>% as.data.frame()
path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% top10$pathway, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>15]=15

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[unique(top10$pathway),]

write.table(path_all_pat_df_plot, "pathway.unique.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

# library(circlize)
f1 = colorRamp2(seq(0, 15, length = 15), c(magma(n=15)))
select_idx = c(1,2,6,9,10,14,16,21,24,36,39,45,46,59,60,64,65,72,74,81,82,90)
ha = rowAnnotation(foo = anno_mark(at = select_idx, 
                                   labels = rownames(path_all_pat_df_plot)[select_idx]))
pdf("cluster.pathway.unique.pdf", height = 6, width = 8)
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           right_annotation = ha,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           # column_names_rot = 45,
           show_row_names = FALSE,
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()


### visualize top 10 heatmap & highlight selected pathway
path_all_pat_df_use = path_all_pat_df[path_all_pat_df$es>2,]
path_all_pat_df_use$count <- table(path_all_pat_df_use$pathway)[path_all_pat_df_use$pathway]
path_all_pat_df_use = path_all_pat_df_use[path_all_pat_df_use$count<3, -4]

top10 <- path_all_pat_df_use %>% group_by(patient) %>% top_n(n = 10, es) %>% as.data.frame()
top20 <- path_all_pat_df_use %>% group_by(patient) %>% top_n(n = 20, es) %>% as.data.frame()
top_list <- top10

path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% top_list$pathway, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>15]=15

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[unique(top_list$pathway),]

write.table(path_all_pat_df_plot, "pathway.unique.filter.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

# library(circlize)
f1 = colorRamp2(seq(0, 15, length = 15), c(magma(n=15)))
pdf("cluster.pathway.unique.filter.pdf")
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()














path_pat_df = data.frame(patient=character(),
                         pathway=character(),
                         es=numeric())

path_all_pat_df = data.frame(patient=character(),
                             pathway=character(),
                             es=numeric())

for(fi in list.files(".", "metascape.zip")){
  # print(paste0("Patient ", pi))
  
  ## for result in monocle
  # meta_zip = paste0("Pat", pi, ".P_pre_P.count.patien2.metascape.zip")
  ## for result in P-M
  pi = gsub(".metascape.zip", "", fi)
  
  unzip(zipfile=fi, files = "metascape_result.xlsx", exdir="./tmp", overwrite = T)
  data <- read.xlsx("tmp/metascape_result.xlsx", 2, header = TRUE) 
  
  data_use = data[grep("Summary", data$GroupID), ]
  category = strsplit(data_use$Category, " ")
  data_use$Category = sapply(category, '[', seq(max(sapply(category,length))))[1,]
  path_pat_df_tmp = data.frame(patient=rep(pi, nrow(data_use)),
                               pathway=paste0(str_to_title(data_use$Description), " (", data_use$Category, ")"),
                               es=-1*data_use$Log.q.value.)
  path_pat_df = rbind(path_pat_df, path_pat_df_tmp)
  
  data_use = data[-grep("Summary", data$GroupID), ]
  category = strsplit(data_use$Category, " ")
  data_use$Category = sapply(category, '[', seq(max(sapply(category,length))))[1,]
  path_pat_all_df_tmp = data.frame(patient=rep(pi, nrow(data_use)),
                                   pathway=paste0(str_to_title(data_use$Description), " (", data_use$Category, ")"),
                                   es=-1*data_use$Log.q.value.)
  path_all_pat_df = rbind(path_all_pat_df, path_pat_all_df_tmp)
}

top10 <- path_all_pat_df[path_all_pat_df$es>2,] %>% group_by(patient) %>% top_n(n = 10) %>% as.data.frame()
top20 <- path_all_pat_df[path_all_pat_df$es>2,] %>% group_by(patient) %>% top_n(n = 20, es) %>% as.data.frame()
top_list <- top10

path_all_pat_df_plot = path_all_pat_df[path_all_pat_df$pathway %in% top_list$pathway, ]
path_all_pat_df_plot = t(reshape(path_all_pat_df_plot, idvar = c("patient"), timevar = "pathway", direction = "wide"))
path_all_pat_df_plot[is.na(path_all_pat_df_plot)]=0
colnames(path_all_pat_df_plot) = path_all_pat_df_plot[1,]
path_all_pat_df_plot = path_all_pat_df_plot[-1,]
path_all_pat_df_plot = apply(path_all_pat_df_plot, 1:2, function(x)as.numeric(as.character(x)))
path_all_pat_df_plot[path_all_pat_df_plot>15]=15

rownames(path_all_pat_df_plot) = gsub("^es.", "", rownames(path_all_pat_df_plot))
path_all_pat_df_plot = path_all_pat_df_plot[unique(top_list$pathway),]

write.table(path_all_pat_df_plot, "pathway.unique.es.matrix.out", col.names = T, row.names = T, quote=F, sep='\t')

library(circlize)
f1 = colorRamp2(seq(0, 15, length = 15), c(magma(n=15)))
select_idx = c(1,2,4,6,9,11,12,17,21,22,24,32,36,44,45,52,56,59,61,69,71,76)
ha = rowAnnotation(foo = anno_mark(at = select_idx, 
                                   labels = rownames(path_all_pat_df_plot)[select_idx]))
pdf("cluster.pathway.unique.pdf", height = 6, width = 8)
ht=Heatmap(path_all_pat_df_plot, 
           cluster_rows = FALSE,
           cluster_columns = FALSE,
           col=f1,
           right_annotation = ha,
           # row_names_max_width = max_text_width(
           #   rownames(path_all_pat_df_plot),
           #   gp = gpar(fontsize = 12)
           # ),
           # show_column_names = FALSE
           # column_names_rot = 45,
           show_row_names = FALSE,
           heatmap_legend_param = list(
             title = "-log(q-value)", direction="horizontal")
)
draw(ht, heatmap_legend_side="bottom")
dev.off()



### calculate gene pathway activity score

### load time-series data
#### latent time

ds_epi = readRDS("/data/itmll/qzqiu/scRNA/analysis/cluster/epithelial/epithelial_harmony_cluster.umap.rds")
ds_epi[["wnn.umap.harmony"]] <- NULL
sample_list = unique(ds_epi$Patien2)

pid_list = c(428145, 437835, 429346, 407537, 415673, 418334, 441036, 446171)
for(pi in 1:8){
  pid = pid_list[pi]
  latent_file = paste0("/data/itmll/qzqiu/scRNA/analysis/scvelo/", pid, "/latent_time.txt")
  pseudo_file = paste0("/data/itmll/qzqiu/scRNA/analysis/monocle2/Pat", pi, ".pseudotime.txt")
  
  latent_dat = read.table(latent_file, header=T)
  pseudo_dat = read.table(pseudo_file, header=T)
  
  es_df = data.frame(latent_time=latent_dat$latent_time,
                     pseudo_time=pseudo_dat$count.celltype)
  
  patient_list = c(sample_list[c(pi, pi+8)])
  Pat1 = subset(ds_epi, Patien2 %in% patient_list)
  
  expr_mat = as.matrix(GetAssayData(object = Pat1, assay = "RNA", slot = "data"))
  path_use = names(sort(table(path_pat_df$Pathway), decreasing = T))
  path_gene_use = path_gene_lst[path_use]
  es = gsva(expr_mat, path_gene_use, method='ssgsea')
  
  es_df = cbind(es_df, t(es))
  es_df = melt(es_df, id.vars = c("latent_time", "pseudo_time"))
  es_df['patient'] = paste0("Pat", pi)
  
  write.table(es_df, paste0("Pat", pi, ".selected_pathway.GSEA.out"), col.names = T, row.names = F, quote = F, sep='\t')
  
}


### merge result together
es_df = data.frame(latent_time=numeric(),
                   pseudo_time=numeric(),
                   pathway=character(),
                   es=numeric(),
                   patient=character(),
                   es.pred.latent=numeric(),
                   es.pred.pseudo=numeric())
for(pi in 1:8){
  es_file = paste0("Pat", pi, ".selected_pathway.GSEA.out")
  es_dat = read.table(es_file, header=T, sep='\t')
  
  es_pred_latent = c(); es_pred_pseudo = c()
  for(i in unique(es_dat$variable)){
    fit1 = vglm(value ~ sm.ns(latent_time, df=3), uninormal, data=es_dat[es_dat$variable==i,])
    pred1 = predict(fit1, type = "response")
    es_pred_latent = c(es_pred_latent, pred1)
    
    fit2 = vglm(value ~ sm.ns(pseudo_time, df=3), uninormal, data=es_dat[es_dat$variable==i,])
    pred2 = predict(fit2, type = "response")
    es_pred_pseudo = c(es_pred_pseudo, pred2)
    
  }
  
  es_dat['es.pred.latent'] = es_pred_latent
  es_dat['es.pred.pseudo'] = es_pred_pseudo
  
  es_df = rbind(es_df, es_dat)
}
colnames(es_df) = c("latent_time", "pseudo_time", "pathway", "es", "patient", "es.pred.latent", "es.pred.pseudo")
write.table(es_df, "merged.selected_pathway.GSEA.out", col.names = T, row.names = F, sep='\t', quote=F)

### visualize by heatmap (per sample)
meta_dat = read.table("/data/itmll/qzqiu/scRNA/analysis/cluster/epithelial/epithelial.metadata.out", 
                      header=T, sep='\t')
library(circlize)
for(pi in 1:8){
  es_file = paste0("Pat", pi, ".selected_pathway.GSEA.out")
  es_dat = read.table(es_file, header=T, sep='\t')
  
  es_dat = reshape(es_dat[,-5], idvar = c("latent_time", "pseudo_time"), timevar = "variable", direction = "wide")
  colnames(es_dat) = gsub("value.", "", colnames(es_dat))
  
  for(i in c("latent_time", "pseudo_time")){
    meta_zip = paste0("Pat", pi, ".P_pre_P.count.patien2.metascape.zip")
    unzip(zipfile=meta_zip, files = "metascape_result.xlsx", exdir="./tmp", overwrite = T)
    data <- read.xlsx("tmp/metascape_result.xlsx", 2, header = TRUE) 
    
    data_use = data[grep("Summary", data$GroupID), ]
    path_use = data_use[1:10, ]$Description
    es_dat_use = t(es_dat[order(es_dat[,i]), path_use])
    
    mat = es_dat_use
    # mat[mat>0.4]=1
    # mat = t(scale(t(es_dat_use)))
    f1 = colorRamp2(seq(min(mat), max(mat), length = 2), c("#EEEEEE", "black"))
    
    anno = meta_dat[grepl(paste0("_", pi), meta_dat$Patien2), ]$LT
    anno = anno[order(es_dat[,i])]
    ha = HeatmapAnnotation(Sample = anno,
                           col = list(Sample = c("L" = "blue", "T" = "red")))
    
    Heatmap(mat, 
            top_annotation = ha, 
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            col=f1,
            # row_names_max_width = max_text_width(
            #   rownames(es_dat_use), 
            #   gp = gpar(fontsize = 12)
            # ),
            show_column_names = FALSE)
    dev.off()
    
    pi = 1; pre_state_list = c(2,5,7)
    patient = paste0("Pat", pi)
    branch = read.table(paste0("/data/itmll/qzqiu/scRNA/analysis/monocle2/",patient, ".branch.txt"), header=T, sep='\t')
    pre_list = branch[branch$count.celltype %in% pre_state_list, 1]
    p_list = meta_dat[meta_dat$Patien2 == paste0("Local_", pi),]$Row.names
    anno2 = rep("Primary", length(p_list))
    anno2[p_list %in% pre_list] = "Precursor"
    ha = HeatmapAnnotation(Sample = anno2,
                           col = list(Sample = c("Primary" = "blue", "Precursor" = "red")))
    mat = es_dat_use[,anno=="L"]
    Heatmap(mat, 
            top_annotation = ha, 
            cluster_rows = FALSE,
            cluster_columns = FALSE,
            col=f1,
            # row_names_max_width = max_text_width(
            #   rownames(es_dat_use), 
            #   gp = gpar(fontsize = 12)
            # ),
            show_column_names = FALSE)
    dev.off()
    
  }
  
  
}



### visualize by curve
library(dplyr)
es_df = read.table("merged.selected_pathway.GSEA.out", header=T, sep='\t')
es_df = es_df %>%
  group_by(patient) %>%
  mutate(pseudo.mod = pseudo_time/max(pseudo_time, na.rm=TRUE)) %>% as.data.frame

# path_use = names(sort(table(path_pat_df$Pathway), decreasing = T)[1:7])
ggplot(aes(pseudo_time, es), data = es_df[es_df$pathway %in% path_use, ]) +
  # geom_point(aes(color = patient), position=position_jitter()) + 
  # geom_smooth(aes(color = patient), se=FALSE) +
  geom_line(aes(x = pseudo_time, y = es.pred.pseudo, color = patient)) + 
  theme_classic() +
  facet_wrap(~pathway, 
             ncol = 2, scales = "free_y") +
  xlab("Pseudo-time") + ylab("Enrichment score")

ggplot(aes(pseudo.mod, es), data = es_df[es_df$pathway %in% path_use, ]) +
  # geom_point(aes(color = patient), position=position_jitter()) + 
  # geom_smooth(aes(color = patient), se=FALSE) +
  geom_line(aes(x = pseudo.mod, y = es.pred.pseudo, color = patient)) + 
  theme_classic() +
  facet_wrap(~pathway, 
             ncol = 2, scales = "free_y") +
  xlab("Pseudotime") + ylab("Pathway activity score")

ggplot(aes(latent_time, es), data = es_df[es_df$pathway %in% path_use, ]) +
  # geom_point(aes(color = patient), position=position_jitter()) + 
  # geom_smooth(aes(color = patient), se=FALSE) +
  geom_line(aes(x = latent_time, y = es.pred.latent, color = patient)) + 
  theme_classic() +
  facet_wrap(~pathway, 
             ncol = 2, scales = "free_y") +
  xlab("Latent time") + ylab("Pathway activity score")

dev.off()
#### pseudotime

### plot gsva in time


plot_genes_in_pseudotime <-function(cds_subset, 
                                    min_expr=NULL, 
                                    cell_size=0.75, 
                                    nrow=NULL, 
                                    ncol=1, 
                                    panel_order=NULL, 
                                    color_by="State",
                                    trend_formula="~ sm.ns(Pseudotime, df=3)",
                                    label_by_short_name=TRUE,
                                    relative_expr=TRUE,
                                    vertical_jitter=NULL,
                                    horizontal_jitter=NULL){
  
  f_id <- NA
  Cell <- NA
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- reshape2::melt(as.matrix(exprs(cds_subset)))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  #cds_exprs$f_id <- as.character(cds_exprs$f_id)
  #cds_exprs$Cell <- as.character(cds_exprs$Cell)
  
  if (integer_expression) {
    cds_exprs$adjusted_expression <- cds_exprs$expression
  }
  else {
    cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  }
  # trend_formula <- paste("adjusted_expression", trend_formula,
  #     sep = "")
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- as.character(cds_exprs$gene_short_name)
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  cds_exprs$f_id <- as.character(cds_exprs$f_id)
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  new_data <- data.frame(Pseudotime = pData(cds_subset)$Pseudotime)
  model_expectation <- genSmoothCurves(cds_subset, cores=1, trend_formula = trend_formula,
                                       relative_expr = T, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- ddply(cds_exprs, .(f_id, Cell), function(x) data.frame("expectation"=model_expectation[x$f_id, x$Cell]))
  cds_exprs <- merge(cds_exprs, expectation)
  #cds_exprs$expectation <- expectation#apply(cds_exprs,1, function(x) model_expectation[x$f_id, x$Cell])
  
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label,
                                      levels = panel_order)
  }
  q <- ggplot(aes(Pseudotime, expression), data = cds_exprs)
  if (is.null(color_by) == FALSE) {
    q <- q + geom_point(aes_string(color = color_by), size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
  }
  else {
    q <- q + geom_point(size = I(cell_size), position=position_jitter(horizontal_jitter, vertical_jitter))
  }
  
  q <- q + geom_line(aes(x = Pseudotime, y = expectation), data = cds_exprs)
  
  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow,
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  if (relative_expr) {
    q <- q + ylab("Relative Expression")
  }
  else {
    q <- q + ylab("Absolute Expression")
  }
  q <- q + xlab("Pseudo-time")
  q <- q + monocle_theme_opts()
  q
}
