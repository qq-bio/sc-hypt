library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(CellChat)
# library(ggsci)
# library(ggupset)
# library(VennDiagram)
# library(ggvenn)
# library(ggtern)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))



setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/")

################################################################################
model = c("AngII", "Salt-sensitive", "Spontaneous"); names(model)=c("mouse", "rat.ss", "rat.sp")

e = readRDS("cellchat.rds")
cc_df = c()
for(i in 1:length(e)){
  
  cellchat = get(names(e)[i], e)
  tissue = strsplit(names(e)[i], "_")[[1]][1]
  strain = gsub("\\.", "/", strsplit(names(e)[i], "_")[[1]][2])
  treatment = strsplit(names(e)[i], "_")[[1]][3]
  
  LR_list = rownames(cellchat@LR$LRsig)
  
  lr_df = c()
  for(lr in LR_list){
    prob_tmp = cellchat@net$prob[, , lr]
    
    lr_df_tmp <- reshape2::melt(prob_tmp, value.name = "value")
    colnames(lr_df_tmp)[1:2] <- c("source", "target")
    lr_df_tmp$LR_pair = lr
    
    lr_df = rbind(lr_df, lr_df_tmp[lr_df_tmp$value>0,])
    
  }
  
  lr_df$model = model[gsub("(.*)\\.([A-Z]+).*", "\\1", tissue, perl = T)]
  lr_df$tissue = gsub("(.*)\\.([A-Z]+).*", "\\2", tissue, perl = T)
  lr_df$strain = strain
  lr_df$treatment = treatment
  lr_df = cbind(lr_df, cellchat@LR$LRsig[lr_df$LR_pair, ])
  
  cc_df = rbind(cc_df, lr_df)
}

write.table(cc_df, "cellchat.result.out", col.names = T, row.names = F, sep = '\t', quote = F)


################################################################################
### merge results
# file_list = list.files(pattern = ".tsv")
# 
# cc_df = c()
# for(i in file_list){
#   dat_tmp = read.table(i, header = T)
#   cc_df = rbind(cc_df, dat_tmp)
# }

cc_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.result.out", header = T, sep = '\t')

cc_df$model = factor(cc_df$model, levels=c("AngII", "Salt-sensitive", "Spontaneous"))
cc_df$strain = factor(cc_df$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
cc_df$tissue = factor(cc_df$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA"))
cc_df$treatment = factor(cc_df$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w"))
cc_df$sxt = paste0(cc_df$strain, "-", cc_df$treatment)

id_vars = setdiff(colnames(cc_df), c("sxt", "value"))
cc_df_reshape = reshape(cc_df, idvar = id_vars, timevar = "sxt", direction = "wide")
cc_df_reshape[is.na(cc_df_reshape)] <- 0

cc_df_reshape$`diff.C57BL/6-AngII 3d` = cc_df_reshape$`value.C57BL/6-AngII 3d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.C57BL/6-AngII 28d` = cc_df_reshape$`value.C57BL/6-AngII 28d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.SS-HS 3d` = cc_df_reshape$`value.SS-HS 3d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SS-HS 21d` = cc_df_reshape$`value.SS-HS 21d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SD-HS 3d` = cc_df_reshape$`value.SD-HS 3d` - cc_df_reshape$`value.SD-LS`
cc_df_reshape$`diff.SHR-26w` = cc_df_reshape$`value.SHR-26w` - cc_df_reshape$`value.SHR-10w`
cc_df_reshape$`diff.WKY-26w` = cc_df_reshape$`value.WKY-26w` - cc_df_reshape$`value.WKY-10w`

diff_cols = colnames(cc_df_reshape)[grepl("diff", colnames(cc_df_reshape))]
dis_cols = c("treatment", colnames(cc_df_reshape)[grepl("value.", colnames(cc_df_reshape))])
id_vars = setdiff(colnames(cc_df_reshape), c(dis_cols, diff_cols))
cc_df_diff = reshape2::melt(cc_df_reshape[, c(id_vars, diff_cols)], id.vars = id_vars, measured.vars=diff_cols,
                  variable.name = "sxt")
cc_df_diff$sxt = gsub("diff.", "", cc_df_diff$sxt)
cc_df_diff$treatment = as.character(lapply(strsplit(cc_df_diff$sxt, "-"), function(x) x[2]))
cc_df_diff$treatment = factor(cc_df_diff$treatment, c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w"))
cc_df_diff = cc_df_diff[cc_df_diff$value!=0, ]


### all changes
cc_df_use = cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
pathway_selected = c("NRXN", "NCAM", "NEGR", "PTN") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & cc_df_use$tissue=="HYP", ]) +
  geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") + # , limits=c(-5,5)
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Differential\ninteraction\nstrength") +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")

# ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.cellchat.heatmap.png", width=494/96, height=227/96, dpi=300)

pathway_selected = c("NRXN", "NCAM", "NEGR", "PTN") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
cc_df_use = cc_df_diff %>% filter(tissue=="HYP" & pathway_name %in% pathway_selected) %>%
  group_by(sxt, model, strain, tissue, treatment, source, target, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

ggplot(cc_df_use) +
  geom_tile(mapping = aes(x = fct_inorder(treatment), y = source, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Communication\nprobability", title = pathway_selected) +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")



### EC relevant communication (sender and receiver)
cc_df_use = cc_df_diff %>% filter(grepl("EC", source) | grepl("EC", source)) %>%
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
ggplot(cc_df_use[abs(cc_df_use$prob_sum)>0.5,]) +
  geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red", limits=c(-5,5)) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Communication\nprobability") +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")



cc_df_use = cc_df_diff %>% filter(grepl("EC", source)) %>%
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

cc_df_use = cc_df_diff %>% filter(grepl("EC", target)) %>%
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
ggplot(cc_df_use[abs(cc_df_use$prob_sum)>0.05,]) +
  geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red", limits=c(-5,5)) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Communication\nprobability") +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")


# pathway_selected = "PTPRM" # SLIT; COLLAGEN
pathway_selected = c("PTPRM") # "SEMA3", "SEMA6", "ADGRL", "APP", "CADM", "PTPRM", "UNC5"
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected,]) +
  geom_tile(mapping = aes(x = fct_inorder(treatment), y = target, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient(low="white", high="red") +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Communication\nprobability", title = pathway_selected) +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")




### EC relevant communication (sender and receiver)
cc_df_use = cc_df_diff %>% 
  filter(tissue=="LK") %>%
  filter(grepl("EC_Mecom", source) | grepl("EC_Mecom", source)) %>%
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

cc_df_use[cc_df_use$prob_sum>5,]$prob_sum <- 5; cc_df_use[cc_df_use$prob_sum< -5]$prob_sum <- -5
cc_df_use %>% 
  # filter(abs(cc_df_use$prob_sum)>0.5) %>%
  ggplot() +
  geom_tile(mapping = aes(x = treatment, y = pathway_name, fill=prob_sum)) +
  # geom_text(mapping = aes(x = comp, y = cell_type, label = annotation), color = "black", vjust = 0.75) +  # Add asterisks for significance
  scale_y_discrete(limits=rev) +
  scale_fill_gradient2(low="blue", high="red", limits=c(-5,5)) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x="", y="", fill="Communication\nprobability") +
  facet_nested(tissue ~ model + strain , scales = "free", space = "free")

