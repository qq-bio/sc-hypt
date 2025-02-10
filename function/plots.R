# .libPaths("/home/qqiu/R/x86_64-pc-linux-gnu-library/4.2")
library(ggh4x)
library(plyr)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(openxlsx)
library(ggsci)
library(cowplot)
library(ggthemes)
library(scales)
library(ggalluvial)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))
options(digits=2)

species_col = pal_jama()(5)
names(species_col) = c("C57BL/6", "SHR", "WKY", "SS", "SD")


# ## color code
# ### type
# type_col = tableau_color_pal("Superfishel Stone")(3) # pal_nejm("default")(3)
# names(type_col) = c("NC", "DN", "HN")
# sample_col = c(colorRampPalette(c("#ffae34", "#FFFFFF"))(8)[1:5],
#                colorRampPalette(c("#ef6f6a", "#FFFFFF"))(8)[1:6],
#                colorRampPalette(c("#6388b4", "#FFFFFF"))(6)[1:4]
# )
# names(sample_col) = c(paste0("D", 1:5), paste0("H", 1:6), paste0("C", 1:3))
# ### cluster
# getPalette = colorRampPalette(brewer.pal(12, "Paired"))
# cell_col = getPalette(35)
# names(cell_col) = c("POD", "PEC", "PT-S1", "PT-S2", "Inj. PT", "TL", "TAL", "MD", "DCT", "CNT", "PC", "IC-A", "IC-B", "EC-PTC", "EC-GC", "EC-AEA", "EC-LYM", "MC", "FIB", "MYOF", "MB-IGHM", "MB-SSPN", "NB", 
#                     "PL", "T-CD4", "Treg", "T-CD8", "NKT", "MAST", "MAC-M2", "cDC", "pDC", "moDC", "ncMON", "NEU")
# ### gene expression: https://ouyanglab.com/singlecell/basic.html
# colGEX = c("grey85", brewer.pal(7, "Oranges"))


study_design_dotplot = function(){
  sample = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_sample_info.txt", header = T, sep = "\t")
  sample = unique(sample[,-1])
  sample = sample[!(sample$tissue=="Middle cerebral artery" & sample$strain %in% c("C57BL/6", "SS")), ]
  sample = sample[!(sample$seqID2 %in% c("RLV6SN ", "RMSA4SN")), ]
  
  cell_count = c("Hypothalamus"="237,815", "Middle cerebral artery"="12,563",
                 "Cardiac left ventricle"="164,839", "Kidney"="179,637", 
                 "3rd mesenteric artery"="18,130"
                 )
  
  # new.tissue <- c("3rd Mesenteric Artery", "Cardiac Left Ventricle", "Hypothalamus", "Middle Cerebral Artery", "Kidney"
  # )
  # names(new.tissue) = c("3rd mesenteric artery", "cardiac left ventricle", "hypothalamus", "middle cerebral artery", "kidney")
  # 
  sample$tissue = factor(sample$tissue, 
                         levels = c("3rd mesenteric artery", "Kidney",
                                    "Cardiac left ventricle", "Middle cerebral artery", "Hypothalamus"
                         ))
  sample$treatment = factor(sample$treatment, 
                            levels = c("Saline 3d", "AngII 3d", "AngII 28d",
                                       "10w", "26w", "LS", "HS 3d", "HS 21d"))
  sample$strain = factor(sample$strain, levels=c("C57BL/6", "SD", "SS", "WKY", "SHR"))
  sample$project = factor(c("AngII", "Spontaneous", "Spontaneous", "Salt-sensitive", "Salt-sensitive")[
    match(sample$strain, c("C57BL/6", "WKY", "SHR", "SD", "SS"))
  ], levels=c("AngII", "Salt-sensitive", "Spontaneous"))
  
  sample = sample %>% group_by(strain, tissue, treatment, project) %>% 
    dplyr::summarise(Count = n()) %>% as.data.frame()
  
  sample$cell_count = factor(cell_count[as.character(sample$tissue)], levels = cell_count)
  sample$Assay = "sn-RNA"
  sample[sample$tissue=="Kidney",]$Assay = "sn-multiome"
  
  # sample[sample$Count>2, ]$Count = 2
  
  min_value=10000;  max_value=250000; palette <- c("white", "#2F79B5")
  color_function <- col_numeric(palette = palette, domain = c(min_value, max_value))
  ridiculous_strips <- strip_nested(
    background_y = elem_list_rect(
      fill = color_function(as.numeric(gsub(",", "", cell_count)))
    ),
    text_y = elem_list_text(angle = c(0)),
    by_layer_y = FALSE
  )
  
  ggplot(data = sample, aes(y = tissue, x = treatment, color = strain, shape = Assay,
                            size = as.factor(Count), group = interaction(treatment, strain))) + 
    geom_point(#aes(color = strain),
               position = position_dodge(width=0.75)) +
    labs(title="", color='Strain', size="Sample size",
         x='Treatment', y='') +
    theme(legend.position = "right",legend.direction = "vertical",
          legend.box="horizontal",
          text = element_text(size=10, colour = 'black'),
          axis.text.y = element_text(size=10, colour = 'black'),
          axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
          axis.title = element_text(size=10, colour = 'black'),
          strip.text = element_text(size=10)) +
    # scale_x_discrete() + 
    # coord_flip() +
    scale_color_manual(values = species_col) + 
    scale_size_manual(values = c("1"=1, "2"=3, "3"=4), labels = c("N=1", "N=2", "N=3")) +
    scale_shape_manual(values = c(19, 17)) + 
    facet_nested(cell_count ~ project + strain, scales = "free", strip = ridiculous_strips)

  ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/sample_comp.png", width=931/96, height=279/96, dpi=300)
  
}


cell_color = function(meta_table, cluster){
  # set the color based on major clusters
  clusters = levels(meta_table[, cluster])
  major_clusters = table(gsub("\\(.+\\)", "", clusters))[unique(gsub("\\(.+\\)", "", clusters))]
  cluster_number = length(major_clusters)
  # col_category = paste0("category20", round(cluster_number/10+0.5), "0")
  getPalette = colorRampPalette(pal_d3("category20")(20))
  
  
  cell_col = c()
  for(i in 1:length(major_clusters)){
    cell_pal = getPalette(length(major_clusters))
    cell_col = c(cell_col,
                 if(major_clusters[i]==1) cell_pal[i] else colorRampPalette(c(cell_pal[i], "#FFFFFF"))(major_clusters[i]+1)[1:major_clusters[i]]
    )
  }
  names(cell_col) = clusters
  
  return(cell_col)
  # show_col(cell_col)
  
}

umap_dotplot = function(seurat_object,
                        marker_list,
                        cluster = "new.cluster.ids_umap",
                        reduction = "UMAP"){
  
  cell_col = cell_color(seurat_object@meta.data, cluster)
  
  Idents(seurat_object) = cluster
  UMAP_plot1 = DimPlot(seurat_object, label = F, reduction = reduction) + theme(legend.position = "none") + 
    scale_color_manual(values = cell_col)
  UMAP_plot2 = DimPlot(seurat_object, label = T, reduction = reduction) + theme(legend.position = "none") + 
    scale_color_manual(values = cell_col)
  Dot_plot = DotPlot(seurat_object, features=marker_list, group.by=cluster) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x="", y="")
  
  print(UMAP_plot1)
  print(UMAP_plot2)
  # return(UMAP_plot2)
  print(Dot_plot)
  # return(Dot_plot)
}


stacked_barplot = function(meta_table, 
                           sample = "seqID2", 
                           cluster = "new.cluster.ids_umap", 
                           species = NULL,
                           project = NULL,
                           treatment = "treatment"){
  
  dat = meta_table[, c(sample, cluster, species, project, treatment)]
  pct = dat %>% add_count(get(sample)) %>% 
    group_by_all() %>% 
    dplyr::summarise(Count = n()) %>% ungroup() %>% 
    mutate(percentage = Count/n) %>% as.data.frame()
  colnames(pct)[1:2] = c("Sample", "Cluster")
  
  cell_col = cell_color(meta_table, cluster)
  
  # pdf("stackedbarplot.pdf", width=4, height = 4)
  pctBar = ggplot(data = pct, aes(y = percentage, x = Sample, fill = Cluster)) + 
    geom_bar(stat="identity", 
             position = position_fill(reverse = TRUE)) +
    labs(title="", fill='',
         x='', y='Cell Proportion') +
    theme(legend.position = "right",
          plot.title = element_text(size=15, colour = 'black'),
          axis.text.y = element_text(size=10, colour = 'black'),
          axis.text.x = element_text(size=10, colour = 'black'),
          axis.title = element_text(size=15, colour = 'black'),
          strip.text = element_text(size=10)) +
    scale_x_discrete(limits = rev(levels(pct$Sample))) + 
    coord_flip() +
    scale_fill_manual(values = cell_col)
  
  if(!(is.null(species) | is.null(project))){
    pctBar = pctBar + facet_nested(project + species + treatment ~., scales = "free_y")
  }else if(! is.null(species)){
    pctBar = pctBar + facet_nested(species + treatment ~., scales = "free_y")
  }else if(! is.null(project)){
    pctBar = pctBar + facet_nested(project + treatment ~., scales = "free_y")
  }else{
    pctBar = pctBar + facet_grid(treatment ~., scales = "free_y")
  }
  
  print(pctBar)
  return(pctBar)
  # ggsave("stackedbarplot.pdf", pctBar, width=4, height = 4)
  
  # totalBar = ggplot(data=total, aes(x=Var1, y=Freq)) +
  #   geom_bar(stat="identity") +
  #   labs(title="", fill='',
  #        x='', y='Cell Number') +
  #   theme(axis.line.y=element_blank(), axis.ticks.y =element_blank(),
  #         axis.text.y=element_blank(),
  #         axis.text.x = element_text(size=15, angle = 45, hjust = 1, colour = 'black'),
  #         axis.title = element_text(size=15, colour = 'black'),
  #         strip.text = element_text(size=15)) + 
  #   # scale_x_discrete(limits = rev(levels(total$Var1))) +
  #   theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  # totalBar
  # 
  # plot_grid(pctBar, NULL, totalBar, align="hv", ncol = 1,
  #           rel_widths=c(2,-0.1, 1))
  # dev.off()
}






stacked_barplot_v2 = function(meta_table, 
                           # sample = "seqID2", 
                           cluster = "subclass_level2", 
                           species = "species",
                           project = "project",
                           treatment = "treatment",
                           tissue = "tissue"){
  
  meta_table = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/multi-HYP.meta.out")
  
  cluster = "subclass_level1"
  strain = "strain"
  project = "project"
  treatment = "treatment"
  tissue = "tissue"
  
  # meta_table[meta_table$subclass_level2=="Pars tuberalis cells", ]$class = "endocrine cells"
  # meta_table[meta_table$subclass_level2=="ABC", ]$class = "endothelial cells"
  
  dat = meta_table[, c(cluster, strain, project, treatment, tissue)]
  pct = dat %>% group_by(strain, treatment, tissue) %>% 
    mutate(n=n()) %>% 
    ungroup() %>% 
    group_by_all() %>% 
    summarise(Count = n()) %>% ungroup() %>% 
    mutate(percentage = Count/n, id = row_number()) %>% as.data.frame()
  
  colnames(pct)[1] = c("Cell_type")
  pct$Cell_type = factor(pct$Cell_type, level=cell_order)
  pct$treatment = factor(pct$treatment, level=treatment_order)
  pct$tissue = factor(pct$tissue, level=tissue_order)
  pct$strain = factor(pct$strain, level=strain_order)
  
  # pdf("stackedbarplot.pdf", width=4, height = 4)
  pctBar = ggplot(data = pct, aes(y = percentage, x = treatment, fill = Cell_type, alluvium = Cell_type)) + 
    geom_bar(stat="identity", width = .7) +
    labs(title="", fill='',
         x='', y='Cell composition') +
    theme(legend.position = "None",
          plot.title = element_text(size=15, colour = 'black'),
          axis.text.y = element_text(size=10, colour = 'black'),
          axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
          axis.title = element_text(size=15, colour = 'black'),
          strip.text = element_text(size=10)) +
    # scale_x_discrete(limits = rev(levels(pct$Sample))) + 
    scale_x_discrete(limits=rev) +
    geom_flow(alpha= 0.3, lty = 1, color = "white",
              curve_type = "linear", width = .7) +
    coord_flip() +
    scale_fill_manual(values = cell_col)
  
  pctBar = pctBar + facet_nested(project + strain ~ tissue, scales = "free_y", space = "free_y")
  print(pctBar)
  
  ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/major_cluster.stacked_barplot.png", width=755/96, height=461/96, dpi=300)
  
  return(pctBar)
}





da_dot_v2 = function(meta_table, 
                              # sample = "seqID2", 
                              cluster = "subclass_level2", 
                              species = "species",
                              project = "project",
                              treatment = "treatment",
                              tissue = "tissue"){
  
  meta_table = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/multi-HYP.meta.out")
  
  percent_df <- meta_table %>% 
    filter(!(subclass_level1=="Microglia" & tissue %in% c("LK", "LV"))) %>% 
    # Calculate counts and initial proportions for each subclass within each unique sample
    group_by(project, strain, treatment, tissue, seqID2) %>%
    dplyr::count(subclass_level1) %>%
    mutate(count = n,
           total = sum(n), 
           proportion = n / total) %>%
    # Create a combined factor column that categorizes data by strain and treatment
    mutate(combined = factor(interaction(strain, treatment, sep = "-"),
                             levels = c("C57BL/6-Saline 3d", "C57BL/6-AngII 3d", "C57BL/6-AngII 28d",
                                        "SS-LS", "SS-HS 3d", "SS-HS 21d", 
                                        "SD-LS", "SD-HS 3d",
                                        "SHR-10w", "SHR-26w", 
                                        "WKY-10w", "WKY-26w"))) %>%
    ungroup() %>%
    group_by(project, strain, treatment, tissue, combined, subclass_level1) %>%
    summarise(mean_proportion = mean(proportion), 
              mean_count = mean(count), .groups = 'drop') %>%
    dplyr::select(project, strain, treatment, tissue, combined, subclass_level1, mean_count, mean_proportion)
  
  percent_df_v2 = merge(percent_df[!(percent_df$combined %in% c("C57BL/6-Saline 3d", "SS-LS", "SD-LS", "SHR-10w", "WKY-10w")), ],
                        percent_df[percent_df$combined %in% c("C57BL/6-Saline 3d", "SS-LS", "SD-LS", "SHR-10w", "WKY-10w"), ],
                        by = c("project", "strain", "tissue", "subclass_level1"), all.x=T
  )
  
  percent_df_v2 <- percent_df_v2 %>%
    mutate(mean_proportion.x = ifelse(is.na(mean_proportion.x), 0.0001, mean_proportion.x),
           mean_proportion.y = ifelse(is.na(mean_proportion.y), 0.0001, mean_proportion.y))
  
  percent_df_v2 <- percent_df_v2 %>%
    mutate(fold_change = log2(mean_proportion.x / mean_proportion.y),
           cell_count = (mean_count.x + mean_count.y)/2,
           cell_type_significance = case_when(
             is.na(mean_proportion.x) & !is.na(mean_proportion.y) ~ "Present in Control Only",
             !is.na(mean_proportion.x) & is.na(mean_proportion.y) ~ "Present in Treatment Only",
             TRUE ~ "Present in Both"
           )) %>%
    mutate(capped_log2FC = pmin(pmax(fold_change, -5), 5))
  
  percent_df_v2$cell_type = factor(percent_df_v2$subclass_level1, levels=unique(cell_order))
  percent_df_v2$tissue = factor(percent_df_v2$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
  
  species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))
  alpha <- setNames(c(1, 0.5, 0.5, 1, 1), unique(percent_df_v2$treatment.x))
  
  combined <- unique(percent_df_v2$combined.x)
  combined_colors <- setNames(species_col[gsub("\\-.*", "", combined)], combined)
  combined_alphas <- setNames(alpha[gsub(".*\\-", "", combined)], combined)
  
  ggplot(percent_df_v2, aes(x = capped_log2FC, y = cell_type)) +
    # geom_segment(aes(yend = cell_type), xend = 0, colour = "grey50") +
    geom_point(aes(colour = combined.x, alpha = combined.x, size=log10(cell_count))) +
    scale_y_discrete(limits=rev) +
    scale_colour_manual(values = combined_colors) +
    scale_alpha_manual(values = combined_alphas) +
    scale_size_continuous(breaks = c(2, 3, 4), labels = c("100", "1,000", "10,000")) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),   # No horizontal grid lines
      # legend.position = c(1, 0.55),           # Put legend inside plot area
      legend.justification = c(1, 0.5),
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      # text = element_text(colour = 'black'),
      strip.text = element_text(colour = 'black'),
      strip.background = element_rect(colour = "black", fill = NA)
    ) +
    labs(x="log2 FC in cell proportion", y="", color="Strain x treatment", alpha="Strain x treatment", size="Number of cells") +
    facet_grid2(cols = vars(project), rows = vars(tissue), 
                scales = "free", space = "free_y") +
    geom_vline(xintercept = c(-1,  1), linetype = "solid", color = "black") +
      geom_vline(xintercept = c(0), linetype = "dashed", color = "black")
  
}












library(tidyverse)
library(ggrepel)
library(emmeans)
library(ggsignif)

DA_plot = function(meta_table, 
                   sample = "seqID2", 
                   cluster = "new.cluster.ids_umap", 
                   treatment = "treatment",
                   species = "species", 
                   project = NULL,
                   outfile){
  
  ### data process
  sampleInfo = data.frame(Sample_ID=meta_table[, sample], 
                          treatment=meta_table[, treatment], 
                          Cluster=meta_table[, cluster], 
                          species=meta_table[, species],
                          Project=meta_table[, project])
  sampleCount = sampleInfo %>% add_count(Sample_ID) %>% 
    group_by_all() %>% dplyr::summarise(Count = n()) %>% 
    ungroup() %>% complete(Sample_ID, Cluster, fill = list(Count = 0)) %>%  fill(n) %>%
    mutate(Cluster = as.character(Cluster)) %>%
    as.data.frame() 
  colnames(sampleCount)[ncol(sampleCount)-1] = 'Total'
  sampleCount["Other"] = sampleCount$Total - sampleCount$Count
  sampleCount["Freq"] = sampleCount$Count/sampleCount$Total
  
  
  ### visualization 
  pdf(outfile, height=3, width=5)
  for(i in unique(sampleCount$Cluster)){
    subtitle = i
    p = ggplot(aes(x = treatment, y = (Count / Total)*100, 
                   color = species), data = sampleCount %>% filter(Cluster == i)) +
      geom_point() + 
      theme_classic() + 
      theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
            axis.text.y = element_text(colour = 'black'),
            legend.text = element_text(size=8),
            legend.position = "right") +
      labs(subtitle = subtitle, x = "", y="Proportion (%)", color="Strain") +
      scale_y_continuous(labels = function(x) sprintf("%.2f", x)) +
      scale_color_manual(values = species_col)
    
    if(length(unique(sampleCount$species))>1){
      p = p + facet_nested(. ~ species, scales = "free_x") #free_x
    }
    
    if(length(unique(sampleCount$Project))>1){
      p = p + facet_nested(. ~ Project + species, scales = "free_x") #free_x
    }
    
    print(p)
  }
  dev.off()
  
}




enrichment_dotplot = function(enrich_res, n=10, outfile){
  
  pathway_list = enrich_res[enrich_res$p.adjust<0.05, ] %>% 
    filter(species %in% c("C57BL/6", "SHR", "SS")) %>% 
    group_by(project, treatment, species) %>% top_n(-p.adjust, n=n) %>%
    as.data.frame()
  dat_use = enrich_res[enrich_res$p.adjust<0.05, ] %>% 
    filter(Description %in% pathway_list$Description) %>%
    as.data.frame()
  dat_use$species = factor(dat_use$species, levels=c("C57BL/6", "WKY", "SHR", "SD", "SS"))
  dat_use$treatment = factor(dat_use$treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d",
                                                           "10w", "26w", "LS", "HS 3d", "HS 21d"))
  dat_use$project = factor(dat_use$project, levels = c("AngII", "Salt-sensitive", "Spontaneous"))
  
  pctDot = ggplot(data = dat_use, aes(y = Description, x = treatment, 
                                      size = -log(p.adjust, 10), color = species,
                                      group = interaction(treatment, species))) + 
    geom_point(position = position_dodge(width=0.75)) +
    labs(title="", x='', y='', 
         size='-log10(p-adjust)',
         color='Strain'
    ) +
    theme(axis.text.y = element_text(size=10, colour = 'black'),
          axis.text.x = element_text(size=10, angle = 45, hjust = 1, colour = 'black'),
          strip.text = element_text(size=10)) +
    scale_color_manual(values = species_col)
  pctDot = pctDot + facet_nested(~ project + paste("vs.", control), scale="free_x")
  
  print(pctDot)
  return(dat_use)
}
