
library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(CellChat)
library(ggh4x)
library(patchwork)


source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")






################################################################################
### GSEA score and visualization (Figure 4a, b)
################################################################################
# functions
make_rank_score <- function(logfc, pval, p_floor = 1e-300) {
  pval2 <- pmax(pval, p_floor)
  -1 * sign(logfc) * (-log10(pval2))
}

get_msigdb_pathways <- function(species = "Mus musculus",
                                category = "H",
                                subcategory = NULL,
                                min_geneset_size = 10,
                                max_geneset_size = 500) {
  
  msig <- msigdbr(species = species, category = category)
  if (!is.null(subcategory)) {
    msig <- msig %>% filter(.data$gs_subcat == subcategory)
  }
  
  pathways <- msig %>%
    distinct(gs_name, gene_symbol) %>%
    group_by(gs_name) %>%
    dplyr::summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
    mutate(size = lengths(genes)) %>%
    filter(size >= min_geneset_size, size <= max_geneset_size) %>%
    select(gs_name, genes)
  
  set_names(pathways$genes, pathways$gs_name)
}

make_ranks_for_group <- function(df_group) {
  df_group %>%
    group_by(gene_name) %>%
    dplyr::summarise(rank_score = rank_score[which.max(abs(rank_score))], .groups = "drop") %>%
    arrange(desc(rank_score)) %>%
    deframe()
}

run_fgsea_one <- function(ranks, pathways, minSize = 10, maxSize = 500) {
  fgsea::fgsea(
    pathways = pathways,
    stats    = ranks,
    minSize  = minSize,
    maxSize  = maxSize,
    nproc    = max(1, parallel::detectCores() - 1)
  ) %>%
    as_tibble() %>%
    arrange(padj, desc(abs(NES)))
}

parse_group_id <- function(group_id, sep = " \\| ") {
  parts <- str_split(group_id, sep, simplify = TRUE)
  as_tibble(parts) %>% setNames(group_vars)
}



deg_l1 <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.rename.out", sep = "\t", header = T)

pathways_H <- get_msigdb_pathways(category = "H")
pathways_reactome <- get_msigdb_pathways(category = "C2", subcategory = "CP:REACTOME")

group_vars <- c("project","strain","tissue","cell_type","control","treatment")

deg_use <- deg_l1 %>%
  dplyr::mutate(
    gene_name = as.character(gene_name),
    rank_score = make_rank_score(avg_log2FC, p_val),
    group_id = do.call(paste, c(across(all_of(group_vars)), sep = " | "))
  ) %>%
  filter(!is.na(rank_score), is.finite(rank_score), !is.na(gene_name), gene_name != "")

ranks_list <- deg_use %>%
  group_by(group_id) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$group_id))) %>%
  map(make_ranks_for_group)

meta_tbl <- tibble(group_id = names(ranks_list)) %>%
  bind_cols(map_dfr(.$group_id, parse_group_id))

fgsea_H_long <- imap_dfr(
  ranks_list,
  ~ run_fgsea_one(.x, pathways_H) %>% mutate(group_id = .y, .before = 1)
)

fgsea_H_long <- fgsea_H_long %>%
  left_join(meta_tbl, by = "group_id") %>%
  mutate(
    leadingEdge = map_chr(leadingEdge, ~ paste(.x, collapse = ";"))
  )

write.table(fgsea_H_long, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/fgsea/DEG.L1.fgsea.hallmark.out", sep = "\t", col.names = T, row.names = F)


htn_strains <- c("C57BL/6", "SHR", "SS")

replication_tbl_grouped <- fgsea_H_long %>%
  group_by(pathway, tissue, cell_type, strain) %>%
  dplyr::summarise(
    NES_strain = median(NES, na.rm = TRUE),
    padj_min = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    strain_group = ifelse(strain %in% htn_strains, "HTN", "NT"),
    dir = case_when(
      NES_strain > 0 ~ "up",
      NES_strain < 0 ~ "down",
      TRUE ~ "zero"
    )
  ) %>%
  group_by(pathway, tissue, cell_type) %>%
  dplyr::summarise(
    n_strain_total = n_distinct(strain),
    
    n_htn = sum(strain_group == "HTN"),
    n_nt  = sum(strain_group == "NT"),
    
    median_NES_htn = ifelse(n_htn > 0, median(NES_strain[strain_group=="HTN"], na.rm = TRUE), NA_real_),
    median_NES_nt  = ifelse(n_nt  > 0, median(NES_strain[strain_group=="NT"],  na.rm = TRUE), NA_real_),
    
    htn_up   = sum(strain_group=="HTN" & NES_strain > 0),
    htn_down = sum(strain_group=="HTN" & NES_strain < 0),
    nt_up    = sum(strain_group=="NT"  & NES_strain > 0),
    nt_down  = sum(strain_group=="NT"  & NES_strain < 0),
    
    frac_major_htn = ifelse(n_htn > 0, pmax(htn_up, htn_down)/n_htn, NA_real_),
    frac_major_nt  = ifelse(n_nt  > 0, pmax(nt_up,  nt_down )/n_nt,  NA_real_),
    
    flip_dir_median = !is.na(median_NES_htn) & !is.na(median_NES_nt) &
      sign(median_NES_htn) == -sign(median_NES_nt),
    
    delta_NES = median_NES_htn - median_NES_nt,
    abs_delta_NES = abs(delta_NES),
    
    any_sig_htn = any(strain_group=="HTN" & padj_min < 0.25, na.rm = TRUE),
    any_sig_nt  = any(strain_group=="NT"  & padj_min < 0.25, na.rm = TRUE),
    
    .groups = "drop"
  )

# concordance barplot (Figure 4a)
replication_tbl_grouped %>% 
  filter(!is.na(frac_major_htn)) %>%
  mutate(label = paste0(frac_major_htn * n_htn, "/", n_htn)) %>%
  dplyr::count(tissue, frac_major_htn, label) %>%
  ggplot(aes(x = tissue, y = n, fill = as.numeric(frac_major_htn))) +
  geom_col(position = "stack", color = "grey30") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black",
    family = "Arial"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "steelblue",
    limits = c(0, 1),
    name = "Hypertensive\nstrain\nconcordance"
  ) +
  labs(
    x = "Tissue",
    y = "Number of pathways"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(family = "Arial"),
    panel.grid.major.x = element_blank()
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fgsea.hallmark.htn.concordance.barplot.png", width=360/96, height=350/96, dpi = 300)



# concordance dotplot (Figure 4b)
replication_tbl_grouped %>% 
  filter(!is.na(frac_major_htn), n_htn>1, any_sig_htn) %>%
  mutate(label = paste0(frac_major_htn * n_htn, "/", n_htn)) %>%
  dplyr::count(tissue, frac_major_htn, label) %>%
  ggplot(aes(x = tissue, y = n, fill = as.numeric(frac_major_htn))) +
  geom_col(position = "stack", color = "grey30") +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "black",
    family = "Arial"
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "steelblue",
    limits = c(0.5, 1),
    name = "Hypertensive\nstrain\nconsistency\nratio"
  ) +
  labs(
    x = "Tissue",
    y = "Number of Hallmark pathway-cell type pairs"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(family = "Arial"),
    panel.grid.major.x = element_blank()
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fgsea.hallmark.htn.concordance.gt_1.barplot.png", width=285/96, height=350/96, dpi = 300)




################################################################################
### Score cell states of neurons and glia (Fig. S9a)
################################################################################
hyp_m = readRDS('/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds')
hyp_m$SxT = factor(hyp_m$SxT, levels = sxt_order)

intrinsic_pro_apoptosis <- c("Bax", "Bak1", "Pmaip1", "Bbc3", "Bid", "Apaf1", "Casp9")
stress_inducers <- c("Trp53", "Ddit3", "Gadd45a", "Gadd45b", "Atf4", "Trib3", "Casp12")
anti_apoptosis_markers <- c("Bcl2", "Bcl2l1", "Xiap", "Mcl1", "Birc5")
proliferation_list = c("Mki67", "Pcna", "Top2a", "Ccnb1", "Ccnb2", 'Cdk1', "Ube2c", "Cenpa", "Cenpf")

gene_score_list = list(
  intrinsic_pro_apoptosis_score = intrinsic_pro_apoptosis,
  stress_score = stress_inducers,
  anti_apoptosis_score = anti_apoptosis_markers,
  proliferation_score = proliferation_list
)

hyp_m_summary = calc_module_scores_and_summary(hyp_m, gene_score_list, c("Project", "Strain",  "Treatment", "SxT", "Cell_type_L1"))$summary_long

hyp_m_score = calc_module_scores_and_summary(hyp_m, gene_score_list, c("Project", "Strain",  "Treatment", "SxT", "Cell_type_L1"), return_cell_scores = T)$cell_scores  

score_merged <- hyp_m_summary %>%
  mutate(SxT = factor(SxT, levels = sxt_order),
         Cell_type_L1 = factor(Cell_type_L1, levels = cell_order))

ct_selected = c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
                "Astrocyte", "OPC", "Myelinating OL", 
                "Tanycyte", "Ependymal cell", "Microglia")

p1 = score_vis(score_merged, score_selected = c("intrinsic_pro_apoptosis_score"), project_selected = c("AngII"), cell_type_selected = ct_selected) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p2 = score_vis(score_merged, score_selected = c("stress_score"), project_selected = c("AngII"), cell_type_selected = ct_selected) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p3 = score_vis(score_merged, score_selected = c("anti_apoptosis_score"), project_selected = c("AngII"), cell_type_selected = ct_selected) + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
p4 = score_vis(score_merged, score_selected = c("proliferation_score"), project_selected = c("AngII"), cell_type_selected = ct_selected)

cowplot::plot_grid(
  p1, p2, p3, p4, 
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1, 1, 2)
)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.angii.score.dotplot.png", width=1600/96, height=550/96, dpi=300)









################################################################################
### Visualize DEGs by Occurrence Rank (Fig. S14b)
################################################################################
deg_merged <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", sep='\t', header=T)
deg_merged <- deg_merged[deg_merged$tissue == "HYP" & deg_merged$p_val_adj < 0.05 & abs(deg_merged$avg_log2FC) > 0.5, ]

deg_merged$cell_type <- factor(deg_merged$cell_type, levels = cell_order)
deg_merged$treatment <- factor(deg_merged$treatment, levels = treatment_order)
deg_merged$strain <- factor(deg_merged$strain, levels = strain_order)


selected_genes <- c("Fth1", "Rora", "Apoe", "Ptgds", "Csmd1", "App", "Ncam1", "Hsp90ab1", "Zbtb16", "Lars2")

p <- deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  group_by(gene_name) %>%
  mutate(deg_num = n()) %>% 
  ungroup() %>%
  select(gene_name, deg_num) %>%
  distinct() %>%
  arrange(desc(deg_num)) %>%
  mutate(
    index = row_number(),
    highlight = ifelse(gene_name %in% selected_genes, "Highlighted", "Normal")
  ) %>%
  arrange(desc(highlight)) 

# Plot DEG occurrences
p <- ggplot(plot_data, aes(x = deg_num, y = index, color = highlight)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("Highlighted" = "red", "Normal" = "black")) + 
  scale_y_continuous(trans = "reverse") +
  ggrepel::geom_label_repel(
    aes(label = ifelse(gene_name %in% selected_genes, gene_name, "")), 
    size = 4, 
    max.overlaps = Inf,
    force = 10,
    box.padding = 1,
    point.padding = 0.5,
    color = "black"
  ) +
  theme_classic() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "none",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Number of DEG Occurrences", y = "Rank of Genes by DEG Occurrence")

print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.reccurent_deg.dot.png", width=353/96, height=287/96, dpi=300)





################################################################################
### top DEGs in each cell type (Fig. S14c)
################################################################################
deg_merged = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", sep='\t', header=T)
deg_merged = deg_merged[deg_merged$tissue=="HYP", ]

deg_merged = deg_merged[deg_merged$p_val_adj<0.05 & abs(deg_merged$avg_log2FC)>0.5, ]
deg_merged$cell_type = factor(deg_merged$cell_type, levels = cell_order)
deg_merged$treatment = factor(deg_merged$treatment, levels = treatment_order)
deg_merged$strain = factor(deg_merged$strain, levels = strain_order)

cell_type=c("Astrocyte", "Microglia", "Activated microglia", "OPC")
for(ci in cell_type){
  
  for(pi in c("Salt-sensitive")){
    
    outfile = sprintf("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/%s.%s.top_gene.dotplot.png", ci, pi)
    outfile = gsub(" ", "_", outfile)
    
    top_genes <- deg_merged %>%
      filter(
        cell_type == ci,
        project == pi,
        strain %in% c("C57BL/6", "SS", "SHR")
      ) %>%
      group_by(gene_name) %>%
      slice_max(order_by = abs(avg_log2FC), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      mutate(direction = factor(ifelse(-avg_log2FC > 0, "up", "down"), levels = c("up", "down"))) %>%
      group_by(direction) %>%
      slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE) %>%
      ungroup() %>%
      arrange(direction, p_val_adj)
    
    gene_levels <- rev(top_genes$gene_name)
    
    p <- deg_merged %>%
      filter(
        cell_type == ci,
        project == pi,
        gene_name %in% top_genes$gene_name
      ) %>%
      mutate(
        gene_name = factor(gene_name, levels = gene_levels),
        p_val = ifelse(p_val == 0, 1.6e-303, p_val)
      ) %>%
      ggplot(aes(x = treatment, y = gene_name)) +
      geom_point(aes(colour = -avg_log2FC, size = -log10(p_val))) +
      scale_color_gradient2(
        midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab"
      ) +
      theme_bw() +
      theme(
        text = element_text(family = "Arial"),
        panel.grid.major.y = element_blank(),
        legend.justification = c(1, 0.5),
        axis.text.y = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
        strip.text = element_text(colour = "black"),
        strip.background = element_rect(colour = "black", fill = NA)
      ) +
      labs(x = "", y = "", color = "log2(FC)", size = "-log10(p-value)") +
      facet_nested(~ cell_type + strain, scales = "free", space = "free")
    
    print(p)
    ggsave(outfile, width=416/96, height=395/96, dpi = 300)
    
  }
}





################################################################################
### reccurrent DEGs (Fig. S14d)
################################################################################
deg_hyp = deg_merged[deg_merged$tissue=="HYP", ]

deg_processed <- deg_hyp %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR")) %>%
  mutate(
    sign_fc = case_when(
      avg_log2FC > 0  ~ 1,
      avg_log2FC < 0  ~ -1,
      TRUE            ~ 0
    )
  )

deg_count <- deg_processed %>%
  distinct(gene_name, cell_type) %>%
  group_by(gene_name) %>%
  summarise(deg_num = n(), .groups = "drop")

final_summary <- deg_count %>%
  left_join(strain_wide, by = "gene_name") %>%
  arrange(desc(deg_num))


selected_consistent <- final_summary %>%
  filter(deg_num >= 7,
         ! is.na(cross_strain_consistent),
         cross_strain_consistent == TRUE,
         `pct_same_sign_C57BL/6` >= 0.7,
         pct_same_sign_SS >= 0.7, 
         pct_same_sign_SHR >= 0.7)
selected_gene <- unique(selected_consistent$gene_name)
title = "Strain-concordant recurrent DEGs"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.recurrent_deg.concordant.dot.png"


selected_inconsistent <- final_summary %>%
  filter(deg_num >= 7,
         ! is.na(cross_strain_consistent),
         cross_strain_consistent == FALSE,
         `pct_same_sign_C57BL/6` >= 0.7,
         pct_same_sign_SS >= 0.7, 
         pct_same_sign_SHR >= 0.7) %>%
  slice_head(n=10)
selected_gene <- unique(selected_inconsistent$gene_name)
title = "Strain-divergent recurrent DEGs"
outfile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.recurrent_deg.divergent.dot.png"








################################################################################
### Visualize Consensus Changes in snRNA-seq & Bulk RNA-seq (Fig. S14e)
################################################################################

bulk <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/bulk/DEG_bulk.csv", header = TRUE, sep = ",")

selected_genes <- c("Snhg11", "Meg3", "Sntg1", "Zbtb16", "Csmd1", "ENSRNOG00000065867",
                    "Cacna1c", "Lars2", "Gpm6b", "Ncam1")

# Plot snRNA-seq differential expression
p1 <- deg_merged %>%
  filter(strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_genes) %>%
  mutate(gene_name = factor(gene_name, levels = selected_genes)) %>%
  ggplot(aes(x = -avg_log2FC, y = gene_name)) +
  geom_point(aes(colour = -avg_log2FC, size = -log10(p_val_adj))) +
  scale_y_discrete(limits = rev) +
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical",
    axis.text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "snRNA-seq", x = "Cell type-wise log2(FC)", y = "", color = "log2(FC)", size = "-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(color = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

# Plot bulk RNA-seq differential expression
p2 <- bulk %>%
  filter(tissue == "HYP", strain %in% c("C57BL/6", "SS", "SHR"),
         gene_name %in% selected_genes) %>%
  mutate(gene_name = factor(gene_name, levels = selected_genes),
         strain = factor(strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))) %>%
  ggplot(aes(x = log2FoldChange, y = gene_name)) +
  geom_vline(xintercept = 0, colour = "black") +
  geom_point(shape = 21, aes(fill = log2FoldChange, size = -log10(padj))) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(midpoint = 0, low = "blue", mid = "white", high = "red", space = "Lab") +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left", 
    legend.box = "vertical",
    axis.text.y = element_blank(),
    axis.text.x = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black'),
    strip.background = element_rect(colour = "black", fill = NA)
  ) +
  labs(title = "Bulk RNA-seq", x = "log2(FC) in bulk", y = "", fill = "log2(FC)", size = "-log10(p-adj)") +
  facet_nested(~ strain, scales = "free") +
  guides(fill = guide_legend(nrow = 1), size = guide_legend(nrow = 1))

p1+p2+ plot_layout(widths = c(1,1))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/HYP.sn_bulk.consist_expr.png", width=832/96, height=356/96, dpi=300)







################################################################################
### Pathway Enrichment Visualization (Figure 2d)
################################################################################
merge_reshape = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/metascape.merge_all.reshape.out", header = T, sep = '\t', quote = "")
merge_reshape = merge_reshape[merge_reshape$tissue=="HYP",]

sxt_rename <- c(
  "Log.q.value..C57BL.6.AngII.3d" = "C57BL/6 - AngII 3d", 
  "Log.q.value..C57BL.6.AngII.28d" = "C57BL/6 - AngII 28d", 
  "Log.q.value..SS.HS.3d" = "SS - HS 3d", 
  "Log.q.value..SS.HS.21d" = "SS - HS 21d",
  "Log.q.value..SHR.26w" = "SHR - 26w"
)

melted_data <- merge_reshape %>%
  filter(category != "Other", summary == "Yes", model_count > 1) %>%
  pivot_longer(cols = starts_with("Log.q.value"), names_to = "condition", values_to = "log_q_value") %>%
  filter(condition %in% names(sxt_rename)) %>%
  mutate(
    condition = sxt_rename[condition],
    log_q_value = ifelse(log_q_value < 1.3, 0, log_q_value),
    pathway = fct_reorder(pathway, cell_type),
    cell_type = fct_inorder(cell_type),
    condition = factor(condition, levels = c("C57BL/6 - AngII 3d", "C57BL/6 - AngII 28d", "SS - HS 3d", "SS - HS 21d", "SHR - 26w"))
  )


# Load group info and merge with data
group_info <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/metascape/HYP.multi_model.pathway.group.out", header = FALSE, sep = '\t')
colnames(group_info) <- c("pathway", "group")
group_info <- unique(group_info)
rownames(group_info) <- group_info$pathway
melted_data$group <- group_info[as.character(melted_data$pathway), "group"]

# Heatmap for pathway enrichment
p_data <- ggplot(melted_data %>% filter(log_q_value > 1.3), aes(x = condition, y = pathway, fill = log_q_value)) +
  geom_tile(color = "black") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient(low = "white", high = "red") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    plot.margin = margin(0, 0, 0, -2, "cm")
  ) +
  labs(x = "", y = "", fill = "-log(q-val)", title = "Shared Pathways in at Least Two Hypertension Models") +
  facet_grid(cols = vars(cell_type), scales = "free", space = "free")

# Bar plot for pathway groups
p_y <- ggplot(melted_data) +
  geom_col(aes(y = pathway, x = 1, fill = group), width = 0.8) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set1")) +
  scale_y_discrete(limits = rev) +
  theme_void() +
  theme(
    axis.text.y = element_text(colour = 'black'),
    legend.position = "right",
    plot.margin = margin(0, -2, 0, 0, "cm")
  ) +
  labs(fill = "Category")

p_y + p_data + plot_layout(widths = c(0.03, 1), guides = "collect")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2d.pathway.heatmap.png", width=1325/96, height=548/96, dpi=300)







################################################################################
### Visualize CellChat Results (Figure 2e and f)
################################################################################
cc_df <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.result.out", 
                    header = TRUE, sep = '\t')

cc_df$model <- factor(cc_df$model, levels = model_order)
cc_df$strain <- factor(cc_df$strain, levels = strain_order)
cc_df$tissue <- factor(cc_df$tissue, levels = tissue_order)
cc_df$treatment <- factor(cc_df$treatment, levels = treatment_order)
cc_df$sxt <- paste0(cc_df$strain, "-", cc_df$treatment)

id_vars <- setdiff(colnames(cc_df), c("sxt", "value"))
cc_df_reshape <- reshape(cc_df, idvar = id_vars, timevar = "sxt", direction = "wide")
cc_df_reshape[is.na(cc_df_reshape)] <- 0

# Calculate differential interaction strength for selected conditions
cc_df_reshape$`diff.C57BL/6-AngII 3d` <- cc_df_reshape$`value.C57BL/6-AngII 3d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.C57BL/6-AngII 28d` <- cc_df_reshape$`value.C57BL/6-AngII 28d` - cc_df_reshape$`value.C57BL/6-Saline 3d`
cc_df_reshape$`diff.SS-HS 3d` <- cc_df_reshape$`value.SS-HS 3d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SS-HS 21d` <- cc_df_reshape$`value.SS-HS 21d` - cc_df_reshape$`value.SS-LS`
cc_df_reshape$`diff.SD-HS 3d` <- cc_df_reshape$`value.SD-HS 3d` - cc_df_reshape$`value.SD-LS`
cc_df_reshape$`diff.SHR-26w` <- cc_df_reshape$`value.SHR-26w` - cc_df_reshape$`value.SHR-10w`
cc_df_reshape$`diff.WKY-26w` <- cc_df_reshape$`value.WKY-26w` - cc_df_reshape$`value.WKY-10w`

# Reshape data to long format for plotting
diff_cols <- grep("diff", colnames(cc_df_reshape), value = TRUE)
cc_df_diff <- reshape2::melt(cc_df_reshape[, c(id_vars, diff_cols)], 
                             id.vars = id_vars, measured.vars = diff_cols,
                             variable.name = "sxt")
cc_df_diff$sxt <- gsub("diff.", "", cc_df_diff$sxt)
cc_df_diff$treatment <- factor(sapply(strsplit(cc_df_diff$sxt, "-"), `[`, 2), levels = treatment_order)
cc_df_diff <- cc_df_diff[cc_df_diff$value != 0, ]

# Summarize the differential values for each pathway interaction
cc_df_use <- cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name, interaction_name_2) %>%
  dplyr::summarise(prob_sum = sum(value)) %>%
  as.data.frame()

# Define selected pathways and interactions for visualization
pathway_selected <- c("NRXN", "NCAM", "NEGR") 
interaction_name_list <- c("Ncam1  - Ncam2", "Negr1  - Negr1", "Nrxn1  - Nlgn1", "Nrxn3  - Nlgn1")
cc_df_use$name_use <- paste0(cc_df_use$pathway_name, ": ", cc_df_use$interaction_name_2)

# Generate heatmap of interaction strengths in selected pathways for HYP tissue
ggplot(cc_df_use[cc_df_use$pathway_name %in% pathway_selected & 
                   cc_df_use$interaction_name_2 %in% interaction_name_list & 
                   cc_df_use$tissue == "HYP", ]) +
  geom_tile(aes(x = treatment, y = name_use, fill = prob_sum)) +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(low = "blue", high = "red") + 
  theme(
    panel.grid.major.y = element_blank(), 
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(x = "", y = "", fill = "Differential\nInteraction\nStrength") +
  facet_nested(~ model + strain, scales = "free", space = "free")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/fig2e.HYP.cellchat.heatmap.png", 
       width = 596/96, height = 241/96, dpi = 300)





e = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.rds")

saline3d=e$`mouse.HYP_C57BL.6_Saline 3d`
angii3d=e$`mouse.HYP_C57BL.6_AngII 3d`
angii28d=e$`mouse.HYP_C57BL.6_AngII 28d`
ss_ls=e$rat.ss.HYP_SS_LS
ss_hs3d=e$`rat.ss.HYP_SS_HS 3d`
ss_hs21d=e$`rat.ss.HYP_SS_HS 21d`
sd_ls=e$rat.ss.HYP_SD_LS
sd_hs3d=e$`rat.ss.HYP_SD_HS 3d`
shr_10w=e$rat.sp.HYP_SHR_10w
shr_26w=e$rat.sp.HYP_SHR_26w
wky_10w=e$rat.sp.HYP_WKY_10w
wky_26w=e$rat.sp.HYP_WKY_26w


title = c("saline3d"="Saline 3d", "angii3d"="Ang II 3d")
for( i in c("saline3d", "angii3d") ){
  cellchat = get(i)
  outfile = paste0("fig2f.", i, "cellchat.png")
  
  png(filename = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/", outfile), 
      width = 786/96*300, height = 291/96*300, res = 300)
  
  par(mfrow=c(1, 3))
  
  for(j in c("NCAM", "NEGR", "NRXN") ){
    pathways.show <- j
    
    netVisual_aggregate(cellchat, signaling = pathways.show, 
                        layout = "circle", title.space = 1,
                        remove.isolate = T,
                        weight.scale = F, 
                        edge.weight.max = 1, 
                        edge.width.max = 8 
    )
    
    title(main = paste0(title[i], " - ", j))
  }
  
  dev.off()
  par(mfrow=c(1, 1))
}















