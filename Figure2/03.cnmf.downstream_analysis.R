suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(tidygraph)
  library(ggraph)
})

source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/Figure2/01.cnmf.functions.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/sc-hypt/utils/00.initial_setting.R")



################################################################################
### Load cNMF results (spectra and usage)
################################################################################

data_dir <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/per_celltype/"

ct_list <- list.files(data_dir)
ct_list <- setdiff(ct_list, c("Immune_cell"))

input_para <- "gene_spectra_score.*dt_0_1"
spectra_list <- list()
for( i in ct_list ){
  
  run_file <- paste(data_dir, i, i, sep="/")
  spectra_score_file <- list.files(run_file, input_para, full.names = T)
  spectra_score <- read.table(spectra_score_file, sep='\t', row.names=1, header=TRUE)
  
  spectra_list[[i]] <- list(spectra_score = spectra_score,
                            spectra_name  = paste(i, rownames(spectra_score), sep = "-"))
  
}


input_para <- "usages.*dt_0_1"
usage_list <- list()
for( i in ct_list ){
  
  run_file <- paste(data_dir, i, i, sep="/")
  usage_file <- list.files(run_file, input_para, full.names = T)
  usage <- read.table(usage_file, sep='\t', row.names=1, header=TRUE)
  
  rownames(usage) <- gsub("\\|.*", "", rownames(usage))
  colnames(usage) <- paste0(i, "-", gsub("X", "", colnames(usage)))
  usage_norm <- as.data.frame(t(apply(usage, 1, function(x) x / sum(x))))
  
  usage_list[[i]] <- usage_norm
  
}

e <- list(spectra_list = spectra_list,
          usage_list = usage_list)

saveRDS(e, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/rCNMF/CNMF.spectra_usage.rds")


### load cell metadata
metadata <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/merged_metadata.tsv", header = T, sep = "\t")
metadata = metadata[!(metadata$Tissue=="MCA" & metadata$Strain %in% c("C57BL/6", "SS")), ]
metadata$Cell_type_L1 = factor(metadata$Cell_type_L1, levels = cell_order)
metadata$Strain = factor(metadata$Strain, levels = strain_order)
metadata$Tissue = factor(metadata$Tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA", "PBMC"))
metadata$Treatment = factor(metadata$Treatment, levels = c("Saline 3d", "AngII 3d", "AngII 28d", "10w", "26w", "LS", "HS 3d", "HS 21d"))

# sample-wise phenotype data
metadata_bp <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/data/Multiomics_sample_bp.txt", header = T, sep = "\t", comment.char = "", na.strings = "#N/A")

sample_bp <- metadata_bp %>% 
  mutate(Seq.ID = gsub(" $", "", Seq.ID)) %>%
  group_by(Seq.ID) %>%
  dplyr::summarise(BP_mean = mean(BP_endpoint, na.rm=TRUE),
                   UACR_mean = mean(UACR, na.rm=TRUE),
                   kidney_fibrosis_KC_mean = mean(kidney_fibrosis_KC, na.rm=TRUE),
                   kidney_fibrosis_OM_mean = mean(kidney_fibrosis_OM, na.rm=TRUE),
                   mesentery_artery_wall_thickness_ratio_mean = mean(mesentery.artery.wall.thickness.ratio, na.rm=TRUE),
                   lv_fibrosis_mean = mean(lv_fibrosis, na.rm=TRUE),
                   cardiomyocyte_cross_sectional_area_mean = mean(cardiomyocyte_cross.sectional_area, na.rm=TRUE))
metadata_merged <- metadata %>% left_join(sample_bp, by = c("Sample_ID"="Seq.ID"))













################################################################################
### Program-phenotype association analysis
################################################################################

program_summary_bp <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "BP_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_uacr <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "UACR_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_kc <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "kidney_fibrosis_KC_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_om <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "kidney_fibrosis_OM_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_msa <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "mesentery_artery_wall_thickness_ratio_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_lv <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "lv_fibrosis_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_summary_cm <- analyze_cnmf_programs_bp(
  usage_list, metadata_merged, 
  bp_col = "cardiomyocyte_cross_sectional_area_mean", bp_transform  = T, bp_group_cols = "Project",
  agg_fun = "mean", min_cells_per_sample = 0
)

program_phenotype_summary <- 
  rbind(program_summary_bp$results_lm %>% filter(term == "BP_num") %>% mutate(term = "BP_mean"),
        program_summary_uacr$results_lm %>% filter(term == "BP_num") %>% mutate(term = "UACR_mean"),
        program_summary_kc$results_lm %>% filter(term == "BP_num") %>% mutate(term = "kidney_fibrosis_KC_mean"),
        program_summary_om$results_lm %>% filter(term == "BP_num") %>% mutate(term = "kidney_fibrosis_OM_mean"),
        program_summary_msa$results_lm %>% filter(term == "BP_num") %>% mutate(term = "mesentery_artery_wall_thickness_ratio_mean"),
        program_summary_lv$results_lm %>% filter(term == "BP_num") %>% mutate(term = "lv_fibrosis_mean"),
        program_summary_cm$results_lm %>% filter(term == "BP_num") %>% mutate(term = "cardiomyocyte_cross_sectional_area_mean"))


write.table(program_phenotype_summary, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/program_phenotype.summary.tsv", col.names = T, row.names = F, sep = "\t")





################################################################################
### Enrichment of program top gene to DEG
################################################################################
DEG_res <- read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/DEG.L1.all.out", header = T)
deg_by_celltype <- DEG_res %>%
  filter(
    strain %in% c("C57BL/6", "SS", "SHR"),
    !is.na(p_val_adj),
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.5,
  ) %>%
  distinct(cell_type, gene_name) %>% 
  mutate(cell_type = gsub("[/ ]", "_", cell_type)) %>% 
  group_by(cell_type) %>%
  dplyr::summarize(genes = list(sort(unique(gene_name))), .groups = "drop") %>%
  tibble::deframe()     
universe <- unique(DEG_res$gene_name)

top_gene_list <- lapply(spectra_list, function(x) {
  get_top_genes(x$spectra_score, k = 200)
})

prog_to_genes <- flatten_top_genes(top_gene_list)

fisher_res <- enrich_prog_vs_deg(prog_to_genes, deg_by_celltype, universe = universe, top_k = 200)
rownames(fisher_res) <- fisher_res$program







################################################################################
### Group programs by cosine similarity using top 200 genes in each program
################################################################################
res <- compare_all_spectra(spectra_list, method = "all")

spectra_name <- colnames(res$cosine)
spectra_ct <- gsub("-\\d+", "", spectra_name)
spectra_ct <- factor(spectra_ct, levels = names(new_cell_col))
spectra_tissue <- ct_to_tissue[as.character(spectra_ct)]

row_ha <- rowAnnotation(
  Cell_type_origin        = spectra_tissue,
  Cell_type_L1  = spectra_ct,
  col = list(
    Cell_type_origin       = new_tissue_col,
    Cell_type_L1 = new_cell_col
  ),
  annotation_legend_param = list(
    Cell_type_origin       = list(title = "Cell type origin"),
    Cell_type_L1 = list(title = "Cell type L1")
  )
)

m <- res$cosine
d  <- as.dist(1-m)
hc <- hclust(d, method = "ward.D2")

h_cut <- 1.4
grp  <- cutree(hc, h = h_cut)

prog_to_grp <- paste0("Group", grp)
names(prog_to_grp) <- names(grp)
program_group_list <- split(names(prog_to_grp), prog_to_grp)

fisher.test(spectra_tissue, prog_to_grp, simulate.p.value=TRUE)
# p-value = 5e-04





################################################################################
### Pathway enrichment analysis for group
################################################################################
# top 1000 genes per group
genes_per_grp_list <- rank_group_genes(program_group_list,
                                       spectra_list,
                                       top_n = 30,        # for visualization
                                       use_abs = FALSE,
                                       drop_na = TRUE)

group_gene_for_meta_1000 <- sapply(seq_along(genes_per_grp_list), function(i) {
  tb <- genes_per_grp_list[[i]]$table
  top_genes <- head(tb$gene, 1000)
  paste(top_genes, collapse = ",") %>% map_genes_to_ensembl("mouse")
})

names(group_gene_for_meta_1000) <- paste0(names(genes_per_grp_list), "_cosine_ward_1000")

group_gene_for_meta_df <- data.frame("#Name" = names(group_gene_for_meta_1000), 
                                     "Genes" = group_gene_for_meta_1000,
                                     check.names = FALSE)

write.table(group_gene_for_meta_df, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/CNMF.per_ct.cosine_ward.1000.for_metascape.ensembl.tsv", col.names = T, row.names = F, sep = "\t", quote = F)



# organize metascape results
input_list <- list.files(
  path = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/metascape/union/",
  pattern = "metascape_result.xlsx$", recursive = TRUE, full.names = TRUE)

group_input_list <- input_list[grepl("Group.*cosine_ward_1000", input_list)]
meta_merged_group_1k <- c()
for(infile in group_input_list){
  
  g <- strsplit(infile, "/")[[1]][12]
  meta <- readxl::read_xlsx(infile, sheet = 2)
  
  if(nrow(meta)>0){
    # meta_sum <- as.data.frame(meta[grepl("Summary", meta$GroupID), ])
    meta_sum <- as.data.frame(meta)
    meta_sum$group <- g
    
    meta_merged_group_1k <- rbind(meta_merged_group_1k, meta_sum)
  }
  
}
write.table(meta_merged_group_1k, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/metascape/group.1k.metascape.tsv", col.names = T, row.names = F, sep = "\t")






################################################################################
### Prioritize programs using BP-association and DEG-enrichment results (Figure 2a-d)
################################################################################

lm_res_bp <- program_summary_bp$results_lm
lm_res_bp <- lm_res_bp[lm_res_bp$term=="BP_num", ]
rownames(lm_res_bp) <- lm_res_bp$program

fisher_res_use <- fisher_res[rownames(lm_res_bp), ]

df <- data.frame(
  program = lm_res_bp$program,
  lm_bp_beta = lm_res_bp$estimate,
  lm_bp_p = lm_res_bp$p.value,
  lm_bp_p_adj = lm_res_bp$p_adj,
  fisher_overlap = fisher_res_use$A,
  fisher_or = fisher_res_use$OR,
  fisher_p = fisher_res_use$p,
  fisher_p_adj = fisher_res_use$p_adj
) %>%
  mutate(fisher_or = ifelse(is.infinite(fisher_or),
                            max(fisher_or[is.finite(fisher_or)]),
                            fisher_or),
         fisher_or_capped = pmin(fisher_or, 50))

highlight_prog_list <- c("T_cells-3", "VSMC-3", "TAL-6", "Myelinating_OL-9", "Fibroblast-13")
df %>%
  mutate(
    both_sig = (fisher_p_adj < 0.05) & (lm_bp_p_adj < 0.05),
    size_pt = ifelse(both_sig, 3, 1.5),
    label = program %in% highlight_prog_list
  ) %>%
  ggplot(aes(-log10(lm_bp_p_adj), -log10(fisher_p_adj), color = both_sig)) +
  geom_point(aes(size = size_pt)) +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "lightgrey")) +
  scale_size_identity() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "darkgrey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.5, color = "darkgrey") +
  geom_text_repel(
    data = ~ subset(.x, label),
    aes(label = program), 
    color = "black",
    size = 4,
    segment.color = "black",
    min.segment.length = 0,
    box.padding = 0.3,
    point.padding = 0.2,
    max.overlaps = Inf
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(legend.position = "none",
        axis.text = element_text(color = "black")) +
  labs(
    x = "Program-BP association\n(-log10 adjusted p-value)",
    y = "Program-DEG enrichment\n(-log10 adjusted p-value)"
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cnmf.program.p_val.dotplot.png", width=320/90, height=311/90, dpi = 300)


program_summary_bp$lm_input %>% 
  filter(program %in% c("VSMC-3")) %>%
  ggplot(aes(BP_num, score, colour = Strain)) +
  geom_point(size=3) +
  annotate("text", x=1.3, y=0.15, hjust = 0, label = "β = 0.088\np = 6.8e-4") +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(axis.text = element_text(color = "black")) +
  scale_color_manual(values = strain_col) +
  labs(x = "Normalized blood pressure", y = "Program usage score (VSMC-3)")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/VSMC-3.bp_usage.dotplot.png", width=360/90, height=311/90, dpi = 300)



program_vis(usage_list, metadata, program="VSMC-3") +
  scale_fill_gradient(limits = c(0, 0.4), low = "#D6604D", high = "#D6604D",
                      oob = scales::squish, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4), limits = c(0, 0.4)) +
  labs(title = "", y = "Program usage score\n(VSMC-3)") +
  theme(panel.spacing = unit(0.3, "lines"))
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/VSMC-3.usage.boxplot.png", width=380/90, height=311/90, dpi = 300)



# candidate genes based on pathway enrichment results
selected_genes <- c("Egfr", "Flt1", "Pparg", 
                    "Fn1", "Itga9", "Thbs1",
                    "Rock2", "Kalrn", "Plekhg2", 
                    "Rgs17", "Tgfb1")

W <- as.matrix(t(spectra_list$VSMC$spectra_score))
Wp <- apply(W, 2, function(x) rank(x, ties.method="average") / length(x))

vsmc_cols <- c(1:ncol(W))
k0 <- 3

w0 <- Wp[, k0]
w_other_max <- apply(Wp[, setdiff(vsmc_cols, k0), drop=FALSE], 1, max)

vsmc_3_gene <- tibble(
  gene      = rownames(W),
  score     = W[, 3], 
  delta_pct = w0 - w_other_max  
) %>%
  mutate(
    score = as.numeric(score),
    delta_pct = as.numeric(delta_pct),
    delta_pos = pmax(delta_pct, 0),
    top_by_score = rank(-score, ties.method = "min") <= 500,
    # highlight = top_by_score & (rank(-delta_pos, ties.method = "min") <= 200)
    highlight = gene %in% selected_genes
  ) %>%
  arrange(highlight)

ggplot(vsmc_3_gene, aes(x = delta_pct, y = score, color = highlight)) +
  geom_point() +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "lightgrey")) +
  geom_text_repel(
    data = subset(vsmc_3_gene, highlight),
    aes(label = gene),
    color = "black",
    size = 4,
    segment.color = "black",
    box.padding = 0.6,
    point.padding = 0.4,
    force = 2,
    force_pull = 0.2,
    max.iter = 20000,
    max.time = 5,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 1
  ) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5)
  ) +
  labs(x = "Gene specificity\n(vs. other VSMC programs)", y = "Gene score (VSMC-3)")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/VSMC-3.gene.specificity.png", width=300/90, height=311/90, dpi = 300)





################################################################################
### Visualize grouping results (Figure 2e)
################################################################################
S <- res$cosine
D <- as.dist(1-S)
umap_res <- umap::umap(as.matrix(D))
umap_coord <- umap_res$layout

umap_df <- data.frame(
  program = rownames(umap_coord),
  group = prog_to_grp[rownames(umap_coord)],
  UMAP1 = umap_coord[,1],
  UMAP2 = umap_coord[,2],
  tissue = spectra_tissue
)

umap_df_rename = umap_df %>%
  mutate(group_rename = gsub("Group", "G", group)) %>%
  mutate(group_rename = case_when(
    group_rename == "G3"  ~ "G3 immune/\ncytokine",
    group_rename == "G20" ~ "G20 stress/\nimmune",
    group_rename == "G24" ~ "G24 small GTPase\n& tube morphogenesis",
    group_rename == "G26" ~ "G26 xenobiotic/\namino acid metab",
    group_rename == "G28" ~ "G28 xenobiotic/\norganic anion",
    group_rename == "G29" ~ "G29 Rho /\nvascular migration",
    group_rename == "G30" ~ "G30 bile acid /\nperoxisome",
    group_rename == "G31" ~ "G31 cell cycle /\nmitosis",
    group_rename == "G33" ~ "G33 ECM / EMT",
    group_rename == "G34" ~ "G34 ECM /\ntube morphogenesis",
    group_rename == "G36" ~ "G36 muscle\ncontraction",
    group_rename == "G39" ~ "G39 OXPHOS /\nmyogenesis",
    TRUE ~ group_rename 
  ))

label_df <- umap_df_rename %>%
  group_by(group_rename) %>%
  dplyr::summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2), .groups = "drop")

ggplot(umap_df_rename, aes(UMAP1, UMAP2, color = group)) +
  geom_point(size = 3, alpha=0.5) +
  ggrepel::geom_text_repel(
    data = label_df,
    aes(label = group_rename),
    color = "black", size = 4, 
    bg.color = "white", bg.r = 0.1,
    box.padding = 0.25, point.padding = 0.25,
    max.overlaps = Inf, show.legend = FALSE,
    lineheight = 0.8, seed = 1
  ) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL)

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cnmf.group.umap.png", width=656/90, height=502/90, dpi = 300)


ggplot(umap_df_rename, aes(UMAP1, UMAP2, color = tissue)) +
  geom_point(size = 3) +
  scale_color_manual(values = new_tissue_col) +
  theme_classic() +
  theme(legend.position = "right") +
  labs(x = NULL, y = NULL, color = "Cell-type origin")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cnmf.group.umap.tissue.png", width=424/90, height=297/90, dpi = 300)







################################################################################
### Visualize program-phenotype network in group 34 and shared leading genes (Figure 2f-g)
################################################################################

selected_grp <- "Group34"

selected_prog <- program_group_list[[selected_grp]]
selected_prog_pheno <- program_phenotype_summary %>%
  filter(program %in% selected_prog & p_adj<0.05) %>%
  dplyr::mutate(term = case_when(
                  term == "BP_mean" ~ "Blood pressure",
                  term == "UACR_mean" ~ "Urine albumin-\ncreatinine ratio",
                  term == "cardiomyocyte_cross_sectional_area_mean" ~ "Cardiomyocyte\ncross-sectional area",
                  term == "kidney_fibrosis_KC_mean" ~ "Kidney fibrosis\n(cortex)",
                  term == "kidney_fibrosis_OM_mean" ~ "Kidney fibrosis\n(outer medulla)",
                  term == "lv_fibrosis_mean" ~ "Left ventricle fibrosis",
                  term == "mesentery_artery_wall_thickness_ratio_mean" ~ "Mesentery artery\nwall thickness ratio",
                  TRUE ~ term
                ))
  

selected_prog_idx <- which(umap_df_rename$group==selected_grp)
edges_idx <- S[selected_prog_idx, selected_prog_idx]
edges_long <- edges_idx %>%
  as.data.frame() %>%
  rownames_to_column(var = "from") %>%
  pivot_longer(
    cols = -from,
    names_to = "to",
    values_to = "Similarity"
  ) %>%
  filter(from != to) %>%
  mutate(
    x    = umap_coord[from, 1],
    y    = umap_coord[from, 2],
    xend = umap_coord[to, 1],
    yend = umap_coord[to, 2]
  ) %>%
  filter(Similarity > 0.1)


nodes <- tibble(
  name = unique(c(edges_long$from, edges_long$to, selected_prog_pheno$program, selected_prog_pheno$term)),
  type = ifelse(name %in% selected_prog_pheno$term, "Phenotype", "Program")
)

edges_program <- edges_long %>%
  transmute(from = from, to = to,
            weight = Similarity,
            edge_type = "program_program")

edges_pheno <- selected_prog_pheno %>%
  transmute(from = program, to = term,
            weight = -log10(p.value),
            edge_type = "program_phenotype")

edges <- bind_rows(edges_program, edges_pheno)

g <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE) %>%
  activate(edges) %>%
  mutate(
    w_prog = if_else(
      edge_type == "program_program",
      weight,  
      0.5 
    ),
    a_pheno = if_else(
      edge_type == "program_phenotype",
      weight, 
      0.5 
    ),
    edge_kind = if_else(edge_type == "program_program",
                        "Program–program", "Program–phenotype")
  )

set.seed(123456)
ggraph(g, layout = "fr") +
  geom_edge_link(
    aes(edge_width  = w_prog,
        edge_alpha  = a_pheno,
        edge_colour = edge_kind),
    show.legend = c(
      edge_width  = TRUE,   # legend for similarity
      edge_alpha  = TRUE,   # legend for -log10(p)
      edge_colour = FALSE   # no separate color legend
    )
  ) +
  scale_edge_width_continuous(
    name  = "Program-program\nsimilarity",
    range = c(0.2, 2)
  ) +
  scale_edge_alpha_continuous(
    name  = "Program-phenotype\n(-log10 p)",
    range = c(0.2, 1)
  ) +
  scale_edge_colour_manual(
    values = c(
      "Program-program"   = "grey60",
      "Program-phenotype" = "firebrick"
    ),
    guide = "none"
  ) +
  geom_node_point(aes(color = type), size = 4) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  scale_color_manual(values = c(
    Program   = "#1f77b4",
    Phenotype = "#d62728"
  ), name = "Node type") +
  theme_void(base_family = "Arial")

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/G34.prog.pheno.network.png", width=373/90, height=239/90, dpi = 300)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/G34.prog.pheno.network.legend.png", width=528/90, height=355/90, dpi = 300)



# figure 2g
selected_genes <- c("Slc34a2", "Hcn1", "Grik2", "Cp", "Nfkbiz", "Casp12", "Anxa1", "Il24", "Tpm1", "Cald1", "Tnik")
df <- genes_per_grp_list$Group34$table %>%
  mutate(gene = fct_reorder(gene, median_score),
         idx = 1:nrow(genes_per_grp_list$Group34$table),
         highlight = gene %in% selected_genes) %>%
  filter(idx<30) %>%
  arrange(highlight)
df %>%
  ggplot(aes(x = median_score, y = idx, color = highlight)) +
  geom_point() +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "lightgrey")) +
  geom_text_repel(
    data = subset(df, highlight),
    aes(label = gene),
    color = "black",
    size = 4,
    segment.color = "black",
    direction = "y",
    hjust = 0,
    # nudge_x = 1000 - subset(df, highlight)$idx,
    max.overlaps = Inf,
    min.segment.length = 0,
    seed = 1
  ) +
  labs(subtitle = "Top 30 ranked genes in G34",
       x = "Median rank across programs\n(scaled 0-1)", y = NULL) +
  theme_classic(base_size = 12, base_family = "Arial") +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5)
  )

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/G34.gene.rank.png", width=230/90, height=351/90, dpi = 300)





