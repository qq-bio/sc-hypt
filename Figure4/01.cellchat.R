
library(CellChat)
library(Seurat)
library(biomaRt)



################################################################################
### Asssembl rat cellchatdb
################################################################################
# Map gene strings like "A+BC" using map_vec (keep unknowns unchanged)
.map_gene_string <- function(x, map_vec) {
  if (is.na(x) || !nzchar(x)) return(x)
  parts <- strsplit(x, "\\+")[[1]]
  parts2 <- vapply(parts, function(g) {
    g2 <- unname(map_vec[g])
    if (is.na(g2) || !nzchar(g2)) g else g2
  }, character(1))
  paste(parts2, collapse = "+")
}

# Map a single gene token if present in map_vec, otherwise keep it
.map_token <- function(tok, map_vec) {
  tok <- trimws(tok)
  if (!nzchar(tok)) return(tok)
  t2 <- unname(map_vec[tok])
  if (is.na(t2) || !nzchar(t2)) tok else t2
}

# Deduplicate while preserving order
.unique_preserve <- function(x) x[!duplicated(x)]

map_interaction_name_2 <- function(x, map_vec, dedup = TRUE) {
  if (is.na(x) || !nzchar(x)) return(x)

  # split on " - " (allow variable spaces)
  parts <- strsplit(x, "\\s*-\\s*", perl = TRUE)[[1]]
  if (length(parts) < 2) return(x)  # unexpected format, keep as-is

  left_raw  <- parts[1]
  right_raw <- paste(parts[-1], collapse = " - ")  # in case " - " appears in names (rare)

  left_mapped <- .map_token(left_raw, map_vec)

  right_raw_trim <- trimws(right_raw)

  # case 1: parentheses form "(B+C+...)"
  if (grepl("^\\(.*\\)$", right_raw_trim)) {
    inner <- sub("^\\((.*)\\)$", "\\1", right_raw_trim)
    toks <- strsplit(inner, "\\+")[[1]]
    toks_m <- vapply(toks, .map_token, character(1), map_vec = map_vec)
    toks_m <- trimws(toks_m)
    toks_m <- toks_m[nzchar(toks_m)]
    if (dedup) toks_m <- .unique_preserve(toks_m)

    # if only one token left, CellChat sometimes uses "(X)" still; keep parentheses to preserve style
    right_mapped <- paste0("(", paste(toks_m, collapse = "+"), ")")
  } else {
    # case 2: simple "B"
    right_mapped <- .map_token(right_raw_trim, map_vec)
  }

  paste0(left_mapped, " - ", right_mapped)
}

map_complex_df <- function(complex, map_vec,
                           subunit_cols = grep("^subunit_", colnames(complex), value = TRUE)) {
  out <- complex
  for (cc in intersect(subunit_cols, colnames(out))) {
    out[[cc]] <- vapply(out[[cc]], function(g) {
      if (is.na(g) || !nzchar(g)) return(g)
      mg <- unname(map_vec[g])
      if (is.na(mg)) g else mg
    }, character(1))
  }
  out
}

map_cofactor_df <- function(cofactor, map_vec,
                            cof_cols = grep("^cofactor", colnames(cofactor), value = TRUE)) {
  out <- cofactor
  for (cc in intersect(cof_cols, colnames(out))) {
    out[[cc]] <- vapply(out[[cc]], function(g) {
      if (is.na(g) || !nzchar(g)) return(g)
      mg <- unname(map_vec[g])
      if (is.na(mg)) g else mg
    }, character(1))
  }
  out
}

map_geneInfo_df <- function(geneInfo, map_vec, symbol_col = "Symbol") {
  out <- geneInfo
  if (symbol_col %in% colnames(out)) {
    out[[symbol_col]] <- vapply(out[[symbol_col]], function(g) {
      if (is.na(g) || !nzchar(g)) return(g)
      mg <- unname(map_vec[g])
      if (is.na(mg)) g else mg
    }, character(1))
  }
  out
}

m2r <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2rat.out.txt", header = T, sep = "\t")
m2r <- m2r[m2r$Rat.orthology.confidence..0.low..1.high.==1, ]

map_df <- unique(m2r[, c("Gene.name", "Rat.gene.name")])
colnames(map_df) <- c("mouse_symbol", "rat_symbol")
map_df <- map_df[!is.na(map_df$mouse_symbol) & !is.na(map_df$rat_symbol) &
                   map_df$mouse_symbol != "" & map_df$rat_symbol != "", , drop = FALSE]

map_vec <- setNames(map_df$rat_symbol, map_df$mouse_symbol)


interaction = CellChatDB.mouse$interaction
complex = CellChatDB.mouse$complex
cofactor = CellChatDB.mouse$cofactor
geneInfo = CellChatDB.mouse$geneInfo

interaction_rat <- interaction
interaction_rat$interaction_name_2 <- vapply(interaction_rat$interaction_name_2, map_interaction_name_2,
                                 character(1), map_vec = map_vec, dedup = TRUE)
interaction_rat$ligand   <- vapply(interaction_rat$ligand,   .map_gene_string, character(1), map_vec = map_vec)
interaction_rat$receptor <- vapply(interaction_rat$receptor, .map_gene_string, character(1), map_vec = map_vec)
interaction_rat$interaction_name <- paste0(interaction_rat$ligand, "_", interaction_rat$receptor)
rownames(interaction_rat) <- make.unique(interaction_rat$interaction_name)

complex_rat     <- map_complex_df(complex, map_vec)
cofactor_rat    <- map_cofactor_df(cofactor, map_vec)
geneInfo_rat    <- map_geneInfo_df(geneInfo, map_vec)

CellChatDB.rat = list()
CellChatDB.rat$interaction = interaction_rat
CellChatDB.rat$complex = complex_rat
CellChatDB.rat$cofactor = cofactor_rat
CellChatDB.rat$geneInfo = geneInfo_rat

save(CellChatDB.rat, file = "/xdisk/mliang1/qqiu/data/cellchatdb/CellChatDB.rat.rda")








################################################################################
### CellChat Analysis 
################################################################################
collapse_rows_by_group_sparse <- function(mat, group) {
  # mat: genes x cells (dgCMatrix ok)
  stopifnot(length(group) == nrow(mat))
  
  f <- factor(group, levels = unique(group))  # preserve order of first appearance
  S <- Matrix::sparseMatrix(
    i = as.integer(f),
    j = seq_along(f),
    x = 1,
    dims = c(nlevels(f), length(f))
  )
  out <- S %*% mat
  rownames(out) <- levels(f)
  out
}

rename_features_by_rebuild_assay <- function(obj, assay = "RNA", map_vec,
                                             remove_ensembl_version = TRUE,
                                             collapse_duplicates = TRUE,
                                             normalize_after = TRUE) {
  
  obj <- Seurat::UpdateSeuratObject(obj)
  
  # Get counts from v4/v5
  counts <- tryCatch(
    GetAssayData(obj, assay = assay, slot = "counts"),
    error = function(e) GetAssayData(obj, assay = assay, layer = "counts")
  )
  
  old <- rownames(counts)
  key <- old
  if (remove_ensembl_version) key <- sub("\\.[0-9]+$", "", key)
  
  new <- unname(map_vec[key])
  new[is.na(new) | new == ""] <- old[is.na(new) | new == ""]
  
  # Collapse duplicates (recommended when Ensembl->symbol creates collisions)
  if (collapse_duplicates) {
    counts2 <- collapse_rows_by_group_sparse(counts, group = new)
    # ensure unique rownames (rare if your group had exact dup labels after factor order)
    rownames(counts2) <- make.unique(rownames(counts2))
  } else {
    counts2 <- counts
    rownames(counts2) <- make.unique(new)
  }
  
  obj[[assay]] <- CreateAssayObject(counts = counts2)
  
  if (normalize_after) {
    obj <- NormalizeData(obj, assay = assay, verbose = FALSE)
  }
  
  obj
}

load("/xdisk/mliang1/qqiu/data/cellchatdb/CellChatDB.rat.rda")

r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")

map_df <- unique(r2h[, c("Gene.name", "Gene.stable.ID")])
colnames(map_df) <- c("rat_symbol", "rat_id")
map_df <- map_df[!is.na(map_df$rat_id) & !is.na(map_df$rat_symbol) &
                   map_df$rat_id != "" & map_df$rat_symbol != "", , drop = FALSE]

rat_map_vec <- setNames(map_df$rat_symbol, map_df$rat_id)

input_file = c(
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
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.rds"
)

e <- new.env()
for(i in input_file){
  
  dataset = gsub("\\.[RNA|multiomics|EC]+.anno.L2.rds", "", basename(i), perl = T)
  
  CellChatDB = CellChatDB.mouse
  if(grepl("mouse", dataset)){
    CellChatDB = CellChatDB.mouse
  }else{
    CellChatDB = CellChatDB.rat
  }
  
  tissue = unlist(lapply(strsplit(dataset, "\\."), function(x) x[length(x)]))
  
  seurat_object = readRDS(i)
  
  if(!grepl("mouse", dataset)){
    seurat_object <- rename_features_by_rebuild_assay(
      seurat_object,
      assay = "RNA",
      map_vec = rat_map_vec,
      remove_ensembl_version = TRUE,
      collapse_duplicates = TRUE
    )
  }
  
  Idents(seurat_object) <- "Cell_type_L1"
  seurat_object@active.ident = factor(seurat_object@meta.data[, "Cell_type_L1"])
  Idents(seurat_object) = "Cell_type_L1"
  
  strain_list = unique(seurat_object$Strain)
  for(si in strain_list){
    
    treatment_list = unique(seurat_object@meta.data[seurat_object$Strain==si, ]$Treatment)
    
    for(ti in treatment_list){
      
      seurat_object_tmp = subset(seurat_object, Strain==si & Treatment==ti)
      cellchat = createCellChat(seurat_object_tmp, group.by = "Cell_type_L1")
      cellchat@DB = CellChatDB
      cellchat = subsetData(cellchat)
      cellchat <- identifyOverExpressedGenes(cellchat)
      cellchat <- identifyOverExpressedInteractions(cellchat)
      cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
      cellchat <- filterCommunication(cellchat, min.cells = 10)
      cellchat <- computeCommunProbPathway(cellchat)
      cellchat <- aggregateNet(cellchat)
      
      cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
      
      output_name = paste0(tissue, "_", gsub("/", ".", si), "_", ti)
      
      with(e, {
        assign(output_name, cellchat)
      })
      
    }
    
  }
  
}

saveRDS(e, "cellchat.rds")



model = c("AngII", "Salt-sensitive", "Salt-sensitive", "Spontaneous", "Spontaneous"); names(model)=c("C57BL/6", "SS", "SD", "SHR", "WKY")

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
  
  lr_df$model = model[strain]
  lr_df$tissue = gsub("(.*)\\.([A-Z]+).*", "\\2", tissue, perl = T)
  lr_df$strain = strain
  lr_df$treatment = treatment
  lr_df = cbind(lr_df, cellchat@LR$LRsig[lr_df$LR_pair, ])
  
  cc_df = rbind(cc_df, lr_df)
}

write.table(cc_df, "cellchat.result.out", col.names = T, row.names = F, sep = '\t', quote = F)


################################################################################
### Visualize CellChat results in hypothalamus (Figure 4c-d) 
################################################################################
cc_df = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellchat/cellchat.result.out", header = T, sep = '\t')

cc_df$model = factor(cc_df$model, levels=c("AngII", "Salt-sensitive", "Spontaneous"))
cc_df$strain = factor(cc_df$strain, levels = c("C57BL/6", "SS", "SD", "SHR", "WKY"))
cc_df$tissue = factor(cc_df$tissue, levels = c("HYP", "MCA", "LV", "LK", "MSA", "PBMC"))
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


### concordance barplot (Figure 4a)
cc_df_use = cc_df_diff %>% 
  group_by(sxt, model, strain, tissue, treatment, pathway_name) %>%
  dplyr::summarise(n_lr = dplyr::n_distinct(interaction_name),
                   prob_sum = sum(value),
                   prob_sum_r = prob_sum/n_lr,
                   prob_abs_sum = sum(abs(value)),
                   prob = abs(prob_sum/prob_abs_sum),
                   prob_score = prob_sum*prob,
                   .groups = "drop") %>%
  as.data.frame()


strain_tbl <- cc_df_use %>% 
  filter(strain %in% c('C57BL/6', 'SS', 'SHR')) %>%
  group_by(pathway_name, tissue, strain) %>%
  dplyr::summarise(
    value_strain = median(prob_sum, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(tissue, pathway_name) %>%
  dplyr::summarise(
    n_pos = sum(.data[["value_strain"]] > 0, na.rm = TRUE),
    n_neg = sum(.data[["value_strain"]] < 0, na.rm = TRUE),
    n_eff = n_pos + n_neg,
    concordance_dir = ifelse(n_eff == 0, NA_real_, pmax(n_pos, n_neg)/n_eff),
    direction = case_when(n_pos > n_neg ~ "up", n_neg > n_pos ~ "down", TRUE ~ "tie"),
    strain_mag = median(abs(.data[["value_strain"]]), na.rm = TRUE),
    strain_score_med = median(.data[["value_strain"]], na.rm = TRUE),
    
    .groups = "drop"
  )

strain_tbl %>% 
  filter(!is.na(concordance_dir), n_eff>1) %>%
  mutate(label = paste0(concordance_dir * n_eff, "/", n_eff)) %>%
  dplyr::count(tissue, concordance_dir, label) %>%
  ggplot(aes(x = tissue, y = n, fill = as.numeric(concordance_dir))) +
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
    y = "Number of cellular communication pathways"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    text = element_text(family = "Arial"),
    panel.grid.major.x = element_blank()
  )
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/cellchat.htn.concordance.gt_1.barplot.png", width=285/96, height=350/96, dpi = 300)



# ranked plot (Figure 4d)
cc_df_ranked <- cc_df_use %>% 
  mutate(strain_type = ifelse(strain %in% c("C57BL/6", "SS", "SHR"), "Hypertensive", "Normotenive")) %>%
  group_by(model, strain, tissue) %>%
  dplyr::mutate(rank = rank(prob_sum, ties.method = "first"),
                rank_norm = (rank - 1) / (n() - 1)) %>%
  arrange(rev(strain_type), rank)

cc_df_ranked_hyp <- cc_df_ranked %>%
  filter(tissue == "HYP", strain_type == "Hypertensive")

label_df <- cc_df_ranked_hyp %>%
  filter(
    tissue == "HYP",
    pathway_name %in% c("NCAM")   # add more pathways here if needed
  )

ggplot(cc_df_ranked_hyp, aes(x = prob_sum, y = rank_norm)) +
  geom_line(linewidth = 1, color = "grey") +
  geom_point(alpha = 0.5, aes(color = strain)) +
  geom_text_repel(
    data = label_df,
    aes(label = pathway_name),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    min.segment.length = 0,
    max.overlaps = Inf        # don't drop labels
  ) +
  facet_nested( ~ strain, scales = "free_x") +
  scale_color_manual(values = strain_col) +
  theme_classic() +
  theme(legend.position = "None",
        text = element_text(family = "Arial"),
        axis.text = element_text(color = "black")) +
  labs(
    x = "Sum of communication probability change",
    y = "Pathway rank (0–1)",
    color = "Strain type"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/hyp.cellchat.rank.png", width=374/96, height=214/96, dpi=300)














