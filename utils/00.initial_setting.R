library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(colorspace)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size) +
            theme(
              text = element_text(family = "Arial")
            ))

################################################################################
### variables
cell_order = c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("EC", "E/P transition cell", "Pericyte", "VSMC", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Immune cell"),
  c("T cells", "NK cells", "B cells", "Plasmablasts", "Plasma cells", 
    "Monocytes", "Neutrophils", "Mast cells", "Erythrocytes", "Megakaryocytes"),
  c("Neuronal")
)

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
cell_col = getPalette(length(cell_order))
names(cell_col) = cell_order

strain_order = c("C57BL/6", "SHR", "WKY", "SS", "SD")
strain_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))

class_order = c("neurons", "glial cells", "muscle cells", "epithelial cells", "endothelial cells", 
                "stromal cells", "immune cells", "adipocytes", "endocrine cells", "blood-related cells")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
class_col = getPalette(length(class_order))
names(class_col) = class_order

tissue_order = c("HYP", "MCA", "LV", "LK", "MSA", "PBMC")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
tissue_col = getPalette(length(tissue_order))
names(tissue_col) = tissue_order

model_order = c("AngII", "Salt-sensitive", "Spontaneous"); names(model_order)=c("mouse", "rat.ss", "rat.sp")

treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")

sxt_order = c("C57BL/6-Saline 3d", "C57BL/6-AngII 3d", "C57BL/6-AngII 28d", 
              "SS-LS", "SS-HS 3d", "SS-HS 21d", "SD-LS", "SD-HS 3d",
              "SHR-10w", "SHR-26w", "WKY-10w", "WKY-26w")

sxt_colors = c(
  lighten(strain_col["C57BL/6"], 0.35), strain_col["C57BL/6"], darken(strain_col["C57BL/6"], 0.25),
  lighten(strain_col["SS"], 0.35), strain_col["SS"], darken(strain_col["SS"], 0.25),
  lighten(strain_col["SD"], 0.35), strain_col["SD"],
  strain_col["SHR"], darken(strain_col["SHR"], 0.25),
  lighten(strain_col["WKY"], 0.35), strain_col["WKY"]
)

names(sxt_colors) = sxt_order





################################################################################
blank_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)






################################################################################
# function to score and visualize
calc_module_scores_and_summary <- function(
    seu,
    gene_sets,
    group_cols, 
    assay = NULL,
    slot = "data",
    nbin = 24,
    ctrl = 100,
    seed = 1,
    score_prefix = "MS",
    pct_cutoff = 0,          # score > pct_cutoff is "positive"
    group_sep = "_",
    na_group_drop = TRUE,
    return_cell_scores = FALSE,
    verbose = FALSE
) {
  stopifnot(inherits(seu, "Seurat"))
  stopifnot(is.list(gene_sets), length(gene_sets) > 0)
  
  if (is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    names(gene_sets) <- paste0("Set", seq_along(gene_sets))
  }
  
  if (is.null(assay)) assay <- Seurat::DefaultAssay(seu)
  if (!assay %in% names(seu@assays)) stop("Assay not found: ", assay)
  Seurat::DefaultAssay(seu) <- assay
  
  present_feats <- rownames(seu[[assay]])
  gene_sets_filt <- lapply(gene_sets, \(gs) unique(gs[gs %in% present_feats]))
  kept_sizes <- vapply(gene_sets_filt, length, integer(1))
  
  if (any(kept_sizes == 0)) {
    empty_sets <- names(gene_sets_filt)[kept_sizes == 0]
    warning("Dropped gene sets with 0 present genes: ", paste(empty_sets, collapse = ", "))
    gene_sets_filt <- gene_sets_filt[kept_sizes > 0]
  }
  if (length(gene_sets_filt) == 0) stop("No gene sets left after filtering to present features.")
  
  set.seed(seed)
  seu2 <- Seurat::AddModuleScore(
    object = seu,
    features = gene_sets_filt,
    assay = assay,
    slot = slot,
    nbin = nbin,
    ctrl = ctrl,
    name = paste0(score_prefix, "_"),
    verbose = verbose
  )
  
  k <- length(gene_sets_filt)
  score_cols  <- paste0(score_prefix, "_", seq_len(k))
  score_names <- names(gene_sets_filt)
  
  if (missing(group_cols) || is.null(group_cols) || length(group_cols) < 1) {
    stop("Please provide group_cols, e.g. group_cols = 'cell_type' or c('treatment','ct').")
  }
  if (!is.character(group_cols)) stop("group_cols must be character column name(s).")
  if (!all(group_cols %in% colnames(seu2@meta.data))) {
    stop("These group_cols are not in seu@meta.data: ",
         paste(setdiff(group_cols, colnames(seu2@meta.data)), collapse = ", "))
  }
  
  group_df <- seu2@meta.data[, group_cols, drop = FALSE]
  group_df$group_id <- apply(group_df, 1, paste, collapse = group_sep)
  
  # fetch cell scores
  df_scores <- Seurat::FetchData(seu2, vars = score_cols)
  df <- cbind(group_df, df_scores)
  
  if (na_group_drop) {
    keep <- stats::complete.cases(df[, c(group_cols, "group_id"), drop = FALSE])
    df <- df[keep, , drop = FALSE]
  }
  
  if (length(pct_cutoff) == 1) {
    cutoff_vec <- rep(pct_cutoff, k); names(cutoff_vec) <- score_names
  } else {
    if (is.null(names(pct_cutoff))) {
      stopifnot(length(pct_cutoff) == k)
      cutoff_vec <- pct_cutoff; names(cutoff_vec) <- score_names
    } else {
      cutoff_vec <- pct_cutoff[score_names]
      if (any(is.na(cutoff_vec))) {
        stop("pct_cutoff missing for: ", paste(score_names[is.na(cutoff_vec)], collapse = ", "))
      }
    }
  }
  
  out_list <- vector("list", k)
  
  for (i in seq_len(k)) {
    sc <- score_cols[i]
    nm <- score_names[i]
    co <- cutoff_vec[[nm]]
    
    # aggregate by group_cols + group_id (group_id is redundant but convenient for wide)
    tmp <- stats::aggregate(
      df[[sc]],
      by = c(df[, c(group_cols, "group_id"), drop = FALSE]),
      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                          mean_pos = mean(x[x > co], na.rm = TRUE),
                          median = median(x, na.rm = TRUE),
                          median_pos = median(x[x > co], na.rm = TRUE),
                          pct  = mean(x > co, na.rm = TRUE),
                          n    = sum(!is.na(x)))
    )
    
    tmp2 <- data.frame(
      tmp[, c(group_cols, "group_id"), drop = FALSE],
      score = nm,
      mean = tmp$x[, "mean"],
      mean_pos = tmp$x[, "mean_pos"],
      median = tmp$x[, "median"],
      median_pos = tmp$x[, "median_pos"],
      pct = tmp$x[, "pct"],
      n_cells = tmp$x[, "n"],
      cutoff = as.numeric(co),
      stringsAsFactors = FALSE
    )
    
    out_list[[i]] <- tmp2
  }
  
  summary_long <- do.call(rbind, out_list)
  
  group_key <- unique(summary_long[, c(group_cols, "group_id"), drop = FALSE])
  
  # AverageExpression-style wide matrices keyed by group_id
  mean_wide <- reshape(
    summary_long[, c("score", "group_id", "mean")],
    idvar = "score", timevar = "group_id", direction = "wide"
  )
  mean_pos_wide <- reshape(
    summary_long[, c("score", "group_id", "mean_pos")],
    idvar = "score", timevar = "group_id", direction = "wide"
  )
  median_wide <- reshape(
    summary_long[, c("score", "group_id", "median")],
    idvar = "score", timevar = "group_id", direction = "wide"
  )
  median_pos_wide <- reshape(
    summary_long[, c("score", "group_id", "median_pos")],
    idvar = "score", timevar = "group_id", direction = "wide"
  )
  pct_wide <- reshape(
    summary_long[, c("score", "group_id", "pct")],
    idvar = "score", timevar = "group_id", direction = "wide"
  )
  colnames(mean_wide) <- sub("^mean\\.", "", colnames(mean_wide))
  colnames(mean_pos_wide) <- sub("^mean_pos\\.", "", colnames(mean_pos_wide))
  colnames(median_wide) <- sub("^median\\.", "", colnames(median_wide))
  colnames(median_pos_wide) <- sub("^median_pos\\.", "", colnames(median_pos_wide))
  colnames(pct_wide)  <- sub("^pct\\.",  "", colnames(pct_wide))
  rownames(mean_wide) <- mean_wide$score
  rownames(mean_pos_wide) <- mean_pos_wide$score
  rownames(median_wide)  <- median_wide$score
  rownames(median_pos_wide)  <- median_pos_wide$score
  rownames(pct_wide)  <- pct_wide$score
  mean_wide$score <- NULL
  mean_pos_wide$score  <- NULL
  median_wide$score <- NULL
  median_pos_wide$score <- NULL
  pct_wide$score  <- NULL
  
  res <- list(
    seurat = seu2,
    summary_long = summary_long,  # includes separate treatment/ct columns + group_id
    mean_wide = mean_wide,        # columns are group_id
    pct_wide = pct_wide,          # columns are group_id
    mean_pos_wide = mean_pos_wide,        # columns are group_id
    median_wide = median_wide,          # columns are group_id
    median_pos_wide = median_pos_wide,          # columns are group_id
    group_key = group_key,        # map group_id -> treatment/ct/etc
    gene_set_sizes_present = kept_sizes[names(gene_sets_filt)]
  )
  
  if (return_cell_scores) {
    cell_scores <- df
    names(cell_scores)[match(score_cols, names(cell_scores))] <- score_names
    res$cell_scores <- cell_scores
  }
  
  res
}


score_vis <- function(score_df, score_selected, project_selected, cell_type_selected){
  dat_use <- score_df %>%
    dplyr::filter(Cell_type_L1 %in% cell_type_selected,
                  score %in% score_selected,
                  Project %in% project_selected) %>% 
    group_by(score) %>%
    ungroup()
  
  p <- ggplot(dat_use, aes(x = SxT, y = score)) +
    geom_point(aes(fill = mean, size=pct*100), shape = 21) +
    scale_y_discrete(limits=rev) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme(
      text = element_text(family = "Arial"), 
      panel.grid.major.y = element_blank(), 
      legend.justification = c(1, 0.5),
      axis.text.y = element_text(colour = 'black'),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
      strip.text = element_text(colour = 'black'),
      legend.position="right",  legend.direction="horizontal"
    ) +
    labs(x="", y="", fill="Mean", size="Percentage") +
    ggh4x::facet_nested( ~ Cell_type_L1, scales = "free", space = "free") +
    guides(
      fill  = guide_colourbar(order = 1),
      size  = guide_legend(order = 2)
    )
  
  return(p)
}




