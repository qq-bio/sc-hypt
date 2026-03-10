suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
})

base_font_size = 12
theme_set(theme_classic(base_size = base_font_size,
                        base_family = "Arial"))



################################################################################
### initial setting
cell_tissue_order <- list(
  "Parenchymal (HYP)" = c("Inhibitory_neuron", "Excitatory_neuron", "Avp+_neuron", 
                          "Astrocyte", "OPC", "NFO", "Premyelinating_OL", "Myelinating_OL", 
                          "Tanycyte", "Ependymal_cell", "Pars_tuberalis_cell"),
  "Parenchymal (LV)" = c("CM", "Neuronal"),
  "Parenchymal (LK)" = c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  "Vascular" = c("EC", "Pericyte", "E_P_transition_cell", "VSMC"),
  "Stromal" = c("Fibroblast", "Adipocyte"),
  "Immune" = c("T_cells", "NK_cells", "NKT", "B_cells", "Plasmablasts", "Plasma_cells", "Macrophages",
               "Monocytes", "Neutrophils", "DC", "Mast_cells", "Microglia", "Activated_microglia"),
  "Blood" = c("Erythrocytes", "Megakaryocytes")
)

ct_to_tissue <- unlist(lapply(names(cell_tissue_order), function(tis) {
  setNames(rep(tis, length(cell_tissue_order[[tis]])), cell_tissue_order[[tis]])
}))

tissue_order = names(cell_tissue_order)
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
new_tissue_col = getPalette(length(tissue_order))
names(new_tissue_col) = tissue_order


cell_order = c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("EC", "E/P transition cell", "Pericyte", "VSMC", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Immune cell"),
  c("T cells", "NK cells", "NKT", "B cells", "Plasmablasts", "Plasma cells", "Macrophages", 
    "Monocytes", "Neutrophils", "DC", "Mast cells", "Erythrocytes", "Megakaryocytes"),
  c("Neuronal")
)

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
new_cell_col = getPalette(length(cell_order))
names(new_cell_col) = gsub(" |/", "_", cell_order)


strain_order = c("C57BL/6", "SHR", "WKY", "SS", "SD")
strain_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))








################################################################################
### organize genes for metascape analysis
map_genes_to_ensembl <- function(genes_str,
                                 species = c("mouse", "rat"),
                                 return_mapping = FALSE) {
  species <- match.arg(species)
  
  if (species == "mouse") {
    OrgDb <- org.Mm.eg.db::org.Mm.eg.db
    ens_pat <- "^ENSMUSG[0-9]+$"
  } else {
    OrgDb <- org.Rn.eg.db::org.Rn.eg.db
    ens_pat <- "^ENSRNOG[0-9]+$"
  }
  
  # parse input
  genes <- unique(trimws(unlist(strsplit(genes_str, ","))))
  genes <- genes[nzchar(genes)]
  if (!length(genes)) return("")
  
  # separate Ensembl vs SYMBOL
  is_ens <- grepl(ens_pat, toupper(genes))
  ens_ids <- genes[is_ens]
  symbols <- genes[!is_ens]
  
  # map SYMBOLs → Ensembl
  mapped <- if (length(symbols)) {
    AnnotationDbi::select(OrgDb,
                          keys = symbols,
                          keytype = "SYMBOL",
                          columns = "ENSEMBL")
  } else data.frame(SYMBOL = character(), ENSEMBL = character())
  
  # combine with original Ensembl IDs
  all_ensembl <- unique(c(ens_ids, mapped$ENSEMBL))
  all_ensembl <- all_ensembl[!is.na(all_ensembl)]
  
  result <- paste(all_ensembl, collapse = ",")
  
  if (return_mapping) {
    return(list(ensembl = result, mapping = mapped))
  } else {
    return(result)
  }
}

get_top_genes <- function(score_mat, k = 100, use_abs = FALSE) {
  apply(score_mat, 1, function(x) {
    ord <- order(if (use_abs) abs(x) else x, decreasing = TRUE)
    names(x)[ord[1:min(k, length(x))]]
  })
}

flatten_top_genes <- function(nested_top_genes) {
  out <- list()
  for (ct in names(nested_top_genes)) {
    sublist <- nested_top_genes[[ct]]
    for (idx in 1:ncol(sublist)) {
      out[[paste0(ct, "-", idx)]] <- sublist[, idx]
    }
  }
  out
}

get_group_union_genes <- function(group_assign, prog_to_genes) {
  groups <- unique(na.omit(group_assign))
  lapply(groups, function(g) {
    progs <- names(group_assign)[group_assign == g]
    genes <- unlist(prog_to_genes[intersect(progs, names(prog_to_genes))])
    # genes_keep <- names(table(genes))[table(genes) >= 2]
    
    genes
  }) |> setNames(groups)
}










################################################################################
### boxplot to visualize cellular program score across animal conditions
program_vis <- function(usage_list, metadata, program, group.by = "Cell_type_L1", tissue = NULL){
  
  ct <- gsub("-.*", "", program)
  df <- usage_list[[ct]] %>%
    tibble::rownames_to_column("Cell_ID") %>%
    dplyr::select(Cell_ID, program = dplyr::all_of(program)) %>%
    dplyr::left_join(metadata, by = "Cell_ID") %>%
    dplyr::filter(!is.na(.data[[group.by]]))
  
  grp <- rlang::sym(group.by)
  df <- df %>% dplyr::mutate(.group = .data[[group.by]])
  
  med_tbl <- df %>%
    dplyr::group_by(Tissue, Strain, Treatment, .group) %>%
    dplyr::summarise(.med = median(.data[["program"]], na.rm = TRUE),
                     .groups = "drop")
  
  df <- df %>%
    dplyr::left_join(med_tbl, by = c("Tissue","Strain","Treatment",".group"))
  
  if(grp=="Cell_type_L1"){
    df$.group <- factor(df$.group, levels = cell_order)
  }else{
    df$.group <- factor(df$.group)
  }
  
  if(!is.null(tissue)){
    df <- df[df$Tissue == tissue, ]
  }
  
  p <- ggplot(df, aes(x = Treatment, y = program, fill = .med)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    guides(color = "none") +
    labs(x = NULL, y = "Program score", title = program, fill = "Median\nscore") +
    # coord_flip() +
    theme_bw(base_size = 22, base_family = "Arial") +
    theme(
      plot.title = ggtext::element_markdown(lineheight = 1.1),
      legend.key.height = unit(0.6, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9),
      text              = element_text(size = 11, color = "black"),
      axis.text.y       = element_text(size = 11, colour = "black"),
      axis.text.x       = element_text(size = 9,  colour = "black", angle = 45, hjust = 1),
      axis.title        = element_text(size = 11, colour = "black"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(colour = "black", fill = NA)
    ) +
    ggh4x::facet_nested(Tissue ~ Project + Strain , scales = "free_x", space = "free_x") +
    scale_fill_gradient(limits = c(0, 1), low = "white", high = "red",
                        oob = scales::squish, name = "Median") +
    scale_color_gradient(limits = c(0, 1), low = "white", high = "red",
                         oob = scales::squish) +
    scale_y_continuous(breaks = c(0, 0.5, 1))
  
  return(p)
  
}







################################################################################
### program-to-phenotype association using linear model
analyze_cnmf_programs_bp <- function(
    usage_list,
    metadata_merged,                 # merged cell-level metadata with BP already joined
    bp_col        = "BP_mean",       # name of numeric BP column in metadata_merged
    bp_transform  = FALSE,           # TRUE -> scale to min=1 within groups
    bp_group_cols = NULL,            # e.g., "Strain", "Project", or c("Project","Strain")
    ct_col        = "Cell_type_L1",
    sample_col    = "Sample_ID",
    tissue_col    = "Tissue",
    strain_col    = "Strain",
    agg_fun       = c("mean","median"),
    min_cells_per_sample = 20,
    verbose       = TRUE
) {
  require(dplyr); require(tidyr); require(purrr); require(broom)
  
  agg_fun <- match.arg(agg_fun)
  stopifnot(is.list(usage_list))
  needed_meta <- c(ct_col, sample_col, tissue_col, strain_col, bp_col)
  mis <- setdiff(needed_meta, colnames(metadata_merged))
  if (length(mis)) stop("Missing columns in metadata_merged: ", paste(mis, collapse = ", "))
  if (is.null(rownames(metadata_merged))) stop("metadata_merged must have rownames = cell barcodes")
  
  .agg <- switch(agg_fun,
                 mean   = function(x) mean(x, na.rm = TRUE),
                 median = function(x) stats::median(x, na.rm = TRUE))
  
  # Build a distinct per-sample BP table from the merged metadata
  bp_tbl <- metadata_merged %>%
    dplyr::select(
      Sample = .data[[sample_col]],
      BP_value = .data[[bp_col]],
      dplyr::all_of(bp_group_cols)
    ) %>%
    dplyr::distinct() %>%
    dplyr::mutate(Sample = as.character(Sample))
  
  # Optional: relative BP scaling within bp_group_cols (min -> 1)
  if (isTRUE(bp_transform) && !is.null(bp_group_cols) && length(bp_group_cols) > 0) {
    bp_tbl <- bp_tbl %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(bp_group_cols))) %>%
      dplyr::mutate(BP_num = BP_value / min(BP_value, na.rm = TRUE)) %>%
      dplyr::ungroup()
  } else {
    bp_tbl <- bp_tbl %>% dplyr::mutate(BP_num = BP_value)
  }
  
  if (verbose) {
    message("BP (per-sample) preview:")
    print(utils::head(bp_tbl, 5))
  }
  
  results_lm <- list()
  all_qc     <- list()
  all_input  <- list()
  
  ct_names <- names(usage_list)
  for (ct in ct_names) {
    if (verbose) message("Processing CT: ", ct)
    
    U <- usage_list[[ct]]
    if (is.null(U) || nrow(as.matrix(U)) == 0 || ncol(as.matrix(U)) == 0) next
    
    # Long-format usage and join cell-level metadata_merged
    df <- as.data.frame(U, check.names = FALSE)
    df$Cell_ID <- rownames(U)
    
    df <- df %>%
      tidyr::pivot_longer(-Cell_ID, names_to = "program", values_to = "score") %>%
      dplyr::left_join(metadata_merged, by = "Cell_ID") %>%
      dplyr::filter(!is.na(.data[[ct_col]])) %>%
      dplyr::mutate(
        CT     = ct,
        Tissue = as.character(.data[[tissue_col]]),
        Strain = as.character(.data[[strain_col]]),
        Sample = as.character(.data[[sample_col]])
      )
    
    # Aggregate to sample level, then join per-sample BP_num
    agg <- df %>%
      dplyr::group_by(CT, program, Sample, Tissue, Strain) %>%
      dplyr::summarise(
        score   = .agg(score),
        n_cells = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(n_cells >= min_cells_per_sample) %>%
      dplyr::left_join(bp_tbl %>% dplyr::select(Sample, BP_value, BP_num), by = "Sample") %>%
      dplyr::filter(!is.na(BP_num)) %>%
      dplyr::ungroup()
    
    # agg$score_int <- qnorm(rank(agg$score) / (nrow(agg) + 1))
    # agg$score <- agg$score_int
    
    agg$CT <- ct
    all_input[[ct]] <- agg  
    
    if (verbose) message("  rows after aggregate + BP join: ", nrow(agg))
    
    qc_ct <- agg %>%
      dplyr::group_by(CT, program) %>%
      dplyr::summarise(
        n_samples = dplyr::n(),
        median_cells_per_sample = stats::median(n_cells),
        .groups = "drop"
      )
    all_qc[[ct]] <- qc_ct
    
    if (nrow(agg) == 0) next
    
    # Fit LM: score ~ BP_num + (Tissue if multi-level) + (Strain if multi-level)
    fit_one <- function(dat, keys) {
      dat <- dat[!is.na(dat$score) & !is.na(dat$BP_num), , drop = FALSE]
      
      n_bpnum  <- dplyr::n_distinct(dat$BP_num)
      n_tissue <- dplyr::n_distinct(dat$Tissue)
      n_strain <- dplyr::n_distinct(dat$Strain)
      
      tissue_value <- if (n_tissue == 1) as.character(dat$Tissue[1]) else "mixed"
      strain_value <- if (n_strain == 1) as.character(dat$Strain[1]) else "mixed"
      
      if (n_bpnum < 2) {
        return(tibble::tibble(
          term        = "insufficient_BP_numeric_levels",
          estimate    = NA_real_, std.error = NA_real_,
          statistic   = NA_real_, p.value   = NA_real_,
          N           = nrow(dat), adj_r2   = NA_real_,
          n_bp_levels = n_bpnum,
          n_tissue_lv = n_tissue, n_strain_lv = n_strain,
          tissue_value = tissue_value, strain_value = strain_value
        ))
      }
      
      if (n_tissue >= 2) dat$Tissue <- factor(dat$Tissue)
      if (n_strain >= 2) dat$Strain <- factor(dat$Strain)
      
      rhs_terms <- c("BP_num",
                     if (n_tissue >= 2) "Tissue",
                     if (n_strain >= 2) "Strain")
      f <- stats::reformulate(rhs_terms, response = "score")
      m <- stats::lm(f, data = dat)
      
      td <- broom::tidy(m)
      gl <- broom::glance(m)
      
      td$N            <- nrow(dat)
      td$adj_r2       <- gl$adj.r.squared
      td$n_bp_levels  <- n_bpnum
      td$n_tissue_lv  <- n_tissue
      td$n_strain_lv  <- n_strain
      td$tissue_value <- tissue_value
      td$strain_value <- strain_value
      
      td
    }
    
    res_ct <- agg %>%
      dplyr::group_by(CT, program) %>%
      dplyr::group_modify(~ fit_one(.x, .y)) %>%
      dplyr::ungroup()
    
    res_ct <- res_ct %>%
      group_by(CT, term) %>% 
      mutate(p_adj = p.adjust(p.value, "BH")) %>%
      ungroup()
    
    results_lm[[ct]] <- res_ct
  }
  
  results_lm <- dplyr::bind_rows(results_lm)
  res_qc     <- dplyr::bind_rows(all_qc)
  input_data <- dplyr::bind_rows(all_input)
  
  list(
    results_lm = results_lm,
    qc         = res_qc,
    lm_input   = input_data
  )
}









################################################################################
### enrichment of program top genes to DEG
fisher_or_sets <- function(topk, deg, universe) {
  topk <- intersect(unique(topk), universe)
  deg  <- intersect(unique(deg),  universe)
  
  A <- length(intersect(topk, deg))                           # in topk & DEG
  B <- length(setdiff(topk, deg))                             # in topk & not DEG
  C <- length(setdiff(deg, topk))                             # in DEG & not topk
  D <- length(setdiff(universe, union(topk, deg)))            # neither
  
  mat <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
  ft  <- tryCatch(fisher.test(mat, alternative = "greater"),
                  error = function(e) list(estimate = NA_real_, p.value = 1))
  
  or_adj <- (A + 0.5) * (D + 0.5) / ((B + 0.5) * (C + 0.5))
  
  list(A = A, B = B, C = C, D = D,
       OR = ifelse(is.null(ft$estimate), NA_real_, unname(ft$estimate)),
       OR_adj = or_adj,
       p = ifelse(is.null(ft$p.value), 1, ft$p.value))
}

enrich_prog_vs_deg <- function(prog_to_genes, deg_by_celltype, universe = NULL, top_k = NULL) {
  
  normalize_ct_from_program <- function(program_name) {
    ct <- gsub("-\\d+$", "", program_name)
    ct
  }
  
  stopifnot(is.list(prog_to_genes), is.list(deg_by_celltype))
  
  # compute universe if not provided
  if (is.null(universe)) {
    universe <- union(unique(unlist(prog_to_genes, use.names = FALSE)),
                      unique(unlist(deg_by_celltype, use.names = FALSE)))
    universe <- unique(universe)
  }
  
  tibble(program = names(prog_to_genes)) %>%
    mutate(
      cell_type_inferred = normalize_ct_from_program(program),
      # get DEG set matched by inferred cell type (NA if absent)
      deg_genes = map(cell_type_inferred, ~ deg_by_celltype[[.x]] %||% character()),
      prog_genes_full = prog_to_genes[program],
      prog_genes = if (is.null(top_k)) prog_genes_full else map(prog_genes_full, ~ head(.x, top_k)),
      res = pmap(list(prog_genes, deg_genes, list(universe)),
                 ~ fisher_or_sets(..1, ..2, ..3))
    ) %>%
    mutate(
      A = map_int(res, "A"),
      B = map_int(res, "B"),
      C = map_int(res, "C"),
      D = map_int(res, "D"),
      OR = map_dbl(res, "OR"),
      OR_adj = map_dbl(res, "OR_adj"),
      p = map_dbl(res, "p")
    ) %>%
    select(program, cell_type_inferred, A, B, C, D, OR, OR_adj, p) %>%
    mutate(p_adj = p.adjust(p, method = "BH"))
}







################################################################################
### calculate pair-wise similarity between programs
flatten_spectra <- function(spectra_list) {
  out <- list()
  for (ct in names(spectra_list)) {
    x   <- spectra_list[[ct]]
    mat <- as.matrix(x$spectra_score)
    rownames(mat) <- x$spectra_name
    
    for (prog in rownames(mat)) {
      out[[prog]] <- mat[prog, ]
    }
  }
  out
}

similarity_pairwise <- function(program_list,
                                method = c("correlation", "cosine", "jaccard", "all"),
                                min_common_genes = 1) {
  method <- match.arg(method)
  progs  <- names(program_list)
  n      <- length(progs)
  
  # precompute presence sets for Jaccard (non-zero genes)
  present_genes <- lapply(program_list, function(v) names(v)[v != 0])
  
  # initialize result matrices only if needed
  Corr <- Cos <- J <- NULL
  
  if (method %in% c("correlation", "all")) {
    Corr <- matrix(0, n, n, dimnames = list(progs, progs))
  }
  if (method %in% c("cosine", "all")) {
    Cos <- matrix(0, n, n, dimnames = list(progs, progs))
  }
  if (method %in% c("jaccard", "all")) {
    J <- matrix(0, n, n, dimnames = list(progs, progs))
  }
  
  for (i in seq_len(n)) {
    if (!is.null(Corr)) Corr[i, i] <- 1
    if (!is.null(Cos))  Cos[i, i]  <- 1
    if (!is.null(J))    J[i, i]    <- 1
    
    if (i < n) {
      for (j in (i + 1):n) {
        xi <- program_list[[i]]
        xj <- program_list[[j]]
        
        # ---------- correlation / cosine: intersect *all* genes ----------
        if (method %in% c("correlation", "cosine", "all")) {
          g <- intersect(names(xi), names(xj))
          
          if (length(g) >= min_common_genes) {
            v1 <- xi[g]
            v2 <- xj[g]
            
            # correlation (Spearman)
            if (!is.null(Corr)) {
              cval <- suppressWarnings(cor(v1, v2, method = "spearman"))
              if (is.na(cval)) cval <- 0
              Corr[i, j] <- Corr[j, i] <- cval
            }
            
            # cosine
            if (!is.null(Cos)) {
              denom <- sqrt(sum(v1^2)) * sqrt(sum(v2^2))
              cval <- if (denom > 0) sum(v1 * v2) / denom else 0
              if (is.na(cval)) cval <- 0
              Cos[i, j] <- Cos[j, i] <- cval
            }
          } else {
            # not enough overlapping genes
            if (!is.null(Corr)) Corr[i, j] <- Corr[j, i] <- 0
            if (!is.null(Cos))  Cos[i, j]  <- Cos[j, i]  <- 0
          }
        }
        
        # ---------- Jaccard: intersect/union of non-zero genes ----------
        if (!is.null(J)) {
          gi <- present_genes[[i]]
          gj <- present_genes[[j]]
          
          if (length(gi) == 0 && length(gj) == 0) {
            jval <- 0
          } else {
            inter <- length(intersect(gi, gj))
            uni   <- length(union(gi, gj))
            jval  <- if (uni > 0) inter / uni else 0
          }
          if (is.na(jval)) jval <- 0
          J[i, j] <- J[j, i] <- jval
        }
      }
    }
  }
  
  out <- list()
  if (!is.null(Corr)) out$correlation <- Corr
  if (!is.null(Cos))  out$cosine      <- Cos
  if (!is.null(J))    out$jaccard     <- J
  out
}

compare_all_spectra <- function(spectra_list,
                                method = c("correlation", "cosine", "jaccard", "all"),
                                min_common_genes = 1,
                                k = NULL) {
  method <- match.arg(method)
  prog_list <- flatten_spectra(spectra_list)
  
  if (!is.null(k)) {
    prog_list <- lapply(prog_list, function(v) {
      o <- order(v, decreasing = TRUE, na.last = NA)
      topk <- head(o, min(k, length(o)))
      v[topk]
    })
  }
  
  similarity_pairwise(prog_list,
                      method = method,
                      min_common_genes = min_common_genes)
}

################################################################################
# rank genes per group
rank_group_genes <- function(program_group_list,
                             spectra_list,
                             top_n = 30,
                             use_abs = FALSE,
                             drop_na = TRUE) {
  library(dplyr); library(tidyr); library(purrr); library(ggplot2)
  
  if (is.atomic(program_group_list) && !is.null(names(program_group_list))) {
    program_to_group <- program_group_list
    program_group <- split(names(program_to_group), program_to_group)
  } else if (is.list(program_group_list)) {
    program_group <- program_group_list
  } else {
    stop("`program_group_list` must be either a named vector (program->group) or a list (group->programs).")
  }
  
  long_tbl <- purrr::map2_dfr(
    .x = names(spectra_list),
    .y = spectra_list,
    .f = function(ct, x) {
      stopifnot(is.list(x), !is.null(x$spectra_score))
      M <- as.matrix(x$spectra_score)
      rn <- x$spectra_name
      if (is.null(rn)) rn <- rownames(M)
      if (is.null(rn)) rn <- paste0(ct, "-", seq_len(nrow(M)))
      rownames(M) <- rn
      
      tibble::as_tibble(M, rownames = "program") %>%
        tidyr::pivot_longer(-program, names_to = "gene", values_to = "score")
    }
  )
  
  if (use_abs) long_tbl <- long_tbl %>% mutate(score = abs(score))
  
  resolver <- imap_dfr(program_group, ~ tibble(program = .x, group = .y))
  
  matched <- long_tbl %>%
    inner_join(resolver, by = "program")
  
  if (nrow(matched) == 0) {
    stop("No matching programs found between `spectra_list` and `program_group_list`.\n",
         "Check if program names in `spectra_name` match the IDs in your program_group_list.")
  }
  
  groups <- unique(matched$group)
  out <- vector("list", length(groups))
  names(out) <- groups
  
  for (g in groups) {
    g_dat <- matched %>% filter(group == g)
    if (drop_na) g_dat <- g_dat %>% filter(!is.na(score))
    
    g_dat <- g_dat %>%
      group_by(program) %>%
      dplyr::mutate(score_rank = 1 - (min_rank(dplyr::desc(score)) - 1) / (dplyr::n() - 1)) %>%
      ungroup()
    
    gene_rank <- g_dat %>%
      group_by(gene) %>%
      dplyr::summarise(
        n_programs   = dplyr::n_distinct(program),
        median_score = median(score_rank, na.rm = TRUE),
        mean_score   = mean(score_rank, na.rm = TRUE),
        mad_score    = mad(score_rank, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(dplyr::desc(median_score)) %>%
      {if (grepl("Group", g)) {
        dplyr::filter(., n_programs > 1)
      } else {
        .
      }}
    
    top_genes <- head(gene_rank$gene, top_n)
    plot_dat <- g_dat %>%
      filter(gene %in% top_genes) %>%
      mutate(gene = factor(gene, levels = top_genes))
    
    p <- ggplot(plot_dat, aes(x = gene, y = score_rank)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, color = "grey") +
      geom_point() +
      coord_flip() +
      labs(
        x = NULL, y = "Within-program rank (0-1)",
        title = sprintf("%s (top %d genes%s)",
                        g, length(top_genes),
                        if (use_abs) ", |score|" else "")
      ) +
      theme_bw(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 9),
        plot.title  = element_text(size = 12, face = "bold")
      )
    
    out[[g]] <- list(
      table = gene_rank,
      plot  = p
    )
  }
  
  out
}

