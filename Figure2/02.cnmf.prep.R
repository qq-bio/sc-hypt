suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(data.table)
})



################################################################################
### Unify mouse and rat gene symbols
################################################################################
# functions
build_key2uni <- function(mapping_df) {
  dt <- as.data.table(mapping_df)[,
                                  .(Gene.name, Rat.gene.name, Gene.stable.ID, Rat.gene.stable.ID, uni_symbol)
  ]
  
  long <- melt(
    dt,
    measure.vars = c("Gene.name", "Rat.gene.name",
                     "Gene.stable.ID", "Rat.gene.stable.ID"),
    value.name = "key"   
  )[, .(key, uni_symbol)]
  
  long <- long[!is.na(key) & key != ""]
  
  setNames(long$uni_symbol, long$key)
}

remap_to_unisymbol <- function(obj, key2uni, assay = "RNA",
                               drop_unmapped = TRUE,
                               agg_mode = c("sum","maxmean")) {
  agg_mode <- match.arg(agg_mode)
  stopifnot(inherits(obj, "Seurat"))
  
  ass <- obj@assays[[assay]]
  is_v5 <- inherits(ass, "Assay5")  # v5 vs v4 flag
  
  # read counts (v4 vs v5)
  M_in <- if (is_v5) {
    tryCatch(GetAssayData(obj, assay = assay, slot = "counts"),
             error = function(e) ass[["counts"]])
  } else {
    ass@counts
  }
  if (is.null(M_in)) stop("Assay has no counts matrix.")
  
  feat   <- rownames(M_in)
  target <- unname(key2uni[feat])  # NA for unmapped
  
  # handle unmapped
  if (drop_unmapped) {
    keep <- !is.na(target)
    if (!any(keep)) stop("No features mapped to uni_symbol.")
  } else {
    keep <- rep(TRUE, length(feat))
    target[is.na(target)] <- feat[is.na(target)]  # retain original names
  }
  
  tgt_kept <- target[keep]
  dup_flag <- any(duplicated(tgt_kept))
  
  M <- M_in[keep, , drop = FALSE]
  if (inherits(M, "lgCMatrix")) M <- as(M, "dgCMatrix")
  orig_cols <- colnames(M_in)
  
  # finalize helpers
  .finalize_v5 <- function(obj, assay, counts_mat, orig_cols) {
    colnames(counts_mat) <- orig_cols
    obj <- SetAssayData(obj, assay = assay, slot = "counts",
                        new.data = as(counts_mat, "dgCMatrix"))
    obj <- SetAssayData(obj, assay = assay, slot = "data",
                        new.data = GetAssayData(obj, assay = assay, slot = "counts"))
    # Assay5 has no scale.data
    suppressWarnings(try(
      obj <- SetAssayData(obj, assay = assay, slot = "scale.data",
                          new.data = matrix(0, 0, 0)), silent = TRUE))
    VariableFeatures(obj, assay = assay) <- character(0)
    obj
  }
  
  .finalize_v4 <- function(obj, assay, counts_mat, orig_cols) {
    A <- obj@assays[[assay]]
    colnames(counts_mat) <- orig_cols
    stopifnot(ncol(counts_mat) == ncol(A@counts))
    
    A@counts        <- as(counts_mat, "dgCMatrix")
    A@data          <- A@counts  
    A@scale.data    <- matrix(0, nrow = 0, ncol = 0)
    A@meta.features <- data.frame(row.names = rownames(A@counts))
    if (!is.null(A@var.features)) A@var.features <- character(0)
    
    obj@assays[[assay]] <- A
    obj
  }
  
  .finalize_assay <- function(obj, assay, counts_mat, orig_cols) {
    if (is_v5) .finalize_v5(obj, assay, counts_mat, orig_cols)
    else       .finalize_v4(obj, assay, counts_mat, orig_cols)
  }
  
  if (!dup_flag) {
    rownames(M) <- tgt_kept
    obj <- .finalize_assay(obj, assay, M, orig_cols)
    
    return(list(
      object = obj,
      report = list(
        mode = "rename_only",
        n_rows_in = length(feat),
        n_mapped_rows = length(tgt_kept),
        n_unmapped_dropped = if (drop_unmapped) sum(is.na(target)) else 0,
        n_uni_symbols_out = nrow(M),
        had_duplicates = FALSE,
        n_collapsed_symbols = 0
      )
    ))
  }
  
  # aggregate duplicates
  f <- factor(tgt_kept, levels = unique(tgt_kept))  # preserves first-seen order
  
  if (agg_mode == "sum") {
    MT <- as(M, "dgTMatrix")
    counts_new <- sparseMatrix(
      i    = as.integer(f)[MT@i + 1],
      j    = MT@j + 1,
      x    = MT@x,
      dims = c(nlevels(f), ncol(M))
    )
    rownames(counts_new) <- levels(f)
    counts_new <- as(counts_new, "dgCMatrix")
  } else {
    # pick row with largest row-mean per group
    ord   <- order(f)
    M_ord <- M[ord, , drop = FALSE]
    f_ord <- f[ord]
    rmv   <- Matrix::rowMeans(M_ord)
    
    r <- rle(as.integer(f_ord))
    starts <- cumsum(c(1, head(r$lengths, -1)))
    ends   <- cumsum(r$lengths)
    
    pick_idx_ord <- integer(length(r$values))
    for (k in seq_along(r$values)) {
      s <- starts[k]; e <- ends[k]
      pick_idx_ord[k] <- s + which.max(rmv[s:e]) - 1
    }
    pick_idx <- ord[pick_idx_ord]
    counts_new <- M[pick_idx, , drop = FALSE]
    rownames(counts_new) <- levels(f)
    counts_new <- as(counts_new, "dgCMatrix")
  }
  
  # finalize
  obj <- .finalize_assay(obj, assay, counts_new, orig_cols)
  
  n_collapsed <- sum(table(tgt_kept) > 1)
  list(
    object = obj,
    report = list(
      mode = paste0("aggregate_", agg_mode),
      n_rows_in = length(feat),
      n_mapped_rows = length(tgt_kept),
      n_unmapped_dropped = if (drop_unmapped) sum(is.na(target)) else 0,
      n_uni_symbols_out = nrow(counts_new),
      had_duplicates = TRUE,
      n_collapsed_symbols = n_collapsed
    )
  )
}


exclude_patterns <- c("^mt-", "^Rps", "^Rpl", "^Mrps", "^Mrpl",
                      "^Igh", "^Igk", "^Igl",
                      "^Tr", "^Hist", "^Hba", "^Hbb", "^Olfr", "^Gm", "^.*Rik$")

m2r_gene_list <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2rat.out.txt", header = T, sep = "\t")
m2r_gene_ortho <- m2r_gene_list[m2r_gene_list$Rat.orthology.confidence..0.low..1.high.==1, ]

all(m2r_gene_ortho$Gene.name!="")
all( !(is.na(m2r_gene_ortho$Gene.name)))
m2r_gene_ortho$uni_symbol <- m2r_gene_ortho$Gene.name

exclude_genes <- m2r_gene_ortho$uni_symbol[grepl(paste(exclude_patterns, collapse="|"),
                                                 m2r_gene_ortho$uni_symbol, ignore.case=FALSE)]

m2r_gene_ortho <- m2r_gene_ortho[!(m2r_gene_ortho$uni_symbol %in% exclude_genes), ]

key2uni <- build_key2uni(m2r_gene_ortho)



################################################################################
### cell type-wise HVG detection
################################################################################
# functions
make_union_tbl <- function(lst) {
  all <- bind_rows(lapply(names(lst), function(ds) {
    inner <- lst[[ds]]
    bind_rows(lapply(names(inner), function(ct) {
      tibble(
        gene = inner[[ct]],
        dataset = ds,
        cell_type = ct
      )
    }))
  }))
  
  all %>%
    dplyr::count(gene, name = "n_hits") %>%
    dplyr::arrange(desc(n_hits), gene)
}


input_file <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LK.multiomics.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.MSA.RNA.anno.L2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LV.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MSA.RNA.anno.L2.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MCA.RNA.anno.L2.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.PBMC.RNA.anno.L2.rds"
)

hvg_list <- list(`2k` = list())
log_df   <- list()

for (i in input_file) {
  message("Processing: ", i)
  so_outfile <- sub("\\.rds$", ".uni_symbol.rds", i)
  
  tryCatch({
    so <- readRDS(i)
    
    res <- remap_to_unisymbol(so, key2uni, assay = "RNA", agg_mode = "sum")
    so  <- res$object
    saveRDS(so, so_outfile)
    
    rm(res)
    
    if (!("Cell_type_L1" %in% colnames(so@meta.data))) {
      stop("Missing 'Cell_type_L1' in meta.data for: ", basename(i))
    }
    
    DefaultAssay(so) <- "RNA"
    so <- NormalizeData(so, verbose = FALSE)
    
    ct_levels <- sort(unique(so$Cell_type_L1))
    
    if(any(c("IMM", "Immune cell") %in% ct_levels)){
      IMM_present = TRUE
      immune_levels <- sort(unique(so@meta.data[so$Cell_type_L1 %in% c("IMM", "Immune cell"), ]$Cell_type_L2))
      ct_levels <- c(ct_levels, immune_levels)
    }else{
      IMM_present = FALSE
      immune_levels = c()
    }
    
    hvg_list[["2k"]][[basename(i)]]     <- list()
    
    for (ct in ct_levels) {
      
      if(ct %in% immune_levels & IMM_present){
        cells_ct <- WhichCells(so, expression = Cell_type_L2 == ct)
        if (length(cells_ct) < 30) { 
          message("  - Skipping small group: ", ct, " (n=", length(cells_ct), ")")
          next
        }
      }else{
        cells_ct <- WhichCells(so, expression = Cell_type_L1 == ct)
        if (length(cells_ct) < 30) { 
          message("  - Skipping small group: ", ct, " (n=", length(cells_ct), ")")
          next
        }
      }
      
      so_ct <- subset(so, cells = cells_ct)
      
      so_ct_vst <- FindVariableFeatures(so_ct, selection.method = "vst",
                                        nfeatures = 2000, verbose = FALSE)
      hvg_vst <- VariableFeatures(so_ct_vst)
      
      hvg_list[["2k"]][[basename(i)]][[ct]] <- hvg_vst
      
      rm(so_ct, so_ct_vst)
    }
    
    rm(so)
    
  }, error = function(e) {
    warning("Failed on: ", i, " | ", conditionMessage(e))
  })
}

union_2k_tbl <- make_union_tbl(hvg_list[["2k"]])
union_2k <- unique(union_2k_tbl$gene)

saveRDS(
  list(
    hvg_per_dataset = hvg_list,
    union_2k_tbl = union_2k_tbl,
    union_2k = union_2k
  ),
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/union_hvg_2k.by_celltype.rds"
)












################################################################################
### Generate per-cell type input for cNMF
################################################################################
# functions
merge_hvg_per_ct <- function(hvg_list) {
  
  hvg_by_ct <- list()
  
  for (tissue in names(hvg_list)) {
    sub <- hvg_list[[tissue]]
    for (ct in names(sub)) {
      genes <- sub[[ct]]
      if (is.null(genes)) next
      hvg_by_ct[[ct]] <- unique(c(hvg_by_ct[[ct]], genes))
    }
  }
  
  return(hvg_by_ct)
  
}

pad_align_sparse <- function(M, target_genes) {
  miss <- setdiff(target_genes, rownames(M))
  if (length(miss)) {
    M0 <- Matrix(0, nrow=length(miss), ncol=ncol(M), sparse=TRUE,
                 dimnames=list(miss, colnames(M)))
    M <- rbind(M, M0)
  }
  M[target_genes, , drop=FALSE]
}


input_files <- c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.HYP.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LV.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/mouse.LK.multiomics.anno.L2.uni_symbol.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.HYP.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LV.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.LK.multiomics.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.MSA.RNA.anno.L2.uni_symbol.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.HYP.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LV.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.LK.multiomics.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MSA.RNA.anno.L2.uni_symbol.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.sp.MCA.RNA.anno.L2.uni_symbol.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/L2_anno/rat.ss.PBMC.RNA.anno.L2.uni_symbol.rds"
)

hvg_list <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/union_hvg_2k.by_celltype.rds")

CT_COL <- "Cell_type_L1"
CT_COL_ALT <- "Cell_type_L2"
SAMPLE_COL <- c("Sample_ID")
MIN_CELLS_TISSUE_CT <- 20 

ct_hvg <- merge_hvg_per_ct(hvg_list$hvg_per_dataset$`2k`)
ct_names <- names(ct_hvg)

per_ct_counts <- setNames(vector("list", length(ct_names)), ct_names)
per_ct_barcodes <- setNames(vector("list", length(ct_names)), ct_names)

for (f in input_files) {
  message("Reading: ", f)
  so <- readRDS(f)
  
  ds_name <- tools::file_path_sans_ext(basename(f))
  
  for (ct in ct_names) {
    
    if (ct %in% as.character(so[[CT_COL]][, 1])) {
      Idents(so) <- CT_COL
    } else if (ct %in% as.character(so[[CT_COL_ALT]][, 1])) {
      Idents(so) <- CT_COL_ALT
    } else {
      next
    }
    
    genes_ct <- ct_hvg[[ct]]
    if (length(genes_ct) < 200) next
    
    cells_ct <- WhichCells(so, idents = ct)
    if (length(cells_ct) < MIN_CELLS_TISSUE_CT) next
    
    if (!length(cells_ct)) next
    
    M <- GetAssayData(so, assay="RNA", slot="counts")[, cells_ct, drop=FALSE]
    M <- M[intersect(rownames(M), genes_ct), , drop=FALSE]
    M <- pad_align_sparse(M, genes_ct)
    
    # unique barcodes include dataset and ct
    bc <- paste0(colnames(M), "|", ds_name, "|", ct)
    colnames(M) <- bc
    
    per_ct_counts[[ct]][[ds_name]] <- if (is.null(per_ct_counts[[ct]][[ds_name]])) M else M
    per_ct_barcodes[[ct]][[ds_name]] <- bc
  }
  rm(so)
}

root_out <- "/xdisk/mliang1/qqiu/project/multiomics-hypertension/CNMF/input/hvg.per_celltype"
dir.create(root_out, recursive=TRUE, showWarnings=FALSE)

for (ct in ct_names) {
  parts <- per_ct_counts[[ct]]
  if (length(parts)==0) next
  counts_ct <- do.call(cbind, parts)
  ct_folder = gsub("[ /]", "_", ct)
  outdir <- file.path(root_out, ct_folder)
  dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
  
  Matrix::writeMM(counts_ct, file.path(outdir, "matrix.mtx"))
  all_bc <- unlist(per_ct_barcodes[[ct]], use.names=FALSE)
  readr::write_tsv(data.frame(all_bc), file.path(outdir, "barcodes.tsv"), col_names=FALSE)
  
  genes_ct <- ct_hvg[[ct]]
  feats <- data.frame(gene_id=genes_ct, gene_name=genes_ct, type="Gene Expression")
  readr::write_tsv(feats, file.path(outdir, "genes.tsv"), col_names=FALSE)
  
  message(sprintf("Wrote %s: genes=%d, cells=%d", ct, nrow(counts_ct), ncol(counts_ct)))
}


































