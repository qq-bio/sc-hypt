library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(sceasy)
library(dplyr)


### overall estimate
seurat = sceasy::convertFormat("/xdisk/mliang1/qqiu/data/HeartCellAtlas/visium-OCT_LV_lognormalised.updated.h5ad", from = "anndata", to = "seurat")

FeaturePlot(seurat, "IL1R1", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "EC6_ven", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "EC7_endocardial", reduction = "spatial", split.by = "sangerID")

FeaturePlot(seurat, "SMAD6", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "EC5_art", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "EC7_endocardial", reduction = "spatial", split.by = "sangerID")

FeaturePlot(seurat, "KCNT2", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "LRRC3B", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "HMCN1", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "UNC5C", reduction = "spatial", split.by = "sangerID")


ec_list = colnames(seurat@meta.data)[grepl("EC", colnames(seurat@meta.data))]
for(i in ec_list){
  tryCatch({
    result <- cor.test(seurat@meta.data[, i], seurat@assays$RNA@data["IL1R1", ])
    print(c(i, result$estimate, result$p.value))
  }, error = function(e) {
  })
}




### marker expression and cell abundance correlation
marker_list = c("FHOD3", "UNC5C", "IL1R1", "HMCN1", "VWF", #m5813
                "NEBL", "ST6GALNAC3", "MAYO10", "LRRC3B", "BTNL9", "KCNT2", "SMAD6", "RGS3" #c1,7,9
)
marker_list_use = intersect(marker_list, rownames(seurat))

expr_matrix = t(as.matrix(seurat@assays$RNA@data[marker_list_use, ]))
cell_matrix = as.matrix(seurat@meta.data[, c(18:89)])
ec_matrix = as.matrix(seurat@meta.data[, c(32:40, 88)])

cor_cell_p <- cor(expr_matrix, cell_matrix, use="pairwise.complete.obs", method = "pearson")
cor_cell_s <- cor(expr_matrix, cell_matrix, use="pairwise.complete.obs", method = "spearman")
cor_ec_p <- cor(expr_matrix, ec_matrix, use="pairwise.complete.obs", method = "pearson")
cor_ec_s <- cor(expr_matrix, ec_matrix, use="pairwise.complete.obs", method = "spearman")

cor_mat = cor_ec_p
cor_matrix_melted <- melt(cor_mat)
ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(min(cor_matrix_melted$value), max(cor_matrix_melted$value))) +
  theme_minimal() +
  labs(title = "Marker expression and cell abundance correlation", 
       x = "", y = "", fill = "Pearson\ncorrelation") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))


cor_mat = cor_ec_s
cor_matrix_melted <- melt(cor_mat)
ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(min(cor_matrix_melted$value), max(cor_matrix_melted$value))) +
  theme_minimal() +
  labs(title = "Marker expression and cell abundance correlation", 
       x = "", y = "", fill = "Spearmam\ncorrelation") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))



# library(lsa)
combined <- cbind(expr_matrix, cell_matrix)
combined <- cbind(expr_matrix, ec_matrix)

sim_all <- cosine(as.matrix(combined))
n1 <- ncol(expr_matrix)
n2 <- ncol(ec_matrix)
sim_expr_vs_cell <- sim_all[1:n1, (n1+1):(n1+n2)]

cor_mat = sim_expr_vs_cell
cor_matrix_melted <- melt(cor_mat)
ggplot(cor_matrix_melted, aes(Var2, Var1, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                       limits = c(min(cor_matrix_melted$value), max(cor_matrix_melted$value))) +
  theme_minimal() +
  labs(title = "Marker expression and cell abundance correlation", 
       x = "", y = "", fill = "Cosine\nsimilarity") + 
  theme(axis.text = element_text(color = "black", size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1))




### cell type specific signature and moran's i
library(spdep)

r2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.rat2human.out.txt", header = T, sep = "\t")
m2h = read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2human.out.txt", header = T, sep = "\t")
r2h_ortho_1 = data.frame(gene_rodent = c(r2h[r2h$Human.orthology.confidence..0.low..1.high.==1, ]$Gene.stable.ID,
                                         r2h[r2h$Human.orthology.confidence..0.low..1.high.==1, ]$Gene.name,
                                         m2h[m2h$Human.orthology.confidence..0.low..1.high.==1, ]$Gene.stable.ID,
                                         m2h[m2h$Human.orthology.confidence..0.low..1.high.==1, ]$Gene.name),
                         gene_human = c(r2h[r2h$Human.orthology.confidence..0.low..1.high.==1, ]$Human.gene.name,
                                        r2h[r2h$Human.orthology.confidence..0.low..1.high.==1, ]$Human.gene.name,
                                        m2h[m2h$Human.orthology.confidence..0.low..1.high.==1, ]$Human.gene.name,
                                        m2h[m2h$Human.orthology.confidence..0.low..1.high.==1, ]$Human.gene.name))
r2h_ortho_1 = unique(r2h_ortho_1[rowSums(is.na(r2h_ortho_1))==0, ])

r2h_ortho_all = data.frame(gene_rodent = c(r2h[r2h$Human.gene.name!="", ]$Gene.stable.ID,
                                         r2h[r2h$Human.gene.name!="", ]$Gene.name,
                                         m2h[m2h$Human.gene.name!="", ]$Gene.stable.ID,
                                         m2h[m2h$Human.gene.name!="", ]$Gene.name),
                         gene_human = c(r2h[r2h$Human.gene.name!="", ]$Human.gene.name,
                                        r2h[r2h$Human.gene.name!="", ]$Human.gene.name,
                                        m2h[m2h$Human.gene.name!="", ]$Human.gene.name,
                                        m2h[m2h$Human.gene.name!="", ]$Human.gene.name))
r2h_ortho_all = unique(r2h_ortho_all[rowSums(is.na(r2h_ortho_all))==0, ])


deg_list = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.merged.DEG_all.out", header = T, sep = " ")
deg_list = deg_list[deg_list$p_val_adj<0.05 & abs(deg_list$avg_log2FC)>0.25, ] %>%
  group_by(cell_type) %>% arrange(desc(pct.diff)) %>% ungroup()

common_ec_list = c("Pecam1", "Egfl7")
lv_ec_list = c("Nav3", "Ablim3", "Ccdc85a", "Ncald", "Nrp1", "Kitlg", "Dach1", "Arhgap18", "Ank3")

top_n = 10
m0610_marker_list = c(deg_list[deg_list$cell_type=="M0610", ]$gene_name[1:top_n], common_ec_list, lv_ec_list)
m0610_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m0610_marker_list, ]$gene_human)

m24_marker_list = c(deg_list[deg_list$cell_type=="M24", ]$gene_name[1:top_n], common_ec_list, lv_ec_list)
m24_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m24_marker_list, ]$gene_human)

m5813_marker_list = c(deg_list[deg_list$cell_type=="M5813", ]$gene_name[1:top_n], common_ec_list, lv_ec_list)
m5813_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m5813_marker_list, ]$gene_human)

c7_marker_list = c(deg_list[deg_list$cell_type=="C7", ]$gene_name[1:top_n], common_ec_list, lv_ec_list)
c7_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% c7_marker_list, ]$gene_human)

marker_list = list("M0610" = m0610_marker_human,
                   "M24" = m24_marker_human,
                   "M5813" = m5813_marker_human,
                   "C7" = c7_marker_human)


m0610_marker_list = c("Hdac9", "Nrp2", "Lamc1", "Myo10")
m0610_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m0610_marker_list, ]$gene_human)

m24_marker_list = c("Hdac9", "Nrp2", "Lamc1", "Ano4")
m24_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m24_marker_list, ]$gene_human)

m5813_marker_list = c("Hdac9", "Fhod3", "Unc5c", "Il1r1", "Hmcn1", "Vwf")
m5813_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% m5813_marker_list, ]$gene_human)

c7_marker_list = c("Myo10", "Lrrc3b", "Btnl9", "Kcnt2", "Nebl", "St6galnac3")
c7_marker_human = unique(r2h_ortho_all[r2h_ortho_all$gene_rodent %in% c7_marker_list, ]$gene_human)

marker_list = list("M0610" = m0610_marker_human,
                   "M24" = m24_marker_human,
                   "M5813" = m5813_marker_human,
                   "C7" = c7_marker_human)

# seurat = AddModuleScore(seurat, features = marker_list, name = names(marker_list), slot = "data")

slide_ids <- unique(seurat$sangerID)
signature_list = paste0(names(marker_list), c(1, 2, 3, 4))

seurat@meta.data[, signature_list] <- 0
for (slide in slide_ids) {
  slide_cells <- WhichCells(seurat, expression = sangerID == slide)
  slide_obj <- subset(seurat, cells = slide_cells)
  
  slide_obj <- AddModuleScore(
    slide_obj,
    features = marker_list,
    name = names(marker_list),
    slot = "data"
  )
  
  # Store back or collect scores
  seurat@meta.data[slide_cells, signature_list] <- slide_obj@meta.data[, signature_list]
}

seurat_meta = seurat@meta.data


FeaturePlot(seurat, "M06101", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "M242", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "M58133", reduction = "spatial", split.by = "sangerID")
FeaturePlot(seurat, "C74", reduction = "spatial", split.by = "sangerID")



all_results <- list()
signature_list = paste0(names(marker_list), c(1, 2, 3, 4))
signature_list = c("EC5_art", "EC6_ven", "EC7_endocardial")

for (slide_id in unique(seurat_meta$sangerID)) {
  slide_results <- list()  # Store results for all signatures for this slide
  
  coords <- seurat_meta[seurat_meta$sangerID == slide_id, c("array_row", "array_col")]
  
  knn <- knearneigh(coords, k = 6)
  nb <- knn2nb(knn)
  listw <- nb2listw(nb, style = "W")
  
  for (signature_id in signature_list) {
    
    # Extract signature score for this slide and signature
    gene_expr <- seurat_meta[seurat_meta$sangerID == slide_id, signature_id]
    
    # Check for enough non-NA values
    if (sum(!is.na(gene_expr)) >= 3) {  # Moran’s I requires at least 3 non-NA values
      moran_i <- moran.test(gene_expr, listw)
      slide_results[[signature_id]] <- moran_i$estimate[["Moran I statistic"]]
    } else {
      slide_results[[signature_id]] <- NA
    }
  }
  
  all_results[[slide_id]] <- slide_results
}
# Combine into a data frame
moran_df <- do.call(cbind, all_results)  # rows = genes, cols = slides
moran_df




all_results <- list()
gene_list = c("IL1R1", "SMAD6")

for (slide_id in unique(seurat_meta$sangerID)) {
  slide_results <- list()  # Store results for all signatures for this slide
  
  coords <- seurat_meta[seurat_meta$sangerID == slide_id, c("array_row", "array_col")]
  
  knn <- knearneigh(coords, k = 6)
  nb <- knn2nb(knn)
  listw <- nb2listw(nb, style = "W")
  
  for (gene in gene_list) {
    
    # Extract signature score for this slide and signature
    gene_expr <- seurat@assays$RNA$data[gene, seurat_meta$sangerID == slide_id]
    
    # Check for enough non-NA values
    if (sum(!is.na(gene_expr)) >= 3) {  # Moran’s I requires at least 3 non-NA values
      moran_i <- moran.test(gene_expr, listw)
      slide_results[[gene]] <- moran_i$estimate[["Moran I statistic"]]
    } else {
      slide_results[[gene]] <- NA
    }
  }
  
  all_results[[slide_id]] <- slide_results
}
# Combine into a data frame
moran_df <- do.call(cbind, all_results)  # rows = genes, cols = slides
moran_df





moran_long <- as.data.frame(moran_df) %>%
  tibble::rownames_to_column("signature") %>%
  tidyr::pivot_longer(
    cols = -signature,
    names_to = "sample",
    values_to = "moran_i"
  ) %>%
  mutate(moran_i = as.numeric(moran_i))

ggplot(moran_long, aes(x = signature, y = moran_i)) +
  geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) +
  stat_summary(fun = mean, geom = "point", color = "red", size = 3, shape = 18) +  # optional mean dot
  labs(
    x = "Signature",
    y = "Moran's I",
    title = "Moran's I per EC Signature"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )

