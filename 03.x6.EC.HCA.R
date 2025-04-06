library(Seurat)
library(SeuratDisk)
library(hdf5r)
library(sceasy)


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











