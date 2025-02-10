library(Matrix)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/scvelo")

lk_ec_mouse = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/mouse.LK.EC.anno.rds")
lk_ec_ss = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.ss.LK.EC.anno.rds")
lk_ec_shr = readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/rat.sp.LK.EC.anno.rds")

lk_ec_mouse = subset(lk_ec_mouse, subclass_level2 %in% setdiff(unique(lk_ec_mouse$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_ss = subset(lk_ec_ss, subclass_level2 %in% setdiff(unique(lk_ec_ss$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))
lk_ec_shr = subset(lk_ec_shr, subclass_level2 %in% setdiff(unique(lk_ec_shr$subclass_level2), c("Neur", "Imm Cell", "PT Cell", "UnID Cell", "VSMC")))

lk_ec_mouse$sxt = paste0(lk_ec_mouse$strain, " - ", lk_ec_mouse$treatment)
lk_ec_ss$sxt = paste0(lk_ec_ss$strain, " - ", lk_ec_ss$treatment)
lk_ec_shr$sxt = paste0(lk_ec_shr$strain, " - ", lk_ec_shr$treatment)


seurat_obj = lk_ec_ss
# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$wnn.umap.harmony@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$wnn.umap.harmony@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='rat.ss.LK.EC.metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='rat.ss.LK.EC.counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$wnn.umap.harmony@cell.embeddings, file='rat.ss.LK.EC.umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='rat.ss.LK.EC.gene_names.csv',
  quote=F,row.names=F,col.names=F
)



seurat_obj = lk_ec_shr
# save metadata table:
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$wnn.umap.harmony@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$wnn.umap.harmony@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file='rat.shr.LK.EC.metadata.csv', quote=F, row.names=F)

# write expression counts matrix
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file='rat.shr.LK.EC.counts.mtx')

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$wnn.umap.harmony@cell.embeddings, file='rat.shr.LK.EC.umap.csv', quote=F, row.names=F)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='rat.shr.LK.EC.gene_names.csv',
  quote=F,row.names=F,col.names=F
)
