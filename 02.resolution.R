library(clustree)
library(ROGUE)

# /xdisk/mliang1/qqiu/project/multiomics/cluster/mouse/mouse.HYP.RNA.cluster.final.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/mouse/mouse.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/mouse/mouse.LK.RNA.cluster.rds


# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.ss.HYP.RNA.cluster.final.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.ss.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.ss.LK.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.sp.HYP.RNA.cluster.final.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.sp.LV.RNA.cluster.rds
# /xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.sp.LK.RNA.cluster.rds


seurat_obj = readRDS("/xdisk/mliang1/qqiu/project/multiomics/cluster/rat/rat.ss.HYP.RNA.cluster.rds")

seurat_expr = seurat_obj@assays$RNA@data
seurat_meta = seurat_obj@meta.data


seurat_expr = seurat_expr
ent.res = SE_fun(seurat_expr)
rogue.value = CalculateRogue(ent.res, platform = "UMI")
rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)

clustree(seurat_meta, prefix = "RNA_snn_res.")





outfolder = paste0(getwd(), "/", resolution, "/")
outfile = paste0(c(species, assay, resolution), collapse = ".")

reso_para = paste0(assay, '_snn_res.', resolution)
seurat_object@active.ident = seurat_object@meta.data[, reso_para]
seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]

### generate plots & tables for first-step cluster result
# cluster_sample_tbl = table(seurat_object$orig.ident, seurat_object$seurat_clusters)
# cluster_sample_sum(cluster_sample_tbl, paste0(outfolder, outfile, ".cluster_sample_sum.xlsx"))

markers = FindAllMarkers(seurat_object, group.by=reso_para, assay = assay, only.pos = TRUE, min.pct = 0.25)
write.table(markers,file=paste0(outfolder, outfile, ".allmarker.0.25.long.txt"), sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                        v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file=paste0(outfolder, outfile, ".allmarker.0.25.wide.txt"), sep="\t", row.names = F)

panglao_anno(marker_tbl, paste0(outfolder, outfile, ".allmarker.0.25.wide.xlsx"))
