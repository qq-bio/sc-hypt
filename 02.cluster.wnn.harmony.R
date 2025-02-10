#!/usr/bin/env Rscript
library(Seurat)
library(ArchR)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)
library(clustree)
library(ROGUE)
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/panglao_anno.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/DE_enrichment.R")
source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/QC_harmony.R")
# source("/xdisk/mliang1/qqiu/project/multiomics-hypertension/src/function/plots.R")

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster")



## process for multi-omics data
angii_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/mouse.LK.multiomics.merged.rds"
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.ss.LK.multiomics.merged.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/QC/rat.sp.LK.multiomics.merged.rds"

infile = sp_infile
seurat_object = readRDS(infile)

### load data & sample QC & remove blacklist genes
DefaultAssay(seurat_object) = "RNA"
seurat_object = subset(seurat_object, nFeature_RNA>200 & 
                         nCount_RNA<20000 & 
                         # doublet_pc.2=="Singlet" &
                         nCount_ATAC < 100000 & 
                         nCount_ATAC > 400 & 
                         TSSEnrichment > 2)


### cell cycle score & to keep the signals separating non-cycling and cycling cells; reference: https://satijalab.org/seurat/articles/cell_cycle_vignette.html
s.genes = intersect(c(str_to_title(cc.genes$s.genes), "Cenpu"), rownames(seurat_object))
g2m.genes = intersect(c(str_to_title(cc.genes$g2m.genes), "Pimreg", "Jpt1"), rownames(seurat_object))
seurat_object = CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
seurat_object$CC.Difference = seurat_object$S.Score - seurat_object$G2M.Score

### process the RNA data
DefaultAssay(seurat_object) = "RNA"
seurat_object = NormalizeData(seurat_object)
seurat_object = FindVariableFeatures(seurat_object)
seurat_object = ScaleData(seurat_object, vars.to.regress = c("percent.mt"))
seurat_object = RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))
seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "pca", 
                           reduction.save = "harmony.rna",assay.use = "RNA", project.dim = F)

## process the ATAC data
DefaultAssay(seurat_object) = "ATAC"
seurat_object@assays$ATAC@counts 
seurat_object= RunTFIDF(seurat_object)
seurat_object = FindTopFeatures(seurat_object, min.cutoff='q0')
seurat_object = RunSVD(seurat_object)
seurat_object = RunHarmony(seurat_object, group.by.vars = "orig.ident", reduction = "lsi", 
                           reduction.save = "harmony.atac", assay.use = "ATAC", project.dim = F)

### integration & cluster
seurat_object = FindMultiModalNeighbors(seurat_object, reduction.list = list("harmony.rna", "harmony.atac"), 
                                        dims.list = list(1:30, 2:30))
seurat_object = RunUMAP(seurat_object, nn.name = "weighted.nn", reduction.name = "wnn.umap.harmony", reduction.key = "wnnUMAPHarmony_")
seurat_object = FindClusters(seurat_object, graph.name = "wsnn", algorithm = 3, 
                             resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1))

# for(assay in c("RNA")){
#   for(resolution in c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)){
#     outfolder = paste0(getwd(), "/multiomics/", resolution, "/")
#     outfile = paste0(c(species, assay, resolution), collapse = ".")
#     
#     reso_para = paste0('wsnn_res.', resolution)
#     seurat_object@active.ident = seurat_object@meta.data[, reso_para]
#     seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]
#     
#     ### generate plots & tables for first-step cluster result
#     cluster_sample_tbl = table(seurat_object$orig.ident, seurat_object$seurat_clusters)
#     cluster_sample_sum(cluster_sample_tbl, paste0(outfolder, outfile, ".cluster_sample_sum.xlsx"))
#     
#     pdf(paste0(outfolder, outfile, ".wnn.umap.pdf"), width=10, height=6)
#     print(DimPlot(seurat_object, reduction = "wnn.umap.harmony", group.by = c("orig.ident", "seurat_clusters"), label=TRUE) & theme(legend.position = "bottom"))
#     # print(DimPlot(seurat_object, reduction = "wnn.umap.harmony", group.by = c("tissue", "treatment"), label=TRUE) & theme(legend.position = "bottom"))
#     print(DimPlot(seurat_object, reduction = "wnn.umap.harmony", group.by = c("Phase"), label=TRUE) & theme(legend.position = "bottom") & 
#             FeaturePlot(seurat_object, features = "percent.mt", label=TRUE) )
#     dev.off()
#     
#   }
# }


reso_para = paste0('wsnn_res.', 0.4)
seurat_object@active.ident = seurat_object@meta.data[, reso_para]
seurat_object$seurat_clusters = seurat_object@meta.data[, reso_para]
saveRDS(seurat_object, paste0(species, ".multiomics.rna.cluster.rds"))






################################################################################
# cluster resolution evaluation by clustree and ROGUE
angii_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.cluster.rds"
ss_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.cluster.rds"
sp_infile = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.cluster.rds"

infile = angii_infile
seurat_object = readRDS(infile)

seurat_meta = seurat_obj@meta.data
clustree(seurat_meta, prefix = "RNA_snn_res.")

for(i in c(0.4)){
  cluster = paste0("wsnn_res.", i)
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/DEG/mouse.LK.reso_", i)
  DE_analysis(seurat_object, cluster = cluster, outfile = outfile)
}



seurat_expr = seurat_obj@assays$RNA@data


seurat_expr = seurat_expr
ent.res = SE_fun(seurat_expr)
rogue.value = CalculateRogue(ent.res, platform = "UMI")
rogue.res <- rogue(expr, labels = meta$ct, samples = meta$Patient, platform = "UMI", span = 0.6)




































DefaultAssay(seurat_object) = "ATAC"
gene.activities = GeneActivity(seurat_object)
seurat_object[['gene.activity']] = CreateAssayObject(counts = gene.activities)
seurat_object = NormalizeData(
  object = seurat_object,
  assay = 'gene.activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(seurat_object$nCount_)
)


seurat_object = readRDS(paste0(species, ".multiomics.rna.cluster.rds"))
markers = FindAllMarkers(seurat_object, assay = "RNA", only.pos = TRUE, min.pct = 0.25)
write.table(markers,file=paste0(outfolder, outfile, ".multiomics.allmarker.0.25.long.txt"), sep="\t")

marker_tbl = markers[markers$p_val_adj<0.05 & markers$avg_log2FC>0.1,] %>% group_by(cluster) %>%
  dplyr::mutate(id = row_number()) %>% as.data.frame() %>%
  select(id, cluster, gene) %>% reshape(., idvar = "id", timevar = "cluster",
                                        v.names="gene", sep=" ", direction = "wide")
colnames(marker_tbl) = sub("gene", "Cluster", colnames(marker_tbl))
write.table(marker_tbl,file=paste0(outfolder, outfile, ".multiomics.allmarker.0.25.wide.txt"), sep="\t", row.names = F)

panglao_anno(marker_tbl, paste0(outfolder, outfile, ".multiomics.allmarker.0.25.wide.xlsx"))






















### generate gene activity score using cicero
library(cicero)
library(SeuratWrappers)
library(rtracklayer)

annotations = import("/scratch/u/qqiu/refdata/cellranger-arc/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz")
gene_anno = unique(annotations[annotations$type=="gene" &
                                seqnames(annotations) %in% standardChromosomes(annotations) , c("gene_name")])

promoter = promoters(gene_anno, upstream=1000, downstream=1000)
promoterDF = data.frame(chromosome=seqnames(promoter),
                        start=start(promoter),
                        end=end(promoter),
                        gene=promoter$gene_name) %>% unique()

mm10 = genomeAnnoMm10$chromSizes
mm10 = as.data.frame(mm10)[,c(1,4)]


seurat_object = readRDS("MLK1_4.multiomics.rna.cluster.rds")
DefaultAssay(seurat_object) = "ATAC"

ds_cds = as.cell_data_set(x = seurat_object, assay="ATAC", default.reduction = "wnn.umap.harmony")
ds_cicero = make_cicero_cds(ds_cds, reduced_coordinates = reducedDims(ds_cds)$WNN.UMAP.HARMONY)

ds_conns = run_cicero(ds_cicero, mm10, sample_num = 100) 



input_cds = annotate_cds_by_site(input_cds, gene_annotation_sub)

# generate unnormalized gene activity matrix
unnorm_ga = build_gene_activity_matrix(input_cds, conns)

# remove any rows/columns with all zeroes
unnorm_ga = unnorm_ga[!Matrix::rowSums(unnorm_ga) == 0, 
                       !Matrix::colSums(unnorm_ga) == 0]

# make a list of num_genes_expressed
num_genes = pData(input_cds)$num_genes_expressed
names(num_genes) = row.names(pData(input_cds))

# normalize
cicero_gene_activities = normalize_gene_activities(unnorm_ga, num_genes)

# if you had two datasets to normalize, you would pass both:
# num_genes should then include all cells from both sets
unnorm_ga2 = unnorm_ga
cicero_gene_activities = normalize_gene_activities(list(unnorm_ga, unnorm_ga2), 
                                                    num_genes)