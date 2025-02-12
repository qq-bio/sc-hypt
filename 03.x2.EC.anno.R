library(Seurat)
library(Signac)
library(sceasy)
library(reticulate)

reticulate::use_condaenv("cell2location", conda = "/opt/ohpc/pub/apps/anaconda/2022.05/bin/conda", required = TRUE)

sc <- import("scanpy", convert = FALSE)
scvi <- import("scvi", convert = FALSE)



seurat_object_org <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.rds")
seurat_object_nb_1k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.rds")
seurat_object_nb_2k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_2k.rds")
seurat_object_nb_3k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_3k.rds")
seurat_object_nb_4k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_4k.rds")
seurat_object_rec <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.rec.rds")

# ###Find clusters, then run UMAP, and visualize

seurat_object_org <- RunUMAP(seurat_object_org, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_org, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_org, split.by = "tissue", reduction = "umap", pt.size = 1)


seurat_object_nb_1k <- RunUMAP(seurat_object_nb_1k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_nb_1k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_nb_1k, split.by = "strain", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_nb_1k, split.by = "tissue", reduction = "umap", pt.size = 1)

seurat_object_nb_2k <- RunUMAP(seurat_object_nb_2k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_nb_2k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_nb_2k, split.by = "tissue", reduction = "umap", pt.size = 1)

seurat_object_nb_3k <- RunUMAP(seurat_object_nb_3k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_nb_3k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_nb_3k, split.by = "tissue", reduction = "umap", pt.size = 1)

# seurat_object_nb_4k <- RunUMAP(seurat_object_nb_4k, dims = 1:10, reduction = "scvi", n.components = 2)
# DimPlot(seurat_object_nb_4k, group.by = "species", reduction = "umap", pt.size = 1)
# DimPlot(seurat_object_nb_4k, split.by = "tissue", reduction = "umap", pt.size = 1)

seurat_object = seurat_object_nb_1k
# seurat_object = seurat_object_nb_2k
# seurat_object = seurat_object_nb_3k
# seurat_object = seurat_object_nb_4k
marker_list = c("Pecam1", "Egfl7", "Vwf", # EC
                "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                # "Flt1",
                "Plvap", # venous EC
                "Rgcc", # capillary EC
                "Ccl21", "Prox1", # lymphatic EC
                "Igfbp5", # kidney capillary EC
                "Ehd3", # glomerular EC
                "Adgrl3", "Slc38a3" # BBB EC
)
FeaturePlot(seurat_object, features = marker_list)

# negative marker list: "Ackr1", "Gpihbp1", "Fcn3", "Ca4", "Lyve1", "Pdpn", "Sost", "Postn", "Selp", "Sele"

seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution =1)

DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T)

DotPlot(seurat_object, features = marker_list, group.by = "RNA_snn_res.1") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


seurat_object <- subset(seurat_object, seurat_clusters %in% c(1:14, 18))
seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 1.5, 2, 2.5, 3))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")



seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T)
DotPlot(seurat_object, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(seurat_object, features = marker_list, order = TRUE)

DimPlot(seurat_object, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object, split.by = "tissue", reduction = "umap", pt.size = 1)
DimPlot(seurat_object, split.by = "strain", reduction = "umap", pt.size = 1)
DimPlot(seurat_object, split.by = "sxtxt", reduction = "umap", pt.size = 1)



seurat_object$sxt <- paste0(seurat_object$strain, "-", seurat_object$treatment)
DimPlot(seurat_object, split.by = "sxt", reduction = "umap", pt.size = 1)

table(seurat_object$seurat_clusters, seurat_object$subclass_level2)
table(seurat_object$seurat_clusters, seurat_object$tissue)
table(seurat_object$seurat_clusters, seurat_object$sxt)

markers = FindAllMarkers(seurat_object, group.by="seurat_clusters", only.pos = TRUE, min.pct = 0.25)
markers$pct.diff = markers$pct.1 - markers$pct.2
markers = markers %>% group_by(cluster) %>%
  dplyr::arrange(desc(pct.diff), .by_group=TRUE) %>%
  dplyr::mutate(rank = row_number()) %>% ungroup()
write.table(markers,file="/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.DEG.out", sep="\t")


treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")
seurat_object$treatment <- factor(seurat_object$treatment, levels=treatment_order)

prop_data <- seurat_object@meta.data %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
  group_by(seurat_clusters, hypt, strain, treatment) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  left_join(
    seurat_object@meta.data %>%
      dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
      count(hypt, strain, treatment, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = c("hypt", "strain", "treatment")  # Ensure correct join
  ) %>%
  group_by(hypt, strain, treatment, seurat_clusters) %>%
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()


ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = treatment, y = seurat_clusters, fill = log_obs_exp)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Tissue", y = "Cluster", fill = "log(obs/exp)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  ) +
  facet_nested( ~ hypt + strain, scales = "free", space = "free_x")


prop_data <- seurat_object@meta.data %>%
  group_by(seurat_clusters, tissue) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  # Compute expected proportion based on total cell distribution per tissue
  left_join(
    seurat_object@meta.data %>%
      count(tissue, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = "tissue"
  ) %>%
  # Compute log(obs/exp) and adjust by variance
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()

ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = tissue, y = seurat_clusters, fill = log_obs_exp)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Tissue", y = "Cluster", fill = "log(obs/exp)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  )



prop_data <- seurat_object@meta.data %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
  group_by(seurat_clusters, hypt, strain, treatment, tissue) %>%
  dplyr::summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  left_join(
    seurat_object@meta.data %>%
      dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
      count(hypt, strain, treatment, tissue, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = c("hypt", "strain", "treatment", "tissue")  # Ensure correct join
  ) %>%
  group_by(hypt, strain, treatment, tissue, seurat_clusters) %>%
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()


ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = treatment, y = seurat_clusters, fill = log_obs_exp)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Tissue", y = "Cluster", fill = "log(obs/exp)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  ) +
  facet_nested( ~ hypt + strain + tissue, scales = "free", space = "free_x")


generic_anno <- c("0"="Cap EC", "1"="Cap EC", "2"="Cap EC", "3"="Cap EC",
                  "4"="Cap EC", "5"="Cap EC", "6"="Cap EC", "7"="Ven Cap EC",
                  "8"="Ven Cap EC", "9"="Art EC", "10"="BBB EC", "11"="Art EC",
                  "12"="Lymph EC", "13"="Art EC", "14"="Cap EC", "15"="Prolif Cap EC",
                  "16"="BBB EC", "17"="Endocardial EC")






DimPlot(seurat_object, group.by = "subclass_level2", split.by = "tissue", reduction = "umap", pt.size = 3)
DimPlot(seurat_object, group.by = "tissue", split.by = "subclass_level2", reduction = "umap", pt.size = 3)

FeaturePlot(seurat_object, features = c("Vwf", "Fabp4", "Pecam1", "Egfl7", "Flt1", "Sema3g", "Postn", "Plvap", "Pdpn", "Prox1", "Slc38a3"))

FeaturePlot(seurat_object, features = c("Tagln", "Acta2", "Pdgfra", "Pdgfrb", "Notch3", "Myh11", "Kcnj8")) # mural cell
FeaturePlot(seurat_object, features = c("Smoc1", "Inhba", "Npr3")) # endocardial_ec

FeaturePlot(seurat_object, features = c("Slco1a2", "Slc38a3", "Slc7a8"))
FeaturePlot(seurat_object, features = c("Cdh5", "Sema3g", "Rgcc", "Ackr1", "Ccl21"))

FeaturePlot(seurat_object, features = c("Ptgds", "Pla1a", "Adgrl3", "Sh3rf3", "Cdh23", "Trpm6")) # brain vein

FeaturePlot(seurat_object, features = c("Mecom"))
