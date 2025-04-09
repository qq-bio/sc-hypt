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

seurat_object_strain_1k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_1k.rds")
seurat_object_strain_2k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_2k.rds")
seurat_object_strain_3k <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.strain_batch.gene_nb.hvg_3k.rds")
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

seurat_object_strain <- RunUMAP(seurat_object_strain, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_strain, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain, split.by = "tissue", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain, split.by = "strain", reduction = "umap", pt.size = 1)

seurat_object_strain_1k <- RunUMAP(seurat_object_strain_1k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_strain_1k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_1k, split.by = "tissue", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_1k, split.by = "strain", reduction = "umap", pt.size = 1)

seurat_object_strain_2k <- RunUMAP(seurat_object_strain_2k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_strain_2k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_2k, split.by = "tissue", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_2k, split.by = "strain", reduction = "umap", pt.size = 1)

seurat_object_strain_3k <- RunUMAP(seurat_object_strain_3k, dims = 1:10, reduction = "scvi", n.components = 2)
DimPlot(seurat_object_strain_3k, group.by = "species", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_3k, split.by = "tissue", reduction = "umap", pt.size = 1)
DimPlot(seurat_object_strain_3k, split.by = "strain", reduction = "umap", pt.size = 1)

# seurat_object_nb_4k <- RunUMAP(seurat_object_nb_4k, dims = 1:10, reduction = "scvi", n.components = 2)
# DimPlot(seurat_object_nb_4k, group.by = "species", reduction = "umap", pt.size = 1)
# DimPlot(seurat_object_nb_4k, split.by = "tissue", reduction = "umap", pt.size = 1)

# seurat_object = seurat_object_strain_1k
# seurat_object = seurat_object_strain_2k
# seurat_object = seurat_object_strain_3k
# seurat_object = seurat_object_nb_2k
# seurat_object = seurat_object_nb_3k
# seurat_object = seurat_object_nb_4k

seurat_object = seurat_object_nb_1k
seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 2))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)

Idents(seurat_object) = "RNA_snn_res.1"

contamination_list = c("Pecam1", "Egfl7", "Vwf", "Tagln", "Acta2", "Pdgfra", "Pdgfrb", "Notch3", "Myh11", "Kcnj8", "Myh6", "Mb", "Lrp2")
# DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T, group.by = "seurat_clusters")
DotPlot(seurat_object, features = contamination_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(seurat_object, features = c("Pecam1", "Mb", "Lrp2"), ncol = 3)
FeaturePlot(seurat_object, features = c("Trpv4", "Kcnd1", "Piezo1"), ncol = 3)


marker_list = c("Pecam1", "Egfl7", "Vwf", # EC
                "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                # "Flt1",
                "Plvap", # venous EC
                "Rgcc", # capillary EC
                "Ccl21", "Prox1", "Lyve1", # lymphatic EC
                "Igfbp5", # kidney capillary EC
                "Ehd3", # "Emcn", "Sost", "Plat", # glomerular EC
                "Adgrl3", "Slc38a3", # BBB EC
                "Npr3"
)
FeaturePlot(seurat_object, features = marker_list)

# negative marker list: "Ackr1", "Gpihbp1", "Fcn3", "Ca4", "Lyve1", "Pdpn", "Sost", "Postn", "Selp", "Sele"


DimPlot(seurat_object, reduction = "umap", pt.size = 1, label = T)

DotPlot(seurat_object, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


seurat_object <- subset(seurat_object, RNA_snn_res.1 %in% c(1:6, 8:14, 18))
seurat_object <- FindNeighbors(seurat_object, dims = 1:10, reduction = "scvi")
seurat_object <- FindClusters(seurat_object, resolution = c(1, 1.5, 2, 2.5, 3))
seurat_object <- RunUMAP(seurat_object, dims = 1:10, reduction = "scvi", n.components = 2)
Idents(seurat_object) = "RNA_snn_res.1.5"
seurat_object$seurat_clusters = seurat_object$RNA_snn_res.1.5

DotPlot(seurat_object, features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
DotPlot(seurat_object, features = contamination_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")


################################################################################

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



################################################################################
### merge clusters
seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.rds")
cluster_order = c(0, 10, 6, 2, 4, 18, 5, 8, 13, 
                  1, 7, 9, 
                  12, 21, 15, 19,
                  20,14, 23, 
                  22, 16, 
                  3, 11, 17)
seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels=cluster_order)
Idents(seurat_object) = "seurat_clusters"

marker_list = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cross-organ_EC/DEG/ec.scvi.gene_nb.hvg_1k.refined.cluster_wise.DEG.out")
marker_list = marker_list %>% arrange(factor(cluster, levels = cluster_order)) %>% filter(rank<=5) %>% select(gene) %>% unique()

marker_list= c("Nav3",'Ablim3',"Ccdc85a",'Ncald',"Nrp1",'Kitlg',"Dach1","Arhgap18",'Ank3',
               "Hdac9","Nrp2","Lamc1",'Ano4',
               "Nav2",'Diaph3',"Mki67","Rad51b","Sdk1",
               'Fhod3',"Unc5c","Il1r1","Hmcn1","Vwf",
               "Myo10","Lrrc3b","Btnl9",'Kcnt2',
               "Nebl",'St6galnac3',
               "Prdm16","Ptprj","Sulf1","Auts2l1",'Eln',"Myof","Pcdh7","Alcam",
               "Slco1a4","Lef1","Gpcpd1",'Ccdc141',
               "Zfp521","Nuak1","Mctp1",
               'Pkhd1l1','Cgnl1',"Cdh11",
               "Meis2","Pbx1","Ldb2",
               "Rbms3","Inpp4b","Fmnl2","Chrm3","Malat1",
               'Pkhd1',"Erbb4",'Ca12')


DotPlot(seurat_object, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pre_marker_list = c("Pecam1", "Egfl7", "Vwf", # EC
                "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                # "Flt1",
                "Plvap", # venous EC
                "Rgcc", # capillary EC
                "Ccl21", "Prox1", "Lyve1", # lymphatic EC
                "Igfbp5", # kidney capillary EC
                "Ehd3", # "Emcn", "Sost", "Plat", # glomerular EC
                "Adgrl3", "Slc38a3", # BBB EC
                "Npr3"
)
DotPlot(seurat_object, features = pre_marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


seurat_object$seurat_clusters <- paste0("C", seurat_object$seurat_clusters)
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C5", "C8", "C13"), ]$seurat_clusters <- "M5813"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C2", "C4"), ]$seurat_clusters <- "M24"
seurat_object@meta.data[seurat_object$seurat_clusters %in% c("C0", "C6", "C10"), ]$seurat_clusters <- "M0610"

saveRDS(seurat_object, "/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.merged.rds")




seurat_object <- readRDS("/xdisk/mliang1/qqiu/project/multiomics-hypertension/subcluster/ec.scvi.gene_nb.hvg_1k.refined.merged.rds")
cluster_order = c("M0610", "M24", "C18", "M5813",
                  "C1", "C7", "C9", 
                  "C12", "C21", "C15", "C19",
                  "C20", "C14", "C23", 
                  "C22", "C16", 
                  "C3", "C11", "C17")

seurat_object$seurat_clusters <- factor(seurat_object$seurat_clusters, levels=cluster_order)
Idents(seurat_object) = "seurat_clusters"

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


marker_list= c("Nav3",'Ablim3',"Ccdc85a",'Ncald',"Nrp1",'Kitlg',"Dach1","Arhgap18",'Ank3',
               "Hdac9","Nrp2","Lamc1",'Ano4',
               "Nav2",'Diaph3',"Mki67","Rad51b","Sdk1",
               'Fhod3',"Unc5c","Il1r1","Hmcn1","Vwf",
               "Myo10","Lrrc3b","Btnl9",'Kcnt2',
               "Nebl",'St6galnac3',
               "Prdm16","Ptprj","Sulf1","Auts2l1",'Eln',"Myof","Pcdh7","Alcam",
               "Slco1a4","Lef1","Gpcpd1",'Ccdc141',
               "Zfp521","Nuak1","Mctp1",
               'Pkhd1l1','Cgnl1',"Cdh11",
               "Meis2","Pbx1","Ldb2",
               "Rbms3","Inpp4b","Fmnl2","Chrm3","Malat1",
               'Pkhd1',"Erbb4",'Ca12')


DotPlot(seurat_object, features = marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pre_marker_list = c("Pecam1", "Egfl7", "Vwf", # EC
                    "Sulf1", "Col8a1", "Eln", "Sema3g", # arterial EC
                    # "Flt1",
                    "Plvap", # venous EC
                    "Rgcc", # capillary EC
                    "Ccl21", "Prox1", "Lyve1", # lymphatic EC
                    "Igfbp5", # kidney capillary EC
                    "Ehd3", # "Emcn", "Sost", "Plat", # glomerular EC
                    "Adgrl3", "Slc38a3", # BBB EC
                    "Npr3"
)
DotPlot(seurat_object, features = pre_marker_list) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


DimPlot(seurat_object, group.by = "seurat_clusters", reduction = "umap", pt.size = 1, label = T)




meta_list <- c(
  # ---- Glycolysis ----
  "Hk1", "Hk2", # Hexokinases
  "Gpi1",       # Glucose-6-phosphate isomerase
  "Pfkp", "Pfkm", "Pfkl", # Phosphofructokinase isoforms
  "Pfkfb3",     # PFKFB3 (regulator of glycolysis)
  "Aldoa",      # Aldolase A
  "Tpi1",       # Triosephosphate isomerase
  "Gapdh",      # Glyceraldehyde-3-phosphate dehydrogenase
  "Pgk1",       # Phosphoglycerate kinase 1
  "Pgam1",      # Phosphoglycerate mutase 1
  "Eno1",       # Enolase 1
  "Pkm",        # Pyruvate kinase, muscle
  "Ldha",       # Lactate dehydrogenase A
  
  # ---- TCA cycle ----
  "Cs",         # Citrate synthase
  "Aco2",       # Aconitase 2
  "Idh3a", "Idh3b", "Idh3g", # Mitochondrial IDH3 subunits
  "Ogdh",       # Oxoglutarate dehydrogenase
  "Sdha", "Sdhb", # Succinate dehydrogenase subunits
  "Fh",         # Fumarate hydratase
  "Mdh1", "Mdh2", # Malate dehydrogenases (cytosolic and mitochondrial)
  
  # ---- Oxidative Phosphorylation (ETC) ----
  # Complex I
  "Ndufa9", "Ndufb5", "Ndufs1",
  # Complex II
  "Sdha", "Sdhb",   # (already included above, but can be repeated if needed)
  # Complex III
  "Uqcrc1", "Uqcrc2",
  # Complex IV
  "Cox4i1", "Cox7a2",
  # Complex V
  "Atp5f1a", "Atp5f1b",
  
  # ---- Fatty Acid β-Oxidation ----
  "Cpt1a", "Cpt2",  # Carnitine palmitoyltransferases
  "Acadl", "Acadm", "Acads",  # Acyl-CoA dehydrogenases (long-, medium-, short-chain)
  "Acox1", "Acox2",  # Acyl-CoA oxidase
  "Acaa2",          # 3-Ketoacyl-CoA thiolase
  "Hadha", "Hadhb", # Hydroxyacyl-CoA dehydrogenase subunits
  
  # ---- PPP (Pentose Phosphate Pathway) ----
  "G6pd",   # Glucose-6-phosphate dehydrogenase
  "Pgls",   # 6-phosphogluconolactonase
  "Pgd",    # 6-phosphogluconate dehydrogenase
  "Tkt",    # Transketolase
  "Taldo1", # Transaldolase
  "Rpia",   # Ribose 5-phosphate isomerase A
  
  # ---- Serine Biosynthesis ----
  "ENSRNOG00000068650", "ENSRNOG00000066787", "ENSRNOG00000069953", "Phgdh",  # Phosphoglycerate dehydrogenase
  "Psat1",  # Phosphoserine aminotransferase
  "Psph",   # Phosphoserine phosphatase
  
  # ---- Glutamine Metabolism (GLS-related) ----
  "Gls",    # Glutaminase (Kidney type)
  "Gls2",   # Glutaminase 2
  "Glud1",  # Glutamate dehydrogenase 1 (GDH1)
  "Slc1a5", # Neutral amino acid transporter B(0)
  "Gpt2",    # Alanine transaminase 2 (ALT2)
  
  # ---- Superoxide Dismutases (SOD) ----
  "Sod1",   # Cytosolic superoxide dismutase
  "Sod2",   # Mitochondrial superoxide dismutase
  "Sod3",   # Extracellular superoxide dismutase
  
  # ---- Catalase ----
  "Cat",    # Catalase
  
  # ---- Glutathione System ----
  "Gclc",   # Glutamate-cysteine ligase, catalytic subunit
  "Gclm",   # Glutamate-cysteine ligase, modifier subunit
  "Gss",    # Glutathione synthetase
  "Gsr",    # Glutathione-disulfide reductase
  "Gpx1",   # Glutathione peroxidase 1
  "Gpx2",   # Glutathione peroxidase 2
  "Gpx4",   # Glutathione peroxidase 4
  "Glrx1",  # Glutaredoxin 1
  "Glrx2",  # Glutaredoxin 2
  
  # ---- Thioredoxin System ----
  "Txn1",   # Thioredoxin 1 (cytosolic)
  "Txn2",   # Thioredoxin 2 (mitochondrial)
  "Txnrd1", # Thioredoxin reductase 1 (cytosolic)
  "Txnrd2", # Thioredoxin reductase 2 (mitochondrial)
  "Txnrd3", # Thioredoxin reductase 3 (testis-specific, but sometimes relevant)
  
  # ---- Peroxiredoxins ----
  "Prdx1",
  "Prdx2",
  "Prdx3",
  "Prdx4",
  "Prdx5",
  "Prdx6",
  
  # ---- Heme Oxygenases ----
  "Hmox1",  # Heme oxygenase 1
  "Hmox2",  # Heme oxygenase 2
  
  # ---- NRF2 Pathway ----
  "Nfe2l2", # NRF2 transcription factor
  "Keap1",   # NRF2 repressor
  
  # ---- NADPH-consuming vasculoprotective genes ----
  # eNOS (nitric oxide synthase 3)
  "Nos3",
  # Prostaglandin G/H synthase 1 (Cyclooxygenase-1)
  "Ptgs1",
  # Glutaredoxin
  "Glrx"
)

laminar_shear_stress_genes <- c(
  "Klf2",   # Kruppel-like factor 2, key shear-responsive transcription factor
  "Klf4",   # Kruppel-like factor 4, also upregulated by laminar flow
  "Nos3",   # eNOS, promoted by laminar flow for NO production
  "Nfe2l2", # NRF2, redox-sensitive TF, upregulated in certain shear contexts
  "Hmox1",  # Heme oxygenase-1, protective enzyme often induced under shear
  "Thbd",   # Thrombomodulin, anti-thrombotic factor upregulated by shear
  "Bmx"     # Non-receptor tyrosine kinase, implicated in flow signaling
)

# ---- Disturbed/Low/Oscillatory Flow-Responsive Genes ----
disturbed_flow_genes <- c(
  "Nfkb1",  # NF-κB subunit, activated under disturbed flow → inflammation
  "Icam1",  # Intercellular adhesion molecule
  "Vcam1",  # Vascular cell adhesion molecule
  "Edn1",   # Endothelin-1, promotes vasoconstriction, often up in disturbed flow
  "Cxcl1",  # Chemokine (C-X-C motif) ligand 1
  "Ccl2",   # Chemokine (C-C motif) ligand 2 (a.k.a. MCP-1)
  "Selp",   # P-selectin, upregulated in activated/endothelium
  "Egr1",    # Early growth response 1, immediate-early gene in stress response
  "Nox4"
)

# ---- Hypoxia-Responsive Genes ----
hypoxia_genes <- c(
  "Hif1a",   # Hypoxia-inducible factor-1 alpha
  "Epas1",   # HIF-2 alpha, also known as EPAS1
  "Egln1",   # PHD2 (prolyl hydroxylase domain protein 2), key HIF regulator
  "Vegfa",   # Vascular endothelial growth factor A
  "Pdk1",    # Pyruvate dehydrogenase kinase 1
  "Slc2a1",  # Glucose transporter 1 (GLUT1)
  "Pgk1",    # Phosphoglycerate kinase 1
  "Ca9",     # Carbonic anhydrase IX
  "Bnip3",   # BCL2 interacting protein 3, involved in autophagy/apoptosis
  "Epo"      # Erythropoietin (may be low in ECs, but a classic hypoxia-inducible gene)
)

# fatty acid transporter
fatty_trans_genes <- c(
  "Slc27a1",
  "Slc27a4",
  "Cd36",
  "Mfsd2a"
)


DotPlot(seurat_object, features = unique(meta_list)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = unique(laminar_shear_stress_genes)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = unique(disturbed_flow_genes)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = unique(hypoxia_genes)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = c("Cpt1a", "Mfsd2a", "Pfkfb3")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

DotPlot(seurat_object, features = unique(fatty_trans_genes)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")
seurat_object$treatment <- factor(seurat_object$treatment, levels=treatment_order)


seurat_object@meta.data %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
  group_by(orig.ident, seurat_clusters, hypt, strain, treatment) %>%
  dplyr::summarise(cell_count = n()) %>%
  ungroup() %>% group_by(orig.ident) %>% 
  dplyr::mutate(proportion = cell_count / sum(cell_count)) %>% 
  group_by(hypt, strain, treatment, seurat_clusters) %>%  # Group for mean and SD calculation
  dplyr::summarise(
    mean_proportion = mean(proportion),
    sd_proportion = sd(proportion)
  ) %>% 
  ggplot(aes(x = treatment, y = mean_proportion, fill = seurat_clusters)) +  # Plot the results
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +  # Create a bar plot with means
  geom_errorbar(aes(ymin = mean_proportion - sd_proportion, ymax = mean_proportion + sd_proportion),
                position = position_dodge(width = 0.8), width = 0.2) +  # Add error bars
  scale_y_continuous(labels = scales::percent) +  # Format the y-axis as percentages
  # scale_fill_manual(values = lk_ec_col) +
  labs(x = "Treatment", y = "Cell proportion") +  # Add labels and title
  theme(    axis.text.y = element_text(colour = 'black'),
            axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
            legend.position = "None"
  ) +
  facet_nested(seurat_clusters ~ hypt + strain, scales = "free") 



prop_data <- seurat_object@meta.data %>%
  dplyr::mutate(hypt = ifelse(strain %in% c("SD", "WKY"), "Normotensive", "Hypertensive")) %>%
  group_by(seurat_clusters, sxtxt, tissue, hypt, strain) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  mutate(proportion = cell_count / sum(cell_count)) %>%
  ungroup() %>%
  # Compute expected proportion based on total cell distribution per tissue
  left_join(
    seurat_object@meta.data %>%
      count(sxtxt, name = "total_cells") %>%
      mutate(expected_prop = total_cells / sum(total_cells)),
    by = "sxtxt"
  ) %>%
  # Compute log(obs/exp) and adjust by variance
  mutate(
    log_obs_exp = log(proportion / expected_prop)
  ) %>%
  ungroup()

ggplot(prop_data[prop_data$log_obs_exp>0, ], aes(x = sxtxt, y = seurat_clusters, fill = log_obs_exp)) +
  geom_tile(color="black") +  # Heatmap-style visualization
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Tissue", y = "Cluster", fill = "log(obs/exp)") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    axis.text.y = element_text(colour = 'black')
  ) +
  facet_nested( ~ tissue + hypt + strain, scales = "free", space = "free_x")




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




################################################################################



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

FeaturePlot(seurat_object, features = c("Tagln", "Acta2", "Pdgfra", "Pdgfrb", "Notch3", "Myh11", "Kcnj8", "Myh6", "Mb")) # mural cell & cardiomyocyte
FeaturePlot(seurat_object, features = c("Smoc1", "Inhba", "Npr3")) # endocardial_ec

FeaturePlot(seurat_object, features = c("Slco1a2", "Slc38a3", "Slc7a8"))
FeaturePlot(seurat_object, features = c("Cdh5", "Sema3g", "Rgcc", "Ackr1", "Ccl21"))

FeaturePlot(seurat_object, features = c("Ptgds", "Pla1a", "Adgrl3", "Sh3rf3", "Cdh23", "Trpm6")) # brain vein

FeaturePlot(seurat_object, features = c("Mecom"))





################################################################################
###

seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 5000)

r2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.rat2human.out.txt", sep = "\t", header = T)
m2h <- read.table("/xdisk/mliang1/qqiu/reference/biomaRt/biomaRt.gene.mouse2human.out.txt", sep = "\t", header = T)

Tip_cell_list = c("ADM", "ANGPT2", "ANKRD37", "APLN", "C1QTNF6", "CD93", "CLDN5", 
                  "COL4A1", "COL4A2", "COTL1", "CXCR4", "DLL4", "EDNRB", "ESM1", 
                  "FSCN1", "GPIHBP1", "HSPG2", "IGFBP3", "INHBB", "ITGA5", "JUP", 
                  "KCNE3", "KCNJ8", "KDR", "LAMA4", "LAMB1", "LAMC1", "LXN", "MARCKS", 
                  "MARCKSL1", "MCAM", "MEST", "MYH9", "MYO1B", "N4BP3", "NID2", 
                  "NOTCH4", "PDGFB", "PGF", "PLOD1", "PLXND1", "PMEPA1", "PTN", 
                  "RAMP3", "RBP1", "RGCC", "RHOC", "SMAD1", "SOX17", "SOX4", "SPARC", 
                  "TCF4", "UNC5B", "VIM")
Stalk_cell_list = c("ACKR1", "AQP1", "C1QTNF9", "CD36", "CSRP2", "EHD4", "FBLN5", 
                    "HSPB1", "LIGP1", "IL6ST", "JAM2", "LGALS3", "LRG1", "MEOX2", 
                    "PLSCR2", "SDPR", "SELP", "SPINT2", "TGFBI", "TGM2", "TMEM176A", 
                    "TMEM176B", "TMEM252", "TSPAN7", "VEGFR1", "VWF")

Tip_cell_list_use = unique(c(r2h[r2h$Human.gene.name %in% Tip_cell_list, ]$Gene.name,
                             m2h[m2h$Human.gene.name %in% Tip_cell_list, ]$Gene.name)) %>%
  setdiff(., c("")) # %>% intersect(., VariableFeatures(seurat_object))

Stalk_cell_list_use = unique(c(r2h[r2h$Human.gene.name %in% Stalk_cell_list, ]$Gene.name,
                             m2h[m2h$Human.gene.name %in% Stalk_cell_list, ]$Gene.name)) %>%
  setdiff(., c("")) # %>% intersect(., VariableFeatures(seurat_object))

seurat_object <- AddModuleScore(seurat_object, features=list(Tip_cell_list_use), name = "tip_cell_score", assay = "RNA")
seurat_object <- AddModuleScore(seurat_object, features=list(Stalk_cell_list_use), name = "stalk_cell_score", assay = "RNA")

FeaturePlot(seurat_object, features = c("tip_cell_score1"))
FeaturePlot(seurat_object, features = c("stalk_cell_score1"))




library(UCell)

gene.sets <- list(tip_cell_score = Tip_cell_list_use,
                  stalk_cell_score = Stalk_cell_list_use)

seurat_object <- AddModuleScore_UCell(seurat_object, features = gene.sets)

FeaturePlot(seurat_object, features = c("tip_cell_score_UCell"))
FeaturePlot(seurat_object, features = c("stalk_cell_score_UCell"))






dflow_genes <- c("Icam1", "Vcam1", "Yap", "Bmp4", "Vegfa")
sflow_genes <- c("Klf2", "Klf4", "Nos3", "Cdh5")

gene.sets <- list(dflow_score = dflow_genes,
                  sflow_score = sflow_genes)

seurat_object <- AddModuleScore_UCell(seurat_object, features = gene.sets)

FeaturePlot(seurat_object, features = c("dflow_score_UCell"))
FeaturePlot(seurat_object, features = c("sflow_score_UCell"))
