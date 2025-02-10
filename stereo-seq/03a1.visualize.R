library(Seurat)
library(sctransform)
library(harmony)
library(dplyr)
library(stringr)
library(ggplot2)

setwd("/scratch/g/mliang/multiomics-CKD/cluster")

# plot_dir = "/scratch/u/qqiu/project/multiomics/plot/UMAP/"/scratch/g/mliang/multiomics-CKD


seurat_object_rna = readRDS("multiomics-CKD.rna.cluster.rds")
# seurat_object_sct = readRDS("multiomics-CKD.sct.cluster.rds")
Idents(seurat_object_rna) <- "RNA_snn_res.0.9"
seurat_object_rna$seurat_clusters = seurat_object_rna$RNA_snn_res.0.9
# Idents(seurat_object_sct) <- "SCT_snn_res.0.6"


disease.id <- c("ADN", "ADN", "ADN", 
                "AHN", "AHN", "AHN", 
                "NC", "NC")

seurat_object_rna@meta.data$"disease.id" = factor(disease.id[match(seurat_object_rna$orig.ident, paste0("patient", c(1:8)))],
                                                  levels= c("NC", "ADN", "AHN"))

DimPlot(seurat_object_rna, reduction = "umap", group.by = c("disease.id")) & theme(legend.position = "bottom")


marker_list = c("NPHS2", "NPHS1", # Podocytes
                "ALDH1A2","RBFOX1", "KCNMB2", # Parietal epithelial cell
                "DACH2", "NCAM1", # Parietal epithelial cells
                "LRP2","CUBN" ,"SLC5A12", # Proximal tubules
                "CDH6", "PTPRD", # Injured proximal tubules
                "SLC12A1", "UMOD", "NOS1",# Thick ascending limbs
                # "ERBB4",
                
                "SLC12A3", # Distal convoluted tubules
                "SLC8A1", # Connecting tubules
                "AQP2", # Collecting ducts
                
                "SLC26A7", # Type A intercalated cells
                "KIT", # Type A intercalated cells
                "SLC26A4", # Type B intercalated cells
                
                "DCN", "CALD1", "PDGFRB", # Fibroblasts
                "PECAM1", "PDK1",# Endothelial cells
                
                "PTPRC", "CD247", "SAMD3", # T cells
                "BANK1", # B cells
                "HLA-DRA", "LYZ",  # Activated macrophages
                "F13A1", # Alternatively activated macrophages
                # "PDK1", "JCHAIN", # Plasma cells
                
                "CNTNAP2", "CD36"
                )

# marker_list = unique(marker_list)

new.cluster.ids <- c("T cells(C0)", "Injured proximal tubules(C1)", "Thick ascending limbs(C2)", 
                     "Collecting ducts(C3)", "Proximal tubules(C4)", "Distal convoluted tubules(C5)", 
                     "Fibroblasts(C6)", "Endothelial cells(C7)", "Thick ascending limbs(C8)", 
                     "B cells(C9)", "Activated macrophages(C10)", "Connecting tubules(C11)", "T cells(C12)", 
                     "Type A intercalated cells(C13)", "Proximal tubules(C14)", "Vascular endothelial cells(C15)", 
                     "Thick ascending limbs(C16)", "Alternatively activated macrophages(C17)",
                     "Proximal tubules(C18)", "Type B intercalated cells(C19)", "(C20)", "Fibroblasts/Mesangial cells(C21)", 
                     "Parietal epithelial cells(C22)", "(C23)", "Parietal epithelial cells(C24)", "Podocytes(C25)", 
                     "Type A intercalated cells(C26)", "Parietal epithelial cells(C27)", "(C28)")

seurat_object_rna@meta.data$"new.cluster.ids" = factor(new.cluster.ids[match(seurat_object_rna$RNA_snn_res.0.9, c(0:28))],
                                                   levels= c("Podocytes(C25)", 
                                                             "Parietal epithelial cells(C22)", "Parietal epithelial cells(C24)", "Parietal epithelial cells(C27)", 
                                                             "Proximal tubules(C4)", "Proximal tubules(C14)", "Proximal tubules(C18)",
                                                             "Injured proximal tubules(C1)", 
                                                             "Thick ascending limbs(C2)", "Thick ascending limbs(C8)", "Thick ascending limbs(C16)", 
                                                             "Distal convoluted tubules(C5)", "Connecting tubules(C11)", "Collecting ducts(C3)", 
                                                             "Type A intercalated cells(C13)", "Type A intercalated cells(C26)", "Type B intercalated cells(C19)",
                                                             "Fibroblasts(C6)", "Fibroblasts/Mesangial cells(C21)", 
                                                             "Endothelial cells(C7)", "Vascular endothelial cells(C15)",
                                                             "T cells(C0)", "T cells(C12)", 
                                                             "B cells(C9)", "Activated macrophages(C10)", "Alternatively activated macrophages(C17)",
                                                             "(C20)", "(C23)", "(C28)"
                                                   ))


Idents(seurat_object_rna) = "new.cluster.ids"
DotPlot(seurat_object_rna, # group.by = "new.cluster.ids", 
        features = marker_list) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))




new.cluster.ids_umap <- c("T cells", "Injured proximal tubules", "Thick\nascending\nlimbs", 
                          "Collecting\nducts", "Proximal tubules", "Distal\nconvoluted\ntubules", 
                          "Fibroblasts", "Endothelial cells", "Thick\nascending\nlimbs", 
                          "B cells", "Activated macrophages", "Connecting\ntubules", "T cells", 
                          "Type A intercalated cells", "Proximal tubules", "Vascular\nendothelial cells", 
                          "Thick\nascending\nlimbs", "Alternatively activated macrophages",
                          "Proximal tubules", "Type B\nintercalated cells", "C20", "Fibroblasts/Mesangial cells", 
                          "Parietal epithelial cells", "C23", "Parietal epithelial cells", "Podocytes", 
                          "Type A intercalated cells", "Parietal epithelial cells", "C28")

seurat_object_rna@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object_rna$RNA_snn_res.0.9, c(0:28))])

# output = paste0(gsub("RNA.cluster.rds", "", infile), reso, ".UMAP.pdf")
# pdf(output, width=6, height=6)
DimPlot(seurat_object_rna, reduction = "umap", group.by = c("new.cluster.ids_umap"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
# dev.off()



seurat_object_rna_nc = subset(seurat_object_rna, disease.id %in% c("NC"))
DimPlot(seurat_object_rna_nc, reduction = "umap", group.by = c("new.cluster.ids"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
DimPlot(seurat_object_rna_nc, reduction = "umap", group.by = c("new.cluster.ids"), label=F) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")

seurat_object_rna_adn = subset(seurat_object_rna, disease.id %in% c("ADN"))
DimPlot(seurat_object_rna_adn, reduction = "umap", group.by = c("new.cluster.ids"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
DimPlot(seurat_object_rna_adn, reduction = "umap", group.by = c("new.cluster.ids"), label=F) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")

seurat_object_rna_ahn = subset(seurat_object_rna, disease.id %in% c("AHN"))
DimPlot(seurat_object_rna_ahn, reduction = "umap", group.by = c("new.cluster.ids"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
DimPlot(seurat_object_rna_ahn, reduction = "umap", group.by = c("new.cluster.ids"), label=F) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")



C23_marker_list = c("SRRM2", "SETD2", "NSD1", "ANKHD1", "RBM6", "TTC3")
c20_marker_list = c("CNTNAP2", "RBFOX1", "PTPRD", "CASC4", "WDPCP", "SIPA1L3", "ZFHX3")
C28_marker_list = c("CCL21", "SEMA3A", "CD36", "MMRN1", "FLRT2")
DotPlot(seurat_object_rna, # group.by = "new.cluster.ids", 
        features = c("PECAM1")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))





infile="mouse.RNA.cluster.rds"
seurat_object = readRDS(infile)

new.cluster.ids <- c("Neurons(C0)", "Endothelial cells(C1)", "Astrocytes(C2)", "Myelinating Oligodendrocytes(C3)", "Cardiomyocytes(C4)", 
                     "Macrophages(C5)", "Neurons(C6)", "Fibroblasts(C7)", "Neurons(C8)", "Vascular mural cells(C9)", "Oligodendrocyte progenitor cells(C10)", 
                     "Neurons(C11)", "Neurons(C12)", "Neurons(C13)", "Tanycytes(C14)", "Cardiomyocytes(C15)", "Microglial cells(C16)", "Neurons(C17)")

seurat_object@meta.data$"new.cluster.ids" = factor(new.cluster.ids[match(seurat_object$RNA_snn_res.0.2, c(0:17))],
                                                   levels= c("Neurons(C0)", "Neurons(C6)", "Neurons(C8)", "Neurons(C11)", "Neurons(C12)", "Neurons(C13)", "Neurons(C17)", 
                                                             "Endothelial cells(C1)", "Astrocytes(C2)", "Myelinating Oligodendrocytes(C3)", "Cardiomyocytes(C4)", "Cardiomyocytes(C15)", 
                                                             "Macrophages(C5)", "Fibroblasts(C7)", "Vascular mural cells(C9)", "Oligodendrocyte progenitor cells(C10)", 
                                                             "Tanycytes(C14)", "Microglial cells(C16)"))

Idents(seurat_object) <- "tissue"
seurat_object_list = SplitObject(seurat_object, split.by = "tissue")
# Idents(seurat_object) = "new.cluster.ids"
DotPlot(seurat_object_list[[3]], group.by = "new.cluster.ids", 
        features = c("Nrg1", "Meg3", 
                     "Flt1", 
                     "Slc1a2",
                     "Mbp", "St18",
                     "Ttn", "Mb",
                     "Mrc1", "Ptprc",
                     "Ranbp3l", "Dcn",
                     "Atp13a5", "Pdgfrb",
                     "Cspg4", "Pdgfra",
                     "Col23a1",
                     "Bnc2", "Foxp2")) + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))











infile="tissue_sep/mouse.HYP.RNA.cluster.rds"
reso="RNA_snn_res.0.3"
seurat_object = readRDS(infile)
Idents(seurat_object) <- reso


gene_list = c(
  "Gad1", "Gad2", "Slc32a1", # inhibitory
  "Slc17a6", "Slc17a8", # excitatory
  ""
)

marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6",
                "Pecam1", 
                "Slc1a2",
                "Mbp", "St18",
                "Ranbp3l", "Dcn", "Ptgds",
                "Atp13a5", "Pdgfrb",
                "Cspg4", "Pdgfra",
                "Col23a1",
                "Bnc2", "Foxp2", "P2ry12")

ependymal = c("Foxj1", "Pifo", "Dynlrb2")

# FeaturePlot(seurat_object, gene_list) & theme(legend.position = "bottom")

DotPlot(seurat_object, 
        features = ependymal)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))



## check about the vascular related genes
gene_list = c("Cldn5", "Adgrf5", "Emcn", "Vtn", "Cspg4", "Atp13a5", "Pth1r", "Kcnj8", "Abcc9", "Apln", "Cd82", "Chst1")
DotPlot(seurat_object, 
        features = gene_list)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
FeaturePlot(seurat_object, gene_list) & theme(legend.position = "bottom")




## merge final information
new.cluster.ids_dot <- c("Inhibitory neurons(C0)", "Excitatory neurons(C1)", "Astrocytes(C2)", "Myelinating Oligodendrocytes(C3)", "Excitatory neurons(C4)", 
                         "Excitatory neurons(C5)", "Oligodendrocyte progenitor cells(C6)", "Tanycytes(C7)", "Inhibitory neurons(C8)", "Microglial cells(C9)", "Endothelial cells(C10)", 
                         "Pericytes(C11)", "Fibroblast(C12)", "Inhibitory neurons(C13)", "Neurons(C14)", "(C15)", "Inhibitory neurons(C16)", "Astrocytes(C17)")

seurat_object@meta.data$"new.cluster.ids_dot" = factor(new.cluster.ids_dot[match(seurat_object$RNA_snn_res.0.3, c(0:17))],
                                                       levels= c("Inhibitory neurons(C0)", "Inhibitory neurons(C8)", "Inhibitory neurons(C13)", "Inhibitory neurons(C16)", 
                                                                 "Excitatory neurons(C1)", "Excitatory neurons(C4)", "Excitatory neurons(C5)", "Neurons(C14)", 
                                                                 "Astrocytes(C2)", "Astrocytes(C17)", "Myelinating Oligodendrocytes(C3)", 
                                                                 "Oligodendrocyte progenitor cells(C6)", "Tanycytes(C7)", "Microglial cells(C9)", "Endothelial cells(C10)", 
                                                                 "Pericytes(C11)", "Fibroblast(C12)", "(C15)"))
Idents(seurat_object) = "new.cluster.ids_dot"
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6",
                "Slc1a2",
                "Mbp", "St18",
                "Cspg4", "Pdgfra",
                "Col23a1",
                "Tgfbr1", 
                "Flt1", 
                "Atp13a5",
                "Ranbp3l",
                "Bnc2")

output = paste0(gsub("RNA.cluster.rds", "", infile), reso, ".dotplot.pdf")
pdf(output, width=10, height=5)
DotPlot(seurat_object, 
        features = marker_list)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()

new.cluster.ids_umap <- c("Inhibitory neurons", "Excitatory neurons", "Astrocytes", "Myelinating Oligodendrocytes", "Excitatory neurons", 
                          "Excitatory neurons", "Oligodendrocyte\nprogenitor cells", "Tanycytes", "Inhibitory neurons", "Microglial cells", "Endothelial cells", 
                          "Pericytes", "Fibroblast", "Inhibitory neurons", "Neurons", "C15", "Inhibitory neurons", "Astrocytes")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$RNA_snn_res.0.3, c(0:17))])

output = paste0(gsub("RNA.cluster.rds", "", infile), reso, ".UMAP.pdf")
pdf(output, width=6, height=6)
DimPlot(seurat_object, reduction = "umap", group.by = c("new.cluster.ids_umap"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
dev.off()






infile="tissue_sep/mouse.LV.RNA.cluster.rds"
seurat_object = readRDS(infile)
Idents(seurat_object) <- "RNA_snn_res.0.4"

FeaturePlot(seurat_object, c("Actn2", "Gja5", "Scn5a", "Irx3", "Cacna2d2")) & theme(legend.position = "bottom")
FeaturePlot(seurat_object, c("Actn2", "Scn5a"))& lims(x=c(6,13), y=c(0,7.5))

seurat_object_tmp = subset(seurat_object, RNA_snn_res.0.4 %in% c(1))
FeatureScatter(seurat_object_tmp, feature1 = "Actn2", feature2 = "Scn5a")
c1=rownames(seurat_object@meta.data[seurat_object@meta.data$RNA_snn_res.0.4==1,])
table(seurat_object@assays$RNA@data['Actn2',c1]>0, 
      seurat_object@assays$RNA@data['Scn5a',c1]>0)

DimPlot(seurat_object, reduction = "umap", group.by = c("RNA_snn_res.0.4"), label=TRUE) & theme(legend.position = "bottom")
FeaturePlot(seurat_object, c("Tspan7", "Ptprc")) & theme(legend.position = "bottom")



new.cluster.ids_dot <- c("Endothelial cells(C0)", "Cardiomyocytes(C1)", "Macrophages(C2)", "Endothelial cells(C3)", 
                         "Fibroblasts(C4)", "Cardiomyocytes(C5)", "Pericytes(C6)", "Smooth muscle cells(C7)", "Fibroblasts(C8)")

seurat_object@meta.data$"new.cluster.ids_dot" = factor(new.cluster.ids_dot[match(seurat_object$RNA_snn_res.0.4, c(0:8))],
                                                       levels= c("Endothelial cells(C0)", "Endothelial cells(C3)", 
                                                                 "Cardiomyocytes(C1)", "Cardiomyocytes(C5)", 
                                                                 "Macrophages(C2)", 
                                                                 "Fibroblasts(C4)", "Fibroblasts(C8)", 
                                                                 "Pericytes(C6)", "Smooth muscle cells(C7)"))

Idents(seurat_object) = "new.cluster.ids_dot"
marker_list = c("Flt1", "Tspan7",
                "Ttn", "Mb",
                "F13a1", "Ptprc",
                "Col8a1", "Matn2",
                "Pdgfrb",
                "Myh11", "Acta2",
                "")
output = paste0(gsub("RNA.cluster.rds", "", infile), reso, ".dotplot.pdf")
pdf(output, width=7, height=3)
DotPlot(seurat_object, 
        features = marker_list)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()


new.cluster.ids_umap <- c("Endothelial cells", "Cardiomyocytes", "Macrophages", "Endothelial cells", 
                          "Fibroblasts", "Cardiomyocytes", "Pericytes", "Smooth muscle cells", "Fibroblasts")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$RNA_snn_res.0.4, c(0:8))])

output = paste0(gsub("RNA.cluster.rds", "", infile), reso, ".UMAP.pdf")
pdf(output, width=6, height=6)
DimPlot(seurat_object, reduction = "umap", group.by = c("new.cluster.ids_umap"), label=TRUE) & theme(legend.position = "bottom") & labs(x="UMAP 1", y="UMAP 2", title="")
dev.off()















infile="tissue_sep/mouse.MCA.RNA.cluster.rds"
seurat_object = readRDS(infile)
Idents(seurat_object) <- "RNA_snn_res.0.3"


FeaturePlot(seurat_object, c("Snhg11")) & theme(legend.position = "bottom")


marker_list = c("Snhg11", "Gad1", "Gad2", "Slc17a6",
                "Mbp", "St18",
                "Slc1a2",
                "Ranbp3l",
                "Bnc2",
                "Cspg4", "Pdgfra",
                "Atp13a5",
                "Tgfbr1", "Ptprc",
                "F13a1",
                "Flt1", 
                "")



new.cluster.ids_dot <- c("Endothelial cells(C0)", "Cardiomyocytes(C1)", "Macrophages(C2)", "Endothelial cells(C3)", 
                         "Fibroblasts(C4)", "Cardiomyocytes(C5)", "(C6)", "Smooth muscle cell(C7)", "(C8)")

seurat_object@meta.data$"new.cluster.ids_dot" = factor(new.cluster.ids_dot[match(seurat_object$RNA_snn_res.0.4, c(0:8))],
                                                       levels= c("Endothelial cells(C0)", "Endothelial cells(C3)", "Cardiomyocytes(C5)", 
                                                                 "Cardiomyocytes(C1)", 
                                                                 "Macrophages(C2)", 
                                                                 "Fibroblasts(C4)", "(C6)", "Smooth muscle cell(C7)", "(C8)"))

Idents(seurat_object) = "new.cluster.ids_dot"
marker_list = c("Flt1", "Pecam1", "Tspan7",
                "Ttn", "Mb",
                "F13a1", "Ptprc",
                "Col8a1", "Dcn",
                "Abcc9",
                "Myh11", "Acta2",
                "Csmd1",
                "")
DotPlot(seurat_object, 
        features = marker_list)+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1))


new.cluster.ids_umap <- c("Endothelial cells", "Cardiomyocytes", "Macrophages", "Endothelial cells", 
                          "Fibroblasts", "Cardiomyocytes", "C6", "Smooth muscle cell", "C8")

seurat_object@meta.data$"new.cluster.ids_umap" = factor(new.cluster.ids_umap[match(seurat_object$RNA_snn_res.0.4, c(0:8))])
DimPlot(seurat_object, reduction = "umap", group.by = c("new.cluster.ids_umap"), label=TRUE) & theme(legend.position = "bottom")

