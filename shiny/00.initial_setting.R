# dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
# dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
# dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

# library(Seurat)
# library(dplyr)
# library(ggplot2)
# library(RColorBrewer)
# library(ggsci)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))

################################################################################
### variables
cell_order = c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", # "Microglia", "Activated microglia", 
    "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("EC", "VSMC", "E/P transition cell", "Pericyte", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
    "NK cells", "NKT", "T cells", "B cells")
)

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
cell_col = getPalette(length(cell_order))
names(cell_col) = cell_order

strain_order = c("C57BL/6", "SHR", "WKY", "SS", "SD")
species_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))

class_order = c("neurons", "glial cells", "muscle cells", "epithelial cells", "endothelial cells", "stromal cells", "immune cells", "adipocytes", "endocrine cells")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
class_col = getPalette(length(class_order))
names(class_col) = class_order

tissue_order = c("HYP", "MCA", "LV", "LK", "MSA")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
tissue_col = getPalette(length(tissue_order))
names(tissue_col) = tissue_order

model_order = c("AngII", "Salt-sensitive", "Spontaneous"); names(model_order)=c("mouse", "rat.ss", "rat.sp")

treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")
sxt_order = c("C57BL/6 - Saline 3d", "C57BL/6 - AngII 3d", "C57BL/6 - AngII 28d", 
              "SS - LS", "SS - HS 3d", "SS - HS 21d", "SD - LS", "SD - HS 3d",
              "SHR - 10w", "SHR - 26w", "WKY - 10w", "WKY - 26w")
################################################################################
### functions
# png <- function(filename, width=480, height=480, dpi=300, ...) {
#   graphics::png(filename, width=width, height=height, units="px", res=dpi, ...)
# }
