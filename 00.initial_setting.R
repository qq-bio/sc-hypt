dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
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

ec_order <- c("C15","C19","C12","C23","C21",
              "C20","C14",
              "M0610","M24","M5813","C1","C7","C9","C18","C22","C16",
              "C3","C11","C17")

ec_colors <- c(
  "C15" = "#D73027",
  "C19" = "#F46D43",
  "C12" = "#FDAE61",
  "C23" = "#FDD49E",
  "C21" = "#FEC44F",
  
  "C20" = "#4575B4",
  "C14" = "#74ADD1",
  
  "M0610" = "#66C2A5",
  "M24"   = "#41B6C4",
  "M5813" = "#1D91C0",
  "C1"    = "#A6DBA0",
  "C7"    = "#74C476",
  "C9"    = "#31A354",
  "C18"   = "#006D2C",
  "C22"   = "#B8E186",
  "C16"   = "#7FC97F",
  
  "C3"  = "#B3B3B3",
  "C11" = "#9E9AC8",
  "C17" = "#807DBA"
)

blank_theme <- theme(
  axis.line = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks = element_blank(),
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  legend.position = "none",
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.background = element_blank()
)

################################################################################
### functions
# png <- function(filename, width=480, height=480, dpi=300, ...) {
#   graphics::png(filename, width=width, height=height, units="px", res=dpi, ...)
# }
