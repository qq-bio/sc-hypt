library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(colorspace)
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size) +
            theme(
              text = element_text(family = "Arial")
            ))

################################################################################
### variables
cell_order = c(
  c("Inhibitory neuron", "Excitatory neuron", "Avp+ neuron", 
    "Astrocyte", "OPC", "NFO", "Premyelinating OL", "Myelinating OL", 
    "Tanycyte", "Ependymal cell", "Pars tuberalis cell"),
  c("CM"),
  c("POD", "PT", "TL", "TAL", "DCT", "CT", "CD", "IC"),
  c("EC", "E/P transition cell", "Pericyte", "VSMC", "Fibroblast", "Adipocyte"),
  c("Microglia", "Activated microglia", "Immune cell"),
  c("T cells", "NK cells", "B cells", "Plasmablasts", "Plasma cells", 
    "Monocytes", "Neutrophils", "Mast cells", "Erythrocytes", "Megakaryocytes"),
  c("Neuronal")
)

getPalette = colorRampPalette(brewer.pal(12, "Paired"))
cell_col = getPalette(length(cell_order))
names(cell_col) = cell_order

strain_order = c("C57BL/6", "SHR", "WKY", "SS", "SD")
strain_col <- setNames(pal_jama()(5), c("C57BL/6", "SHR", "WKY", "SS", "SD"))

class_order = c("neurons", "glial cells", "muscle cells", "epithelial cells", "endothelial cells", 
                "stromal cells", "immune cells", "adipocytes", "endocrine cells", "blood-related cells")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
class_col = getPalette(length(class_order))
names(class_col) = class_order

tissue_order = c("HYP", "MCA", "LV", "LK", "MSA", "PBMC")
getPalette = colorRampPalette(brewer.pal(12, "Paired"))
tissue_col = getPalette(length(tissue_order))
names(tissue_col) = tissue_order

model_order = c("AngII", "Salt-sensitive", "Spontaneous"); names(model_order)=c("mouse", "rat.ss", "rat.sp")

treatment_order = c("Saline 3d", "AngII 3d", "AngII 28d", "LS", "HS 3d", "HS 21d", "10w", "26w")

sxt_order = c("C57BL/6-Saline 3d", "C57BL/6-AngII 3d", "C57BL/6-AngII 28d", 
              "SS-LS", "SS-HS 3d", "SS-HS 21d", "SD-LS", "SD-HS 3d",
              "SHR-10w", "SHR-26w", "WKY-10w", "WKY-26w")

sxt_colors = c(
  lighten(strain_col["C57BL/6"], 0.35), strain_col["C57BL/6"], darken(strain_col["C57BL/6"], 0.25),
  lighten(strain_col["SS"], 0.35), strain_col["SS"], darken(strain_col["SS"], 0.25),
  lighten(strain_col["SD"], 0.35), strain_col["SD"],
  strain_col["SHR"], darken(strain_col["SHR"], 0.25),
  lighten(strain_col["WKY"], 0.35), strain_col["WKY"]
)

names(sxt_colors) = sxt_order



################################################################################
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

