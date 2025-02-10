


################################################################################
#### fine-tuned marker list

### mouse-HYP
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", "Avp",
                "Slc1a2", "Cx3cr1", "P2ry12", "Tgfbr1", "Cspg4", "Pdgfra","Mbp", "St18", 
                "Col23a1", "Tmem212", "Flt1", "Pecam1", "Ebf1",
                "Atp13a5", "Ptgds")


### ss-HYP
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Tgfbr1", "Arhgap15", "Lyn", "Cd74",
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Lama2")



### sp-HYP
marker_list = c("Syt1", "Gad1", "Gad2", "Slc17a6", 
                "Slc1a2", 
                "Tgfbr1", "Arhgap15", "Lyn", 
                "Cspg4", "Pdgfra", "Fyn", "Enpp6", "Mbp", "Plp1", 
                "Col23a1", 
                "Cga", 
                "Mecom", "Ebf1", "Atp13a5", "Lama2", "Ptgds")



### mouse-LV
marker_list = c("Flt1", "Pecam1", "Fhl2", "Ttn", "Ptprc", "Dcn", "Pdgfrb", "Rgs5")


### ss-LV
marker_list = c("Vwf", "Pecam1", "Dcn", "Fhl2", "Ttn", "Pdgfrb", "Rgs5", "Ptprc", "Flt4", "Prox1", "Nrxn1")


### sp-LV
marker_list = c("Vwf", "Pecam1", "Dcn", "Fhl2", "Ttn", "Pdgfrb", "Rgs5", "Ptprc", "Nrxn1")




marker_list = c("Il3ra", "Thbd", "Cd1c", 
                "Fcgr3a", "Cd14",
                "Pecam1", "Dcn", "Fhl2", "Ttn", "Pdgfrb", "Rgs5", "Ptprc", "Nrxn1")







### immune cells
# c("Microglia", "Activated microglia", "Monocytes", "Macrophages", "DC", "Neutrophils",
#   "NK cells", "NKT", "T cells", "B cells")

markers <- c(
  c("Cx3cr1", "Tmem119", "P2ry12"), # Microglia
  c("Aif1", "Cd68", "Ccl3"), # Activated microglia
  c("Ly6c2", "Ccr2", "Sell"), # Monocytes
  c("Adgre1", "Mrc1", "Cd163"), # Macrophages
  c("Itgax", "Flt3", "Cd74"), # DC
  c("S100a8", "Ly6g", "Mpo"), # Neutrophils
  c("Ncr1", "Klrb1c", "Gzmb"), # NK cells
  c("Cd3e", "Klrb1c", "Ccl5"), # NKT
  c("Cd3d", "Cd4", "Cd8a"), # T cells
  c("Cd19", "Ms4a1", "Cd79a") # B cells
)




################################################################################
#### other backups
### mouse-LV-immune
gene_list = c("Ptprc", "Mitf", "Mctp1", "Ophn1", "Diaph2",
              "Ptprg", "Cd36", "Myl2", "Myh6", "Fgr", "Gpr141", "Ebf1", "Pax5", "Gm2682", "Skap1")
### ss-LV-immune
gene_list = c("Ptprc", "Alcam", "Rbpj", "F13a1", "Cd163", "Skap1", "Ets1", "Ebf1", "Pax5",
              "Bank1", "Cd74", "Mki67", "Top2a", "Cxcr2", "Ttn")









#### from reference
### HYP
## cell-2018-endo
marker_list = c("Cldn5", "Adgrf5", "Emcn", "Vtn", "Cspg4", "Atp13a5", "Pth1r", "Kcnj8", 
                "Abcc9", "Apln", "Cd82", "Chst1", "Tagln", "Pln", "Bmx", "Gkn3", "Lgfbp2",
                "Dcn", "Lum", "Pdgfra", "Il33", "Ptgds", "Nnat", "Rspo3", "Nov", "Slc47a1",
                "Sox10", "Foxd3", "Aldh1a3", "Anxa11", "Slc18a2", "Klhl30", "Gfra3")


# reference: https://www.frontiersin.org/articles/10.3389/fcell.2021.714169/full
marker_list = c("Cspg4", "Pdgfra", "Bmp4", "Gpr17", "Bcas1", "Enpp6", "Plp1",
                "Tmem2", "Prom1", "Cx47", "Mag", "Mobp")




### EC
marker_list = c(# pan-organ
                "Pecam1", "Kdr", "Cldn5", "Emcn", "Cdh5", "Tie1", "Egfl7",
                "Apj", "Nr2f2", "Vwf",
                "Gja4", "Mecom",
                # heart-adult
                "Wt1", "Slc28a2", "Eepd1", "Kcna5", "Car8", "Fbln1", "Meox2", "Rftn1", "Lamb1", "Myadm",
                "Apln", "Dll4", "Notch1", "Mecom", "Igfbp3", "Cxcr4", "Cx40",
                # brain
                "Slco1c1", "Slco1a4", "Slc22a8", "Mfsd2a", "Slc38a3", "Spock2", "Foxf2", "Edn3", "Stra6", "Slc38ak",
                "Gkn3", "Sema3g", "Efnb2",
                "Slc38a5",
                "Slc16a1",
                # kidney
                "Egfl7", "Dram1", "Dkk2", "Esm1", "Igfbp5", "Pbx1", "Boc", "Igfbp3", "Irx3", "Tnfaip2", "Ptpru",
                "Pi16", "Plat", "Ehd3", "Cyb4b1", "Tspan7", "Lpl",
                "Igfbp3", "Plvap", "Npr3",
                "Igf1", "Cryab", "Igfbp7", "Cd36", "Aqp1", "Ifi27l2a"
                )

# cell-2020 heart EC
marker_list = c("Fbln5", "Hey1", "Mecom", 
                "Cxcl12", "Rbp7", 
                "Kdr", "Endou",
                "Vcam1", "Fmo1", "Fmo2",
                "Mgp", "Cfh", "Bgn", "Vwf",
                "Isg15", "Ifit3", "Ifi203", "Ifit1", "Ifi3b", 
                "Col4a2", "Apln", "Sparcl1", "Aplnr", "Trp53i11",
                "Ccl21a", "Prss23", "Lyve1", "Fxyd6", "Cp",
                "Arhgap18", "Mb", "Tnnt2", "Nrp2", "Notch1", "Usp18", "Mki67")


# EC phenotype; https://academic.oup.com/eurheartj/article/40/30/2507/5510786
marker_list = c('Arhgap18', 'Adm', 'Hspb1', 'CD36',
  'Ifit1', 'Ifit2', 'Ifit3', 'Ifit3b', 'Usp18', 'Cxcl10',
  'Myl2', 'Mb', 'Myl3', 'Tnnt2', 'Tnni3', 'Actc1',
  'Klra3', 'Klra9', 'Klra10',
  'Dll4', 'Notch1', 'Hey1', 'Jag1', 'Gja4',
  'Plvap', 'Lrg1', 'Pbp1', 'Bgn', 'Vwf',
  'Ackr1', 'Ehd4', 'Tmem176a', 'Tmem252', 'Tmem176b', 'Selp',
  'Fbln2', 'Anxa2', 'Col5a2', 'Emilin1', 'Hmcn1', 'Bgn', 'Mgp',
  'Serpina1b', 'Serpina1d', 'Serpina1e',
  'Mki67', 'Top2a', 'Cenpf', 'Cks2', 'Birc5', 'Cenpa', 'Ube2c', 'Cdc20')

marker_list = c("Adam17", "Aph1a", "Crebbp", "Ctbp1", "Ctbp2", "Dll1", "Dll3", "Dll4",
                "Dtx1", "Dtx2", "Dtx3", "Dtx3l", "Dtx4", "Dvl1", "Dvl2", "Dvl3",
                "Ep300", "Hdac1", "Hdac2", "Hes1")

marker_list = c("Pdgfb", "Itga4", "Nf1", "Ddit3", "Grb10", "Enpp1")

marker_list = c("Mecom", "Fbln5", "Slc6a6", "Vegfc", "Mast4", "Nebl", "Cdk19", # artery
                "Mgll", "Cxcl12", "Btnl9", "Unc5b", "Mcf2l", "Slc26a10", "Mcc", "Dlc1", # cap.art
                "Fam117b", "Kdr", "Ncald", "Tnnt2", "Epb41l4a", "Sorbs1",
                "Hmcn1", "Fmo1", "Kcnb1", "Calcr1", 
                "Mgp", "Vwf", 
                "Reln", "Nrp2", "Myl2", 'Mb', "Tnnt2", "Mki67",
                "Flt1", "Vegfa", "Vegfb")

### KPMP-EC
marker_list = c("CD34", 'PECAM1', 'PTPRB', "MEIS2", "FLT1", "EMCN", # EC 
                "EMCN", "HECW2", "PLAT", "ITGA8", "EHD3", "KDR", "SOST", # EC-GC
                "BTNL9", "ADAMTS6", "PALMD", 'AQP1', "TM4SF1", "VEGFC", "CCDC3", "CDH5", "SERPINE2", "FBLN5", "CXCL12", "SOX17", # EC-AEA 
                "BTNL9", "ADAMTS6", "PALMD", "AQP1", "TM4SF1", "MCTP1", "SLC14A1", "ENPP2", 'LYPD6B', # EC-DVR
                "CEACAM1", "DNASE1L3", "PLVAP", "PITPNC1", "GRB10", "SLCO2A1", "RAPGEF4", # EC-PTC
                "CEACAM1", "DNASE1L3", "PLVAP", "GPM6A", "EDIL3", "TLL1", "ZNF385D", "NR2F2", # EC-AVR
                "MMRN1", "CD36", "TBX1", "PKHD1L1", "PROX1") # EC-LYM

marker_list = unique(r2h[r2h$Human.gene.name %in% marker_list, ]$Gene.name)





### score list
##science-2023
cytotoxicity_NK = c("GZMA", "GZMB", "GZMH", "GZMM", "GZMK", "GNLY", "PRF1", "CTSW")
inflammatory_NK = c("CCL2", "CCL3", "CCL4", "CCL5", "CXCL10", "CXCL9", "IL1B", "IL6", "IL7", "IL15", "IL18")
stress_NK = c("BAG3", "CALU", "DNAJB1", "DUSP1", "EGR1", "FOS", "FOSB", "HIF1A", "HSP90AA1", "HSP90AB1", "HSP90B1", "HSPA1A", "HSPA1B", 
              "HSPA6", "HSPB1", "HSPH1", "IER2", "JUN", "JUNB", "NFKBIA", "NFKBIZ", "RGS2", "SLC2A3", "SOCS3", "UBC", "ZFAND2A", "ZFP36", "ZFP36L1")

