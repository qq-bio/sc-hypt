.libPaths("/xdisk/mliang1/qqiu/R/library_Rv4.2")

dyn.load("/opt/ohpc/pub/apps/gdal/3.3.2/lib/libgdal.so.29")
dyn.load("/opt/ohpc/pub/apps/proj/7.2.1/lib/libproj.so.19")
dyn.load("/opt/ohpc/pub/libs/gnu8/hdf5/1.10.5/lib/libhdf5_hl.so.100")

library(dplyr)
library(Seurat)
library(MAGMA.Celltyping)
library(ewceData)
library(EWCE)
library(SingleCellExperiment)

# gh::gh_token()

storage_dir = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA"


# 
# # munge GWAS
# input_file = c(
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-b-12493.elsworth.ht.vcf",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.vcf",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018952_buildGRCh37.tsv.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000064_buildGRCh37.tsv",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90081638_buildGRCh38.tsv.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.tsv",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.tsv.gz"
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.tsv",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.tsv",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.txt.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.tsv.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.tsv"
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/BP-ICE_PA_HTN_15-04-2020.txt",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079962_buildGRCh38.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.tsv",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018812_buildGRCh37.tsv.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90080024_buildGRCh38.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST011365_buildGRCh37.tsv",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079971_buildGRCh38.tsv.gz",
#   # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079970_buildGRCh38.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018759_buildGRCh37.tsv.gz"
#                )
# 
# 
# for(i in input_file){
#   
#   # gwas_sumstats_path = input_file[length(input_file)]
#   gwas_sumstats_path = i
#   save_path = gsub("vcf|txt.gz|tsv.gz|tsv|zip|txt", "formatted.tsv.gz", gwas_sumstats_path)
#   path_formatted = MungeSumstats::format_sumstats(path=gwas_sumstats_path,
#                                                   save_path = save_path,
#                                                   ref_genome ="GRCh37")
#   
# }
# 
# 
# # map snps to genes
# input_file = c(#"/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-b-12493.elsworth.ht.formatted.tsv.gz", # 463010
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.formatted.tsv.gz", # 337199
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.formatted.tsv.gz",
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018952_buildGRCh37.formatted.tsv.gz",
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000064_buildGRCh37.formatted.tsv.gz",
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90081638_buildGRCh38.formatted.tsv.gz",
# #                "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv.gz"
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.formatted.tsv.gz", # need more resources
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.formatted.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.formatted.tsv.gz",
#   "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.formatted.tsv.gz"
#                )
# # gwas_sumstats_path = input_file[length(input_file)]
# gwas_sumstats_path = i
# genesOutPath_intelligence = MAGMA.Celltyping::map_snps_to_genes(
#   path_formatted = gwas_sumstats_path,
#   population = c("eur", "eas"),
#   genome_build = "GRCh37")
# 
# 
# 


################################################################################
# cell type dataset
input_file = c(
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  # "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.ss.MCA.RNA.anno.rds",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MSA.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/rat.sp.MCA.RNA.anno.rds"
)

for(i in input_file){
  
  dataset = gsub("\\.(RNA|multiomics)+.anno(.v2)*.rds", "", basename(i), perl = T)
  seurat_object = readRDS(i)
  sce = as.SingleCellExperiment(seurat_object, assay="RNA")
  
  ### convert gene symbol from rat to mouse
  
  sce = fix_bad_mgi_symbols(sce)
  
  sce$subclass_level2 = paste(sce$strain, sce$subclass_level1, sep="-")
  
  input_species = strsplit(basename(i), "\\.")[[1]][1]
  
  sce_dropped = drop_uninformative_genes(exp=sce,
                                         input_species = input_species,
                                         convert_orths = T,
                                         level2annot=sce$subclass_level2)
  
  annotLevels = list(# level1class=as.character(sce$subclass_level1),
                     level2class=as.character(sce$subclass_level2))
  
  ctd = generate_celltype_data(exp=sce_dropped,
                               annotLevels=annotLevels,
                               groupName=dataset,
                               savePath=storage_dir) 
  
}








# https://neurogenomics.github.io/MAGMA_Celltyping/articles/full_workflow.html

input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_mouse.HYP.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_mouse.LK.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_mouse.LV.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_mouse.MCA.rda",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.sp.HYP.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.sp.LK.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.sp.LV.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.sp.MCA.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.sp.MSA.rda",
  
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.ss.HYP.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.ss.LK.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.ss.LV.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.ss.MCA.rda",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ctd_rat.ss.MSA.rda"
)


magma_dirs = "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.formatted.tsv.gz"
magma_dirs = "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-b-12493.elsworth.ht.formatted.tsv.gz"
magma_dirs = "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.formatted.tsv.gz"
magma_dirs = "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.formatted.tsv.gz"
magma_dirs = "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv.gz"

magma_dirs_list = c("/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
"/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.formatted.tsv.gz",
"/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.formatted.tsv.gz",
"/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.formatted.tsv.gz",
"/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.formatted.tsv.gz")

for(magma_dirs in magma_dirs_list){
  
  result_merge = c()
  for( i in input_file ){
    
    dataset = gsub("ctd_(.+).rda", "\\1", basename(i), perl = T)
    ctd = load_rdata(i)
    
    MAGMA_results <- MAGMA.Celltyping::celltype_associations_pipeline(
      magma_dirs = magma_dirs,
      ctd = ctd,
      force_new = T,
      ctd_name = dataset,
      run_linear = TRUE,
      run_top10 = TRUE,
      save_dir = storage_dir)
    
    merged_results <- MAGMA.Celltyping::merge_results(
      MAGMA_results = MAGMA_results)
    
    result_merge = rbind(result_merge, merged_results)
    
  }
  
  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/", gsub("formatted.tsv.gz", "", basename(magma_dirs)), "res.out")
  write.table(result_merge, outfile, row.names=F, col.names=T, sep='\t')
  
  
}

# 
# 
# # other visualization method: https://star-protocols.cell.com/protocols/1392
# base_font_size = 12
# theme_set(theme_classic(base_size = base_font_size))
# 
# model = c("AngII", "Salt-sensitive", "Spontaneous")
# names(model) = c("mouse", "rat.ss", "rat.sp")
# 
# result_merge = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.res.out", header = T, sep = '\t')
# result_merge = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ukb-b-12493.elsworth.ht.res.out", header = T, sep = '\t')
# 
# result_merge$model = model[gsub("(.*)\\.([A-Z]+)_.*", "\\1", result_merge$analysis_name, perl = T)]
# result_merge$model = factor(result_merge$model, levels=model)
# result_merge$tissue = gsub("(.*)\\.([A-Z]+)_.*", "\\2", result_merge$analysis_name, perl = T)
# result_merge$tissue = factor(result_merge$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
# result_merge$Celltype_id = gsub("_", " ", result_merge$Celltype_id)
# result_merge$Celltype_id = factor(result_merge$Celltype_id, levels=c(unique(cell_order), "IMM"))
# 
# ggplot(result_merge[result_merge$EnrichmentMode=="Top 10%" & result_merge$FDR<0.1, ], aes(x = model, y = Celltype_id,
#                    size = BETA, fill = -1 * log10p)) +
#   scale_fill_gradient(low = "white", limits = c(0, 5), high = "darkred") +
#   scale_y_discrete(limits=rev) +
#   scale_size_area(limits = c(0, 0.2)) +
#   geom_point(shape = 21) +
#   theme(
#     panel.grid.major.y = element_blank(),   # No horizontal grid lines
#     # legend.position = c(1, 0.55),           # Put legend inside plot area
#     legend.justification = c(1, 0.5),
#     axis.text.y = element_text(colour = 'black'),
#     axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
#     # text = element_text(colour = 'black'),
#     strip.text = element_text(colour = 'black')
#   ) +
#   labs(fill="-log10(p-value)", y="", x="") +
#   facet_grid(rows = vars(tissue), 
#              scales = "free", space = "free_y")
# 
# 
# 
# 
# heat <- MAGMA.Celltyping::results_heatmap(
#   merged_results = merged_results, 
#   title = "Alzheimer's Disease (ieu-a-298) vs. nervous system cell-types (Zeisel2015)",
#   fdr_thresh = 1)
# 
# top_enrich <- merged_results %>% 
#   dplyr::group_by(EnrichmentMode, GWAS) %>%
#   dplyr::slice_min(FDR, n = 2)
# knitr::kable(top_enrich) 
# 
# 







