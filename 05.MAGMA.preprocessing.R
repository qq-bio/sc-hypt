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

### .Renviron
# gh::gh_token()

storage_dir = "/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA"



# munge GWAS
input_file = c(
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-b-12493.elsworth.ht.vcf",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.vcf",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.txt.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018952_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000064_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90081638_buildGRCh38.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.txt.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.tsv"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BP-ICE_PA_HTN_15-04-2020.txt",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079962_buildGRCh38.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018812_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90080024_buildGRCh38.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST011365_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079971_buildGRCh38.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079970_buildGRCh38.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018759_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90100220_buildGRCh37.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29212778-GCST005194-EFO_0000378-build37.f.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90278639.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132314_buildGRCh37.tsv",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.tsv.gz"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435706.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435714.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018832_buildGRCh37.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079903_buildGRCh38.tsv.gz"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018890_buildGRCh37.tsv.gz"
               )


for(i in input_file){

  # gwas_sumstats_path = input_file[length(input_file)]
  gwas_sumstats_path = i
  save_path = gsub("vcf|txt.gz|tsv.gz|tsv|zip|txt", "formatted.tsv.gz", gwas_sumstats_path)

  tryCatch({
    path_formatted = MungeSumstats::format_sumstats(path = gwas_sumstats_path,
                                                    save_path = save_path,
                                                    ref_genome = "GRCh38")
  }, error = function(e) {
    message(paste("Error processing file:", gwas_sumstats_path))
    message("Error message:", e$message)
  })

}


# map snps to genes
input_file = c(
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-b-12493.elsworth.ht.formatted.tsv.gz", # 463010
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/ukb-a-531.neale.ht.formatted.tsv.gz", # 337199
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018952_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000064_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90081638_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv.gz"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.formatted.tsv.gz", # need more resources
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.formatted.tsv.gz"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BP-ICE_PA_HTN_15-04-2020.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018812_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST011365_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018759_buildGRCh37.formatted.tsv.gz",

  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079962_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90080024_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079971_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079970_buildGRCh38.formatted.tsv.gz"
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.formatted.tsv.gz",
  
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132314_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29212778-GCST005194-EFO_0000378-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90278639.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.formatted.tsv.gz"
  
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018832_buildGRCh37.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435706.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079903_buildGRCh38.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90435714.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018890_buildGRCh37.formatted.tsv.gz"
               )

input_df = data.frame(input_file = input_file,
                   genome_ref = c("eur", "eur", "eur", "eur", "eur"),
                   N = c(585264, 388955, 387930, 398198, 660791))

for(i in 1:nrow(input_df)){
  
  gwas_sumstats_path = input_df$input_file[i]
  genome_ref = input_df$genome_ref[i]
  genome_ref_path = paste0("/xdisk/mliang1/qqiu/reference/MAGMA.Celltyping/g1000_", genome_ref,"/g1000_", genome_ref)
  N = input_df$N[i]
  ### need github token, set in .Renviron
  genesOutPath_intelligence = MAGMA.Celltyping::map_snps_to_genes(
    path_formatted = gwas_sumstats_path,
    genome_ref_path = genome_ref_path,
    population = genome_ref,
    N = N,
    genome_build = "GRCh38")
  
}





# cell type dataset
input_file = c(
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.HYP.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LV.RNA.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.LK.multiomics.anno.rds",
  "/xdisk/mliang1/qqiu/project/multiomics-hypertension/cluster/mouse.MCA.RNA.anno.rds",

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
  sce = fix_bad_mgi_symbols(sce)

  # if(! "subclass_level2" %in% colnames(sce@colData)){
  #   sce$subclass_level2 = sce$subclass_level1
  # }

  input_species = strsplit(basename(i), "\\.")[[1]][1]

  sce_dropped = drop_uninformative_genes(exp=sce,
                                         input_species = input_species,
                                         convert_orths = T,
                                         level2annot=sce$subclass_level1)

  annotLevels = list(level1class=as.character(sce$subclass_level1)
                     # level2class=as.character(sce$subclass_level2)
                     )

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


magma_dirs_list = c(
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132315_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.formatted.tsv.gz",

  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BP-ICE_PA_HTN_15-04-2020.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018812_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST011365_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018759_buildGRCh37.formatted.tsv.gz",
  #
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079962_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90080024_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079971_buildGRCh38.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90079970_buildGRCh38.formatted.tsv.gz"

  # "/xdisk/mliang1/qqiu/data/HT-GWAS/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90104537_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/BUN_overall_ALL_YL_20171017_METAL1_nstud_33.dbgap.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90018948_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/UKB.v2.albuminuria.n382500.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000065_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000066_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90000063_buildGRCh37.formatted.tsv.gz",
  
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008029-EFO_0006336-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008044-EFO_0006335-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/31217584-GCST008036-EFO_0000537-build37.f.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29531354-GCST005841-EFO_1001504.h.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90132314_buildGRCh37.formatted.tsv.gz",
  # "/xdisk/mliang1/qqiu/data/HT-GWAS/29212778-GCST005194-EFO_0000378-build37.f.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/formatted_20180205-MA_overall-ALL-nstud_18-SumMac_400.tbl.rsid.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90278639.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310294.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310295.formatted.tsv.gz",
  "/xdisk/mliang1/qqiu/data/HT-GWAS/GCST90310296.formatted.tsv.gz"
  
)

for(magma_dirs in magma_dirs_list){

  result_merge = c()
  for( i in input_file ){

    dataset = gsub("ctd_(.+).rda", "\\1", basename(i), perl = T)
    ctd = load_rdata(i)

    tryCatch({
      
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
      
    }, error = function(e) {
      message(paste("Error processing file:", i))
    })
    

  }

  outfile = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/", gsub("formatted.tsv.gz", "", basename(magma_dirs)), "res.out")
  write.table(result_merge, outfile, row.names=F, col.names=T, sep='\t')


}



# other visualization method: https://star-protocols.cell.com/protocols/1392
base_font_size = 12
theme_set(theme_classic(base_size = base_font_size))
setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA")

model = c("AngII", "Salt-sensitive", "Spontaneous")
names(model) = c("mouse", "rat.ss", "rat.sp")

# result_file = list.files("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA", "*res.out")
trait_use = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/magma.trait_use.txt", header=T, sep='\t')
### remove the nature-2019 paper & albuminuria paper
trait_use = trait_use[! grepl("Wojcik|Haas", trait_use$Name), ]
result_file = paste0("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/", trait_use$Files, ".res.out")

result_merge = c()
for(i in result_file){

  result_tmp = read.table(i, header = T, sep = '\t')
  
  file = gsub(".res.out", "", basename(i))
  trait = trait_use[trait_use$Files==file, ]$Name

  result_tmp$model = model[gsub("(.*)\\.([A-Z]+)_.*", "\\1", result_tmp$analysis_name, perl = T)]
  result_tmp$model = factor(result_tmp$model, levels=model)
  result_tmp$tissue = gsub("(.*)\\.([A-Z]+)_.*", "\\2", result_tmp$analysis_name, perl = T)
  result_tmp$tissue = factor(result_tmp$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
  result_tmp$Celltype_id = ifelse(result_tmp$Celltype==" E _ P _ t r a n s i t i o n _ c e l l ", 
                                  "E/P transition cell", gsub("_", " ", result_tmp$Celltype_id))
  result_tmp$Celltype_id = factor(result_tmp$Celltype_id, levels=cell_order)
  result_tmp$trait = factor(trait, levels=trait_use$Name)
  
  result_merge = rbind(result_merge, result_tmp)

}

result_merge$study = gsub(", et al.", "", gsub(".*- ", "", result_merge$trait))
result_merge$trait = factor(gsub(" -.*", "", result_merge$trait), 
levels = c("Diastolic BP","Systolic BP","Pulse pressure","Hypertension","Stroke",
           "CAD","eGFR","BUN","Microalbuminuria","Albuminuria"))
p = ggplot(result_merge[result_merge$EnrichmentMode=="Top 10%" & result_merge$FDR<0.05, ], 
           aes(x = model, y = Celltype_id,
               size = BETA, fill = -1 * log10p)) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits=rev) +
  # scale_size_area(limits = c(0, 0.2)) +
  geom_point(shape = 21) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    panel.spacing.y=unit(0.2, "lines"),
    panel.spacing.x=unit(0.1, "lines"),
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(fill="-log10(p-value)", y="", x="") +
  facet_nested(tissue ~ trait + study,
    # rows = vars(tissue),
             # cols = vars(trait),
             scales = "free", space = "free")
print(p)


result_merge_v2 <- result_merge %>%
  filter(EnrichmentMode == "Top 10%" & FDR<0.05) %>%
  group_by(model, tissue, trait, Celltype_id) %>%
  slice_min(order_by = FDR, n = 1) %>%
  ungroup() %>% as.data.frame()
p = ggplot(result_merge_v2, 
           aes(x = model, y = Celltype_id,
               size = BETA, fill = -1 * log10p)) +
  scale_fill_gradient(low = "white", high = "darkred") +
  scale_y_discrete(limits=rev) +
  # scale_size_area(limits = c(0, 0.2)) +
  geom_point(shape = 21) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    panel.spacing.y=unit(0.2, "lines"),
    panel.spacing.x=unit(0.1, "lines"),
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 10),
    legend.text = element_text(size=8),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text.x = element_text(colour = 'black', angle = 90, hjust = 0, margin=margin(l=0))
  ) +
  labs(fill="-log10(p-value)", y="", x="") +
  facet_nested(tissue ~ trait,
               scales = "free", space = "free")
print(p)
ggsave("/xdisk/mliang1/qqiu/project/multiomics-hypertension/figure/magma.dotplot.png", width=769/96, height=743/96, dpi=300)

### traits-wise enrichment analysis
trait_celltypes <- result_merge %>%
  filter(EnrichmentMode == "Top 10%") %>%
  group_by(model, tissue, trait, Celltype_id) %>%
  slice_min(order_by = FDR, n = 1) %>%
  ungroup() %>%
  group_by(model, trait) %>%
  summarise(
    all_celltypes = list(unique(paste(Celltype_id, tissue, sep = "_"))),
    sig_celltypes = list(paste(Celltype_id[FDR < 0.05], tissue[FDR < 0.05], sep = "_"))
  ) %>%
  ungroup() %>%
  mutate(trait = as.character(trait))

calculate_p_value <- function(overlap_count, total_count, n_total) {
  return(phyper(overlap_count - 1, total_count, n_total - total_count, total_count, lower.tail = FALSE))
}

overlap_results <- trait_celltypes %>%
  group_by(model) %>%
  summarise(trait_pairs = list(as.data.frame(t(combn(unique(trait), 2)))), .groups = 'drop') %>%
  unnest(trait_pairs) %>%
  rename(trait1 = V1, trait2 = V2) %>%
  rowwise() %>%
  mutate(
    overlap = list(intersect(
      unlist(trait_celltypes$sig_celltypes[trait_celltypes$trait == trait1 & trait_celltypes$model == model]),
      unlist(trait_celltypes$sig_celltypes[trait_celltypes$trait == trait2 & trait_celltypes$model == model])
    )),
    overlap_count = length(overlap),
    total = list(unique(c(
      unlist(trait_celltypes$sig_celltypes[trait_celltypes$trait == trait1 & trait_celltypes$model == model]),
      unlist(trait_celltypes$sig_celltypes[trait_celltypes$trait == trait2 & trait_celltypes$model == model])
    ))),
    total_count = length(total),
    n_total = length(unique(c(
      unlist(trait_celltypes$all_celltypes[trait_celltypes$trait == trait1 & trait_celltypes$model == model]),
      unlist(trait_celltypes$all_celltypes[trait_celltypes$trait == trait2 & trait_celltypes$model == model])
    ))),
  ) %>%
  mutate(p_value = calculate_p_value(overlap_count, total_count, n_total)) %>%
  dplyr::select(model, trait1, trait2, overlap_count, total_count, p_value)

print(overlap_results)


# Prepare the data
desired_order <- c("Diastolic BP vs Systolic BP", "Diastolic BP vs Pulse pressure", 
                   "Systolic BP vs Pulse pressure", "Diastolic BP vs CAD", 
                   "Systolic BP vs CAD", "Pulse pressure vs CAD")
heatmap_data <- overlap_results %>%
  filter(p_value < 0.05) %>%
  mutate(pair = paste(trait1, trait2, sep = " vs "),
         log_p_value = -log10(p_value)) %>%
  dplyr::select(model, pair, log_p_value) %>%
  spread(model, log_p_value, fill = 0) %>%
  gather(model, log_p_value, -pair) %>%
  mutate(log_p_value = ifelse(log_p_value == 0, NA, log_p_value)) %>%  # Optional: Convert zeros to NA for better visualization
  mutate(pair = factor(pair, levels = desired_order))

# Create the heatmap
ggplot(heatmap_data, aes(x = model, y = pair, fill = log_p_value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "grey90") +
  theme_minimal() +
  theme(axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_discrete(limits = rev(levels(heatmap_data$pair))) + # Reverse the y-axis order
  labs(# title = "Hypergeometric test of trait pairs",
       x = "Model",
       y = "Trait Pairs",
       fill = "-log10(p-value)")






result_merge = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/20171016_MW_eGFR_overall_ALL_nstud61.dbgap.res.out", header = T, sep = '\t')
result_merge = read.table("/xdisk/mliang1/qqiu/project/multiomics-hypertension/MAGMA/ukb-b-12493.elsworth.ht.res.out", header = T, sep = '\t')

result_merge$model = model[gsub("(.*)\\.([A-Z]+)_.*", "\\1", result_merge$analysis_name, perl = T)]
result_merge$model = factor(result_merge$model, levels=model)
result_merge$tissue = gsub("(.*)\\.([A-Z]+)_.*", "\\2", result_merge$analysis_name, perl = T)
result_merge$tissue = factor(result_merge$tissue, levels=c("HYP", "MCA", "LV", "LK", "MSA"))
result_merge$Celltype_id = gsub("_", " ", result_merge$Celltype_id)
result_merge$Celltype_id = factor(result_merge$Celltype_id, levels=c(unique(cell_order), "IMM"))

ggplot(result_merge[result_merge$EnrichmentMode=="Top 10%" & result_merge$FDR<0.1, ], aes(x = model, y = Celltype_id,
                   size = BETA, fill = -1 * log10p)) +
  scale_fill_gradient(low = "white", limits = c(0, 5), high = "darkred") +
  scale_y_discrete(limits=rev) +
  scale_size_area(limits = c(0, 0.2)) +
  geom_point(shape = 21) +
  theme(
    panel.grid.major.y = element_blank(),   # No horizontal grid lines
    # legend.position = c(1, 0.55),           # Put legend inside plot area
    legend.justification = c(1, 0.5),
    axis.text.y = element_text(colour = 'black'),
    axis.text.x = element_text(angle = 45, hjust = 1, colour = 'black'),
    # text = element_text(colour = 'black'),
    strip.text = element_text(colour = 'black')
  ) +
  labs(fill="-log10(p-value)", y="", x="") +
  facet_grid(rows = vars(tissue),
             scales = "free", space = "free_y")




heat <- MAGMA.Celltyping::results_heatmap(
  merged_results = merged_results,
  title = "Alzheimer's Disease (ieu-a-298) vs. nervous system cell-types (Zeisel2015)",
  fdr_thresh = 1)

top_enrich <- merged_results %>%
  dplyr::group_by(EnrichmentMode, GWAS) %>%
  dplyr::slice_min(FDR, n = 2)
knitr::kable(top_enrich)









