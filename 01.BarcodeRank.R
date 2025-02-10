library(DropletUtils)

setwd("/xdisk/mliang1/qqiu/project/multiomics-hypertension/cellranger/")


sample_list = c( paste0("MLK", c(1:6)),  paste0("RLK", c(1:7, 10:13)), "RLK82", "RLK92",
                 paste0("RLKS", c(1:4)), paste0("RLKW", c(1:4)))


for(sample in sample_list){
  file = paste0(sample, "/outs/", "raw_feature_bc_matrix.h5")
  dat = read10xCounts(file, sample)
  
  barcode_rank = barcodeRanks(counts(dat))
  
  # p = plot(barcode_rank$rank, barcode_rank$total, log="xy", xlab="Rank", ylab="Total", main=sample)
  
  assign(paste0(sample, "_rank"), barcode_rank)
  print(sample)
}

save.image(file='LK.barcode_rank.RData')



sample="MLK1"; cr=13344; xlim=c(10000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="MLK2"; cr=11922; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 18000)
sample="MLK3"; cr=14107; xlim=c(10000, 50000); ylim=c(1, 10000); v=c(cr, 19000)
sample="MLK4"; cr=17966; xlim=c(15000, 50000); ylim=c(1, 10000); v=c(cr, 25000)
sample="MLK5"; cr=13808; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 22000)
sample="MLK6"; cr=14226; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 23000)
sample="RLK1"; cr=9068; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="RLK2"; cr=11060; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="RLK3"; cr=9041; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 16000)
sample="RLK4"; cr=7798; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 22000)
sample="RLK5"; cr=2284; xlim=c(2000, 50000); ylim=c(1, 10000); v=c(cr, 8000)
sample="RLK6"; cr=13677; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="RLK7"; cr=16265; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 25000)
sample="RLK82"; cr=18525; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 33000)
sample="RLK92"; cr=17832; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 32000)
sample="RLK10"; cr=13151; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 26000)
sample="RLK11"; cr=14466; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 22000)
sample="RLK12"; cr=15119; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 24000)
sample="RLK13"; cr=19088; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 32000)
sample="RLKS1"; cr=6432; xlim=c(4000, 50000); ylim=c(1, 10000); v=c(cr, 12000)
sample="RLKS2"; cr=10125; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 28000)
sample="RLKS3"; cr=20000; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 33000)
sample="RLKS4"; cr=20000; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 39000)
sample="RLKW1"; cr=6894; xlim=c(5000, 50000); ylim=c(1, 10000); v=c(cr, 16000)
sample="RLKW2"; cr=7734; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="RLKW3"; cr=13518; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 20000)
sample="RLKW4"; cr=9956; xlim=c(8000, 50000); ylim=c(1, 10000); v=c(cr, 16000)

barcode_rank = get(paste0(sample, "_rank"))
# plot(barcode_rank$rank, barcode_rank$total, log="xy", xlab="Rank", ylab="Total", main=sample)
plot(barcode_rank$rank, barcode_rank$total, log="xy", xlab="Rank", ylab="Total", main=sample,
     xlim=xlim, ylim=ylim) + abline(v=v, col="red")
