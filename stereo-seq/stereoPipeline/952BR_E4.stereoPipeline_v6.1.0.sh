#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=mliang1
#SBATCH --job-name=952BR_E4.stereoPipeline
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=4gb
#SBATCH --time=24:00:00
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/952BR_E4.stereoPipeline\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

# module load singularity
module load python
# module load hdf5/openmpi/gcc/64

dataDir=/xdisk/mliang1/qqiu/data/multiomics-hypertension/STO011/952BR_E4/00.Rawdata

sh /xdisk/mliang1/qqiu/software/stereoPipeline_v6.1.mod.sh \
    -genomeSize 2.75 \
    -splitCount 16 \
    -maskFile /xdisk/mliang1/qqiu/data/multiomics-hypertension/STO011/952BR_E4/SS200000952BR_E4.barcodeToPos.h5 \
    -fq1 $dataDir/barcode_101/FP200007286_L01_101_10.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_11.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_12.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_13.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_14.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_15.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_16.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_1.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_2.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_3.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_4.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_5.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_6.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_7.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_8.fq.gz,$dataDir/barcode_101/FP200007286_L01_101_9.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_10.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_11.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_12.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_13.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_14.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_15.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_16.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_1.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_2.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_3.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_4.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_5.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_6.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_7.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_8.fq.gz,$dataDir/barcode_102/FP200007286_L01_102_9.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_10.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_11.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_12.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_13.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_14.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_15.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_16.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_1.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_2.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_3.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_4.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_5.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_6.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_7.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_8.fq.gz,$dataDir/barcode_98/FP200007286_L01_98_9.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_10.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_11.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_12.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_13.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_14.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_15.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_16.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_1.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_2.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_3.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_4.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_5.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_6.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_7.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_8.fq.gz,$dataDir/barcode_99/FP200007286_L01_99_9.fq.gz \
    -refIndex /xdisk/mliang1/qqiu/reference/SAW/rn7/STAR_SJ100/ \
    -speciesName rat \
    -tissueType kidney \
    -imageRecordFile /xdisk/mliang1/qqiu/data/multiomics-hypertension/STO011/952BR_E4/SS200000952BR_E4__SC_20230118_120516_1.1.4.ipr \
    -imageCompressedFile /xdisk/mliang1/qqiu/data/multiomics-hypertension/STO011/952BR_E4/SS200000952BR_E4__SC_20230118_120516_1.1.4.tar.gz \
    -annotationFile /xdisk/mliang1/qqiu/reference/SAW/rn7/genes/genes.gtf \
    -outDir /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/952BR_E4 \
    -sif /xdisk/mliang1/qqiu/software/SAW_6.1.sif \
    -doCellBin Y \
    -threads 16
