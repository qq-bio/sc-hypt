#!/bin/bash
#SBATCH --partition=standard
#SBATCH --account=mliang1
#SBATCH --time=96:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem-per-cpu=5gb
#SBATCH --job-name=952BR_B4.stereoPipeline_6.1
#SBATCH --output=/xdisk/mliang1/qqiu/project/multiomics-hypertension/_log/952BR_B4.stereoPipeline\-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=qqiu@arizona.edu

# module load singularity
module load python
# module load hdf5/openmpi/gcc/64

dataDir=/xdisk/mliang1/qqiu/data/multiomics-hypertension/stereo_seq_062623/A02095A3/00.Rawdata

sh /home/u1/qqiu/software/stereoPipeline_v6.1.mod.sh \
    -genomeSize 2.75 \
    -splitCount 1 \
    -maskFile $dataDir/mask/A02095A3.barcodeToPos.h5 \
    -fq1 $dataDir/reads/E100067636_L01_104_1.fq.gz,$dataDir/reads/E100067636_L01_97_1.fq.gz,$dataDir/reads/E100067636_L01_98_1.fq.gz,$dataDir/reads/E100067636_L01_99_1.fq.gz \
    -fq2 $dataDir/reads/E100067636_L01_104_2.fq.gz,$dataDir/reads/E100067636_L01_97_2.fq.gz,$dataDir/reads/E100067636_L01_98_2.fq.gz,$dataDir/reads/E100067636_L01_99_2.fq.gz \
    -refIndex /xdisk/mliang1/qqiu/reference/SAW/rn7/STAR_SJ100/ \
    -speciesName rat \
    -tissueType kidney \
    -imageRecordFile $dataDir/image/A02095A3_SC_20230523_130654_1.1.4.ipr \
    -imageCompressedFile $dataDir/image/A02095A3_SC_20230523_130654_1.1.4.tar.gz \
    -annotationFile /xdisk/mliang1/qqiu/reference/SAW/rn7/genes/genes.gtf \
    -outDir /xdisk/mliang1/qqiu/project/multiomics-hypertension/stereo_seq/result/952BR_B4 \
    -sif /xdisk/mliang1/qqiu/software/SAW_6.1.sif \
    -doCellBin Y \
    -threads 16
